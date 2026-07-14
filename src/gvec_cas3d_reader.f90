module gvec_cas3d_reader
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_netcdf, only: read_dimension, read_harmonic_component, &
        read_integer_scalar, read_integer_vector, read_optional_real_vector, &
        read_real_scalar, read_real_vector, reader_coordinate_error, &
        reader_data_error, reader_ok, reader_open_error, reader_schema_error
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_full, radial_grid_half
    use netcdf_c_api, only: nc_close_file, nc_enotatt, nc_get_global_text, &
        nc_noerr, nc_open_read
    implicit none
    private

    real(dp), parameter :: coordinate_tolerance = 1.0e-10_dp

    public :: read_gvec_cas3d_file
    public :: reader_coordinate_error
    public :: reader_data_error
    public :: reader_ok
    public :: reader_open_error
    public :: reader_schema_error

contains

    subroutine read_gvec_cas3d_file(filename, equilibrium, info)
        character(len=*), intent(in) :: filename
        type(gvec_cas3d_equilibrium_t), intent(out) :: equilibrium
        integer, intent(out) :: info
        integer :: ncid, status

        info = reader_open_error
        status = nc_open_read(filename, ncid)
        if (status /= nc_noerr) return

        call read_file_contents(ncid, equilibrium, info)
        status = nc_close_file(ncid)
        if (status /= nc_noerr) then
            if (info == reader_ok) info = reader_data_error
        end if
    end subroutine read_gvec_cas3d_file

    subroutine read_file_contents(ncid, equilibrium, info)
        integer, intent(in) :: ncid
        type(gvec_cas3d_equilibrium_t), intent(out) :: equilibrium
        integer, intent(out) :: info
        integer :: radial_count, poloidal_count, toroidal_count

        call read_schema_metadata(ncid, equilibrium%schema_version, info)
        if (info /= reader_ok) return
        call read_dimension(ncid, "s", radial_count, info)
        if (info /= reader_ok) return
        call read_dimension(ncid, "m", poloidal_count, info)
        if (info /= reader_ok) return
        call read_dimension(ncid, "n", toroidal_count, info)
        if (info /= reader_ok) return
        call read_metadata(ncid, equilibrium, info)
        if (info /= reader_ok) return
        call read_coordinates(ncid, radial_count, poloidal_count, &
            toroidal_count, equilibrium, info)
        if (info /= reader_ok) return
        call read_profiles(ncid, radial_count, equilibrium, info)
        if (info /= reader_ok) return
        call read_harmonics(ncid, radial_count, poloidal_count, &
            toroidal_count, equilibrium, info)
    end subroutine read_file_contents

    subroutine read_schema_metadata(ncid, schema_version, info)
        integer, intent(in) :: ncid
        integer, intent(out) :: schema_version, info
        character(len=32) :: schema_name, version_text
        integer :: name_status, version_status

        schema_version = 0
        name_status = nc_get_global_text(ncid, "gliss_schema", schema_name)
        version_status = nc_get_global_text(ncid, "gliss_schema_version", &
            version_text)
        if (name_status == nc_enotatt .and. version_status == nc_enotatt) then
            info = reader_ok
            return
        end if
        info = reader_schema_error
        if (name_status /= nc_noerr .or. version_status /= nc_noerr) return
        if (trim(schema_name) /= "gvec-cas3d-export") return
        if (trim(version_text) /= "1") return
        schema_version = 1
        info = reader_ok
    end subroutine read_schema_metadata

    subroutine read_metadata(ncid, equilibrium, info)
        integer, intent(in) :: ncid
        type(gvec_cas3d_equilibrium_t), intent(inout) :: equilibrium
        integer, intent(out) :: info
        character(len=16) :: symmetry
        integer :: status

        symmetry = ""
        call read_integer_scalar(ncid, "N_FP", equilibrium%field_periods, info)
        if (info /= reader_ok) return
        if (equilibrium%field_periods < 1) then
            info = reader_coordinate_error
            return
        end if
        call read_integer_scalar(ncid, "winding", equilibrium%winding, info)
        if (info /= reader_ok) return
        call read_position_frame(ncid, equilibrium%has_boozer_position_frame, &
            info)
        if (info /= reader_ok) return
        call read_real_scalar(ncid, "beta_avg", equilibrium%beta_average, info)
        if (info /= reader_ok) return

        status = nc_get_global_text(ncid, "stellarator_symmetry", symmetry)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        select case (trim(adjustl(symmetry)))
        case ("True")
            equilibrium%stellarator_symmetric = .true.
        case ("False")
            equilibrium%stellarator_symmetric = .false.
        case default
            info = reader_schema_error
            return
        end select
        info = reader_ok
    end subroutine read_metadata

    subroutine read_position_frame(ncid, compatible, info)
        integer, intent(in) :: ncid
        logical, intent(out) :: compatible
        integer, intent(out) :: info
        character(len=64) :: frame
        integer :: status

        compatible = .false.
        frame = ""
        status = nc_get_global_text(ncid, "position_frame", frame)
        if (status == nc_enotatt) then
            info = reader_ok
            return
        end if
        info = reader_schema_error
        if (status /= nc_noerr) return
        if (trim(frame) /= "xhat,yhat rotated by winding*zeta_B") return
        compatible = .true.
        info = reader_ok
    end subroutine read_position_frame

    subroutine read_coordinates(ncid, radial_count, poloidal_count, &
            toroidal_count, equilibrium, info)
        integer, intent(in) :: ncid, radial_count, poloidal_count
        integer, intent(in) :: toroidal_count
        type(gvec_cas3d_equilibrium_t), intent(inout) :: equilibrium
        integer, intent(out) :: info
        real(dp), allocatable :: file_s(:)
        logical :: has_s

        call read_integer_vector(ncid, "m", "m", poloidal_count, &
            equilibrium%poloidal_modes, info)
        if (info /= reader_ok) return
        call read_integer_vector(ncid, "n", "n", toroidal_count, &
            equilibrium%toroidal_modes, info)
        if (info /= reader_ok) return
        call validate_mode_order(equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, info)
        if (info /= reader_ok) return
        call read_real_vector(ncid, "rho", "s", radial_count, &
            equilibrium%rho, info)
        if (info /= reader_ok) return
        call read_optional_real_vector(ncid, "s", "s", radial_count, &
            file_s, has_s, info)
        if (info /= reader_ok) return
        call classify_radial_coordinates(equilibrium%rho, file_s, has_s, &
            equilibrium%s, equilibrium%radial_grid, info)
    end subroutine read_coordinates

    subroutine read_profiles(ncid, radial_count, equilibrium, info)
        integer, intent(in) :: ncid, radial_count
        type(gvec_cas3d_equilibrium_t), intent(inout) :: equilibrium
        integer, intent(out) :: info

        call read_real_vector(ncid, "p", "s", radial_count, &
            equilibrium%pressure, info)
        if (info /= reader_ok) return
        call read_real_vector(ncid, "B_theta_avg", "s", radial_count, &
            equilibrium%b_theta_average, info)
        if (info /= reader_ok) return
        call read_real_vector(ncid, "B_zeta_avg", "s", radial_count, &
            equilibrium%b_zeta_average, info)
        if (info /= reader_ok) return
        call read_real_vector(ncid, "Phi", "s", radial_count, &
            equilibrium%toroidal_flux, info)
        if (info /= reader_ok) return
        call read_real_vector(ncid, "chi", "s", radial_count, &
            equilibrium%poloidal_flux, info)
        if (info /= reader_ok) return
        call read_real_vector(ncid, "iota", "s", radial_count, &
            equilibrium%rotational_transform, info)
    end subroutine read_profiles

    subroutine read_harmonics(ncid, ns, nm, nn, equilibrium, info)
        integer, intent(in) :: ncid, ns, nm, nn
        type(gvec_cas3d_equilibrium_t), intent(inout) :: equilibrium
        integer, intent(out) :: info
        logical :: require_both, found_st, found_sz

        require_both = .not. equilibrium%stellarator_symmetric
        call read_even_pair(ncid, "mod_B", ns, nm, nn, require_both, &
            equilibrium%mod_b, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "xhat", ns, nm, nn, require_both, &
            equilibrium%xhat, info)
        if (info /= reader_ok) return
        call read_odd_pair(ncid, "yhat", ns, nm, nn, require_both, &
            equilibrium%yhat, info)
        if (info /= reader_ok) return
        call read_odd_pair(ncid, "zhat", ns, nm, nn, require_both, &
            equilibrium%zhat, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "Jac", ns, nm, nn, require_both, &
            equilibrium%jacobian, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "g_tt", ns, nm, nn, require_both, &
            equilibrium%g_tt, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "g_tz", ns, nm, nn, require_both, &
            equilibrium%g_tz, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "g_zz", ns, nm, nn, require_both, &
            equilibrium%g_zz, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "II_tt", ns, nm, nn, require_both, &
            equilibrium%second_form_tt, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "II_tz", ns, nm, nn, require_both, &
            equilibrium%second_form_tz, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "II_zz", ns, nm, nn, require_both, &
            equilibrium%second_form_zz, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "B_contra_t", ns, nm, nn, require_both, &
            equilibrium%b_contravariant_theta, info)
        if (info /= reader_ok) return
        call read_even_pair(ncid, "B_contra_z", ns, nm, nn, require_both, &
            equilibrium%b_contravariant_zeta, info)
        if (info /= reader_ok) return
        call read_optional_pair(ncid, "g_st", ns, nm, nn, &
            equilibrium%g_st, found_st, info)
        if (info /= reader_ok) return
        call read_optional_pair(ncid, "g_sz", ns, nm, nn, &
            equilibrium%g_sz, found_sz, info)
        if (info /= reader_ok) return
        equilibrium%has_chart_metric = found_st .and. found_sz
    end subroutine read_harmonics

    subroutine read_optional_pair(ncid, name, ns, nm, nn, pair, found, info)
        integer, intent(in) :: ncid, ns, nm, nn
        character(len=*), intent(in) :: name
        type(harmonic_pair_t), intent(out) :: pair
        logical, intent(out) :: found
        integer, intent(out) :: info
        logical :: has_cosine, has_sine

        found = .false.
        call read_harmonic_component(ncid, trim(name) // "_mnc", ns, nm, nn, &
            pair%cosine, has_cosine, info)
        if (info /= reader_ok) return
        call read_harmonic_component(ncid, trim(name) // "_mns", ns, nm, nn, &
            pair%sine, has_sine, info)
        if (info /= reader_ok) return
        found = has_cosine .or. has_sine
    end subroutine read_optional_pair

    subroutine read_even_pair(ncid, name, ns, nm, nn, require_both, pair, info)
        integer, intent(in) :: ncid, ns, nm, nn
        character(len=*), intent(in) :: name
        logical, intent(in) :: require_both
        type(harmonic_pair_t), intent(out) :: pair
        integer, intent(out) :: info

        call read_harmonic_pair(ncid, name, ns, nm, nn, .true., &
            require_both, pair, info)
    end subroutine read_even_pair

    subroutine read_odd_pair(ncid, name, ns, nm, nn, require_both, pair, info)
        integer, intent(in) :: ncid, ns, nm, nn
        character(len=*), intent(in) :: name
        logical, intent(in) :: require_both
        type(harmonic_pair_t), intent(out) :: pair
        integer, intent(out) :: info

        call read_harmonic_pair(ncid, name, ns, nm, nn, require_both, &
            .true., pair, info)
    end subroutine read_odd_pair

    subroutine read_harmonic_pair(ncid, name, ns, nm, nn, require_cosine, &
            require_sine, pair, info)
        integer, intent(in) :: ncid, ns, nm, nn
        character(len=*), intent(in) :: name
        logical, intent(in) :: require_cosine, require_sine
        type(harmonic_pair_t), intent(out) :: pair
        integer, intent(out) :: info
        logical :: has_cosine, has_sine

        call read_harmonic_component(ncid, trim(name) // "_mnc", ns, nm, nn, &
            pair%cosine, has_cosine, info)
        if (info /= reader_ok) return
        call read_harmonic_component(ncid, trim(name) // "_mns", ns, nm, nn, &
            pair%sine, has_sine, info)
        if (info /= reader_ok) return
        if (require_cosine .and. .not. has_cosine) info = reader_schema_error
        if (require_sine .and. .not. has_sine) info = reader_schema_error
    end subroutine read_harmonic_pair

    pure subroutine validate_mode_order(poloidal_modes, toroidal_modes, info)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        integer, intent(out) :: info
        integer :: index, maximum_toroidal_mode

        info = reader_coordinate_error
        do index = 1, size(poloidal_modes)
            if (poloidal_modes(index) /= index - 1) return
        end do
        if (mod(size(toroidal_modes), 2) /= 1) return
        maximum_toroidal_mode = (size(toroidal_modes) - 1) / 2
        do index = 1, maximum_toroidal_mode + 1
            if (toroidal_modes(index) /= index - 1) return
        end do
        do index = 1, maximum_toroidal_mode
            if (toroidal_modes(maximum_toroidal_mode + 1 + index) /= &
                index - maximum_toroidal_mode - 1) return
        end do
        info = reader_ok
    end subroutine validate_mode_order

    pure subroutine classify_radial_coordinates(rho, file_s, has_s, s, &
            radial_grid, info)
        real(dp), intent(in) :: rho(:), file_s(:)
        logical, intent(in) :: has_s
        real(dp), allocatable, intent(out) :: s(:)
        integer, intent(out) :: radial_grid, info
        real(dp) :: expected
        integer :: index

        radial_grid = 0
        info = reader_coordinate_error
        if (.not. all(ieee_is_finite(rho))) return
        if (any(rho <= 0.0_dp)) return
        if (any(rho > 1.0_dp)) return
        allocate (s(size(rho)))
        s = rho**2
        if (has_s) then
            if (maxval(abs(file_s - s)) > coordinate_tolerance) return
        end if
        do index = 1, size(s)
            expected = (real(index, dp) - 0.5_dp) / real(size(s), dp)
            if (abs(s(index) - expected) > coordinate_tolerance) exit
        end do
        if (index > size(s)) then
            radial_grid = radial_grid_half
            info = reader_ok
            return
        end if
        if (size(s) < 2) return
        if (abs(s(1) - 1.0e-8_dp) > coordinate_tolerance) return
        do index = 2, size(s)
            expected = real(index - 1, dp) / real(size(s) - 1, dp)
            if (abs(s(index) - expected) > coordinate_tolerance) return
        end do
        radial_grid = radial_grid_full
        info = reader_ok
    end subroutine classify_radial_coordinates

end module gvec_cas3d_reader
