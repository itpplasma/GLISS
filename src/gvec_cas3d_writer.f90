module gvec_cas3d_writer
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_full, radial_grid_half
    use netcdf_c_api, only: nc_close_file, nc_create_netcdf4_exclusive, &
        nc_def_dimension, nc_def_scalar, nc_def_variable, nc_double, &
        nc_end_definitions, nc_int64, nc_noerr, nc_put_global_text, &
        nc_put_integer, nc_put_real
    implicit none
    private

    integer, parameter, public :: writer_ok = 0
    integer, parameter, public :: writer_invalid = 1
    integer, parameter, public :: writer_open_error = 2
    integer, parameter, public :: writer_netcdf_error = 3
    integer, parameter :: profile_count = 6, pair_count = 15
    character(len=11), parameter :: profile_names(profile_count) = &
        [character(len=11) :: "p", "B_theta_avg", "B_zeta_avg", "Phi", &
        "chi", "iota"]
    character(len=11), parameter :: pair_names(pair_count) = &
        [character(len=11) :: "mod_B", "xhat", "yhat", "zhat", "Jac", &
        "g_tt", "g_tz", "g_zz", "II_tt", "II_tz", "II_zz", &
        "B_contra_t", "B_contra_z", "g_st", "g_sz"]
    real(dp), parameter :: coordinate_tolerance = 1.0e-10_dp

    type :: writer_ids_t
        integer :: ncid = -1
        integer :: field_periods = -1, beta_average = -1, winding = -1
        integer :: poloidal_modes = -1, toroidal_modes = -1
        integer :: rho = -1, s = -1
        integer :: profiles(profile_count) = -1
        integer :: pairs(2, pair_count) = -1
    end type writer_ids_t

    public :: write_gvec_cas3d_file

contains

    subroutine write_gvec_cas3d_file(filename, equilibrium, info)
        character(len=*), intent(in) :: filename
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(out) :: info
        type(writer_ids_t) :: ids
        integer :: close_status

        info = writer_invalid
        if (.not. valid_equilibrium(equilibrium)) return
        if (nc_create_netcdf4_exclusive(filename, ids%ncid) /= nc_noerr) then
            info = writer_open_error
            return
        end if
        call define_file(equilibrium, ids, info)
        if (info == writer_ok) call write_file(equilibrium, ids, info)
        close_status = nc_close_file(ids%ncid)
        if (info == writer_ok .and. close_status /= nc_noerr) &
            info = writer_netcdf_error
        if (info /= writer_ok) call delete_partial_file(filename)
    end subroutine write_gvec_cas3d_file

    subroutine define_file(equilibrium, ids, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(writer_ids_t), intent(inout) :: ids
        integer, intent(out) :: info
        integer :: dimensions(3), dim_s, dim_m, dim_n, index
        integer :: vector_dimension(1)

        info = writer_netcdf_error
        if (nc_def_dimension(ids%ncid, "s", size(equilibrium%s), dim_s) &
            /= nc_noerr) return
        if (nc_def_dimension(ids%ncid, "m", &
            size(equilibrium%poloidal_modes), dim_m) /= nc_noerr) return
        if (nc_def_dimension(ids%ncid, "n", &
            size(equilibrium%toroidal_modes), dim_n) /= nc_noerr) return
        if (nc_def_scalar(ids%ncid, "N_FP", nc_int64, ids%field_periods) &
            /= nc_noerr) return
        if (nc_def_scalar(ids%ncid, "beta_avg", nc_double, ids%beta_average) &
            /= nc_noerr) return
        if (nc_def_scalar(ids%ncid, "winding", nc_int64, ids%winding) &
            /= nc_noerr) return
        vector_dimension(1) = dim_m
        if (nc_def_variable(ids%ncid, "m", nc_int64, vector_dimension, &
            ids%poloidal_modes) /= nc_noerr) return
        vector_dimension(1) = dim_n
        if (nc_def_variable(ids%ncid, "n", nc_int64, vector_dimension, &
            ids%toroidal_modes) /= nc_noerr) return
        vector_dimension(1) = dim_s
        if (nc_def_variable(ids%ncid, "rho", nc_double, vector_dimension, &
            ids%rho) /= nc_noerr) return
        if (nc_def_variable(ids%ncid, "s", nc_double, vector_dimension, ids%s) &
            /= nc_noerr) return
        do index = 1, profile_count
            if (nc_def_variable(ids%ncid, trim(profile_names(index)), &
                nc_double, vector_dimension, ids%profiles(index)) &
                /= nc_noerr) return
        end do
        dimensions(1) = dim_s
        dimensions(2) = dim_m
        dimensions(3) = dim_n
        do index = 1, pair_count
            if (index > 13 .and. .not. equilibrium%has_chart_metric) cycle
            if (nc_def_variable(ids%ncid, trim(pair_names(index)) // "_mnc", &
                nc_double, dimensions, ids%pairs(1, index)) /= nc_noerr) return
            if (nc_def_variable(ids%ncid, trim(pair_names(index)) // "_mns", &
                nc_double, dimensions, ids%pairs(2, index)) /= nc_noerr) return
        end do
        if (nc_put_global_text(ids%ncid, "gliss_schema", &
            "gvec-cas3d-export") /= nc_noerr) return
        if (nc_put_global_text(ids%ncid, "gliss_schema_version", "1") &
            /= nc_noerr) return
        if (nc_put_global_text(ids%ncid, "stellarator_symmetry", &
            merge("True ", "False", equilibrium%stellarator_symmetric)) &
            /= nc_noerr) return
        if (equilibrium%has_boozer_position_frame) then
            if (nc_put_global_text(ids%ncid, "position_frame", &
                "xhat,yhat rotated by winding*zeta_B") /= nc_noerr) return
        end if
        if (nc_end_definitions(ids%ncid) /= nc_noerr) return
        info = writer_ok
    end subroutine define_file

    subroutine write_file(equilibrium, ids, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(writer_ids_t), intent(in) :: ids
        integer, intent(out) :: info
        integer :: index

        info = writer_netcdf_error
        if (nc_put_integer(ids%ncid, ids%field_periods, &
            equilibrium%field_periods) /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%beta_average, equilibrium%beta_average) &
            /= nc_noerr) return
        if (nc_put_integer(ids%ncid, ids%winding, equilibrium%winding) &
            /= nc_noerr) return
        if (nc_put_integer(ids%ncid, ids%poloidal_modes, &
            equilibrium%poloidal_modes) /= nc_noerr) return
        if (nc_put_integer(ids%ncid, ids%toroidal_modes, &
            equilibrium%toroidal_modes) /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%rho, equilibrium%rho) /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%s, equilibrium%s) /= nc_noerr) return
        call write_profiles(equilibrium, ids, info)
        if (info /= writer_ok) return
        do index = 1, pair_count
            if (ids%pairs(1, index) < 0) cycle
            call write_named_pair(equilibrium, index, ids%ncid, &
                ids%pairs(:, index), info)
            if (info /= writer_ok) return
        end do
        info = writer_ok
    end subroutine write_file

    subroutine write_profiles(equilibrium, ids, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(writer_ids_t), intent(in) :: ids
        integer, intent(out) :: info

        info = writer_netcdf_error
        if (nc_put_real(ids%ncid, ids%profiles(1), equilibrium%pressure) &
            /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%profiles(2), &
            equilibrium%b_theta_average) /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%profiles(3), &
            equilibrium%b_zeta_average) /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%profiles(4), equilibrium%toroidal_flux) &
            /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%profiles(5), equilibrium%poloidal_flux) &
            /= nc_noerr) return
        if (nc_put_real(ids%ncid, ids%profiles(6), &
            equilibrium%rotational_transform) /= nc_noerr) return
        info = writer_ok
    end subroutine write_profiles

    subroutine write_named_pair(equilibrium, index, ncid, ids, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: index, ncid, ids(2)
        integer, intent(out) :: info

        select case (index)
        case (1); call write_pair(ncid, ids, equilibrium%mod_b, info)
        case (2); call write_pair(ncid, ids, equilibrium%xhat, info)
        case (3); call write_pair(ncid, ids, equilibrium%yhat, info)
        case (4); call write_pair(ncid, ids, equilibrium%zhat, info)
        case (5); call write_pair(ncid, ids, equilibrium%jacobian, info)
        case (6); call write_pair(ncid, ids, equilibrium%g_tt, info)
        case (7); call write_pair(ncid, ids, equilibrium%g_tz, info)
        case (8); call write_pair(ncid, ids, equilibrium%g_zz, info)
        case (9); call write_pair(ncid, ids, equilibrium%second_form_tt, info)
        case (10); call write_pair(ncid, ids, equilibrium%second_form_tz, info)
        case (11); call write_pair(ncid, ids, equilibrium%second_form_zz, info)
        case (12)
            call write_pair(ncid, ids, equilibrium%b_contravariant_theta, info)
        case (13)
            call write_pair(ncid, ids, equilibrium%b_contravariant_zeta, info)
        case (14); call write_pair(ncid, ids, equilibrium%g_st, info)
        case (15); call write_pair(ncid, ids, equilibrium%g_sz, info)
        case default; info = writer_invalid
        end select
    end subroutine write_named_pair

    subroutine write_pair(ncid, ids, pair, info)
        integer, intent(in) :: ncid, ids(2)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(out) :: info
        real(dp), allocatable :: file_values(:, :, :)

        call file_order(pair%cosine, file_values)
        if (nc_put_real(ncid, ids(1), file_values) /= nc_noerr) then
            info = writer_netcdf_error
            return
        end if
        call file_order(pair%sine, file_values)
        if (nc_put_real(ncid, ids(2), file_values) /= nc_noerr) then
            info = writer_netcdf_error
            return
        end if
        info = writer_ok
    end subroutine write_pair

    subroutine file_order(values, file_values)
        real(dp), intent(in) :: values(:, :, :)
        real(dp), allocatable, intent(out) :: file_values(:, :, :)
        integer :: radial, poloidal, toroidal

        allocate (file_values(size(values, 3), size(values, 2), size(values, 1)))
        do toroidal = 1, size(values, 3)
            do poloidal = 1, size(values, 2)
                do radial = 1, size(values, 1)
                    file_values(toroidal, poloidal, radial) = &
                        values(radial, poloidal, toroidal)
                end do
            end do
        end do
    end subroutine file_order

    logical function valid_equilibrium(equilibrium) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer :: expected_shape(3), index

        valid = .false.
        if (equilibrium%schema_version < 0 .or. equilibrium%schema_version > 1) return
        if (equilibrium%field_periods < 1) return
        if (.not. ieee_is_finite(equilibrium%beta_average)) return
        if (.not. valid_coordinates(equilibrium)) return
        if (.not. valid_modes(equilibrium)) return
        if (.not. valid_profiles(equilibrium)) return
        expected_shape(1) = size(equilibrium%s)
        expected_shape(2) = size(equilibrium%poloidal_modes)
        expected_shape(3) = size(equilibrium%toroidal_modes)
        do index = 1, 13
            if (.not. valid_named_pair(equilibrium, index, expected_shape)) return
        end do
        if (equilibrium%has_chart_metric) then
            if (.not. valid_pair(equilibrium%g_st, expected_shape)) return
            if (.not. valid_pair(equilibrium%g_sz, expected_shape)) return
        end if
        valid = .true.
    end function valid_equilibrium

    logical function valid_coordinates(equilibrium) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer :: index, count
        real(dp) :: expected

        valid = .false.
        if (.not. allocated(equilibrium%rho) .or. &
            .not. allocated(equilibrium%s)) return
        count = size(equilibrium%s)
        if (count < 1 .or. size(equilibrium%rho) /= count) return
        if (.not. all(ieee_is_finite(equilibrium%rho))) return
        if (.not. all(ieee_is_finite(equilibrium%s))) return
        if (any(equilibrium%rho <= 0.0_dp) .or. &
            any(equilibrium%rho > 1.0_dp)) return
        if (any(equilibrium%rho**2 /= equilibrium%s)) return
        if (equilibrium%radial_grid == radial_grid_half) then
            do index = 1, count
                expected = (real(index, dp) - 0.5_dp) / real(count, dp)
                if (abs(equilibrium%s(index) - expected) &
                    > coordinate_tolerance) return
            end do
        else if (equilibrium%radial_grid == radial_grid_full) then
            if (count < 2) return
            if (abs(equilibrium%s(1) - 1.0e-8_dp) &
                > coordinate_tolerance) return
            do index = 2, count
                expected = real(index - 1, dp) / real(count - 1, dp)
                if (abs(equilibrium%s(index) - expected) &
                    > coordinate_tolerance) return
            end do
        else
            return
        end if
        valid = .true.
    end function valid_coordinates

    logical function valid_modes(equilibrium) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer :: index, maximum

        valid = .false.
        if (.not. allocated(equilibrium%poloidal_modes)) return
        if (.not. allocated(equilibrium%toroidal_modes)) return
        if (size(equilibrium%poloidal_modes) < 1) return
        if (mod(size(equilibrium%toroidal_modes), 2) /= 1) return
        do index = 1, size(equilibrium%poloidal_modes)
            if (equilibrium%poloidal_modes(index) /= index - 1) return
        end do
        maximum = (size(equilibrium%toroidal_modes) - 1) / 2
        do index = 1, maximum + 1
            if (equilibrium%toroidal_modes(index) /= index - 1) return
        end do
        do index = 1, maximum
            if (equilibrium%toroidal_modes(maximum + 1 + index) /= &
                index - maximum - 1) return
        end do
        valid = .true.
    end function valid_modes

    logical function valid_profiles(equilibrium) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer :: count

        valid = .false.
        count = size(equilibrium%s)
        if (.not. valid_vector(equilibrium%pressure, count)) return
        if (.not. valid_vector(equilibrium%b_theta_average, count)) return
        if (.not. valid_vector(equilibrium%b_zeta_average, count)) return
        if (.not. valid_vector(equilibrium%toroidal_flux, count)) return
        if (.not. valid_vector(equilibrium%poloidal_flux, count)) return
        if (.not. valid_vector(equilibrium%rotational_transform, count)) return
        valid = .true.
    end function valid_profiles

    logical function valid_vector(values, expected_size) result(valid)
        real(dp), allocatable, intent(in) :: values(:)
        integer, intent(in) :: expected_size

        valid = allocated(values)
        if (.not. valid) return
        valid = size(values) == expected_size .and. all(ieee_is_finite(values))
    end function valid_vector

    logical function valid_named_pair(equilibrium, index, expected_shape) &
            result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: index, expected_shape(3)

        select case (index)
        case (1); valid = valid_pair(equilibrium%mod_b, expected_shape)
        case (2); valid = valid_pair(equilibrium%xhat, expected_shape)
        case (3); valid = valid_pair(equilibrium%yhat, expected_shape)
        case (4); valid = valid_pair(equilibrium%zhat, expected_shape)
        case (5); valid = valid_pair(equilibrium%jacobian, expected_shape)
        case (6); valid = valid_pair(equilibrium%g_tt, expected_shape)
        case (7); valid = valid_pair(equilibrium%g_tz, expected_shape)
        case (8); valid = valid_pair(equilibrium%g_zz, expected_shape)
        case (9); valid = valid_pair(equilibrium%second_form_tt, expected_shape)
        case (10); valid = valid_pair(equilibrium%second_form_tz, expected_shape)
        case (11); valid = valid_pair(equilibrium%second_form_zz, expected_shape)
        case (12)
            valid = valid_pair(equilibrium%b_contravariant_theta, expected_shape)
        case (13)
            valid = valid_pair(equilibrium%b_contravariant_zeta, expected_shape)
        case default; valid = .false.
        end select
    end function valid_named_pair

    logical function valid_pair(pair, expected_shape) result(valid)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: expected_shape(3)

        valid = allocated(pair%cosine) .and. allocated(pair%sine)
        if (.not. valid) return
        valid = all(shape(pair%cosine) == expected_shape) &
            .and. all(shape(pair%sine) == expected_shape)
        if (.not. valid) return
        valid = all(ieee_is_finite(pair%cosine)) &
            .and. all(ieee_is_finite(pair%sine))
    end function valid_pair

    subroutine delete_partial_file(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, status

        open (newunit=unit, file=filename, status="old", iostat=status)
        if (status == 0) close (unit, status="delete")
    end subroutine delete_partial_file

end module gvec_cas3d_writer
