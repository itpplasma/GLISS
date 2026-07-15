program test_gvec_cas3d_reader
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, &
        reader_coordinate_error, reader_data_error, reader_ok, &
        reader_open_error, reader_schema_error
    use gvec_cas3d_types, only: equilibrium_is_axisymmetric, &
        gvec_cas3d_equilibrium_t, radial_grid_full, radial_grid_half
    use netcdf_c_api, only: nc_close_file, nc_create_netcdf4, &
        nc_def_dimension, nc_def_scalar, nc_def_variable, nc_double, &
        nc_end_definitions, nc_inquire_variable_id, nc_int64, nc_noerr, &
        nc_open_write, nc_put_global_text, nc_put_integer, nc_put_real
    implicit none

    integer, parameter :: field_count = 13
    integer, parameter :: profile_count = 6
    integer, parameter :: ns = 3, nm = 2, nn = 3
    character(len=11), parameter :: field_names(field_count) = &
        [character(len=11) :: "mod_B", "xhat", "yhat", "zhat", "Jac", &
        "g_tt", "g_tz", "g_zz", "II_tt", "II_tz", "II_zz", &
        "B_contra_t", "B_contra_z"]
    character(len=11), parameter :: profile_names(profile_count) = &
        [character(len=11) :: "p", "B_theta_avg", "B_zeta_avg", &
        "Phi", "chi", "iota"]
    character(len=64) :: half_file, full_file, symmetric_file
    character(len=64) :: corrupt_file, nan_file, mode_file, missing_file
    character(len=64) :: schema_file, future_file, incomplete_file
    character(len=64) :: malformed_file, frame_file, bad_frame_file

    interface
        function get_process_id() result(process_id) bind(c, name="getpid")
            import :: c_int
            integer(c_int) :: process_id
        end function get_process_id
    end interface

    type :: fixture_ids_t
        integer :: ncid = -1
        integer :: nfp = -1, beta = -1, winding = -1
        integer :: m = -1, n = -1, rho = -1, s = -1
        integer :: profiles(profile_count) = -1
        integer :: pairs(2, field_count) = -1
    end type fixture_ids_t

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    integer :: info

    write (half_file, '("reader_half_",i0,".nc")') get_process_id()
    write (full_file, '("reader_full_",i0,".nc")') get_process_id()
    write (symmetric_file, '("reader_symmetric_",i0,".nc")') get_process_id()
    write (corrupt_file, '("reader_corrupt_",i0,".nc")') get_process_id()
    write (nan_file, '("reader_nan_",i0,".nc")') get_process_id()
    write (mode_file, '("reader_mode_",i0,".nc")') get_process_id()
    write (missing_file, '("reader_missing_",i0,".nc")') get_process_id()
    write (schema_file, '("reader_schema_",i0,".nc")') get_process_id()
    write (future_file, '("reader_future_",i0,".nc")') get_process_id()
    write (incomplete_file, '("reader_incomplete_",i0,".nc")') get_process_id()
    write (malformed_file, '("reader_malformed_",i0,".nc")') get_process_id()
    write (frame_file, '("reader_frame_",i0,".nc")') get_process_id()
    write (bad_frame_file, '("reader_bad_frame_",i0,".nc")') get_process_id()

    call create_fixture(half_file, radial_grid_half, .false., .false.)
    call read_gvec_cas3d_file(half_file, equilibrium, info)
    call require(info == reader_ok, "half-mesh fixture was rejected")
    call verify_half_fixture(equilibrium)
    equilibrium%toroidal_modes = 0
    call require(equilibrium_is_axisymmetric(equilibrium), &
        "axisymmetric harmonic table was rejected")
    equilibrium%toroidal_modes(2) = 1
    call require(.not. equilibrium_is_axisymmetric(equilibrium), &
        "nonaxisymmetric harmonic table was accepted")

    call create_fixture(full_file, radial_grid_full, .false., .false.)
    call read_gvec_cas3d_file(full_file, equilibrium, info)
    call require(info == reader_ok, "full-mesh fixture was rejected")
    call require(equilibrium%radial_grid == radial_grid_full, &
        "full mesh was not classified")

    call create_fixture(symmetric_file, radial_grid_half, .true., .false.)
    call read_gvec_cas3d_file(symmetric_file, equilibrium, info)
    call require(info == reader_ok, "stellarator-symmetric fixture was rejected")
    call require(maxval(abs(equilibrium%yhat%cosine)) == 0.0_dp, &
        "missing odd cosine component was not zero-filled")
    call require(maxval(abs(equilibrium%xhat%sine)) == 0.0_dp, &
        "missing even sine component was not zero-filled")

    call create_fixture(corrupt_file, radial_grid_half, .false., .true.)
    call read_gvec_cas3d_file(corrupt_file, equilibrium, info)
    call require(info == reader_schema_error, &
        "missing required harmonic component was accepted")
    call create_fixture(nan_file, radial_grid_half, .false., .false.)
    call overwrite_pressure_with_nan(nan_file)
    call read_gvec_cas3d_file(nan_file, equilibrium, info)
    call require(info == reader_data_error, "nonfinite profile was accepted")
    call create_fixture(mode_file, radial_grid_half, .false., .false.)
    call overwrite_toroidal_modes(mode_file)
    call read_gvec_cas3d_file(mode_file, equilibrium, info)
    call require(info == reader_coordinate_error, &
        "invalid toroidal mode order was accepted")
    call read_gvec_cas3d_file(missing_file, equilibrium, info)
    call require(info == reader_open_error, "missing file was accepted")
    call create_fixture(schema_file, radial_grid_half, .false., .false., 1)
    call read_gvec_cas3d_file(schema_file, equilibrium, info)
    call require(info == reader_ok .and. equilibrium%schema_version == 1, &
        "schema version 1 was not identified")
    call create_fixture(future_file, radial_grid_half, .false., .false., 2)
    call read_gvec_cas3d_file(future_file, equilibrium, info)
    call require(info == reader_schema_error, &
        "unsupported schema version was accepted")
    call create_fixture(incomplete_file, radial_grid_half, .false., .false., -1)
    call read_gvec_cas3d_file(incomplete_file, equilibrium, info)
    call require(info == reader_schema_error, &
        "incomplete schema metadata was accepted")
    call create_fixture(malformed_file, radial_grid_half, .false., .false., -2)
    call read_gvec_cas3d_file(malformed_file, equilibrium, info)
    call require(info == reader_schema_error, &
        "malformed schema metadata was accepted as legacy")
    call create_fixture(frame_file, radial_grid_half, .false., .false., &
        position_frame="xhat,yhat rotated by winding*zeta_B")
    call read_gvec_cas3d_file(frame_file, equilibrium, info)
    call require(info == reader_ok &
        .and. equilibrium%has_boozer_position_frame, &
        "Boozer position frame was not identified")
    call create_fixture(bad_frame_file, radial_grid_half, .false., .false., &
        position_frame="computational zeta")
    call read_gvec_cas3d_file(bad_frame_file, equilibrium, info)
    call require(info == reader_schema_error, &
        "unknown position frame was accepted")

    call delete_fixture(half_file)
    call delete_fixture(full_file)
    call delete_fixture(symmetric_file)
    call delete_fixture(corrupt_file)
    call delete_fixture(nan_file)
    call delete_fixture(mode_file)
    call delete_fixture(schema_file)
    call delete_fixture(future_file)
    call delete_fixture(incomplete_file)
    call delete_fixture(malformed_file)
    call delete_fixture(frame_file)
    call delete_fixture(bad_frame_file)
    write (*, "(a)") "PASS"

contains

    subroutine verify_half_fixture(state)
        type(gvec_cas3d_equilibrium_t), intent(in) :: state
        real(dp), parameter :: expected_s(ns) = [1.0_dp / 6.0_dp, &
            1.0_dp / 2.0_dp, 5.0_dp / 6.0_dp]

        call require(state%field_periods == 5, "field periods are wrong")
        call require(state%winding == 1, "winding is wrong")
        call require(abs(state%beta_average - 0.025_dp) < 1.0e-14_dp, &
            "average beta is wrong")
        call require(state%radial_grid == radial_grid_half, &
            "half mesh was not classified")
        call require(all(state%poloidal_modes == [0, 1]), &
            "poloidal mode order is wrong")
        call require(all(state%toroidal_modes == [0, 1, -1]), &
            "toroidal mode order is wrong")
        call require(maxval(abs(state%s - expected_s)) < 1.0e-14_dp, &
            "s was not recovered from rho")
        call require(abs(state%mod_b%cosine(2, 2, 3) - 1223.0_dp) < &
            1.0e-14_dp, "harmonic storage order is wrong")
        call require(abs(state%mod_b%sine(2, 2, 3) + 1223.0_dp) < &
            1.0e-14_dp, "sine harmonic value is wrong")
    end subroutine verify_half_fixture

    subroutine create_fixture(filename, grid_kind, symmetric, corrupt, &
            schema_version, position_frame)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: grid_kind
        logical, intent(in) :: symmetric, corrupt
        integer, intent(in), optional :: schema_version
        character(len=*), intent(in), optional :: position_frame
        type(fixture_ids_t) :: ids

        call define_fixture(filename, grid_kind, symmetric, corrupt, ids, &
            schema_version, position_frame)
        call write_fixture(grid_kind, ids)
        call require_netcdf(nc_close_file(ids%ncid))
    end subroutine create_fixture

    subroutine define_fixture(filename, grid_kind, symmetric, corrupt, ids, &
            schema_version, position_frame)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: grid_kind
        logical, intent(in) :: symmetric, corrupt
        type(fixture_ids_t), intent(out) :: ids
        integer, intent(in), optional :: schema_version
        character(len=*), intent(in), optional :: position_frame
        integer :: dimensions(3), dim_s, dim_m, dim_n
        character(len=16) :: version_text

        call require_netcdf(nc_create_netcdf4(filename, ids%ncid))
        call require_netcdf(nc_def_dimension(ids%ncid, "s", ns, dim_s))
        call require_netcdf(nc_def_dimension(ids%ncid, "m", nm, dim_m))
        call require_netcdf(nc_def_dimension(ids%ncid, "n", nn, dim_n))
        call define_metadata_variables(ids)
        call define_coordinate_variables(ids, grid_kind, dim_s, dim_m, dim_n)
        call define_profile_variables(ids, dim_s)
        dimensions(1) = dim_s
        dimensions(2) = dim_m
        dimensions(3) = dim_n
        call define_harmonic_variables(ids, symmetric, corrupt, dimensions)
        call require_netcdf(nc_put_global_text(ids%ncid, &
            "stellarator_symmetry", merge("True ", "False", symmetric)))
        if (present(schema_version)) then
            if (schema_version == -2) then
                call require_netcdf(nc_put_global_text(ids%ncid, &
                    "gliss_schema", repeat("x", 40)))
                call require_netcdf(nc_put_global_text(ids%ncid, &
                    "gliss_schema_version", repeat("2", 40)))
            else
                if (schema_version >= 0) call require_netcdf( &
                    nc_put_global_text(ids%ncid, "gliss_schema", &
                    "gvec-cas3d-export"))
                write (version_text, "(i0)") max(1, schema_version)
                call require_netcdf(nc_put_global_text(ids%ncid, &
                    "gliss_schema_version", trim(version_text)))
            end if
        end if
        if (present(position_frame)) call require_netcdf(nc_put_global_text( &
            ids%ncid, "position_frame", position_frame))
        call require_netcdf(nc_end_definitions(ids%ncid))
    end subroutine define_fixture

    subroutine define_metadata_variables(ids)
        type(fixture_ids_t), intent(inout) :: ids

        call require_netcdf(nc_def_scalar(ids%ncid, "N_FP", nc_int64, ids%nfp))
        call require_netcdf(nc_def_scalar(ids%ncid, "beta_avg", nc_double, &
            ids%beta))
        call require_netcdf(nc_def_scalar(ids%ncid, "winding", nc_int64, &
            ids%winding))
    end subroutine define_metadata_variables

    subroutine define_coordinate_variables(ids, grid_kind, dim_s, dim_m, dim_n)
        type(fixture_ids_t), intent(inout) :: ids
        integer, intent(in) :: grid_kind, dim_s, dim_m, dim_n
        integer :: dimension(1)

        dimension(1) = dim_m
        call require_netcdf(nc_def_variable(ids%ncid, "m", nc_int64, &
            dimension, ids%m))
        dimension(1) = dim_n
        call require_netcdf(nc_def_variable(ids%ncid, "n", nc_int64, &
            dimension, ids%n))
        dimension(1) = dim_s
        call require_netcdf(nc_def_variable(ids%ncid, "rho", nc_double, &
            dimension, ids%rho))
        if (grid_kind == radial_grid_full) then
            call require_netcdf(nc_def_variable(ids%ncid, "s", nc_double, &
                dimension, ids%s))
        end if
    end subroutine define_coordinate_variables

    subroutine define_profile_variables(ids, dim_s)
        type(fixture_ids_t), intent(inout) :: ids
        integer, intent(in) :: dim_s
        integer :: dimension(1), profile

        dimension(1) = dim_s
        do profile = 1, profile_count
            call require_netcdf(nc_def_variable(ids%ncid, &
                trim(profile_names(profile)), nc_double, dimension, &
                ids%profiles(profile)))
        end do
    end subroutine define_profile_variables

    subroutine define_harmonic_variables(ids, symmetric, corrupt, dimensions)
        type(fixture_ids_t), intent(inout) :: ids
        logical, intent(in) :: symmetric, corrupt
        integer, intent(in) :: dimensions(3)
        integer :: field
        logical :: odd

        do field = 1, field_count
            odd = trim(field_names(field)) == "yhat"
            if (trim(field_names(field)) == "zhat") odd = .true.
            if (.not. symmetric .or. .not. odd) then
                call require_netcdf(nc_def_variable(ids%ncid, &
                    trim(field_names(field)) // "_mnc", nc_double, &
                    dimensions, ids%pairs(1, field)))
            end if
            if (.not. symmetric .or. odd) then
                if (corrupt .and. trim(field_names(field)) == "Jac") cycle
                call require_netcdf(nc_def_variable(ids%ncid, &
                    trim(field_names(field)) // "_mns", nc_double, &
                    dimensions, ids%pairs(2, field)))
            end if
        end do
    end subroutine define_harmonic_variables

    subroutine write_fixture(grid_kind, ids)
        integer, intent(in) :: grid_kind
        type(fixture_ids_t), intent(in) :: ids
        real(dp) :: s_values(ns), rho_values(ns)
        real(dp) :: profile_values(ns)
        integer :: profile

        call radial_coordinates(grid_kind, s_values, rho_values)
        call require_netcdf(nc_put_integer(ids%ncid, ids%nfp, 5))
        call require_netcdf(nc_put_real(ids%ncid, ids%beta, 0.025_dp))
        call require_netcdf(nc_put_integer(ids%ncid, ids%winding, 1))
        call require_netcdf(nc_put_integer(ids%ncid, ids%m, [0, 1]))
        call require_netcdf(nc_put_integer(ids%ncid, ids%n, [0, 1, -1]))
        call require_netcdf(nc_put_real(ids%ncid, ids%rho, rho_values))
        if (ids%s > 0) call require_netcdf(nc_put_real(ids%ncid, ids%s, s_values))
        do profile = 1, profile_count
            profile_values = real(profile, dp) + s_values
            call require_netcdf(nc_put_real(ids%ncid, ids%profiles(profile), &
                profile_values))
        end do
        call write_harmonics(ids)
    end subroutine write_fixture

    subroutine radial_coordinates(grid_kind, s_values, rho_values)
        integer, intent(in) :: grid_kind
        real(dp), intent(out) :: s_values(ns), rho_values(ns)

        if (grid_kind == radial_grid_half) then
            s_values = [1.0_dp / 6.0_dp, 1.0_dp / 2.0_dp, 5.0_dp / 6.0_dp]
        else
            s_values = [1.0e-8_dp, 0.5_dp, 1.0_dp]
        end if
        rho_values = sqrt(s_values)
    end subroutine radial_coordinates

    subroutine write_harmonics(ids)
        type(fixture_ids_t), intent(in) :: ids
        real(dp) :: negative_values(nn, nm, ns), values(nn, nm, ns)
        integer :: field, radial, poloidal, toroidal

        do field = 1, field_count
            do radial = 1, ns
                do poloidal = 1, nm
                    do toroidal = 1, nn
                        values(toroidal, poloidal, radial) = &
                            1000.0_dp * real(field, dp) + &
                            100.0_dp * real(radial, dp) + &
                            10.0_dp * real(poloidal, dp) + real(toroidal, dp)
                    end do
                end do
            end do
            if (ids%pairs(1, field) > 0) then
                call require_netcdf(nc_put_real(ids%ncid, &
                    ids%pairs(1, field), values))
            end if
            if (ids%pairs(2, field) > 0) then
                do radial = 1, ns
                    do poloidal = 1, nm
                        do toroidal = 1, nn
                            negative_values(toroidal, poloidal, radial) = &
                                -values(toroidal, poloidal, radial)
                        end do
                    end do
                end do
                call require_netcdf(nc_put_real(ids%ncid, &
                    ids%pairs(2, field), negative_values))
            end if
        end do
    end subroutine write_harmonics

    subroutine overwrite_pressure_with_nan(filename)
        character(len=*), intent(in) :: filename
        real(dp) :: values(ns)
        integer :: ncid, varid

        values(1) = 1.0_dp
        values(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        values(3) = 3.0_dp
        call require_netcdf(nc_open_write(filename, ncid))
        call require_netcdf(nc_inquire_variable_id(ncid, "p", varid))
        call require_netcdf(nc_put_real(ncid, varid, values))
        call require_netcdf(nc_close_file(ncid))
    end subroutine overwrite_pressure_with_nan

    subroutine overwrite_toroidal_modes(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid, varid

        call require_netcdf(nc_open_write(filename, ncid))
        call require_netcdf(nc_inquire_variable_id(ncid, "n", varid))
        call require_netcdf(nc_put_integer(ncid, varid, [0, -1, 1]))
        call require_netcdf(nc_close_file(ncid))
    end subroutine overwrite_toroidal_modes

    subroutine delete_fixture(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, status

        open (newunit=unit, file=filename, status="old", action="read", &
            iostat=status)
        if (status == 0) close (unit, status="delete")
    end subroutine delete_fixture

    subroutine require_netcdf(status)
        integer, intent(in) :: status

        call require(status == nc_noerr, "NetCDF fixture operation failed")
    end subroutine require_netcdf

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_gvec_cas3d_reader
