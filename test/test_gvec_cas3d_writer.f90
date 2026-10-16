program test_gvec_cas3d_writer
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_boundary_t, &
        gvec_cas3d_equilibrium_t, harmonic_pair_t, radial_grid_full
    use gvec_cas3d_writer, only: write_gvec_cas3d_file, writer_invalid, &
        writer_ok, writer_open_error
    implicit none

    character(len=64) :: source_file, first_file, second_file
    character(len=64) :: plain_source_file, plain_file, full_file, invalid_file
    type(gvec_cas3d_equilibrium_t) :: source, first, second, invalid
    type(gvec_cas3d_equilibrium_t) :: plain, full, roundtrip
    integer :: info, index

    interface
        function get_process_id() result(process_id) bind(c, name="getpid")
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int) :: process_id
        end function get_process_id
    end interface

    write (source_file, '("writer_source_",i0,".nc")') get_process_id()
    write (first_file, '("writer_first_",i0,".nc")') get_process_id()
    write (second_file, '("writer_second_",i0,".nc")') get_process_id()
    write (plain_source_file, '("writer_plain_source_",i0,".nc")') &
        get_process_id()
    write (plain_file, '("writer_plain_",i0,".nc")') get_process_id()
    write (full_file, '("writer_full_",i0,".nc")') get_process_id()
    write (invalid_file, '("writer_invalid_",i0,".nc")') get_process_id()

    call create_cylinder_fixture(source_file, chart_shift=0.2_dp, surfaces=9)
    call read_gvec_cas3d_file(source_file, source, info)
    call require(info == reader_ok, "source read failed")
    call require(source%schema_version == 0, &
        "unversioned source was not identified as legacy")
    call require(all(source%rho**2 == source%s), &
        "reader did not canonicalize radial coordinates")

    call write_gvec_cas3d_file(first_file, source, info)
    call require(info == writer_ok, "first write failed")
    call read_gvec_cas3d_file(first_file, first, info)
    call require(info == reader_ok, "first round-trip read failed")
    call require(first%schema_version == 1, &
        "writer did not emit schema version 1")
    call require_equal(source, first)

    call write_gvec_cas3d_file(second_file, first, info)
    call require(info == writer_ok, "second write failed")
    call read_gvec_cas3d_file(second_file, second, info)
    call require(info == reader_ok, "second round-trip read failed")
    call require_equal(first, second)

    call create_cylinder_fixture(plain_source_file, surfaces=9)
    call read_gvec_cas3d_file(plain_source_file, plain, info)
    call require(info == reader_ok .and. .not. plain%has_chart_metric, &
        "plain source did not omit the chart metric")
    call write_gvec_cas3d_file(plain_file, plain, info)
    call require(info == writer_ok, "plain write failed")
    call read_gvec_cas3d_file(plain_file, roundtrip, info)
    call require(info == reader_ok .and. .not. roundtrip%has_chart_metric, &
        "plain round trip gained a chart metric")
    call require_equal(plain, roundtrip)

    full = source
    full%radial_grid = radial_grid_full
    full%s(1) = 1.0e-8_dp
    do index = 2, size(full%s)
        full%s(index) = real(index - 1, dp) / real(size(full%s) - 1, dp)
    end do
    full%rho = sqrt(full%s)
    full%s = full%rho**2
    call write_gvec_cas3d_file(full_file, full, info)
    call require(info == writer_ok, "full-grid write failed")
    call read_gvec_cas3d_file(full_file, roundtrip, info)
    call require(info == reader_ok, "full-grid round-trip read failed")
    call require_equal(full, roundtrip)

    call test_boundary_roundtrip(source, plain)
    call test_invalid_boundary(source, plain, invalid_file)

    invalid = source
    invalid%pressure(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call write_gvec_cas3d_file(invalid_file, invalid, info)
    call require(info == writer_invalid, "nonfinite state was accepted")
    call require(.not. file_exists(invalid_file), &
        "invalid write left an output file")

    invalid = source
    invalid%s(1) = nearest(invalid%s(1), 1.0_dp)
    call write_gvec_cas3d_file(invalid_file, invalid, info)
    call require(info == writer_invalid, &
        "noncanonical radial coordinates were accepted")
    call require(.not. file_exists(invalid_file), &
        "noncanonical write left an output file")

    call write_gvec_cas3d_file(first_file, source, info)
    call require(info == writer_open_error, "existing output was overwritten")
    call read_gvec_cas3d_file(first_file, first, info)
    call require(info == reader_ok, "rejected overwrite damaged output")

    call delete_file(source_file)
    call delete_file(first_file)
    call delete_file(second_file)
    call delete_file(plain_source_file)
    call delete_file(plain_file)
    call delete_file(full_file)
    write (*, "(a)") "PASS"

contains

    subroutine require_equal(expected, actual)
        type(gvec_cas3d_equilibrium_t), intent(in) :: expected, actual

        call require(actual%field_periods == expected%field_periods, &
            "field periods changed")
        call require(actual%winding == expected%winding, "winding changed")
        call require(actual%radial_grid == expected%radial_grid, &
            "radial grid changed")
        call require(actual%stellarator_symmetric .eqv. &
            expected%stellarator_symmetric, "symmetry changed")
        call require(actual%has_chart_metric .eqv. expected%has_chart_metric, &
            "chart-metric presence changed")
        call require(actual%has_boozer_position_frame .eqv. &
            expected%has_boozer_position_frame, "position frame changed")
        call require(actual%beta_average == expected%beta_average, &
            "average beta changed")
        call require(all(actual%poloidal_modes == expected%poloidal_modes), &
            "poloidal modes changed")
        call require(all(actual%toroidal_modes == expected%toroidal_modes), &
            "toroidal modes changed")
        call require(all(actual%rho == expected%rho), "rho changed")
        call require(all(actual%s == expected%s), "s changed")
        call require(all(actual%pressure == expected%pressure), &
            "pressure changed")
        call require(all(actual%b_theta_average == expected%b_theta_average), &
            "B_theta_average changed")
        call require(all(actual%b_zeta_average == expected%b_zeta_average), &
            "B_zeta_average changed")
        call require(all(actual%toroidal_flux == expected%toroidal_flux), &
            "toroidal flux changed")
        call require(all(actual%poloidal_flux == expected%poloidal_flux), &
            "poloidal flux changed")
        call require(all(actual%rotational_transform == &
            expected%rotational_transform), "rotational transform changed")
        call require_pair(expected%mod_b, actual%mod_b, "mod_B")
        call require_pair(expected%xhat, actual%xhat, "xhat")
        call require_pair(expected%yhat, actual%yhat, "yhat")
        call require_pair(expected%zhat, actual%zhat, "zhat")
        call require_pair(expected%jacobian, actual%jacobian, "Jac")
        call require_pair(expected%g_tt, actual%g_tt, "g_tt")
        call require_pair(expected%g_tz, actual%g_tz, "g_tz")
        call require_pair(expected%g_zz, actual%g_zz, "g_zz")
        call require_pair(expected%g_st, actual%g_st, "g_st")
        call require_pair(expected%g_sz, actual%g_sz, "g_sz")
        call require_pair(expected%second_form_tt, actual%second_form_tt, &
            "II_tt")
        call require_pair(expected%second_form_tz, actual%second_form_tz, &
            "II_tz")
        call require_pair(expected%second_form_zz, actual%second_form_zz, &
            "II_zz")
        call require_pair(expected%b_contravariant_theta, &
            actual%b_contravariant_theta, "B_contra_t")
        call require_pair(expected%b_contravariant_zeta, &
            actual%b_contravariant_zeta, "B_contra_z")
        call require_boundary_equal(expected%boundary, actual%boundary)
    end subroutine require_equal

    subroutine test_boundary_roundtrip(asymmetric_source, symmetric_source)
        type(gvec_cas3d_equilibrium_t), intent(in) :: asymmetric_source
        type(gvec_cas3d_equilibrium_t), intent(in) :: symmetric_source
        type(gvec_cas3d_equilibrium_t) :: expected, actual
        character(len=64) :: filename
        integer :: status, symmetry

        do symmetry = 0, 1
            write (filename, '("writer_boundary_",i0,"_",i0,".nc")') &
                symmetry, get_process_id()
            if (symmetry == 0) then
                expected = asymmetric_source
            else
                expected = symmetric_source
                expected%stellarator_symmetric = .true.
            end if
            call set_boundary(expected)
            if (symmetry == 1) expected%boundary%poloidal_flux_slope = 0.0_dp
            call write_gvec_cas3d_file(filename, expected, status)
            call require(status == writer_ok, "boundary write failed")
            call read_gvec_cas3d_file(filename, actual, status)
            call require(status == reader_ok, "boundary round-trip read failed")
            call require_equal(expected, actual)
            call delete_file(filename)

            expected = actual
            call write_gvec_cas3d_file(filename, expected, status)
            call require(status == writer_ok, "boundary rewrite failed")
            call read_gvec_cas3d_file(filename, actual, status)
            call require(status == reader_ok, "boundary rewrite read failed")
            call require_equal(expected, actual)
            call delete_file(filename)
        end do
    end subroutine test_boundary_roundtrip

    subroutine set_boundary(equilibrium)
        type(gvec_cas3d_equilibrium_t), intent(inout) :: equilibrium
        integer :: nm, nn, poloidal, toroidal
        real(dp) :: coefficient

        nm = size(equilibrium%poloidal_modes)
        nn = size(equilibrium%toroidal_modes)
        equilibrium%boundary%schema_version = 1
        equilibrium%boundary%toroidal_flux_slope = -1.25_dp
        equilibrium%boundary%poloidal_flux_slope = 0.375_dp
        allocate (equilibrium%boundary%xhat%cosine(1, nm, nn), &
            equilibrium%boundary%xhat%sine(1, nm, nn), &
            equilibrium%boundary%yhat%cosine(1, nm, nn), &
            equilibrium%boundary%yhat%sine(1, nm, nn), &
            equilibrium%boundary%zhat%cosine(1, nm, nn), &
            equilibrium%boundary%zhat%sine(1, nm, nn))
        do toroidal = 1, nn
            do poloidal = 1, nm
                coefficient = real(poloidal + 8 * toroidal, dp) / 16.0_dp
                equilibrium%boundary%xhat%cosine(1, poloidal, toroidal) = &
                    coefficient
                equilibrium%boundary%xhat%sine(1, poloidal, toroidal) = &
                    -coefficient / 2.0_dp
                equilibrium%boundary%yhat%cosine(1, poloidal, toroidal) = &
                    coefficient / 4.0_dp
                equilibrium%boundary%yhat%sine(1, poloidal, toroidal) = &
                    -coefficient / 8.0_dp
                equilibrium%boundary%zhat%cosine(1, poloidal, toroidal) = &
                    coefficient / 16.0_dp
                equilibrium%boundary%zhat%sine(1, poloidal, toroidal) = &
                    -coefficient / 32.0_dp
            end do
        end do
        if (equilibrium%stellarator_symmetric) then
            equilibrium%boundary%xhat%sine = 0.0_dp
            equilibrium%boundary%yhat%cosine = 0.0_dp
            equilibrium%boundary%zhat%cosine = 0.0_dp
        end if
    end subroutine set_boundary

    subroutine test_invalid_boundary(asymmetric_source, symmetric_source, filename)
        type(gvec_cas3d_equilibrium_t), intent(in) :: asymmetric_source
        type(gvec_cas3d_equilibrium_t), intent(in) :: symmetric_source
        character(len=*), intent(in) :: filename
        type(gvec_cas3d_equilibrium_t) :: valid, bad
        integer :: nm, nn, coordinate

        valid = asymmetric_source
        call set_boundary(valid)
        nm = size(valid%poloidal_modes)
        nn = size(valid%toroidal_modes)

        bad = valid
        bad%boundary%schema_version = 2
        call require_invalid_boundary(bad, filename, "unsupported boundary schema")
        bad = valid
        bad%boundary%toroidal_flux_slope = 0.0_dp
        call require_invalid_boundary(bad, filename, "zero boundary toroidal slope")
        bad = valid
        bad%boundary%toroidal_flux_slope = ieee_value(0.0_dp, ieee_quiet_nan)
        call require_invalid_boundary(bad, filename, &
            "nonfinite boundary toroidal slope")
        bad = valid
        bad%boundary%poloidal_flux_slope = ieee_value(0.0_dp, ieee_quiet_nan)
        call require_invalid_boundary(bad, filename, &
            "nonfinite boundary poloidal slope")
        bad = valid
        deallocate (bad%boundary%xhat%cosine)
        allocate (bad%boundary%xhat%cosine(2, nm, nn), source=0.0_dp)
        call require_invalid_boundary(bad, filename, "boundary radial dimension")
        bad = valid
        deallocate (bad%boundary%yhat%sine)
        call require_invalid_boundary(bad, filename, &
            "missing boundary sine coefficients")
        bad = valid
        deallocate (bad%boundary%zhat%cosine)
        allocate (bad%boundary%zhat%cosine(1, nm + 1, nn), source=0.0_dp)
        call require_invalid_boundary(bad, filename, "boundary mode dimension")
        bad = valid
        bad%boundary%yhat%sine(1, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call require_invalid_boundary(bad, filename, "nonfinite boundary coefficient")

        valid = symmetric_source
        valid%stellarator_symmetric = .true.
        call set_boundary(valid)
        do coordinate = 1, 3
            bad = valid
            select case (coordinate)
            case (1); bad%boundary%xhat%sine(1, 1, 1) = 0.25_dp
            case (2); bad%boundary%yhat%cosine(1, 1, 1) = 0.25_dp
            case (3); bad%boundary%zhat%cosine(1, 1, 1) = 0.25_dp
            end select
            call require_invalid_boundary(bad, filename, "boundary symmetry violation")
        end do

        bad = asymmetric_source
        allocate (bad%boundary%xhat%cosine(1, nm, nn), source=0.0_dp)
        call require_invalid_boundary(bad, filename, &
            "unversioned boundary coefficients")
        bad = asymmetric_source
        bad%boundary%toroidal_flux_slope = 1.0_dp
        call require_invalid_boundary(bad, filename, &
            "unversioned boundary toroidal slope")
        bad = asymmetric_source
        bad%boundary%poloidal_flux_slope = 1.0_dp
        call require_invalid_boundary(bad, filename, &
            "unversioned boundary poloidal slope")
    end subroutine test_invalid_boundary

    subroutine require_invalid_boundary(equilibrium, filename, name)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        character(len=*), intent(in) :: filename, name
        integer :: status

        call write_gvec_cas3d_file(filename, equilibrium, status)
        call require(status == writer_invalid, name // " was accepted")
        call require(.not. file_exists(filename), name // " created output")
    end subroutine require_invalid_boundary

    subroutine require_boundary_equal(expected, actual)
        type(gvec_cas3d_boundary_t), intent(in) :: expected, actual

        call require(actual%schema_version == expected%schema_version, &
            "boundary schema changed")
        call require(actual%toroidal_flux_slope == expected%toroidal_flux_slope, &
            "boundary toroidal flux slope changed")
        call require(actual%poloidal_flux_slope == expected%poloidal_flux_slope, &
            "boundary poloidal flux slope changed")
        if (expected%schema_version == 0) then
            call require(.not. allocated(actual%xhat%cosine), &
                "round trip gained boundary x cosine coefficients")
            call require(.not. allocated(actual%xhat%sine), &
                "round trip gained boundary x sine coefficients")
            call require(.not. allocated(actual%yhat%cosine), &
                "round trip gained boundary y cosine coefficients")
            call require(.not. allocated(actual%yhat%sine), &
                "round trip gained boundary y sine coefficients")
            call require(.not. allocated(actual%zhat%cosine), &
                "round trip gained boundary z cosine coefficients")
            call require(.not. allocated(actual%zhat%sine), &
                "round trip gained boundary z sine coefficients")
            return
        end if
        call require_pair(expected%xhat, actual%xhat, "boundary xhat")
        call require_pair(expected%yhat, actual%yhat, "boundary yhat")
        call require_pair(expected%zhat, actual%zhat, "boundary zhat")
    end subroutine require_boundary_equal

    subroutine require_pair(expected, actual, name)
        type(harmonic_pair_t), intent(in) :: expected, actual
        character(len=*), intent(in) :: name

        call require(all(shape(actual%cosine) == shape(expected%cosine)), &
            name // " cosine shape changed")
        call require(all(shape(actual%sine) == shape(expected%sine)), &
            name // " sine shape changed")
        call require(all(actual%cosine == expected%cosine), &
            name // " cosine changed")
        call require(all(actual%sine == expected%sine), &
            name // " sine changed")
    end subroutine require_pair

    logical function file_exists(filename)
        character(len=*), intent(in) :: filename

        inquire (file=filename, exist=file_exists)
    end function file_exists

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, status

        open (newunit=unit, file=filename, status="old", iostat=status)
        call require(status == 0, "failed to open output for deletion")
        close (unit, status="delete", iostat=status)
        call require(status == 0, "failed to delete output")
    end subroutine delete_file

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_gvec_cas3d_writer
