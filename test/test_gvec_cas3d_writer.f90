program test_gvec_cas3d_writer
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_full
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
    end subroutine require_equal

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
