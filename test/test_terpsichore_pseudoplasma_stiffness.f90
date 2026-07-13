program test_terpsichore_pseudoplasma_stiffness
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use terpsichore_pseudoplasma_fixture, only: &
        pseudoplasma_fixture_invalid, pseudoplasma_fixture_ok, &
        read_terpsichore_pseudoplasma_fixture, &
        terpsichore_pseudoplasma_fixture_is_valid, &
        terpsichore_pseudoplasma_fixture_t
    use terpsichore_pseudoplasma_stiffness, only: &
        assemble_terpsichore_pseudoplasma_stiffness, &
        pseudoplasma_stiffness_ok
    use vacuum_schur, only: eliminate_vacuum, vacuum_schur_ok
    implicit none

    call test_analytic_matrix()
    call test_fixture_reader()
    call test_invalid_fixture()
    write (*, "(a)") "PASS"

contains

    subroutine test_analytic_matrix()
        type(terpsichore_pseudoplasma_fixture_t) :: fixture
        real(dp), allocatable :: stiffness(:, :), expected(:, :)
        real(dp), allocatable :: effective(:, :), response(:, :)
        integer :: info

        call build_fixture(fixture)
        call assemble_terpsichore_pseudoplasma_stiffness(fixture, &
            stiffness, info)
        call require(info == pseudoplasma_stiffness_ok, &
            "analytic pseudoplasma fixture was rejected")
        allocate (expected(4, 4), source=0.0_dp)
        expected(1, 1) = 16.5_dp
        expected(1, 2) = -15.5_dp
        expected(2, 1) = -15.5_dp
        expected(2, 2) = 33.0_dp
        expected(3, 3) = 2.0_dp
        expected(4, 4) = 2.0_dp
        call require(maxval(abs(stiffness - expected)) < 1.0e-14_dp, &
            "analytic pseudoplasma matrix is wrong")
        call require(maxval(abs(stiffness - transpose(stiffness))) < &
            1.0e-14_dp, "pseudoplasma matrix is not symmetric")
        call eliminate_vacuum(stiffness(:1, :1), stiffness(2:, 2:), &
            stiffness(:1, 2:), effective, response, info)
        call require(info == vacuum_schur_ok, &
            "pseudoplasma Schur reduction failed")
        call require(abs(effective(1, 1) &
            - (16.5_dp - 15.5_dp**2 / 33.0_dp)) < 1.0e-14_dp, &
            "pseudoplasma interface energy is wrong")
    end subroutine test_analytic_matrix

    subroutine test_fixture_reader()
        type(terpsichore_pseudoplasma_fixture_t) :: source, restored
        integer, parameter :: magic = int(z'47565031'), schema = 1
        integer :: info, unit

        call build_fixture(source)
        open (newunit=unit, status="scratch", access="sequential", &
            form="unformatted", action="readwrite")
        write (unit) magic, schema, source%plasma_intervals, &
            source%vacuum_intervals, source%modes
        write (unit) source%s, source%flux_t_slope, source%flux_p_slope
        write (unit) source%mode_m, source%mode_n
        write (unit) source%coefficient(1, :, :, :), &
            source%coefficient(2, :, :, :), &
            source%coefficient(3, :, :, :), &
            source%coefficient(4, :, :, :), &
            source%coefficient(5, :, :, :), &
            source%coefficient(6, :, :, :)
        rewind (unit)
        call read_terpsichore_pseudoplasma_fixture(unit, restored, info)
        close (unit)
        call require(info == pseudoplasma_fixture_ok, &
            "valid pseudoplasma diagnostic was rejected")
        call require(all(restored%s == source%s), &
            "pseudoplasma radial grid changed during ingestion")
        call require(all(restored%coefficient == source%coefficient), &
            "pseudoplasma coefficients changed during ingestion")
    end subroutine test_fixture_reader

    subroutine test_invalid_fixture()
        type(terpsichore_pseudoplasma_fixture_t) :: fixture, restored
        integer, parameter :: magic = int(z'47565031')
        integer :: info, unit

        call build_fixture(fixture)
        fixture%s(2) = fixture%s(1)
        call require(.not. terpsichore_pseudoplasma_fixture_is_valid(fixture), &
            "nonmonotone vacuum grid was accepted")
        open (newunit=unit, status="scratch", access="sequential", &
            form="unformatted", action="readwrite")
        write (unit) magic, 2, 4, 2, 1
        rewind (unit)
        call read_terpsichore_pseudoplasma_fixture(unit, restored, info)
        close (unit)
        call require(info == pseudoplasma_fixture_invalid, &
            "unknown pseudoplasma schema was accepted")
        call build_fixture(fixture)
        fixture%mode_m(1) = 1000001
        call require(.not. terpsichore_pseudoplasma_fixture_is_valid(fixture), &
            "unsafe mode number was accepted")
        open (newunit=unit, status="scratch", access="sequential", &
            form="unformatted", action="readwrite")
        write (unit) magic, 1, 4, 2, 1
        rewind (unit)
        call read_terpsichore_pseudoplasma_fixture(unit, restored, info)
        close (unit)
        call require(info == pseudoplasma_fixture_invalid, &
            "truncated pseudoplasma fixture was accepted")
    end subroutine test_invalid_fixture

    subroutine build_fixture(fixture)
        type(terpsichore_pseudoplasma_fixture_t), intent(out) :: fixture

        fixture%plasma_intervals = 4
        fixture%vacuum_intervals = 2
        fixture%modes = 1
        fixture%flux_p_slope = 1.0_dp
        fixture%flux_t_slope = 0.0_dp
        allocate (fixture%s(0:2), source=[1.0_dp, 1.5_dp, 2.0_dp])
        allocate (fixture%mode_m(1), source=1)
        allocate (fixture%mode_n(1), source=0)
        allocate (fixture%coefficient(6, 1, 1, 2), source=0.0_dp)
        fixture%coefficient(1, 1, 1, :) = -1.0_dp
        fixture%coefficient(3, 1, 1, :) = -1.0_dp
        fixture%coefficient(6, 1, 1, :) = -2.0_dp
    end subroutine build_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_terpsichore_pseudoplasma_stiffness
