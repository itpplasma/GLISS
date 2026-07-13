program test_terpsichore_solution_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: dynamic_family_layout_t
    use fourier_phase_kind, only: phase_sine
    use terpsichore_reduced_layout, only: &
        build_terpsichore_reduced_free_boundary_layout
    use terpsichore_solution_fixture, only: &
        build_terpsichore_plasma_solution, &
        read_terpsichore_solution_fixture, terpsichore_solution_fixture_t, &
        terpsichore_solution_invalid, terpsichore_solution_ok
    implicit none

    type(terpsichore_solution_fixture_t) :: fixture
    type(dynamic_family_layout_t) :: layout
    real(dp), allocatable :: vector(:)
    integer, allocatable :: element_map(:, :)
    integer :: info, unit

    open (newunit=unit, status="scratch", access="sequential", &
        form="unformatted", action="readwrite")
    call write_fixture(unit)
    rewind (unit)
    call read_terpsichore_solution_fixture(unit, 1, fixture, info)
    close (unit)
    call require(info == terpsichore_solution_ok, &
        "valid TERPSICHORE solution fixture was rejected")
    call require(fixture%potential_energy == -6.0_dp &
        .and. fixture%kinetic_energy == 3.0_dp, &
        "TERPSICHORE reference energies are wrong")
    call build_terpsichore_reduced_free_boundary_layout([1], [0], &
        [phase_sine], 2, layout, element_map, info)
    call build_terpsichore_plasma_solution(fixture, layout, vector, info)
    call require(info == terpsichore_solution_ok, &
        "TERPSICHORE plasma solution mapping failed")
    call require(all(vector == [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp]), &
        "TERPSICHORE plasma solution map is wrong")
    open (newunit=unit, status="scratch", access="sequential", &
        form="unformatted", action="readwrite")
    call write_fixture(unit)
    rewind (unit)
    call read_terpsichore_solution_fixture(unit, 2, fixture, info)
    close (unit)
    call require(info == terpsichore_solution_invalid, &
        "inconsistent IVAC was accepted by the solution reader")
    write (*, "(a)") "PASS"

contains

    subroutine write_fixture(output_unit)
        integer, intent(in) :: output_unit
        integer, parameter :: mode_m(1) = 1, mode_n(1) = 0
        real(dp), parameter :: normal(1, 0:3) = &
            reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [1, 4])
        real(dp), parameter :: tangential(1, 3) = &
            reshape([4.0_dp, 5.0_dp, 6.0_dp], [1, 3])
        real(dp) :: cell(2), mode_surface(1, 0:2), point_surface(1, 0:2)
        real(dp) :: surface(3)

        cell = 0.0_dp
        surface = 0.0_dp
        mode_surface = 0.0_dp
        point_surface = 0.0_dp
        write (output_unit) 2, 1, 1, 1, 1, 1
        write (output_unit) 0.0_dp
        write (output_unit) mode_m, mode_n, [0.0_dp], normal, tangential, &
            reshape([0.0_dp, 0.0_dp, 0.0_dp], [1, 3])
        write (output_unit) 0.0_dp
        write (output_unit) 0.0_dp
        write (output_unit) 1, mode_m, mode_n, mode_surface, point_surface, &
            surface, surface, cell, cell, -6.0_dp, 3.0_dp
    end subroutine write_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_terpsichore_solution_fixture
