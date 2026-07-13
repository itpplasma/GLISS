program test_fixed_boundary_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_invalid, fixed_boundary_ok, &
        fixed_boundary_problem_t, fixed_boundary_spectrum_result_t, &
        solve_fixed_boundary_class
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none

    character(len=*), parameter :: fixture = "fixed_boundary_spectrum.nc"
    real(dp), parameter :: legacy_lowest(2) = &
        [-1.2278816610508129e2_dp, -1.2278816623062953e2_dp]
    real(dp), parameter :: legacy_relative_tolerance = 1.0e-9_dp
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(fixed_boundary_problem_t) :: problem
    type(fixed_boundary_spectrum_result_t) :: first, second, repeated
    integer :: info

    call create_cylinder_fixture(fixture)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "spectrum fixture read failed")

    call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
        2.0_dp, 1.0_dp, [1, 2], [1, 1], 1, problem, info)
    call require(info == fixed_boundary_ok, "problem construction failed")
    call solve_fixed_boundary_class(problem, 1, first, info)
    call require(info == fixed_boundary_ok, "class-one solve failed")
    call solve_fixed_boundary_class(problem, 2, second, info)
    call require(info == fixed_boundary_ok, "class-two solve failed")
    call check_result(first, 1)
    call check_result(second, 2)
    call require(abs(first%lowest_eigenvalue - legacy_lowest(1)) &
        <= legacy_relative_tolerance * abs(legacy_lowest(1)), &
        "class-one result changed from the legacy CLI")
    call require(abs(second%lowest_eigenvalue - legacy_lowest(2)) &
        <= legacy_relative_tolerance * abs(legacy_lowest(2)), &
        "class-two result changed from the legacy CLI")

    call solve_fixed_boundary_class(problem, 1, repeated, info)
    call require(info == fixed_boundary_ok, "repeated solve failed")
    call require(repeated%lowest_eigenvalue == first%lowest_eigenvalue, &
        "repeated solve changed the eigenvalue")
    call require(all(repeated%eigenvector == first%eigenvector), &
        "repeated solve changed the eigenvector")

    call check_invalid_inputs(equilibrium)
    call delete_fixture()
    write (*, "(a)") "PASS"

contains

    subroutine check_result(result, parity_class)
        type(fixed_boundary_spectrum_result_t), intent(in) :: result
        integer, intent(in) :: parity_class

        call require(result%parity_class == parity_class, &
            "result parity class is wrong")
        call require(result%field_periods == 1, &
            "result field-period count is wrong")
        call require(result%mode_count == 2, "result mode count is wrong")
        call require(result%unknowns == 196, "result unknown count is wrong")
        call require(result%normal_unknowns == 64, &
            "normal unknown count is wrong")
        call require(result%eta_unknowns == 66, &
            "eta unknown count is wrong")
        call require(result%mu_unknowns == 66, &
            "mu unknown count is wrong")
        call require(size(result%eigenvector) == result%unknowns, &
            "eigenvector size is wrong")
        call require(all(ieee_is_finite(result%eigenvector)), &
            "eigenvector contains nonfinite values")
        call require(ieee_is_finite(result%lowest_eigenvalue), &
            "eigenvalue is nonfinite")
        call require(result%certificate == result%inertia_interval &
            + result%eigenpair_residual + result%eigenpair_resolution, &
            "certificate components do not close")
        call require(result%negative_count >= 0, "negative count is invalid")
        call require(result%floor_count >= 0, "floor count is invalid")
    end subroutine check_result

    subroutine check_invalid_inputs(local_equilibrium)
        type(gvec_cas3d_equilibrium_t), intent(in) :: local_equilibrium
        type(fixed_boundary_problem_t) :: invalid_problem
        type(fixed_boundary_spectrum_result_t) :: invalid_result
        real(dp) :: nan
        integer :: status

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call build_fixed_boundary_problem(local_equilibrium, nan, 2.0_dp, &
            1.0_dp, [1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "nonfinite adiabatic index was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            0.0_dp, 1.0_dp, [1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "zero density was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 0.0_dp, [1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "zero floor was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [1], 3, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "unknown radial quadrature was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1, 1], [1, 1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "duplicate mode was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [-1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "negative poloidal mode was accepted")
        call solve_fixed_boundary_class(problem, 0, invalid_result, status)
        call require(status == fixed_boundary_invalid, &
            "invalid parity class was accepted")
    end subroutine check_invalid_inputs

    subroutine delete_fixture()
        integer :: unit, status

        open (newunit=unit, file=fixture, status="old", iostat=status)
        call require(status == 0, "failed to open spectrum fixture")
        close (unit, status="delete", iostat=status)
        call require(status == 0, "failed to delete spectrum fixture")
    end subroutine delete_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_fixed_boundary_spectrum
