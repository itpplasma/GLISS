program test_marginality_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_two_component_problem, only: &
        compatible_two_component_problem_t
    use marginality_spectrum, only: marginality_spectrum_invalid, &
        marginality_spectrum_ok, marginality_spectrum_result_t, &
        solve_compatible_marginality_problem
    implicit none

    type(compatible_two_component_problem_t) :: invalid, problem
    type(marginality_spectrum_result_t) :: result
    real(dp), allocatable :: vector(:)
    real(dp) :: mass_image(2)
    character(len=128) :: message
    integer :: info

    call build_problem(problem)
    call solve_compatible_marginality_problem(problem, .true., result, info, &
        message, vector)
    if (info /= marginality_spectrum_ok) &
        write (error_unit, "(a,i0,2a)") "solver status ", info, ": ", &
        trim(message)
    call require(info == marginality_spectrum_ok, &
        "valid compatible problem solve failed")
    call require(size(vector) == 2, "refined eigenvector size differs")
    call require(result%has_eigenpair, "eigenpair solve was not identified")
    call require(result%negative_count == 1, "negative inertia count differs")
    call require(result%lowest_eigenvalue < 0.0_dp, &
        "lowest eigenvalue has the wrong sign")
    mass_image = matmul(problem%mass, vector)
    call require(abs(dot_product(vector, mass_image) - 1.0_dp) < 1.0e-13_dp, &
        "lowest eigenvector is not mass normalized")
    call require(result%eigenpair_residual < 1.0e-13_dp, &
        "exact eigenpair residual is too large")

    call solve_compatible_marginality_problem(problem, .false., result, info, &
        message, vector)
    call require(info == marginality_spectrum_ok, "inertia-only solve failed")
    call require(size(vector) == 0, "inertia-only solve returned a vector")
    call require(.not. result%has_eigenpair, &
        "inertia-only solve was identified as an eigenpair")
    call require(result%negative_count == 1, &
        "inertia-only negative count differs")

    call solve_compatible_marginality_problem(invalid, .true., result, info, &
        message, vector)
    call require(info == marginality_spectrum_invalid, &
        "unallocated problem was accepted")
    call require(index(message, "invalid") > 0, &
        "invalid problem did not return a useful message")
    problem%mass(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call solve_compatible_marginality_problem(problem, .true., result, info, &
        message, vector)
    call require(info == marginality_spectrum_invalid, &
        "nonfinite mass matrix was accepted")
    write (*, "(a)") "PASS"

contains

    subroutine build_problem(value)
        type(compatible_two_component_problem_t), intent(out) :: value

        allocate (value%stiffness(2, 2), source=0.0_dp)
        allocate (value%mass(2, 2), source=0.0_dp)
        allocate (value%stiffness_terms(2, 2, 4), source=0.0_dp)
        value%stiffness(1, 1) = -2.0_dp
        value%stiffness(2, 2) = 3.0_dp
        value%stiffness(1, 2) = 0.25_dp
        value%stiffness(2, 1) = 0.25_dp
        value%mass(1, 1) = 1.0_dp
        value%mass(2, 2) = 2.0_dp
        value%stiffness_terms(:, :, 1) = value%stiffness
        value%normal_unknowns = 1
        value%eta_unknowns = 1
    end subroutine build_problem

    subroutine require(condition, message_text)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message_text

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message_text
        error stop 1
    end subroutine require

end program test_marginality_spectrum
