program test_compatible_vector_evaluation
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_two_component_problem, only: &
        compatible_problem_invalid, compatible_problem_ok, &
        evaluate_compatible_two_component_vector
    implicit none

    integer, parameter :: mode_m(2) = [1, 2], mode_n(2) = [1, 1]
    real(dp), parameter :: stored_power(2) = [0.5_dp, 0.0_dp]
    real(dp), parameter :: eta_power(2) = [0.5_dp, 0.0_dp]
    real(dp), parameter :: coordinates(2) = [0.25_dp, 0.75_dp]
    real(dp), parameter :: vector(6) = [2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, &
        11.0_dp, 13.0_dp]
    real(dp), allocatable :: eta(:, :), normal(:, :)
    real(dp) :: bad_coordinates(2)
    integer :: info

    call evaluate_compatible_two_component_vector(2, mode_m, mode_n, &
        stored_power, eta_power, 1, 1, coordinates, vector, normal, eta, info)
    call require(info == compatible_problem_ok, "valid vector was rejected")
    call require(all(shape(normal) == [2, 2]), "normal shape differs")
    call require(all(shape(eta) == [2, 2]), "eta shape differs")
    call require(abs(normal(1, 1) - 2.0_dp) < 1.0e-14_dp, &
        "weighted normal first point differs")
    call require(abs(normal(1, 2) - 1.5_dp) < 1.0e-14_dp, &
        "unweighted normal first point differs")
    call require(abs(normal(2, 1) - 1.0_dp / sqrt(0.75_dp)) &
        < 1.0e-14_dp, "weighted normal second point differs")
    call require(abs(normal(2, 2) - 1.5_dp) < 1.0e-14_dp, &
        "unweighted normal second point differs")
    call require(abs(eta(1, 1) - 10.0_dp) < 1.0e-14_dp, &
        "weighted eta first point differs")
    call require(abs(eta(1, 2) - 7.0_dp) < 1.0e-14_dp, &
        "unweighted eta first point differs")
    call require(abs(eta(2, 1) - 11.0_dp / sqrt(0.75_dp)) &
        < 1.0e-14_dp, "weighted eta second point differs")
    call require(abs(eta(2, 2) - 13.0_dp) < 1.0e-14_dp, &
        "unweighted eta second point differs")

    bad_coordinates = coordinates
    bad_coordinates(1) = 0.0_dp
    call evaluate_compatible_two_component_vector(2, mode_m, mode_n, &
        stored_power, eta_power, 1, 1, bad_coordinates, vector, normal, eta, &
        info)
    call require(info == compatible_problem_invalid, &
        "axis coordinate was accepted")
    bad_coordinates = coordinates
    bad_coordinates(2) = ieee_value(0.0_dp, ieee_quiet_nan)
    call evaluate_compatible_two_component_vector(2, mode_m, mode_n, &
        stored_power, eta_power, 1, 1, bad_coordinates, vector, normal, eta, &
        info)
    call require(info == compatible_problem_invalid, &
        "nonfinite coordinate was accepted")
    call evaluate_compatible_two_component_vector(2, mode_m, mode_n, &
        stored_power, eta_power, 1, 1, coordinates, vector(:5), normal, eta, &
        info)
    call require(info == compatible_problem_invalid, &
        "wrong vector length was accepted")

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_compatible_vector_evaluation
