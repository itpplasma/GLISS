program test_helical_cylinder_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use helical_cylinder_limit, only: elliptical_wall_radius_ratio, &
        helical_limit_invalid_elongation, helical_limit_ok, &
        helical_vertical_margin, helical_vertical_threshold
    implicit none

    real(dp), parameter :: tolerance = 1.0e-14_dp
    real(dp) :: critical_fraction, wall_ratio
    integer :: info

    critical_fraction = helical_vertical_threshold(2.0_dp)
    call require(abs(critical_fraction - 0.4_dp) < tolerance, &
        "elongation-two vertical threshold is wrong")
    call require(abs(helical_vertical_margin(2.0_dp, critical_fraction)) &
        < tolerance, "vertical energy is not marginal at the threshold")
    call require(helical_vertical_margin(2.0_dp, 0.5_dp) > 0.0_dp, &
        "stable side of the vertical threshold is wrong")
    call require(helical_vertical_margin(2.0_dp, 0.3_dp) < 0.0_dp, &
        "unstable side of the vertical threshold is wrong")
    call require(helical_vertical_threshold(1.0_dp) == 0.0_dp, &
        "circular helical-cylinder threshold is wrong")
    call elliptical_wall_radius_ratio(2.0_dp, wall_ratio, info)
    call require(info == helical_limit_ok, &
        "valid elliptical conducting-wall input was rejected")
    call require(abs(wall_ratio - sqrt(3.0_dp)) &
        < tolerance, "elliptical conducting-wall threshold is wrong")
    call require(abs(wall_ratio &
        - sqrt((2.0_dp + 1.0_dp) / (2.0_dp + 1.0_dp))) > 0.5_dp, &
        "corrupted elliptical-wall denominator was not detected")
    call elliptical_wall_radius_ratio(1.0_dp, wall_ratio, info)
    call require(info == helical_limit_invalid_elongation, &
        "singular elliptical conducting-wall input was accepted")
    call elliptical_wall_radius_ratio(ieee_value(0.0_dp, ieee_quiet_nan), &
        wall_ratio, info)
    call require(info == helical_limit_invalid_elongation, &
        "nonfinite elliptical conducting-wall input was accepted")
    call require(abs(critical_fraction - (2.0_dp**2 - 2.0_dp) &
        / (2.0_dp**2 - 1.0_dp)) > 0.1_dp, &
        "corrupted helical-cylinder denominator was not detected")

    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_helical_cylinder_limit
