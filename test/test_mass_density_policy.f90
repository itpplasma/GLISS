program test_mass_density_policy
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use mass_density_policy, only: evaluate_mass_density, &
        mass_density_invalid, mass_density_ok, mass_density_outside, &
        mass_density_profile_t, validate_mass_density_profile
    implicit none

    type(mass_density_profile_t) :: profile
    real(dp) :: density
    integer :: info

    profile%s = [0.0_dp, 0.4_dp, 1.0_dp]
    profile%kilograms_per_cubic_metre = [2.0_dp, 3.0_dp, 6.0_dp]
    call validate_mass_density_profile(profile, info)
    call require(info == mass_density_ok, "valid density profile was rejected")
    call evaluate_mass_density(profile, 0.7_dp, density, info)
    call require(info == mass_density_ok, "density interpolation failed")
    call require(abs(density - 4.5_dp) < 1.0e-14_dp, &
        "density interpolation is wrong")
    call evaluate_mass_density(profile, 1.0_dp, density, info)
    call require(info == mass_density_ok .and. density == 6.0_dp, &
        "density endpoint evaluation is wrong")

    call evaluate_mass_density(profile, -0.1_dp, density, info)
    call require(info == mass_density_outside, &
        "out-of-domain density coordinate was accepted")
    profile%kilograms_per_cubic_metre(2) = 0.0_dp
    call validate_mass_density_profile(profile, info)
    call require(info == mass_density_invalid, &
        "nonpositive mass density was accepted")
    profile%kilograms_per_cubic_metre(2) = ieee_value(density, ieee_quiet_nan)
    call validate_mass_density_profile(profile, info)
    call require(info == mass_density_invalid, &
        "nonfinite mass density was accepted")
    profile%kilograms_per_cubic_metre = [2.0_dp, 3.0_dp, 6.0_dp]
    profile%s = [0.0_dp, 0.7_dp, 0.6_dp]
    call validate_mass_density_profile(profile, info)
    call require(info == mass_density_invalid, &
        "nonmonotone density grid was accepted")

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

end program test_mass_density_policy
