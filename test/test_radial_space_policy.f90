program test_radial_space_policy
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use radial_space_policy, only: evaluate_normal_basis, form_s_power_edge, &
        radial_space_config_t, radial_space_invalid, radial_space_ok
    implicit none

    type(radial_space_config_t) :: config
    real(dp) :: values(2), derivatives(2), plus(2), minus(2), step
    real(dp) :: identity_values(2), identity_derivatives(2)
    integer :: info

    call evaluate_normal_basis(config, 2, 0.4_dp, 0.1_dp, 0.25_dp, &
        values, derivatives, info)
    call require(info == radial_space_ok, "default radial space failed")
    call require(all(values == [0.75_dp, 0.25_dp]), &
        "P1 normal values are wrong")
    call require(all(derivatives == [-10.0_dp, 10.0_dp]), &
        "P1 normal derivatives are wrong")

    config%form_policy = form_s_power_edge
    call evaluate_normal_basis(config, 2, 0.4_dp, 0.1_dp, 0.25_dp, &
        values, derivatives, info)
    call require(info == radial_space_ok, "form-function radial space failed")
    step = 1.0e-7_dp
    call evaluate_normal_basis(config, 2, 0.4_dp + step, 0.1_dp, &
        0.25_dp + step / 0.1_dp, plus, derivatives, info)
    call require(info == radial_space_ok, "positive difference point failed")
    call evaluate_normal_basis(config, 2, 0.4_dp - step, 0.1_dp, &
        0.25_dp - step / 0.1_dp, minus, derivatives, info)
    call require(info == radial_space_ok, "negative difference point failed")
    call evaluate_normal_basis(config, 2, 0.4_dp, 0.1_dp, 0.25_dp, &
        values, derivatives, info)
    call require(maxval(abs(derivatives - (plus - minus) / (2.0_dp * step))) &
        < 1.0e-8_dp, "form-function product rule is wrong")

    config = radial_space_config_t()
    call evaluate_normal_basis(config, 1, 0.4_dp, 0.1_dp, 0.25_dp, &
        identity_values, identity_derivatives, info)
    call evaluate_normal_basis(config, 1, 0.4_dp, 0.1_dp, 0.25_dp, &
        values, derivatives, info, 0.0_dp)
    call require(info == radial_space_ok, "zero stored power was rejected")
    call require(all(values == identity_values), &
        "zero stored power changes normal values")
    call require(all(derivatives == identity_derivatives), &
        "zero stored power changes normal derivatives")
    call evaluate_normal_basis(config, 1, 0.4_dp + step, 0.1_dp, &
        0.25_dp + step / 0.1_dp, plus, derivatives, info, 0.25_dp)
    call require(info == radial_space_ok, "positive stored-power point failed")
    call evaluate_normal_basis(config, 1, 0.4_dp - step, 0.1_dp, &
        0.25_dp - step / 0.1_dp, minus, derivatives, info, 0.25_dp)
    call require(info == radial_space_ok, "negative stored-power point failed")
    call evaluate_normal_basis(config, 1, 0.4_dp, 0.1_dp, 0.25_dp, &
        values, derivatives, info, 0.25_dp)
    call require(info == radial_space_ok, "stored-power radial space failed")
    call require(maxval(abs(values - 0.4_dp**(-0.25_dp) &
        * identity_values)) < 1.0e-15_dp, &
        "stored normal variable has the wrong power direction")
    call require(maxval(abs(derivatives - (plus - minus) &
        / (2.0_dp * step))) < 1.0e-8_dp, &
        "stored normal variable product rule is wrong")
    call evaluate_normal_basis(config, 1, 0.4_dp + step, 0.1_dp, &
        0.25_dp + step / 0.1_dp, plus, derivatives, info, -0.25_dp)
    call require(info == radial_space_ok, &
        "positive negative-power point failed")
    call evaluate_normal_basis(config, 1, 0.4_dp - step, 0.1_dp, &
        0.25_dp - step / 0.1_dp, minus, derivatives, info, -0.25_dp)
    call require(info == radial_space_ok, &
        "negative negative-power point failed")
    call evaluate_normal_basis(config, 1, 0.4_dp, 0.1_dp, 0.25_dp, &
        values, derivatives, info, -0.25_dp)
    call require(info == radial_space_ok, "negative stored power failed")
    call require(maxval(abs(values - 0.4_dp**0.25_dp * identity_values)) &
        < 1.0e-15_dp, "negative stored power has the wrong direction")
    call require(maxval(abs(derivatives - (plus - minus) &
        / (2.0_dp * step))) < 1.0e-8_dp, &
        "negative stored power product rule is wrong")
    call evaluate_normal_basis(config, 1, 0.0_dp, 0.1_dp, 0.0_dp, &
        values, derivatives, info, 0.25_dp)
    call require(info == radial_space_invalid, &
        "nonzero stored power was evaluated on the axis")

    config%normal_degree = 2
    call evaluate_normal_basis(config, 2, 0.4_dp, 0.1_dp, 0.25_dp, &
        values, derivatives, info)
    call require(info == radial_space_invalid, &
        "unsupported radial degree was accepted")
    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_radial_space_policy
