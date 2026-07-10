program test_gvec_cas3d_integrals
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_integrals, only: integrate_half_mesh_volume, &
        integration_invalid_grid, integration_invalid_orientation, &
        integration_nonfinite_input, integration_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, radial_grid_full, &
        radial_grid_half
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    real(dp) :: full_device_volume, signed_period_volume
    integer :: info

    equilibrium%field_periods = 5
    equilibrium%radial_grid = radial_grid_half
    equilibrium%s = [0.125_dp, 0.375_dp, 0.625_dp, 0.875_dp]
    equilibrium%poloidal_modes = [0, 1]
    equilibrium%toroidal_modes = [0, 1, -1]
    allocate (equilibrium%jacobian%cosine(4, 2, 3))
    equilibrium%jacobian%cosine = 0.0_dp
    equilibrium%jacobian%cosine(:, 1, 1) = &
        [-2.0_dp, -4.0_dp, -6.0_dp, -8.0_dp]
    call integrate_half_mesh_volume(equilibrium, signed_period_volume, &
        full_device_volume, info)
    call require(info == integration_ok, "half-mesh volume failed")
    call require(abs(signed_period_volume + 5.0_dp) < 1.0e-14_dp, &
        "signed one-period volume is wrong")
    call require(abs(full_device_volume - 25.0_dp) < 1.0e-14_dp, &
        "positive full-device volume is wrong")

    equilibrium%radial_grid = radial_grid_full
    call integrate_half_mesh_volume(equilibrium, signed_period_volume, &
        full_device_volume, info)
    call require(info == integration_invalid_grid, "full grid was accepted")
    equilibrium%radial_grid = radial_grid_half
    equilibrium%s(1) = 0.0_dp
    call integrate_half_mesh_volume(equilibrium, signed_period_volume, &
        full_device_volume, info)
    call require(info == integration_invalid_grid, &
        "invalid half mesh was accepted")
    equilibrium%s(1) = 0.125_dp
    equilibrium%jacobian%cosine(:, 1, 1) = 2.0_dp
    call integrate_half_mesh_volume(equilibrium, signed_period_volume, &
        full_device_volume, info)
    call require(info == integration_invalid_orientation, &
        "right-handed volume was accepted")
    equilibrium%jacobian%cosine(1, 1, 1) = &
        ieee_value(0.0_dp, ieee_quiet_nan)
    call integrate_half_mesh_volume(equilibrium, signed_period_volume, &
        full_device_volume, info)
    call require(info == integration_nonfinite_input, &
        "nonfinite Jacobian was accepted")
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

end program test_gvec_cas3d_integrals
