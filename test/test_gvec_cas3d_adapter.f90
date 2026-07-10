program test_gvec_cas3d_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use gvec_cas3d_adapter, only: adapter_invalid_field_periods, &
        adapter_invalid_orientation, adapter_invalid_radius, &
        adapter_nonfinite_input, adapter_ok, adapter_shape_mismatch, &
        cas3d_coordinate_jacobian, export_period_to_full_device_zeta, &
        map_gvec_to_cas3d, reconstruct_fourier_scalar
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tolerance = 2.0e-13_dp
    real(dp) :: s, theta, zeta_full, value, derivative_theta
    real(dp) :: derivative_zeta_full
    real(dp) :: shifted_value, unused_theta, unused_zeta_full
    integer :: info

    call map_gvec_to_cas3d(0.5_dp, pi / 3.0_dp, -pi / 5.0_dp, 1, &
        s, theta, zeta_full, info)
    call require(info == adapter_ok, "positive-orientation map failed")
    call require(abs(s - 0.25_dp) < tolerance, "radial flux map is wrong")
    call require(abs(theta + 1.0_dp / 6.0_dp) < tolerance, &
        "poloidal-angle map is wrong")
    call require(abs(zeta_full + 1.0_dp / 10.0_dp) < tolerance, &
        "toroidal-angle map is wrong")
    call require(cas3d_coordinate_jacobian(0.5_dp) < 0.0_dp, &
        "coordinate map did not reverse handedness")
    call require(abs(cas3d_coordinate_jacobian(0.5_dp) + &
        0.5_dp / (2.0_dp * pi**2)) < tolerance, &
        "coordinate-map Jacobian is wrong")

    call map_gvec_to_cas3d(0.5_dp, pi / 3.0_dp, -pi / 5.0_dp, -1, &
        s, theta, zeta_full, info)
    call require(info == adapter_ok, "negative-orientation map failed")
    call require(abs(theta - 1.0_dp / 6.0_dp) < tolerance, &
        "negative-orientation poloidal map is wrong")
    call require(abs(zeta_full - 1.0_dp / 10.0_dp) < tolerance, &
        "negative-orientation toroidal map is wrong")

    call test_invalid_coordinates()
    call export_period_to_full_device_zeta(0.5_dp, 5, zeta_full, info)
    call require(info == adapter_ok, "exported-zeta conversion failed")
    call require(abs(zeta_full - 0.1_dp) < tolerance, &
        "exported one-period zeta was not normalized")
    call reconstruct_fourier_scalar(1.0_dp / 8.0_dp, 1.0_dp / 40.0_dp, 5, &
        [0, 1, 2], [0, 1, -1], [2.0_dp, 3.0_dp, 4.0_dp], &
        [0.0_dp, 5.0_dp, 6.0_dp], value, derivative_theta, &
        derivative_zeta_full, info)
    call require(info == adapter_ok, "Fourier reconstruction failed")
    call require(abs(value - (5.0_dp + sqrt(2.0_dp))) < tolerance, &
        "Fourier value is wrong")
    call require(abs(derivative_theta - &
        (10.0_dp * pi - 20.0_dp * pi * sqrt(2.0_dp))) < tolerance, &
        "poloidal derivative is wrong")
    call require(abs(derivative_zeta_full + &
        50.0_dp * pi * (1.0_dp + sqrt(2.0_dp))) < tolerance, &
        "toroidal derivative is wrong")

    call reconstruct_fourier_scalar(1.0_dp / 8.0_dp, &
        1.0_dp / 40.0_dp + 1.0_dp / 5.0_dp, 5, [0, 1, 2], [0, 1, -1], &
        [2.0_dp, 3.0_dp, 4.0_dp], [0.0_dp, 5.0_dp, 6.0_dp], shifted_value, &
        unused_theta, unused_zeta_full, info)
    call require(info == adapter_ok, "period-shifted reconstruction failed")
    call require(abs(shifted_value - value) < tolerance, &
        "field-period shift changed the Fourier value")

    call test_invalid_fourier_inputs()
    write (*, "(a)") "PASS"

contains

    subroutine test_invalid_coordinates()
        real(dp) :: nan

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call map_gvec_to_cas3d(-0.1_dp, 0.0_dp, 0.0_dp, 1, &
            s, theta, zeta_full, info)
        call require(info == adapter_invalid_radius, &
            "negative radius was accepted")
        call map_gvec_to_cas3d(0.5_dp, 0.0_dp, 0.0_dp, 0, &
            s, theta, zeta_full, info)
        call require(info == adapter_invalid_orientation, &
            "invalid orientation was accepted")
        call export_period_to_full_device_zeta(0.0_dp, 0, zeta_full, info)
        call require(info == adapter_invalid_field_periods, &
            "invalid exported-zeta field periods were accepted")
        call map_gvec_to_cas3d(nan, 0.0_dp, 0.0_dp, 1, &
            s, theta, zeta_full, info)
        call require(info == adapter_nonfinite_input, &
            "nonfinite coordinate input was accepted")
    end subroutine test_invalid_coordinates

    subroutine test_invalid_fourier_inputs()
        real(dp) :: nan

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call reconstruct_fourier_scalar(0.0_dp, 0.0_dp, 0, [0], [0], &
            [1.0_dp], [0.0_dp], value, derivative_theta, &
            derivative_zeta_full, info)
        call require(info == adapter_invalid_field_periods, &
            "zero field periods were accepted")
        call reconstruct_fourier_scalar(0.0_dp, 0.0_dp, 1, [0, 1], [0], &
            [1.0_dp], [0.0_dp], value, derivative_theta, &
            derivative_zeta_full, info)
        call require(info == adapter_shape_mismatch, &
            "mismatched Fourier arrays were accepted")
        call reconstruct_fourier_scalar(0.0_dp, 0.0_dp, 1, [0], [0], &
            [nan], [0.0_dp], value, derivative_theta, &
            derivative_zeta_full, info)
        call require(info == adapter_nonfinite_input, &
            "nonfinite Fourier coefficient was accepted")
    end subroutine test_invalid_fourier_inputs

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_gvec_cas3d_adapter
