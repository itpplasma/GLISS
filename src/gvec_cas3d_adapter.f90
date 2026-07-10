module gvec_cas3d_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter, public :: adapter_ok = 0
    integer, parameter, public :: adapter_invalid_radius = 1
    integer, parameter, public :: adapter_invalid_orientation = 2
    integer, parameter, public :: adapter_invalid_field_periods = 3
    integer, parameter, public :: adapter_shape_mismatch = 4
    integer, parameter, public :: adapter_nonfinite_input = 5

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: two_pi = 2.0_dp * pi

    public :: cas3d_coordinate_jacobian
    public :: export_period_to_full_device_zeta
    public :: map_gvec_to_cas3d
    public :: reconstruct_fourier_scalar

contains

    pure subroutine map_gvec_to_cas3d(rho, theta_b, zeta_b, orientation, &
            s, theta, zeta_full, info)
        real(dp), intent(in) :: rho, theta_b, zeta_b
        integer, intent(in) :: orientation
        real(dp), intent(out) :: s, theta, zeta_full
        integer, intent(out) :: info

        s = 0.0_dp
        theta = 0.0_dp
        zeta_full = 0.0_dp
        info = adapter_nonfinite_input
        if (.not. ieee_is_finite(rho)) return
        if (.not. ieee_is_finite(theta_b)) return
        if (.not. ieee_is_finite(zeta_b)) return

        info = adapter_invalid_radius
        if (rho < 0.0_dp) return
        if (rho > 1.0_dp) return

        info = adapter_invalid_orientation
        if (orientation /= -1) then
            if (orientation /= 1) return
        end if

        s = rho**2
        theta = -real(orientation, dp) * theta_b / two_pi
        zeta_full = real(orientation, dp) * zeta_b / two_pi
        info = adapter_ok
    end subroutine map_gvec_to_cas3d

    pure subroutine export_period_to_full_device_zeta(zeta_period, &
            field_periods, zeta_full, info)
        real(dp), intent(in) :: zeta_period
        integer, intent(in) :: field_periods
        real(dp), intent(out) :: zeta_full
        integer, intent(out) :: info

        zeta_full = 0.0_dp
        info = adapter_nonfinite_input
        if (.not. ieee_is_finite(zeta_period)) return

        info = adapter_invalid_field_periods
        if (field_periods < 1) return

        zeta_full = zeta_period / real(field_periods, dp)
        info = adapter_ok
    end subroutine export_period_to_full_device_zeta

    pure function cas3d_coordinate_jacobian(rho) result(jacobian)
        real(dp), intent(in) :: rho
        real(dp) :: jacobian

        jacobian = -rho / (2.0_dp * pi**2)
    end function cas3d_coordinate_jacobian

    pure subroutine reconstruct_fourier_scalar(theta, zeta_full, field_periods, &
            poloidal_modes, toroidal_modes, cosine_coefficients, &
            sine_coefficients, value, derivative_theta, derivative_zeta_full, &
            info)
        real(dp), intent(in) :: theta, zeta_full
        integer, intent(in) :: field_periods
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        real(dp), intent(in) :: cosine_coefficients(:), sine_coefficients(:)
        real(dp), intent(out) :: value, derivative_theta, derivative_zeta_full
        integer, intent(out) :: info
        real(dp) :: phase, phase_theta, phase_zeta, phase_derivative
        integer :: mode

        value = 0.0_dp
        derivative_theta = 0.0_dp
        derivative_zeta_full = 0.0_dp
        info = adapter_nonfinite_input
        if (.not. ieee_is_finite(theta)) return
        if (.not. ieee_is_finite(zeta_full)) return
        if (.not. all(ieee_is_finite(cosine_coefficients))) return
        if (.not. all(ieee_is_finite(sine_coefficients))) return

        info = adapter_invalid_field_periods
        if (field_periods < 1) return

        info = adapter_shape_mismatch
        if (size(toroidal_modes) /= size(poloidal_modes)) return
        if (size(cosine_coefficients) /= size(poloidal_modes)) return
        if (size(sine_coefficients) /= size(poloidal_modes)) return

        do mode = 1, size(poloidal_modes)
            phase_theta = two_pi * real(poloidal_modes(mode), dp)
            phase_zeta = -two_pi * real(toroidal_modes(mode), dp) * &
                real(field_periods, dp)
            phase = phase_theta * theta + phase_zeta * zeta_full
            phase_derivative = -cosine_coefficients(mode) * sin(phase) + &
                sine_coefficients(mode) * cos(phase)
            value = value + cosine_coefficients(mode) * cos(phase) + &
                sine_coefficients(mode) * sin(phase)
            derivative_theta = derivative_theta + &
                phase_theta * phase_derivative
            derivative_zeta_full = derivative_zeta_full + &
                phase_zeta * phase_derivative
        end do
        info = adapter_ok
    end subroutine reconstruct_fourier_scalar

end module gvec_cas3d_adapter
