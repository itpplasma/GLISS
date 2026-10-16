module two_component_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gliss_energy_generated, only: generated_bending, generated_shear, &
        generated_compression, generated_components, generated_component_coefficients, &
        generated_two_energy, generated_two_jvp, generated_two_vjp
    implicit none
    private

    public :: bending_component_value
    public :: compression_component_value
    public :: shear_component_value
    public :: two_component_components
    public :: two_component_component_values
    public :: two_component_density

    public :: two_component_coefficients, two_component_energy_jvp, two_component_energy_vjp
    integer, parameter, public :: two_component_input_count = 20

contains

    pure function bending_component_value(flux_t_slope, flux_p_slope, &
            sqrtg, grad_s2, xi_s_theta, xi_s_zeta) result(component)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope, sqrtg, grad_s2
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta
        real(dp) :: component, values(1)

        call generated_bending(flux_t_slope, flux_p_slope, sqrtg, grad_s2, &
            xi_s_theta, xi_s_zeta, values)
        component = values(1)
    end function bending_component_value

    pure function shear_component_value(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, sqrtg, bmag, grad_s2, j_dot_b, &
            sigma_tilde, xi_s, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta) &
            result(component)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, sigma_tilde, xi_s
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta, eta_theta, eta_zeta
        real(dp) :: component, values(1)

        call generated_shear(flux_t_slope, flux_p_slope, flux_t_curve, flux_p_curve, &
            sqrtg, bmag, grad_s2, j_dot_b, sigma_tilde, xi_s, xi_s_theta, xi_s_zeta, &
            eta_theta, eta_zeta, values)
        component = values(1)
    end function shear_component_value

    pure function compression_component_value(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, sqrtg, bmag, &
            pressure_slope, beta_tilde, xi_s, xi_s_s, xi_s_theta, &
            xi_s_zeta, eta_theta, eta_zeta) result(component)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: pressure_slope, beta_tilde, xi_s, xi_s_s
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta, eta_theta, eta_zeta
        real(dp) :: component, values(1)

        call generated_compression(flux_t_slope, flux_p_slope, flux_t_curve, &
            flux_p_curve, current_i, current_j, sqrtg, bmag, pressure_slope, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta, &
            values)
        component = values(1)
    end function compression_component_value

    pure function two_component_component_values(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            sqrtg, bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta) result(components)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: sigma_tilde, beta_tilde
        real(dp), intent(in) :: xi_s, xi_s_s, xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp) :: components(3)

        call generated_components(flux_t_slope, flux_p_slope, flux_t_curve, &
            flux_p_curve, current_i, current_j, sqrtg, bmag, grad_s2, j_dot_b, &
            pressure_slope, sigma_tilde, beta_tilde, xi_s, xi_s_s, xi_s_theta, &
            xi_s_zeta, eta_theta, eta_zeta, components)
    end function two_component_component_values

    pure subroutine two_component_components(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, sqrtg, &
            bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta, c_bending, c_shear, c_compression)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: sigma_tilde, beta_tilde
        real(dp), intent(in) :: xi_s, xi_s_s, xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp), intent(out) :: c_bending, c_shear, c_compression
        real(dp) :: components(3)

        components = two_component_component_values(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            sqrtg, bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta)
        c_bending = components(1)
        c_shear = components(2)
        c_compression = components(3)
    end subroutine two_component_components

    pure subroutine two_component_density(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, sqrtg, &
            bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, drive_a, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, &
            eta_theta, eta_zeta, density)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: sigma_tilde, beta_tilde, drive_a
        real(dp), intent(in) :: xi_s, xi_s_s, xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp), intent(out) :: density
        real(dp) :: values(1)

        call generated_two_energy(flux_t_slope, flux_p_slope, flux_t_curve, &
            flux_p_curve, current_i, current_j, sqrtg, bmag, grad_s2, j_dot_b, &
            pressure_slope, sigma_tilde, beta_tilde, drive_a, xi_s, xi_s_s, &
            xi_s_theta, xi_s_zeta, eta_theta, eta_zeta, values)
        density = values(1)
    end subroutine two_component_density

    ! Fields follow two_component_components up to beta_tilde. Columns act on
    ! (xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta).
    pure subroutine two_component_coefficients(fields, coefficients)
        real(dp), intent(in) :: fields(13)
        real(dp), intent(out) :: coefficients(3, 6)

        call generated_component_coefficients(fields(1), fields(2), fields(3), &
            fields(4), fields(5), fields(6), fields(7), fields(8), fields(9), &
            fields(10), fields(11), fields(12), fields(13), coefficients)
    end subroutine two_component_coefficients

    ! Packed inputs follow the corresponding scalar energy/density argument order.
    ! Domain: nonzero signed Jacobian, positive B and grad_s2; flux norm nonzero.
    pure subroutine two_component_energy_jvp(inputs, tangent, energy, derivative)
        real(dp), intent(in) :: inputs(20), tangent(20)
        real(dp), intent(out) :: energy, derivative
        real(dp) :: values(2)

        call generated_two_jvp(inputs(1), inputs(2), inputs(3), inputs(4), &
            inputs(5), inputs(6), inputs(7), inputs(8), inputs(9), inputs(10), &
            inputs(11), inputs(12), inputs(13), inputs(14), inputs(15), inputs(16), &
            inputs(17), inputs(18), inputs(19), inputs(20), tangent(1), tangent(2), &
            tangent(3), tangent(4), tangent(5), tangent(6), tangent(7), tangent(8), &
            tangent(9), tangent(10), tangent(11), tangent(12), tangent(13), &
            tangent(14), tangent(15), tangent(16), tangent(17), tangent(18), &
            tangent(19), tangent(20), values)
        energy = values(1)
        derivative = values(2)
    end subroutine two_component_energy_jvp

    pure subroutine two_component_energy_vjp(inputs, cotangent, energy, gradient)
        real(dp), intent(in) :: inputs(20), cotangent
        real(dp), intent(out) :: energy, gradient(20)
        real(dp) :: values(21)

        call generated_two_vjp(inputs(1), inputs(2), inputs(3), inputs(4), &
            inputs(5), inputs(6), inputs(7), inputs(8), inputs(9), inputs(10), &
            inputs(11), inputs(12), inputs(13), inputs(14), inputs(15), inputs(16), &
            inputs(17), inputs(18), inputs(19), inputs(20), cotangent, values)
        energy = values(1)
        gradient = values(2:)
    end subroutine two_component_energy_vjp

end module two_component_kernel
