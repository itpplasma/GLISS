program test_three_component_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use physical_constants, only: vacuum_permeability
    use three_component_kernel, only: compressible_divergence, &
        three_component_density
    use two_component_kernel, only: two_component_density
    implicit none

    real(dp), parameter :: ft = 1.1_dp, fp = -0.7_dp, signed_sqrtg = -1.3_dp
    real(dp), parameter :: xi_s_radial = 0.17_dp
    real(dp), parameter :: eta_theta = -0.23_dp, eta_zeta = 0.19_dp
    real(dp), parameter :: mu_theta = 0.11_dp, mu_zeta = -0.09_dp
    real(dp) :: divergence, direct, full_density, reduced_density
    real(dp) :: gamma_pressure, expected_difference

    call compressible_divergence(ft, fp, signed_sqrtg, &
        signed_sqrtg * xi_s_radial, signed_sqrtg * eta_theta, &
        signed_sqrtg * eta_zeta, mu_theta, mu_zeta, divergence)
    direct = direct_constant_metric_divergence(ft, fp, signed_sqrtg, &
        xi_s_radial, eta_theta, eta_zeta, mu_theta, mu_zeta)
    call require(abs(divergence - direct) < 1.0e-14_dp, &
        "compressible divergence does not reconstruct the vector divergence")

    call evaluate_densities(0.0_dp, 0.0_dp, full_density, reduced_density)
    call require(abs(full_density - reduced_density / vacuum_permeability) &
        < 1.0e-10_dp, &
        "zero-gamma physical energy has the wrong vacuum-permeability scale")
    gamma_pressure = 0.9_dp
    call evaluate_densities(gamma_pressure, 1.0_dp, full_density, &
        reduced_density)
    expected_difference = gamma_pressure * divergence**2 * abs(signed_sqrtg)
    call require(abs(full_density - reduced_density / vacuum_permeability &
        - expected_difference) < 1.0e-10_dp, &
        "fluid-compression energy has the wrong SI scale")

    write (*, "(a)") "PASS"

contains

    pure function direct_constant_metric_divergence(ft, fp, jacobian, &
            xi_radial, eta_theta, eta_zeta, mu_theta, mu_zeta) result(value)
        real(dp), intent(in) :: ft, fp, jacobian, xi_radial
        real(dp), intent(in) :: eta_theta, eta_zeta, mu_theta, mu_zeta
        real(dp) :: value, flux_norm_squared
        real(dp) :: xi_theta_theta, xi_zeta_zeta

        flux_norm_squared = ft**2 + fp**2
        xi_theta_theta = (ft * eta_theta + fp * mu_theta / jacobian) &
            / flux_norm_squared
        xi_zeta_zeta = (-fp * eta_zeta + ft * mu_zeta / jacobian) &
            / flux_norm_squared
        value = xi_radial + xi_theta_theta + xi_zeta_zeta
    end function direct_constant_metric_divergence

    subroutine evaluate_densities(gamma_pressure, mu_scale, full, reduced)
        real(dp), intent(in) :: gamma_pressure, mu_scale
        real(dp), intent(out) :: full, reduced
        real(dp) :: sqrtg_xi_s_radial, sqrtg_eta_theta, sqrtg_eta_zeta

        sqrtg_xi_s_radial = signed_sqrtg * xi_s_radial
        sqrtg_eta_theta = signed_sqrtg * eta_theta
        sqrtg_eta_zeta = signed_sqrtg * eta_zeta
        call three_component_density(ft, fp, 0.04_dp, -0.03_dp, 0.8_dp, &
            0.6_dp, signed_sqrtg, 1.4_dp, 1.3_dp, 0.2_dp, -0.15_dp, &
            0.1_dp, -0.12_dp, 0.05_dp, gamma_pressure, 0.3_dp, &
            xi_s_radial, -0.2_dp, 0.1_dp, eta_theta, eta_zeta, &
            sqrtg_xi_s_radial, sqrtg_eta_theta, sqrtg_eta_zeta, &
            mu_scale * mu_theta, mu_scale * mu_zeta, full)
        call two_component_density(ft, fp, 0.04_dp, -0.03_dp, 0.8_dp, &
            0.6_dp, signed_sqrtg, 1.4_dp, 1.3_dp, 0.2_dp, -0.15_dp, &
            0.1_dp, -0.12_dp, 0.05_dp, 0.3_dp, xi_s_radial, -0.2_dp, &
            0.1_dp, eta_theta, eta_zeta, reduced)
    end subroutine evaluate_densities

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_three_component_kernel
