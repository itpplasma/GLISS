program test_three_component_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use physical_constants, only: vacuum_permeability
    use three_component_kernel, only: compressible_divergence, &
        three_component_density, three_component_energy_jvp, &
        three_component_energy_vjp, three_component_input_count
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

    call check_energy_actions()

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

    subroutine check_energy_actions()
        real(dp) :: points(three_component_input_count, 2), inputs(26)
        real(dp) :: tangent(26), energy, derivative, magnetic, fluid, scale
        integer :: point, orientation

        points(:, 1) = [1.2_dp, 0.7_dp, 0.4_dp, -0.3_dp, 0.8_dp, &
            0.6_dp, 1.1_dp, 1.4_dp, 1.3_dp, 0.2_dp, -0.15_dp, &
            0.3_dp, -0.4_dp, 0.5_dp, 0.9_dp, 0.6_dp, -0.2_dp, &
            0.3_dp, -0.4_dp, 0.5_dp, -0.3_dp, 0.2_dp, -0.4_dp, &
            0.3_dp, -0.2_dp, 0.5_dp]
        points(:, 2) = [-0.9_dp, 1.3_dp, -0.2_dp, 0.5_dp, -0.7_dp, &
            0.4_dp, 1.6_dp, 0.9_dp, 0.8_dp, -0.3_dp, 0.25_dp, &
            -0.5_dp, 0.2_dp, -0.4_dp, 0.7_dp, -0.7_dp, 0.4_dp, &
            -0.5_dp, 0.2_dp, -0.6_dp, 0.7_dp, -0.5_dp, 0.6_dp, &
            -0.3_dp, 0.4_dp, -0.2_dp]
        do point = 1, 2
            do orientation = -1, 1, 2
                inputs = points(:, point)
                inputs(7) = real(orientation, dp) * inputs(7)
                call check_energy_point(inputs)

                ! Cancel the radial and tangential divergence to near zero.
                inputs(22) = -(inputs(1) * inputs(23) &
                    - inputs(2) * inputs(24) + inputs(2) * inputs(25) &
                    + inputs(1) * inputs(26)) / sum(inputs(1:2)**2) &
                    + 1.0e-8_dp
                tangent = 0.0_dp
                tangent(15) = 1.0_dp
                call three_component_energy_jvp(inputs, tangent, energy, derivative)
                call direct_energy_parts(inputs, magnetic, fluid)
                scale = abs(inputs(7)) * (1.0e-8_dp / inputs(7))**2
                call require(abs(derivative - scale) < 1.0e-22_dp, &
                    "near-solenoidal gamma derivative loses its small square")
                call require(abs(energy - magnetic - fluid) &
                    < 64.0_dp * epsilon(energy) * max(1.0_dp, abs(magnetic)), &
                    "near-solenoidal energy disagrees with the direct form")
            end do
        end do
    end subroutine check_energy_actions

    subroutine check_energy_point(inputs)
        real(dp), intent(in) :: inputs(three_component_input_count)
        real(dp), parameter :: steps(3) = [2.0e-4_dp, 2.0e-5_dp, 2.0e-6_dp]
        real(dp), parameter :: cotangent = -0.7_dp
        real(dp) :: gradient(26), weighted_gradient(26), tangent(26)
        real(dp) :: plus(26), minus(26), energy, reverse_energy, value
        real(dp) :: magnetic, fluid, magnetic_plus, magnetic_minus
        real(dp) :: fluid_plus, fluid_minus, derivative, estimate, step
        real(dp) :: tolerance, error, best_error, difference_scale
        integer :: coordinate, level

        call three_component_energy_vjp(inputs, 1.0_dp, energy, gradient)
        call direct_energy_parts(inputs, magnetic, fluid)
        call require(abs(energy - magnetic - fluid) &
            < 3.0e-14_dp * max(1.0_dp, abs(energy)), &
            "energy action value disagrees with the direct physical form")
        do coordinate = 1, three_component_input_count
            tangent = 0.0_dp
            tangent(coordinate) = 1.0_dp
            call three_component_energy_jvp(inputs, tangent, value, derivative)
            call require(abs(value - energy) < 3.0e-14_dp &
                * max(1.0_dp, abs(energy)), "JVP and VJP values disagree")
            call require(abs(derivative - gradient(coordinate)) < 2.0e-12_dp &
                * max(1.0_dp, abs(derivative)), "coordinate JVP/VJP disagree")
            best_error = huge(1.0_dp)
            do level = 1, size(steps)
                step = steps(level) * max(1.0_dp, abs(inputs(coordinate)))
                plus = inputs
                minus = inputs
                plus(coordinate) = plus(coordinate) + step
                minus(coordinate) = minus(coordinate) - step
                call direct_energy_parts(plus, magnetic_plus, fluid_plus)
                call direct_energy_parts(minus, magnetic_minus, fluid_minus)
                ! Keep the small fluid response independent of the SI magnetic scale.
                estimate = ((magnetic_plus - magnetic_minus) &
                    + (fluid_plus - fluid_minus)) / (2.0_dp * step)
                difference_scale = max(1.0_dp, abs(fluid_plus), abs(fluid_minus))
                if (coordinate /= 15 .and. coordinate < 22) then
                    difference_scale = max(difference_scale, &
                        abs(magnetic_plus), abs(magnetic_minus))
                end if
                error = abs(estimate - derivative)
                best_error = min(best_error, error)
                tolerance = 3.0e-7_dp * max(1.0_dp, abs(derivative)) &
                    + 32.0_dp * epsilon(energy) * difference_scale / step
                if (error > tolerance) then
                    write (error_unit, "(a,2i4,3es20.10)") &
                        "coordinate FD failure: ", coordinate, level, &
                        derivative, estimate, tolerance
                    error stop 1
                end if
            end do
            call require(best_error < 3.0e-8_dp * max(1.0_dp, abs(derivative)), &
                "coordinate derivative has no finite-difference plateau")
        end do
        do coordinate = 1, three_component_input_count
            tangent(coordinate) = (-1.0_dp)**coordinate &
                * real(coordinate + 1, dp) / real(three_component_input_count, dp)
        end do
        call three_component_energy_jvp(inputs, tangent, value, derivative)
        call three_component_energy_vjp(inputs, cotangent, reverse_energy, &
            weighted_gradient)
        call require(abs(cotangent * derivative &
            - dot_product(tangent, weighted_gradient)) < 3.0e-12_dp &
            * max(1.0_dp, abs(derivative)), "mixed JVP/VJP duality fails")
        call require(maxval(abs(weighted_gradient - cotangent * gradient)) &
            < 3.0e-12_dp * max(1.0_dp, maxval(abs(gradient))), &
            "VJP does not scale with the energy cotangent")
    end subroutine check_energy_point

    pure subroutine direct_energy_parts(inputs, magnetic, fluid)
        real(dp), intent(in) :: inputs(three_component_input_count)
        real(dp), intent(out) :: magnetic, fluid
        real(dp) :: angular_xi, angular_eta, shear, compression, bending
        real(dp) :: divergence, flux_norm_squared

        angular_xi = inputs(2) * inputs(18) + inputs(1) * inputs(19)
        angular_eta = inputs(2) * inputs(20) + inputs(1) * inputs(21)
        shear = angular_eta - (inputs(1) * inputs(4) &
            - inputs(3) * inputs(2)) * inputs(16) &
            + inputs(10) * inputs(7) * inputs(16) / inputs(9) &
            + inputs(12) * inputs(8) * angular_xi / inputs(9)
        compression = inputs(6) * inputs(21) - inputs(5) * inputs(20) &
            - (inputs(1) * inputs(5) + inputs(2) * inputs(6)) * inputs(17) &
            - (inputs(6) * inputs(4) + inputs(5) * inputs(3)) * inputs(16) &
            - inputs(11) * inputs(7) * inputs(16) + inputs(13) * angular_xi
        bending = angular_xi / (inputs(7) * sqrt(inputs(9)))
        shear = -sqrt(inputs(9)) * shear / (inputs(8) * inputs(7))
        compression = compression / (inputs(8) * inputs(7))
        flux_norm_squared = sum(inputs(1:2)**2)
        divergence = (inputs(22) + (inputs(1) * inputs(23) &
            - inputs(2) * inputs(24) + inputs(2) * inputs(25) &
            + inputs(1) * inputs(26)) / flux_norm_squared) / inputs(7)
        magnetic = (bending**2 + shear**2 + compression**2 &
            - inputs(14) * inputs(16)**2) &
            * abs(inputs(7)) / vacuum_permeability
        fluid = inputs(15) * divergence**2 * abs(inputs(7))
    end subroutine direct_energy_parts

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_three_component_kernel
