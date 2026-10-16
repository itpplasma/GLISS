program test_two_component_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use two_component_kernel, only: two_component_components, &
        two_component_density, two_component_coefficients, &
        two_component_energy_jvp, two_component_energy_vjp, &
        two_component_input_count
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: len = 6.0_dp * pi
    real(dp), parameter :: b1 = 0.3_dp, b3 = 0.4_dp, bz = 1.0_dp

    call check_reference(0.3_dp, 0.7_dp, -0.2_dp, 0.5_dp, -0.3_dp, &
        0.4_dp, 0.1_dp, 1.0822536130248883e-2_dp, &
        -4.9458926773973732e-1_dp, -2.3075713033433682_dp)
    call check_reference(0.5_dp, -0.4_dp, 0.3_dp, -0.25_dp, 0.2_dp, &
        -0.5_dp, 0.3_dp, -5.3051647697298445e-3_dp, &
        3.2939224238475452e-1_dp, 6.0633033498439228e-1_dp)
    call check_gradient()
    call check_cartesian_chart()
    call check_energy_actions()
    write (*, "(a)") "PASS"

contains

    pure function b_theta(radius) result(value)
        real(dp), intent(in) :: radius
        real(dp) :: value

        value = b1 * radius + b3 * radius**3
    end function b_theta

    subroutine kernel_inputs(radius, inputs)
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: inputs(12)
        real(dp) :: r

        r = radius
        inputs(1) = 2.0_dp * pi * r * bz
        inputs(2) = len * b_theta(r)
        inputs(3) = 2.0_dp * pi * bz
        inputs(4) = len * (b1 + 3.0_dp * b3 * r**2)
        inputs(5) = len * bz
        inputs(6) = 2.0_dp * pi * r * b_theta(r)
        inputs(7) = 2.0_dp * pi * len * r
        inputs(8) = sqrt(b_theta(r)**2 + bz**2)
        inputs(9) = 1.0_dp
        inputs(10) = (2.0_dp * b1 * r + 4.0_dp * b3 * r**3) * bz / r
        inputs(11) = -b_theta(r) * (2.0_dp * b1 * r &
            + 4.0_dp * b3 * r**3) / r
        inputs(12) = 0.0_dp
    end subroutine kernel_inputs

    subroutine check_reference(radius, xi_s, xi_s_s, xi_s_theta, &
            xi_s_zeta, eta_theta, eta_zeta, c1_ref, c2_ref, c3_ref)
        real(dp), intent(in) :: radius, xi_s, xi_s_s, xi_s_theta
        real(dp), intent(in) :: xi_s_zeta, eta_theta, eta_zeta
        real(dp), intent(in) :: c1_ref, c2_ref, c3_ref
        real(dp) :: inputs(12), c1, c2, c3

        call kernel_inputs(radius, inputs)
        call two_component_components(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta, &
            c1, c2, c3)
        call require(abs(c1 - c1_ref) < 1.0e-14_dp * max(1.0_dp, &
            abs(c1_ref)), "bending component is wrong")
        call require(abs(c2 - c2_ref) < 1.0e-14_dp * max(1.0_dp, &
            abs(c2_ref)), "shear component is wrong")
        call require(abs(c3 - c3_ref) < 1.0e-13_dp * max(1.0_dp, &
            abs(c3_ref)), "compression component is wrong")
    end subroutine check_reference

    subroutine check_gradient()
        real(dp) :: inputs(12), density_plus, density_minus, density
        real(dp) :: analytic, step, c1, c2, c3, drive

        call kernel_inputs(0.4_dp, inputs)
        drive = 2.0_dp * b_theta(0.4_dp) * (2.0_dp * b1 * 0.4_dp &
            + 4.0_dp * b3 * 0.4_dp**3) / 0.4_dp**2
        step = 1.0e-6_dp
        call two_component_density(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            drive, 0.6_dp + step, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, &
            -0.15_dp, density_plus)
        call two_component_density(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            drive, 0.6_dp - step, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, &
            -0.15_dp, density_minus)
        call two_component_density(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            drive, 0.6_dp, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, -0.15_dp, &
            density)
        call two_component_components(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            0.6_dp, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, -0.15_dp, &
            c1, c2, c3)
        analytic = (2.0_dp * c2 * shear_slope(inputs) &
            + 2.0_dp * c3 * compression_slope(inputs) &
            - 2.0_dp * drive * 0.6_dp) * abs(inputs(7))
        call require(abs((density_plus - density_minus) &
            / (2.0_dp * step) - analytic) < 1.0e-5_dp * max(1.0_dp, &
            abs(analytic)), "kernel gradient does not match")
    end subroutine check_gradient

    pure function shear_slope(inputs) result(value)
        real(dp), intent(in) :: inputs(12)
        real(dp) :: value

        value = -(sqrt(inputs(9)) / (inputs(8) * inputs(7))) &
            * (-(inputs(1) * inputs(4) - inputs(3) * inputs(2)) &
            + inputs(10) * inputs(7) / inputs(9))
    end function shear_slope

    pure function compression_slope(inputs) result(value)
        real(dp), intent(in) :: inputs(12)
        real(dp) :: value

        value = (1.0_dp / (inputs(8) * inputs(7))) &
            * (-(inputs(6) * inputs(4) + inputs(5) * inputs(3)) &
            - inputs(11) * inputs(7))
    end function compression_slope

    subroutine check_cartesian_chart()
        real(dp), parameter :: toroidal = 1.2_dp, poloidal = 0.7_dp
        real(dp) :: inputs(20), components(3), expected(3), curl(3)
        real(dp) :: field(3), binormal(3), tangent(20), energy, derivative
        real(dp) :: bnorm, jacobian, scale
        integer :: orientation, fixture

        bnorm = sqrt(toroidal**2 + poloidal**2)
        do orientation = -1, 1, 2
            jacobian = real(orientation, dp)
            inputs = 0.0_dp
            inputs(1:2) = [toroidal, poloidal]
            inputs(5) = toroidal * jacobian
            inputs(6) = poloidal / jacobian
            inputs(7) = jacobian
            inputs(8) = bnorm
            inputs(9) = 1.0_dp
            field(1) = 0.0_dp
            field(2) = poloidal / (jacobian * bnorm)
            field(3) = toroidal / bnorm
            binormal(1) = 0.0_dp
            binormal(2) = -toroidal / bnorm
            binormal(3) = poloidal / (jacobian * bnorm)
            do fixture = 1, 2
                inputs(15:20) = [0.6_dp, -0.2_dp, 0.3_dp, -0.4_dp, &
                    0.5_dp, -0.3_dp]
                if (fixture == 2) then
                    inputs(17) = toroidal
                    inputs(18) = -poloidal + 1.0e-8_dp
                    inputs(19) = toroidal
                    inputs(20) = -poloidal + 1.0e-8_dp
                end if
                ! In the Cartesian chart r=(s,theta,j*zeta), compute curl(xi x B)
                ! before projecting onto the normal, binormal and magnetic field.
                curl(1) = (poloidal * inputs(17) + toroidal * inputs(18)) / jacobian
                curl(2) = (inputs(20) - poloidal * inputs(16)) / jacobian
                curl(3) = -inputs(19) - toroidal * inputs(16)
                expected(1) = curl(1)
                expected(2) = dot_product(curl, binormal)
                expected(3) = dot_product(curl, field)
                call direct_components(inputs, components)
                call require(maxval(abs(components - expected)) < 3.0e-15_dp, &
                    "direct components disagree with the Cartesian curl oracle")
                call two_component_components(inputs(1), inputs(2), inputs(3), &
                    inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
                    inputs(9), inputs(10), inputs(11), inputs(12), inputs(13), &
                    inputs(15), inputs(16), inputs(17), inputs(18), inputs(19), &
                    inputs(20), components(1), components(2), components(3))
                call require(maxval(abs(components - expected)) < 3.0e-15_dp, &
                    "generated components disagree with the Cartesian curl oracle")
                tangent = 0.0_dp
                tangent(18) = 1.0_dp
                call two_component_energy_jvp(inputs, tangent, energy, derivative)
                scale = max(1.0_dp, dot_product(curl, curl))
                call require(abs(energy - dot_product(curl, curl)) &
                    < 3.0e-14_dp * scale, &
                    "constant-chart energy disagrees with the Cartesian curl norm")
                call require(abs(derivative - 2.0_dp * curl(1) &
                    * toroidal / jacobian) < 3.0e-14_dp * scale, &
                    "constant-chart angular derivative is wrong")
            end do
        end do
    end subroutine check_cartesian_chart

    subroutine check_energy_actions()
        real(dp) :: points(two_component_input_count, 2), inputs(20)
        real(dp) :: components(3), coefficients(3, 6), reconstructed(3), energy, derivative
        real(dp) :: tangent(20), scale
        integer :: point, orientation

        points(:, 1) = [1.2_dp, 0.7_dp, 0.4_dp, -0.3_dp, 0.8_dp, &
            0.6_dp, 1.1_dp, 1.4_dp, 1.3_dp, 0.2_dp, -0.15_dp, &
            0.3_dp, -0.4_dp, 0.5_dp, 0.6_dp, -0.2_dp, 0.3_dp, &
            -0.4_dp, 0.5_dp, -0.3_dp]
        points(:, 2) = [-0.9_dp, 1.3_dp, -0.2_dp, 0.5_dp, -0.7_dp, &
            0.4_dp, 1.6_dp, 0.9_dp, 0.8_dp, -0.3_dp, 0.25_dp, &
            -0.5_dp, 0.2_dp, -0.4_dp, -0.7_dp, 0.4_dp, -0.5_dp, &
            0.2_dp, -0.6_dp, 0.7_dp]
        do point = 1, 2
            do orientation = -1, 1, 2
                inputs = points(:, point)
                inputs(7) = real(orientation, dp) * inputs(7)
                call check_energy_point(inputs)
                call direct_components(inputs, components)
                call two_component_coefficients(inputs(1:13), coefficients)
                reconstructed = matmul(coefficients, inputs(15:20))
                call require(maxval(abs(reconstructed - components)) < 3.0e-14_dp &
                    * max(1.0_dp, maxval(abs(components))), &
                    "component coefficients do not reconstruct the direct fields")

                ! Match the drive and positive squares nearly to marginality.
                scale = sum(components**2)
                inputs(14) = scale * (1.0_dp - 1.0e-9_dp) / inputs(15)**2
                tangent = 0.0_dp
                tangent(14) = 1.0_dp
                call two_component_energy_jvp(inputs, tangent, energy, derivative)
                call require(abs(energy - direct_energy(inputs)) &
                    < 128.0_dp * epsilon(energy) * max(1.0_dp, &
                    abs(inputs(7)) * scale), &
                    "near-marginal energy loses its cancellation scale")
                call require(abs(derivative + inputs(15)**2 * abs(inputs(7))) &
                    < 2.0e-14_dp, "near-marginal drive derivative is wrong")
            end do
        end do
    end subroutine check_energy_actions

    subroutine check_energy_point(inputs)
        real(dp), intent(in) :: inputs(two_component_input_count)
        real(dp), parameter :: steps(3) = [2.0e-4_dp, 2.0e-5_dp, 2.0e-6_dp]
        real(dp), parameter :: cotangent = -0.7_dp
        real(dp) :: gradient(20), weighted_gradient(20), tangent(20)
        real(dp) :: plus(20), minus(20), energy, reverse_energy, value
        real(dp) :: derivative, estimate, step, tolerance, error, best_error
        integer :: coordinate, level

        call two_component_energy_vjp(inputs, 1.0_dp, energy, gradient)
        call require(abs(energy - direct_energy(inputs)) &
            < 3.0e-14_dp * max(1.0_dp, abs(energy)), &
            "energy action value disagrees with the direct physical form")
        do coordinate = 1, two_component_input_count
            tangent = 0.0_dp
            tangent(coordinate) = 1.0_dp
            call two_component_energy_jvp(inputs, tangent, value, derivative)
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
                estimate = (direct_energy(plus) - direct_energy(minus)) &
                    / (2.0_dp * step)
                error = abs(estimate - derivative)
                best_error = min(best_error, error)
                tolerance = 3.0e-7_dp * max(1.0_dp, abs(derivative)) &
                    + 32.0_dp * epsilon(energy) * max(1.0_dp, abs(energy)) / step
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
        do coordinate = 1, two_component_input_count
            tangent(coordinate) = (-1.0_dp)**coordinate &
                * real(coordinate + 1, dp) / real(two_component_input_count, dp)
        end do
        call two_component_energy_jvp(inputs, tangent, value, derivative)
        call two_component_energy_vjp(inputs, cotangent, reverse_energy, &
            weighted_gradient)
        call require(abs(cotangent * derivative &
            - dot_product(tangent, weighted_gradient)) < 3.0e-12_dp &
            * max(1.0_dp, abs(derivative)), "mixed JVP/VJP duality fails")
        call require(maxval(abs(weighted_gradient - cotangent * gradient)) &
            < 3.0e-12_dp * max(1.0_dp, maxval(abs(gradient))), &
            "VJP does not scale with the energy cotangent")
    end subroutine check_energy_point

    pure subroutine direct_components(inputs, components)
        real(dp), intent(in) :: inputs(two_component_input_count)
        real(dp), intent(out) :: components(3)
        real(dp) :: angular_xi, angular_eta, shear, compression

        angular_xi = inputs(2) * inputs(17) + inputs(1) * inputs(18)
        angular_eta = inputs(2) * inputs(19) + inputs(1) * inputs(20)
        shear = angular_eta - (inputs(1) * inputs(4) &
            - inputs(3) * inputs(2)) * inputs(15) &
            + inputs(10) * inputs(7) * inputs(15) / inputs(9) &
            + inputs(12) * inputs(8) * angular_xi / inputs(9)
        compression = inputs(6) * inputs(20) - inputs(5) * inputs(19) &
            - (inputs(1) * inputs(5) + inputs(2) * inputs(6)) * inputs(16) &
            - (inputs(6) * inputs(4) + inputs(5) * inputs(3)) * inputs(15) &
            - inputs(11) * inputs(7) * inputs(15) + inputs(13) * angular_xi
        components(1) = angular_xi / (inputs(7) * sqrt(inputs(9)))
        components(2) = -sqrt(inputs(9)) * shear / (inputs(8) * inputs(7))
        components(3) = compression / (inputs(8) * inputs(7))
    end subroutine direct_components

    pure function direct_energy(inputs) result(value)
        real(dp), intent(in) :: inputs(two_component_input_count)
        real(dp) :: value, components(3)

        call direct_components(inputs, components)
        value = (sum(components**2) - inputs(14) * inputs(15)**2) * abs(inputs(7))
    end function direct_energy

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_two_component_kernel
