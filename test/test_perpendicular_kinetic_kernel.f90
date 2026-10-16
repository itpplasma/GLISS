program test_perpendicular_kinetic_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use physical_mass_kernel, only: physical_mass_matrix
    use perpendicular_kinetic_kernel, only: perpendicular_kinetic_energy, &
        perpendicular_kinetic_matrix, perpendicular_kinetic_energy_jvp, &
        perpendicular_kinetic_energy_vjp, perpendicular_kinetic_input_count
    implicit none

    real(dp), parameter :: flux_t_slope = 1.2_dp
    real(dp), parameter :: flux_p_slope = 0.7_dp
    real(dp), parameter :: current_i = 0.8_dp
    real(dp), parameter :: current_j = 0.6_dp
    real(dp), parameter :: signed_sqrtg = -1.1_dp
    real(dp), parameter :: bmag = 1.4_dp
    real(dp), parameter :: grad_s2 = 1.3_dp
    real(dp), parameter :: sigma_tilde = 0.2_dp
    real(dp), parameter :: beta_tilde = -0.15_dp
    real(dp), parameter :: density = 2.3_dp
    real(dp) :: perpendicular(2, 2), physical(3, 3), expected(2, 2)
    real(dp) :: reflected(2, 2)
    real(dp) :: parallel(2), displacement(2), energy, image(2), scale
    real(dp) :: delta, flux_norm_squared
    integer :: column, row

    call perpendicular_kinetic_matrix(signed_sqrtg, bmag, grad_s2, &
        sigma_tilde, density, perpendicular)
    scale = density * abs(signed_sqrtg)
    expected(1, 1) = scale * (1.0_dp + sigma_tilde**2) / grad_s2
    expected(1, 2) = scale * sigma_tilde / bmag
    expected(2, 1) = expected(1, 2)
    expected(2, 2) = scale * grad_s2 / bmag**2
    call require_close(perpendicular, expected, 2.0e-14_dp, &
        "perpendicular kinetic Gram matrix is wrong")

    call physical_mass_matrix(flux_t_slope, flux_p_slope, current_i, &
        current_j, signed_sqrtg, bmag, grad_s2, sigma_tilde, beta_tilde, &
        density, physical)
    flux_norm_squared = flux_t_slope**2 + flux_p_slope**2
    delta = current_i * flux_p_slope - current_j * flux_t_slope
    parallel(1) = beta_tilde / bmag
    parallel(2) = -delta / (bmag * flux_norm_squared)
    do column = 1, 2
        do row = 1, 2
            expected(row, column) = physical(row, column) &
                - scale * parallel(row) * parallel(column)
        end do
    end do
    call require_close(perpendicular, expected, 2.0e-14_dp, &
        "perpendicular norm does not remove only the parallel square")

    call require(abs(perpendicular(1, 1) * perpendicular(2, 2) &
        - perpendicular(1, 2)**2 &
        - scale**2 / bmag**2) < 2.0e-14_dp * scale**2, &
        "perpendicular kinetic determinant is wrong")
    call perpendicular_kinetic_matrix(-signed_sqrtg, bmag, grad_s2, &
        sigma_tilde, density, reflected)
    call require_close(reflected, perpendicular, 2.0e-14_dp, &
        "perpendicular norm depends on orientation")
    call perpendicular_kinetic_matrix(signed_sqrtg, bmag, grad_s2, &
        -sigma_tilde, density, reflected)
    call require(reflected(1, 1) == perpendicular(1, 1) .and. &
        reflected(2, 2) == perpendicular(2, 2) .and. &
        reflected(1, 2) == -perpendicular(1, 2) .and. &
        reflected(2, 1) == -perpendicular(2, 1), &
        "sigma sign does not change only the kinetic coupling")
    displacement = [0.4_dp, -0.3_dp]
    energy = perpendicular_kinetic_energy(signed_sqrtg, bmag, grad_s2, &
        sigma_tilde, density, displacement)
    image(1) = perpendicular(1, 1) * displacement(1) &
        + perpendicular(1, 2) * displacement(2)
    image(2) = perpendicular(2, 1) * displacement(1) &
        + perpendicular(2, 2) * displacement(2)
    call require(abs(energy - 0.5_dp * dot_product(displacement, image)) &
        < 2.0e-14_dp * max(1.0_dp, energy), &
        "perpendicular kinetic energy is inconsistent with its matrix")

    call check_derivatives()
    call check_tangential_cancellation()

    write (*, "(a)") "PASS"

contains

    subroutine check_derivatives()
        real(dp), parameter :: points(7, 2) = reshape([ &
            1.1_dp, 1.4_dp, 1.3_dp, 0.2_dp, 2.3_dp, 0.4_dp, -0.3_dp, &
            0.85_dp, 0.9_dp, 2.1_dp, -0.45_dp, 0.75_dp, -0.65_dp, 0.55_dp], &
            [7, 2])
        real(dp), parameter :: steps(3) = [3.0e-4_dp, 1.0e-4_dp, 3.0e-5_dp]
        real(dp), parameter :: mixed_direction(7) = &
            [-0.13_dp, 0.17_dp, -0.09_dp, 0.11_dp, -0.19_dp, 0.07_dp, 0.21_dp]
        real(dp) :: inputs(perpendicular_kinetic_input_count)
        real(dp) :: tangent(perpendicular_kinetic_input_count)
        real(dp) :: gradient(perpendicular_kinetic_input_count)
        real(dp) :: plus(perpendicular_kinetic_input_count)
        real(dp) :: minus(perpendicular_kinetic_input_count)
        real(dp) :: value, reverse_value, derivative, cotangent, h, finite_difference
        real(dp) :: errors(3), derivative_scale, radial, tangential, gram_energy
        integer :: point, orientation, argument, step

        call require(perpendicular_kinetic_input_count == 7, &
            "perpendicular kinetic packed input count is wrong")
        do point = 1, size(points, 2)
            do orientation = -1, 1, 2
                inputs = points(:, point)
                inputs(1) = real(orientation, dp) * inputs(1)
                radial = inputs(6) / sqrt(inputs(3))
                tangential = inputs(4) * radial &
                    + sqrt(inputs(3)) * inputs(7) / inputs(2)
                gram_energy = 0.5_dp * inputs(5) * abs(inputs(1)) &
                    * (radial**2 + tangential**2)
                cotangent = -0.7_dp + 0.2_dp * real(point, dp)
                call perpendicular_kinetic_energy_vjp(inputs, cotangent, &
                    reverse_value, gradient)
                call require_scalar_close(reverse_value, original_energy(inputs), &
                    2.0e-14_dp, "perpendicular VJP value differs from original energy")
                call require_scalar_close(reverse_value, gram_energy, 2.0e-14_dp, &
                    "perpendicular VJP value differs from orthonormal Gram energy")
                do argument = 1, perpendicular_kinetic_input_count
                    tangent = 0.0_dp
                    tangent(argument) = 1.0_dp
                    call perpendicular_kinetic_energy_jvp(inputs, tangent, &
                        value, derivative)
                    call require_scalar_close(value, gram_energy, 2.0e-14_dp, &
                        "perpendicular coordinate JVP value is wrong")
                    call require_scalar_close(gradient(argument), cotangent * derivative, &
                        4.0e-13_dp, "perpendicular coordinate JVP/VJP duality failed")
                    derivative_scale = max(1.0_dp, abs(derivative))
                    do step = 1, size(steps)
                        h = steps(step) * max(1.0_dp, abs(inputs(argument)))
                        plus = inputs
                        minus = inputs
                        plus(argument) = plus(argument) + h
                        minus(argument) = minus(argument) - h
                        finite_difference = (original_energy(plus) &
                            - original_energy(minus)) / (2.0_dp * h)
                        errors(step) = abs(finite_difference - derivative)
                    end do
                    call require(maxval(errors) < 4.0e-7_dp * derivative_scale, &
                        "perpendicular coordinate derivative failed a centered FD step")
                    call require(errors(size(steps)) < 2.0e-8_dp * derivative_scale, &
                        "perpendicular coordinate derivative lacks a fine-step FD match")
                end do
                call perpendicular_kinetic_energy_jvp(inputs, mixed_direction, &
                    value, derivative)
                call require_scalar_close(value, original_energy(inputs), 2.0e-14_dp, &
                    "perpendicular mixed JVP value differs from original energy")
                call require_scalar_close(dot_product(gradient, mixed_direction), &
                    cotangent * derivative, 4.0e-13_dp, &
                    "perpendicular mixed JVP/VJP duality failed")
                do step = 1, size(steps)
                    h = steps(step)
                    plus = inputs + h * mixed_direction
                    minus = inputs - h * mixed_direction
                    finite_difference = (original_energy(plus) &
                        - original_energy(minus)) / (2.0_dp * h)
                    call require_scalar_close(derivative, finite_difference, 4.0e-8_dp, &
                        "perpendicular mixed derivative failed a centered FD step")
                end do
                call require_scalar_close(gradient(1), &
                    cotangent * gram_energy / inputs(1), 4.0e-13_dp, &
                    "perpendicular Jacobian derivative has the wrong orientation sign")
            end do
        end do
    end subroutine check_derivatives

    subroutine check_tangential_cancellation()
        real(dp), parameter :: offset = 2.0e-6_dp
        real(dp) :: inputs(perpendicular_kinetic_input_count)
        real(dp) :: tangent(perpendicular_kinetic_input_count)
        real(dp) :: gradient(perpendicular_kinetic_input_count)
        real(dp) :: value, reverse_value, derivative, expected_derivative, radial
        real(dp) :: expected_energy
        integer :: orientation

        do orientation = -1, 1, 2
            inputs(1) = real(orientation, dp) * 1.1_dp
            inputs(2) = 1.4_dp
            inputs(3) = 1.3_dp
            inputs(4) = -0.7_dp
            inputs(5) = 2.3_dp
            inputs(6) = 0.8_dp
            inputs(7) = 0.0_dp
            inputs(7) = -inputs(4) * inputs(6) * inputs(2) / inputs(3) + offset
            radial = inputs(6) / sqrt(inputs(3))
            expected_energy = 0.5_dp * inputs(5) * abs(inputs(1)) &
                * (radial**2 + inputs(3) * offset**2 / inputs(2)**2)
            expected_derivative = inputs(5) * abs(inputs(1)) &
                * inputs(3) * offset / inputs(2)**2
            tangent = 0.0_dp
            tangent(7) = 1.0_dp
            call perpendicular_kinetic_energy_jvp(inputs, tangent, value, derivative)
            call perpendicular_kinetic_energy_vjp(inputs, 1.0_dp, reverse_value, gradient)
            call require_scalar_close(value, expected_energy, 2.0e-14_dp, &
                "perpendicular energy is wrong near tangential cancellation")
            call require_scalar_close(reverse_value, original_energy(inputs), &
                2.0e-14_dp, "perpendicular cancellation VJP value is wrong")
            call require_scalar_close(derivative, expected_derivative, 2.0e-13_dp, &
                "perpendicular JVP loses the nearly cancelled tangential derivative")
            call require_scalar_close(gradient(7), expected_derivative, 2.0e-13_dp, &
                "perpendicular VJP loses the nearly cancelled tangential derivative")
        end do
    end subroutine check_tangential_cancellation

    function original_energy(inputs) result(value)
        real(dp), intent(in) :: inputs(perpendicular_kinetic_input_count)
        real(dp) :: value

        value = perpendicular_kinetic_energy(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6:7))
    end function original_energy

    subroutine require_scalar_close(first, second, tolerance, message)
        real(dp), intent(in) :: first, second, tolerance
        character(len=*), intent(in) :: message

        call require(abs(first - second) <= tolerance * max(1.0_dp, abs(second)), message)
    end subroutine require_scalar_close

    subroutine require_close(first, second, tolerance, message)
        real(dp), intent(in) :: first(:, :), second(:, :), tolerance
        character(len=*), intent(in) :: message

        call require(maxval(abs(first - second)) <= tolerance &
            * max(1.0_dp, maxval(abs(second))), message)
    end subroutine require_close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_perpendicular_kinetic_kernel
