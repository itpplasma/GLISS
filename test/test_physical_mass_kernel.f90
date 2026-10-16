program test_physical_mass_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use physical_mass_kernel, only: physical_mass_energy, physical_mass_matrix, &
        physical_mass_energy_jvp, physical_mass_energy_vjp, physical_mass_input_count
    implicit none

    real(dp), parameter :: metric(3, 3) = reshape([ &
        1.6_dp, 0.2_dp, -0.1_dp, &
        0.2_dp, 1.3_dp, 0.15_dp, &
        -0.1_dp, 0.15_dp, 1.1_dp], [3, 3])
    real(dp), parameter :: b_contravariant(2) = [0.4_dp, -0.7_dp]
    real(dp), parameter :: variables(3) = [0.3_dp, -0.2_dp, 0.15_dp]
    real(dp), parameter :: density_kg_m3 = 2.4_dp
    real(dp) :: mass(3, 3), scaled_mass(3, 3), direct_energy, scalar_energy
    real(dp) :: metric_image(3), xi_contravariant(3)
    real(dp) :: flux_t_slope, flux_p_slope
    real(dp) :: current_i, current_j, signed_sqrtg, bmag, grad_s2
    real(dp) :: signed_sigma_tilde, beta_tilde, determinant
    real(dp) :: principal_minor

    determinant = determinant_three(metric)
    signed_sqrtg = -sqrt(determinant)
    flux_p_slope = signed_sqrtg * b_contravariant(1)
    flux_t_slope = signed_sqrtg * b_contravariant(2)
    current_j = metric(2, 2) * b_contravariant(1) &
        + metric(2, 3) * b_contravariant(2)
    current_i = metric(3, 2) * b_contravariant(1) &
        + metric(3, 3) * b_contravariant(2)
    beta_tilde = metric(1, 2) * b_contravariant(1) &
        + metric(1, 3) * b_contravariant(2)
    bmag = sqrt(b_contravariant(1) * current_j &
        + b_contravariant(2) * current_i)
    principal_minor = metric(2, 2) * metric(3, 3) - metric(2, 3)**2
    grad_s2 = principal_minor / determinant
    signed_sigma_tilde = (current_i * metric(1, 2) &
        - current_j * metric(1, 3)) / (signed_sqrtg * bmag)

    call physical_mass_matrix(flux_t_slope, flux_p_slope, current_i, &
        current_j, signed_sqrtg, bmag, grad_s2, signed_sigma_tilde, &
        beta_tilde, &
        density_kg_m3, mass)
    call reconstruct_displacement(flux_t_slope, flux_p_slope, signed_sqrtg, &
        variables, xi_contravariant)
    call require(abs(flux_t_slope * xi_contravariant(2) &
        - flux_p_slope * xi_contravariant(3) - variables(2)) < 1.0e-14_dp, &
        "reconstruction uses the wrong production eta convention")
    call require(abs(signed_sqrtg * (flux_t_slope * xi_contravariant(3) &
        + flux_p_slope * xi_contravariant(2)) - variables(3)) < 1.0e-14_dp, &
        "reconstruction uses the wrong signed mu convention")
    call matrix_vector_product(metric, xi_contravariant, metric_image)
    direct_energy = 0.5_dp * density_kg_m3 * abs(signed_sqrtg) &
        * dot_product(xi_contravariant, metric_image)
    scalar_energy = physical_mass_energy(flux_t_slope, flux_p_slope, &
        current_i, current_j, signed_sqrtg, bmag, grad_s2, &
        signed_sigma_tilde, &
        beta_tilde, density_kg_m3, variables)

    call require(abs(scalar_energy - direct_energy) < 1.0e-13_dp, &
        "scalar mass form does not equal the vector kinetic energy")
    call require(maxval(abs(mass - transpose(mass))) < 1.0e-14_dp, &
        "physical mass matrix is not symmetric")
    call require(positive_definite_three(mass), &
        "physical mass matrix is not positive definite")
    call physical_mass_matrix(flux_t_slope, flux_p_slope, current_i, &
        current_j, signed_sqrtg, bmag, grad_s2, signed_sigma_tilde, &
        beta_tilde, &
        3.0_dp * density_kg_m3, scaled_mass)
    call require(maxval(abs(scaled_mass - 3.0_dp * mass)) < 1.0e-13_dp, &
        "physical mass matrix is not linear in SI mass density")

    call check_metric_orientations()
    call check_derivatives()

    write (*, "(a)") "PASS"

contains

    subroutine check_metric_orientations()
        real(dp), parameter :: alternate_metric(3, 3) = reshape([ &
            2.2_dp, -0.3_dp, 0.25_dp, &
            -0.3_dp, 1.1_dp, -0.17_dp, &
            0.25_dp, -0.17_dp, 1.7_dp], [3, 3])
        real(dp), parameter :: alternate_field(2) = [-0.8_dp, 0.25_dp]
        integer :: orientation

        do orientation = -1, 1, 2
            call check_metric_state(metric, b_contravariant, orientation)
            call check_metric_state(alternate_metric, alternate_field, orientation)
        end do
    end subroutine check_metric_orientations

    subroutine check_metric_state(g, field, orientation)
        real(dp), intent(in) :: g(3, 3), field(2)
        integer, intent(in) :: orientation
        real(dp) :: ft, fp, ci, cj, jacobian, field_norm, gs2, sigma, beta
        real(dp) :: xi(3), image(3), local_mass(3, 3), actual, expected

        jacobian = real(orientation, dp)*sqrt(determinant_three(g))
        fp = jacobian*field(1)
        ft = jacobian*field(2)
        cj = g(2, 2)*field(1) + g(2, 3)*field(2)
        ci = g(3, 2)*field(1) + g(3, 3)*field(2)
        beta = g(1, 2)*field(1) + g(1, 3)*field(2)
        field_norm = sqrt(field(1)*cj + field(2)*ci)
        gs2 = (g(2, 2)*g(3, 3) - g(2, 3)**2)/determinant_three(g)
        sigma = (ci*g(1, 2) - cj*g(1, 3))/(jacobian*field_norm)
        call reconstruct_displacement(ft, fp, jacobian, variables, xi)
        call matrix_vector_product(g, xi, image)
        expected = 0.5_dp*density_kg_m3*abs(jacobian)*dot_product(xi, image)
        actual = physical_mass_energy(ft, fp, ci, cj, jacobian, field_norm, &
            gs2, sigma, beta, density_kg_m3, variables)
        call physical_mass_matrix(ft, fp, ci, cj, jacobian, field_norm, &
            gs2, sigma, beta, density_kg_m3, local_mass)
        call require(abs(actual - expected) < 2.0e-13_dp*expected, &
            'generated energy violates the independent metric oracle')
        call matrix_vector_product(local_mass, variables, image)
        call require(abs(0.5_dp*dot_product(variables, image) - expected) &
            < 2.0e-13_dp*expected, 'generated matrix violates the metric oracle')
        call require(positive_definite_three(local_mass), &
            'generated mass lost positivity on an admissible metric')
    end subroutine check_metric_state

    subroutine check_derivatives()
        real(dp), parameter :: steps(3) = [1.0e-3_dp, 1.0e-4_dp, 1.0e-5_dp]
        real(dp) :: inputs(physical_mass_input_count), tangent(physical_mass_input_count)
        real(dp) :: perturbed(physical_mass_input_count)
        real(dp) :: gradient(physical_mass_input_count), opposite(physical_mass_input_count)
        real(dp) :: value, forward, reverse_value, difference, plus, minus, step, scale
        real(dp) :: cotangent, reference, opposite_value
        integer :: point, orientation, entry, index

        do point = 1, 2
            do orientation = -1, 1, 2
                inputs = [1.2_dp, 0.7_dp, 0.8_dp, 0.6_dp, 1.1_dp, 1.4_dp, &
                    1.3_dp, 0.2_dp, -0.15_dp, 2.3_dp, 0.4_dp, -0.3_dp, 0.2_dp]
                do entry = 1, size(inputs)
                    inputs(entry) = inputs(entry)*(1.0_dp &
                        + 0.17_dp*real(point - 1, dp)*sin(real(entry, dp)))
                    tangent(entry) = 0.13_dp*cos(real(entry*point, dp))
                end do
                inputs(5) = real(orientation, dp)*inputs(5)
                reference = packed_energy(inputs)
                cotangent = -0.37_dp*real(point, dp)
                call physical_mass_energy_jvp(inputs, tangent, value, forward)
                call physical_mass_energy_vjp(inputs, cotangent, reverse_value, gradient)
                call require(abs(value - reference) < 2.0e-13_dp*reference, &
                    'JVP primal differs from physical energy')
                call require(abs(reverse_value - reference) < 2.0e-13_dp*reference, &
                    'VJP primal differs from physical energy')
                scale = max(abs(cotangent*forward), abs(dot_product(gradient, tangent)))
                call require(abs(dot_product(gradient, tangent) - cotangent*forward) &
                    < 2.0e-13_dp*max(scale, reference), 'mass JVP/VJP duality failed')
                call physical_mass_energy_vjp(inputs, 1.0_dp, value, gradient)
                perturbed = inputs
                perturbed(5) = -inputs(5)
                call physical_mass_energy_vjp(perturbed, 1.0_dp, opposite_value, opposite)
                call require(abs(opposite_value - value) < 2.0e-13_dp*value, &
                    'mass energy depends on the Jacobian orientation')
                opposite(5) = -opposite(5)
                call require(maxval(abs(opposite - gradient)) &
                    < 2.0e-13_dp*maxval(abs(gradient)), &
                    'mass derivatives use the wrong abs(J) branch')
                do entry = 1, size(inputs)
                    do index = 1, size(steps)
                        step = steps(index)*max(1.0_dp, abs(inputs(entry)))
                        perturbed = inputs
                        perturbed(entry) = inputs(entry) + step
                        plus = packed_energy(perturbed)
                        perturbed(entry) = inputs(entry) - step
                        minus = packed_energy(perturbed)
                        difference = (plus - minus)/(2.0_dp*step)
                        scale = max(abs(gradient(entry)), abs(difference), &
                            1.0e-3_dp*reference/(1.0_dp + abs(inputs(entry))))
                        call require(abs(difference - gradient(entry)) &
                            < 32.0_dp*steps(index)**2*scale &
                            + 128.0_dp*epsilon(1.0_dp)*reference/step, &
                            'mass coordinate derivative step study failed')
                    end do
                end do
                do index = 1, size(steps)
                    step = steps(index)
                    perturbed = inputs + step*tangent
                    plus = packed_energy(perturbed)
                    perturbed = inputs - step*tangent
                    minus = packed_energy(perturbed)
                    difference = (plus - minus)/(2.0_dp*step)
                    scale = max(abs(forward), reference)
                    call require(abs(difference - forward) &
                        < 32.0_dp*step**2*scale &
                        + 128.0_dp*epsilon(1.0_dp)*reference/step, &
                        'mass directional JVP step study failed')
                end do
            end do
        end do
    end subroutine check_derivatives

    pure function packed_energy(inputs) result(energy)
        real(dp), intent(in) :: inputs(physical_mass_input_count)
        real(dp) :: energy

        energy = physical_mass_energy(inputs(1), inputs(2), inputs(3), inputs(4), &
            inputs(5), inputs(6), inputs(7), inputs(8), inputs(9), inputs(10), inputs(11:13))
    end function packed_energy

    pure subroutine matrix_vector_product(matrix, vector, image)
        real(dp), intent(in) :: matrix(3, 3), vector(3)
        real(dp), intent(out) :: image(3)
        integer :: column, row

        do row = 1, 3
            image(row) = 0.0_dp
            do column = 1, 3
                image(row) = image(row) + matrix(row, column) * vector(column)
            end do
        end do
    end subroutine matrix_vector_product

    pure subroutine reconstruct_displacement(ft, fp, jacobian, values, xi)
        real(dp), intent(in) :: ft, fp, jacobian, values(3)
        real(dp), intent(out) :: xi(3)
        real(dp) :: flux_norm_squared

        flux_norm_squared = ft**2 + fp**2
        xi(1) = values(1)
        xi(2) = (ft * values(2) + fp * values(3) / jacobian) &
            / flux_norm_squared
        xi(3) = (-fp * values(2) + ft * values(3) / jacobian) &
            / flux_norm_squared
    end subroutine reconstruct_displacement

    pure function determinant_three(matrix) result(determinant)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp) :: determinant

        determinant = matrix(1, 1) * (matrix(2, 2) * matrix(3, 3) &
            - matrix(2, 3) * matrix(3, 2)) &
            - matrix(1, 2) * (matrix(2, 1) * matrix(3, 3) &
            - matrix(2, 3) * matrix(3, 1)) &
            + matrix(1, 3) * (matrix(2, 1) * matrix(3, 2) &
            - matrix(2, 2) * matrix(3, 1))
    end function determinant_three

    pure function positive_definite_three(matrix) result(positive)
        real(dp), intent(in) :: matrix(3, 3)
        logical :: positive

        positive = matrix(1, 1) > 0.0_dp
        positive = positive .and. matrix(1, 1) * matrix(2, 2) &
            - matrix(1, 2)**2 > 0.0_dp
        positive = positive .and. determinant_three(matrix) > 0.0_dp
    end function positive_definite_three

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_physical_mass_kernel
