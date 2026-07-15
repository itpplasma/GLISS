program test_physical_mass_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use physical_mass_kernel, only: physical_mass_energy, physical_mass_matrix
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

    write (*, "(a)") "PASS"

contains

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
