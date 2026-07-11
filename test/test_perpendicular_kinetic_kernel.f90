program test_perpendicular_kinetic_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use physical_mass_kernel, only: physical_mass_matrix
    use perpendicular_kinetic_kernel, only: perpendicular_kinetic_energy, &
        perpendicular_kinetic_matrix
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
    real(dp) :: parallel(2), displacement(2), energy, scale
    real(dp) :: delta, flux_norm_squared

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
    parallel = [beta_tilde / bmag, &
        -delta / (bmag * flux_norm_squared)]
    expected = physical(1:2, 1:2) &
        - scale * spread(parallel, 2, 2) * spread(parallel, 1, 2)
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
    call require(abs(energy - 0.5_dp * dot_product(displacement, &
        matmul(perpendicular, displacement))) &
        < 2.0e-14_dp * max(1.0_dp, energy), &
        "perpendicular kinetic energy is inconsistent with its matrix")

    write (*, "(a)") "PASS"

contains

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
