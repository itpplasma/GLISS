program test_dense_generalized_inverse_iteration
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dense_generalized_inverse_iteration, only: &
        dense_inverse_invalid, dense_inverse_mass_not_spd, dense_inverse_ok, &
        solve_dense_generalized_near_shift, &
        solve_dense_generalized_subspace_near_shift
    implicit none

    real(dp) :: backward_error, eigenvalue, equilibrated_residual, mass(3, 3)
    real(dp) :: mass_reciprocal_condition, residual, standard_residual_norm
    real(dp) :: stiffness(3, 3)
    real(dp) :: block_initial(6, 3), block_mass(6, 6)
    real(dp) :: block_stiffness(6, 6)
    real(dp), allocatable :: block_eigenvalues(:), block_eigenvalues_repeat(:)
    real(dp), allocatable :: block_residuals(:), block_residuals_repeat(:)
    real(dp), allocatable :: block_vectors(:, :), block_vectors_repeat(:, :)
    real(dp), allocatable :: vector(:)
    integer :: block_iterations, block_iterations_repeat, info, iterations

    stiffness = 0.0_dp
    mass = 0.0_dp
    mass(1, 1) = 1.0e-12_dp
    mass(2, 2) = 1.0_dp
    mass(3, 3) = 1.0e12_dp
    stiffness(1, 1) = 3.0e-15_dp
    stiffness(2, 2) = 2.0e-2_dp
    stiffness(3, 3) = 5.0e11_dp
    call solve_dense_generalized_near_shift(stiffness, mass, 2.9e-3_dp, &
        eigenvalue, vector, residual, equilibrated_residual, backward_error, &
        standard_residual_norm, mass_reciprocal_condition, iterations, info)
    call require(info == dense_inverse_ok, "ill-scaled solve failed")
    call require(abs(eigenvalue - 3.0e-3_dp) < 1.0e-14_dp, &
        "ill-scaled eigenvalue differs")
    call require(ieee_is_finite(residual), &
        "ill-scaled raw residual is not finite")
    call require(ieee_is_finite(equilibrated_residual), &
        "ill-scaled equilibrated residual is not finite")
    call require(backward_error <= 1.0e-12_dp, &
        "ill-scaled backward error is too large")
    call require(iterations >= 2 .and. iterations < 20, &
        "ill-scaled iteration count differs")
    call require(mass_reciprocal_condition > 0.99_dp, &
        "diagonally equilibrated mass condition differs")

    mass = 0.0_dp
    stiffness = 0.0_dp
    mass(1, 1) = 4.0_dp
    mass(1, 2) = 2.0_dp
    mass(2, 1) = 2.0_dp
    mass(2, 2) = 10.0_dp
    mass(2, 3) = 1.5_dp
    mass(3, 2) = 1.5_dp
    mass(3, 3) = 16.25_dp
    stiffness(1, 1) = 0.012_dp
    stiffness(1, 2) = 0.006_dp
    stiffness(2, 1) = 0.006_dp
    stiffness(2, 2) = 0.183_dp
    stiffness(2, 3) = 0.03_dp
    stiffness(3, 2) = 0.03_dp
    stiffness(3, 3) = 8.005_dp
    call solve_dense_generalized_near_shift(stiffness, mass, 2.9e-3_dp, &
        eigenvalue, vector, residual, equilibrated_residual, backward_error, &
        standard_residual_norm, mass_reciprocal_condition, iterations, info)
    call require(info == dense_inverse_ok, "coupled SPD solve failed")
    call require(abs(eigenvalue - 3.0e-3_dp) < 1.0e-14_dp, &
        "coupled SPD eigenvalue differs")
    call require(ieee_is_finite(residual), &
        "coupled SPD raw residual is not finite")
    call require(backward_error <= 1.0e-12_dp, &
        "coupled SPD backward error is too large")

    call solve_dense_generalized_near_shift(stiffness(:, :2), mass, &
        2.9e-3_dp, eigenvalue, vector, residual, equilibrated_residual, &
        backward_error, standard_residual_norm, mass_reciprocal_condition, &
        iterations, info)
    call require(info == dense_inverse_invalid, "nonsquare problem was accepted")
    call solve_dense_generalized_near_shift(stiffness, mass, &
        ieee_value(0.0_dp, ieee_quiet_nan), eigenvalue, vector, residual, &
        equilibrated_residual, backward_error, standard_residual_norm, &
        mass_reciprocal_condition, iterations, info)
    call require(info == dense_inverse_invalid, "nonfinite shift was accepted")
    stiffness(1, 2) = stiffness(1, 2) + 1.0e-3_dp
    call solve_dense_generalized_near_shift(stiffness, mass, 2.9e-3_dp, &
        eigenvalue, vector, residual, equilibrated_residual, backward_error, &
        standard_residual_norm, mass_reciprocal_condition, iterations, info)
    call require(info == dense_inverse_invalid, &
        "nonsymmetric stiffness was accepted")
    stiffness(1, 2) = stiffness(1, 2) - 1.0e-3_dp
    stiffness(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call solve_dense_generalized_near_shift(stiffness, mass, 2.9e-3_dp, &
        eigenvalue, vector, residual, equilibrated_residual, backward_error, &
        standard_residual_norm, mass_reciprocal_condition, iterations, info)
    call require(info == dense_inverse_invalid, &
        "nonfinite stiffness was accepted")
    stiffness(1, 1) = 0.012_dp
    mass(1, 1) = -1.0_dp
    call solve_dense_generalized_near_shift(stiffness, mass, 2.9e-3_dp, &
        eigenvalue, vector, residual, equilibrated_residual, backward_error, &
        standard_residual_norm, mass_reciprocal_condition, iterations, info)
    call require(info == dense_inverse_mass_not_spd, &
        "negative mass diagonal was accepted")
    mass(1, 1) = 1.0_dp
    mass(1, 2) = 2.0_dp
    mass(2, 1) = 2.0_dp
    mass(2, 2) = 1.0_dp
    mass(2, 3) = 0.0_dp
    mass(3, 2) = 0.0_dp
    mass(3, 3) = 1.0_dp
    call solve_dense_generalized_near_shift(stiffness, mass, 2.9e-3_dp, &
        eigenvalue, vector, residual, equilibrated_residual, backward_error, &
        standard_residual_norm, mass_reciprocal_condition, iterations, info)
    call require(info == dense_inverse_mass_not_spd, &
        "indefinite mass matrix was accepted")

    block_mass = 0.0_dp
    block_stiffness = 0.0_dp
    block_mass(1, 1) = 1.0_dp
    block_mass(2, 2) = 2.0_dp
    block_mass(3, 3) = 4.0_dp
    block_mass(4, 4) = 8.0_dp
    block_mass(5, 5) = 16.0_dp
    block_mass(6, 6) = 32.0_dp
    block_stiffness(1, 1) = 1.0_dp
    block_stiffness(2, 2) = 4.0_dp
    block_stiffness(3, 3) = 12.0_dp
    block_stiffness(4, 4) = 32.0_dp
    block_stiffness(5, 5) = 80.0_dp
    block_stiffness(6, 6) = 192.0_dp
    block_initial = 0.0_dp
    block_initial(1, 1) = 1.0_dp
    block_initial(4, 1) = 0.3_dp
    block_initial(5, 1) = 0.2_dp
    block_initial(2, 2) = 1.0_dp
    block_initial(4, 2) = 0.2_dp
    block_initial(6, 2) = -0.1_dp
    block_initial(3, 3) = 1.0_dp
    block_initial(5, 3) = 0.1_dp
    block_initial(6, 3) = 0.2_dp
    call solve_dense_generalized_subspace_near_shift(block_stiffness, &
        block_mass, 0.0_dp, block_initial, 80, block_eigenvalues, &
        block_vectors, block_residuals, block_iterations, info)
    call require(info == dense_inverse_ok, "block inverse iteration failed")
    call require(block_iterations == 80, "block iteration count differs")
    call require(maxval(abs(block_eigenvalues - [1.0_dp, 2.0_dp, 3.0_dp])) &
        < 1.0e-12_dp, "block eigenvalues differ")
    call require(maxval(block_residuals) < 1.0e-9_dp, &
        "block residuals are too large")
    call solve_dense_generalized_subspace_near_shift(block_stiffness, &
        block_mass, 0.0_dp, block_initial, 80, block_eigenvalues_repeat, &
        block_vectors_repeat, block_residuals_repeat, &
        block_iterations_repeat, info)
    call require(info == dense_inverse_ok, "repeated block iteration failed")
    call require(block_iterations_repeat == block_iterations, &
        "repeated block iteration count differs")
    call require(all(block_eigenvalues_repeat == block_eigenvalues), &
        "repeated block eigenvalues are not bitwise identical")
    call require(all(block_vectors_repeat == block_vectors), &
        "repeated block vectors are not bitwise identical")
    call require(all(block_residuals_repeat == block_residuals), &
        "repeated block residuals are not bitwise identical")

    block_initial(:, 2) = block_initial(:, 1)
    call solve_dense_generalized_subspace_near_shift(block_stiffness, &
        block_mass, 0.0_dp, block_initial, 80, block_eigenvalues, &
        block_vectors, block_residuals, block_iterations, info)
    call require(info == dense_inverse_mass_not_spd, &
        "rank-deficient initial subspace was accepted")
    call solve_dense_generalized_subspace_near_shift(block_stiffness, &
        block_mass, 0.0_dp, block_initial(:, :2), 0, block_eigenvalues, &
        block_vectors, block_residuals, block_iterations, info)
    call require(info == dense_inverse_invalid, &
        "invalid block iteration limit was accepted")

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_dense_generalized_inverse_iteration
