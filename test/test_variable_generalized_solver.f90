program test_variable_generalized_solver
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dense_spectrum_support, only: certify_dense_spectrum_inertia, &
        certify_dense_spectrum_orthogonality, dense_spectrum_ok, &
        diagnose_dense_spectrum, refine_dense_spectrum
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use stable_reduction, only: stable_norm2
    use symmetric_eigensolver, only: solve_symmetric_generalized, &
        solve_symmetric_generalized_allocated
    use variable_block_tridiagonal, only: &
        apply_variable_block_tridiagonal, pack_variable_blocks, &
        variable_block_to_dense, variable_block_tridiagonal_t
    use variable_generalized_equilibration, only: &
        equilibrate_variable_generalized, undo_variable_congruence, &
        variable_equilibration_ok
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_generalized_diagnostics, &
        variable_generalized_inertia, variable_generalized_invalid, &
        variable_generalized_mass_not_spd, variable_generalized_ok
    implicit none

    integer, parameter :: widths(3) = [2, 3, 1]
    type(variable_block_tridiagonal_t) :: stiffness, mass, corrupt
    real(dp) :: dense_k(6, 6), dense_m(6, 6), eigenvalue, quotient, residual
    real(dp) :: resolution, diagnostic_residual, diagnostic_resolution
    real(dp), allocatable :: bad_initial(:), eigenvalues(:)
    real(dp), allocatable :: eigenvectors(:, :), vector(:)
    real(dp) :: shift
    integer :: count, i, info

    call build_fixture(dense_k, dense_m)
    call solve_symmetric_generalized(dense_k, dense_m, eigenvalues, &
        eigenvectors, info)
    call require(info == 0, "dense variable generalized oracle failed")
    call pack_variable_blocks(dense_k, widths, stiffness, info)
    call require(info == 0, "variable stiffness packing failed")
    call pack_variable_blocks(dense_m, widths, mass, info)
    call require(info == 0, "variable mass packing failed")
    call check_equilibration(stiffness, mass, eigenvalues)
    call check_indexed_dense_refinement(stiffness, mass, eigenvalues, &
        eigenvectors)

    call variable_generalized_inertia(stiffness, mass, eigenvalues(1) - 1.0_dp, &
        count, info)
    call require(info == variable_generalized_ok .and. count == 0, &
        "variable generalized inertia below spectrum is wrong")
    do i = 1, size(eigenvalues) - 1
        shift = 0.5_dp * (eigenvalues(i) + eigenvalues(i + 1))
        call variable_generalized_inertia(stiffness, mass, shift, count, info)
        call require(info == variable_generalized_ok .and. count == i, &
            "variable generalized inertia between eigenvalues is wrong")
    end do
    call variable_generalized_inertia(stiffness, mass, eigenvalues(6) + 1.0_dp, &
        count, info)
    call require(info == variable_generalized_ok .and. count == 6, &
        "variable generalized inertia above spectrum is wrong")

    shift = eigenvalues(1) - 0.25_dp * (eigenvalues(2) - eigenvalues(1))
    call iterate_variable_generalized_eigenvalue(stiffness, mass, shift, &
        eigenvalue, vector, residual, resolution, info)
    call require(info == variable_generalized_ok, &
        "variable generalized iteration failed")
    call require(abs(eigenvalue - eigenvalues(1)) < 1.0e-11_dp, &
        "variable generalized eigenvalue disagrees with dense oracle")
    call require(abs(mass_norm(vector, mass) - 1.0_dp) < 2.0e-14_dp, &
        "variable generalized vector is not mass normalized")
    call require(residual < 1.0e-11_dp, &
        "variable generalized residual is too large")
    call variable_generalized_diagnostics(stiffness, mass, vector, eigenvalue, &
        quotient, diagnostic_residual, diagnostic_resolution, info)
    call require(info == variable_generalized_ok, &
        "variable generalized diagnostics failed")
    call require(abs(quotient - eigenvalue) < 1.0e-14_dp, &
        "variable generalized quotient is inconsistent")
    call require(abs(diagnostic_residual - residual) < 1.0e-14_dp, &
        "variable generalized diagnostic residual is inconsistent")
    call require(abs(diagnostic_resolution - resolution) < 1.0e-14_dp, &
        "variable generalized resolution is inconsistent")
    call iterate_variable_generalized_eigenvalue(stiffness, mass, shift, &
        eigenvalue, vector, residual, resolution, info, initial=[1.0_dp])
    call require(info == variable_generalized_invalid, &
        "wrong-sized inverse-iteration seed was accepted")
    allocate (bad_initial(6), source=ieee_value(0.0_dp, ieee_quiet_nan))
    call iterate_variable_generalized_eigenvalue(stiffness, mass, shift, &
        eigenvalue, vector, residual, resolution, info, initial=bad_initial)
    call require(info == variable_generalized_invalid, &
        "nonfinite inverse-iteration seed was accepted")
    call check_long_chain_resolution()
    call check_scaled_compensated_norm()
    call check_ill_scaled_dense_pencil()

    corrupt = stiffness
    corrupt%diagonal(1)%values(1, 2) = &
        corrupt%diagonal(1)%values(1, 2) + 0.1_dp
    call variable_generalized_inertia(corrupt, mass, 0.0_dp, count, info)
    call require(info == variable_generalized_invalid, &
        "nonsymmetric variable stiffness was accepted")
    dense_m(1, 1) = -1.0_dp
    call pack_variable_blocks(dense_m, widths, corrupt, info)
    call variable_generalized_inertia(stiffness, corrupt, 0.0_dp, count, info)
    call require(info == variable_generalized_mass_not_spd, &
        "indefinite variable mass was accepted")
    vector = ieee_value(0.0_dp, ieee_quiet_nan)
    call variable_generalized_diagnostics(stiffness, mass, vector, eigenvalue, &
        quotient, residual, resolution, info)
    call require(info == variable_generalized_invalid, &
        "nonfinite variable generalized vector was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine check_equilibration(stiffness, mass, expected)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: expected(:)
        type(variable_block_tridiagonal_t) :: balanced_k, balanced_m
        real(dp), allocatable :: dense_k(:, :), dense_m(:, :)
        real(dp), allocatable :: scales(:), values(:), vectors(:, :), vector(:)
        integer :: balanced_count, info, original_count

        call equilibrate_variable_generalized(stiffness, mass, balanced_k, &
            balanced_m, scales, info)
        call require(info == variable_equilibration_ok, &
            "variable generalized equilibration failed")
        call variable_block_to_dense(balanced_k, dense_k, info)
        call require(info == 0, "balanced stiffness unpacking failed")
        call variable_block_to_dense(balanced_m, dense_m, info)
        call require(info == 0, "balanced mass unpacking failed")
        call solve_symmetric_generalized(dense_k, dense_m, values, vectors, &
            info)
        call require(info == 0, "balanced dense oracle failed")
        call require(maxval(abs(values - expected)) < 1.0e-11_dp, &
            "congruence equilibration changed the spectrum")
        call variable_generalized_inertia(stiffness, mass, 0.0_dp, &
            original_count, info)
        call require(info == variable_generalized_ok, &
            "original inertia check failed")
        call variable_generalized_inertia(balanced_k, balanced_m, 0.0_dp, &
            balanced_count, info)
        call require(info == variable_generalized_ok &
            .and. balanced_count == original_count, &
            "congruence equilibration changed inertia")
        call undo_variable_congruence(scales, vectors(:, 1), vector, info)
        call require(info == variable_equilibration_ok, &
            "balanced eigenvector back-transform failed")
        call require(abs(mass_norm(vector, mass) - 1.0_dp) < 2.0e-14_dp, &
            "back-transformed vector is not mass normalized")
    end subroutine check_equilibration

    subroutine check_indexed_dense_refinement(stiffness, mass, values, vectors)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: values(:), vectors(:, :)
        type(fixed_boundary_solver_controls_t) :: controls
        real(dp), allocatable :: refined(:), refined_vectors(:, :)
        real(dp), allocatable :: rayleigh(:), residuals(:), resolutions(:)
        integer :: info

        allocate (refined, source=values)
        allocate (refined_vectors, source=vectors)
        refined(3) = refined(2)
        refined_vectors(:, 3) = refined_vectors(:, 2)
        call refine_dense_spectrum(stiffness, mass, controls, refined, &
            refined_vectors, info)
        call require(info == dense_spectrum_ok, &
            "indexed dense refinement failed")
        call require(maxval(abs(refined - values)) < 1.0e-11_dp, &
            "indexed dense refinement skipped an eigenvalue")
        call certify_dense_spectrum_orthogonality(mass, refined_vectors, info)
        call require(info == dense_spectrum_ok, &
            "refined dense vectors are not mass orthonormal")
        call diagnose_dense_spectrum(stiffness, mass, refined, &
            refined_vectors, rayleigh, residuals, resolutions, info)
        call require(info == dense_spectrum_ok, &
            "refined dense diagnostics failed")
        call certify_dense_spectrum_inertia(stiffness, mass, refined, &
            residuals, resolutions, info)
        call require(info == dense_spectrum_ok, &
            "refined dense inertia certificate failed")
    end subroutine check_indexed_dense_refinement

    subroutine build_fixture(stiffness, mass)
        real(dp), intent(out) :: stiffness(:, :), mass(:, :)
        integer :: i

        stiffness = 0.0_dp
        mass = 0.0_dp
        do i = 1, 6
            stiffness(i, i) = 3.0_dp + 0.4_dp * i
            mass(i, i) = 1.5_dp + 0.1_dp * i
        end do
        call set_symmetric(stiffness, 1, 2, 0.2_dp)
        call set_symmetric(stiffness, 2, 3, -0.3_dp)
        call set_symmetric(stiffness, 1, 4, 0.1_dp)
        call set_symmetric(stiffness, 3, 4, 0.15_dp)
        call set_symmetric(stiffness, 4, 5, -0.2_dp)
        call set_symmetric(stiffness, 5, 6, 0.25_dp)
        call set_symmetric(mass, 1, 2, 0.05_dp)
        call set_symmetric(mass, 2, 3, 0.04_dp)
        call set_symmetric(mass, 3, 4, 0.03_dp)
        call set_symmetric(mass, 5, 6, 0.02_dp)
    end subroutine build_fixture

    pure subroutine set_symmetric(matrix, i, j, value)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(in) :: i, j
        real(dp), intent(in) :: value

        matrix(i, j) = value
        matrix(j, i) = value
    end subroutine set_symmetric

    subroutine check_long_chain_resolution()
        integer, parameter :: n = 128
        type(variable_block_tridiagonal_t) :: chain_k, chain_m
        real(dp), allocatable :: dense_k(:, :), dense_m(:, :), vector(:)
        integer, allocatable :: chain_widths(:)
        real(dp) :: expected, factor, quotient, residual, resolution
        integer :: i, info

        allocate (dense_k(n, n), source=0.0_dp)
        allocate (dense_m(n, n), source=0.0_dp)
        allocate (vector(n), source=1.0_dp / sqrt(real(n, dp)))
        allocate (chain_widths(n), source=1)
        do i = 1, n
            dense_k(i, i) = 2.0_dp
            dense_m(i, i) = 1.0_dp
        end do
        call pack_variable_blocks(dense_k, chain_widths, chain_k, info)
        call require(info == 0, "long-chain stiffness packing failed")
        call pack_variable_blocks(dense_m, chain_widths, chain_m, info)
        call require(info == 0, "long-chain mass packing failed")
        call variable_generalized_diagnostics(chain_k, chain_m, vector, &
            2.0_dp, quotient, residual, resolution, info)
        call require(info == variable_generalized_ok, &
            "long-chain diagnostics failed")
        factor = 66.0_dp * epsilon(1.0_dp) + 16.0_dp * real(n, dp) &
            * epsilon(1.0_dp)**2
        factor = factor / (1.0_dp - factor)
        expected = 4.0_dp * factor
        call require(abs(resolution / expected - 1.0_dp) &
            < 32.0_dp * epsilon(1.0_dp), &
            "long-chain compensated resolution bound changed")
    end subroutine check_long_chain_resolution

    subroutine check_scaled_compensated_norm()
        real(dp) :: values(4)

        values = [3.0e200_dp, -4.0e200_dp, 0.0_dp, 0.0_dp]
        call require(abs(stable_norm2(values) / 5.0e200_dp - 1.0_dp) &
            < 4.0_dp * epsilon(1.0_dp), &
            "scaled compensated norm overflowed")
        values = [3.0e-200_dp, -4.0e-200_dp, 0.0_dp, 0.0_dp]
        call require(abs(stable_norm2(values) / 5.0e-200_dp - 1.0_dp) &
            < 4.0_dp * epsilon(1.0_dp), &
            "scaled compensated norm underflowed")
        call require(stable_norm2([real(dp) ::]) == 0.0_dp, &
            "empty compensated norm is nonzero")
    end subroutine check_scaled_compensated_norm

    subroutine check_ill_scaled_dense_pencil()
        integer, parameter :: n = 32
        integer, parameter :: packed_widths(1) = [n]
        type(variable_block_tridiagonal_t) :: packed_k, packed_m
        real(dp) :: dense_k(n, n), dense_m(n, n), mass_image(n)
        real(dp), allocatable :: stiffness_copy(:, :), mass_copy(:, :)
        real(dp), allocatable :: values(:), vectors(:, :)
        real(dp), allocatable :: seed_rayleigh(:)
        real(dp), allocatable :: seed_residuals(:), seed_resolutions(:)
        integer :: info

        call build_ill_scaled_pencil(dense_k, dense_m)
        allocate (stiffness_copy, source=dense_k)
        allocate (mass_copy, source=dense_m)
        call solve_symmetric_generalized_allocated(stiffness_copy, mass_copy, &
            values, vectors, info, equilibrate=.true.)
        call require(info == 0, "ill-scaled dense eigensolve failed")
        call require(abs(values(1) - (2.0_dp - 0.5_dp * cos( &
            acos(-1.0_dp) / real(n + 1, dp)))) < 1.0e-12_dp, &
            "ill-scaled dense eigenvalue changed")
        mass_image = matmul(dense_m, vectors(:, 1))
        call require(abs(dot_product(vectors(:, 1), mass_image) - 1.0_dp) &
            < 1.0e-8_dp, "ill-scaled eigenvector is not mass normalized")
        call pack_variable_blocks(dense_k, packed_widths, packed_k, info)
        call require(info == 0, "ill-scaled stiffness packing failed")
        call pack_variable_blocks(dense_m, packed_widths, packed_m, info)
        call require(info == 0, "ill-scaled mass packing failed")
        call diagnose_dense_spectrum(packed_k, packed_m, values, vectors, &
            seed_rayleigh, seed_residuals, seed_resolutions, info)
        call require(info == dense_spectrum_ok, &
            "ill-scaled seed diagnostics failed")
        call require(all(values(2:) - values(:n - 1) &
            > seed_residuals(2:) + seed_resolutions(2:) &
            + seed_residuals(:n - 1) + seed_resolutions(:n - 1)), &
            "ill-scaled seed gaps are unresolved")
        call certify_dense_spectrum_inertia(packed_k, packed_m, values, &
            seed_residuals, seed_resolutions, info)
        call require(info == dense_spectrum_ok, &
            "ill-scaled seed inertia certificate failed")
        call certify_dense_spectrum_orthogonality(packed_m, vectors, info)
        call require(info == dense_spectrum_ok, &
            "ill-scaled seed vectors are not mass orthonormal")
    end subroutine check_ill_scaled_dense_pencil

    subroutine build_ill_scaled_pencil(stiffness, mass)
        real(dp), intent(out) :: stiffness(:, :), mass(:, :)
        integer :: i, n
        real(dp) :: scale(size(mass, 1))

        n = size(mass, 1)
        stiffness = 0.0_dp
        mass = 0.0_dp
        do i = 1, n
            scale(i) = 10.0_dp**(-5.0_dp + 10.0_dp * real(i - 1, dp) &
                / real(n - 1, dp))
            mass(i, i) = scale(i)**2
            stiffness(i, i) = 2.0_dp * mass(i, i)
        end do
        do i = 1, n - 1
            call set_symmetric(stiffness, i, i + 1, &
                0.25_dp * scale(i) * scale(i + 1))
        end do
    end subroutine build_ill_scaled_pencil

    function mass_norm(vector, mass) result(squared_norm)
        real(dp), intent(in) :: vector(:)
        type(variable_block_tridiagonal_t), intent(in) :: mass
        real(dp) :: squared_norm, image(size(vector))
        integer :: info

        call apply_variable_block_tridiagonal(mass, vector, image, info)
        call require(info == 0, "variable mass apply failed")
        squared_norm = dot_product(vector, image)
    end function mass_norm

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_variable_generalized_solver
