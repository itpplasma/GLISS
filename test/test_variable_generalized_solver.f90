program test_variable_generalized_solver
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use symmetric_eigensolver, only: solve_symmetric_generalized
    use variable_block_tridiagonal, only: &
        apply_variable_block_tridiagonal, pack_variable_blocks, &
        variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_generalized_diagnostics, &
        variable_generalized_inertia, variable_generalized_invalid, &
        variable_generalized_mass_not_spd, variable_generalized_ok
    implicit none

    integer, parameter :: widths(3) = [2, 3, 1]
    type(variable_block_tridiagonal_t) :: stiffness, mass, corrupt
    real(dp) :: dense_k(6, 6), dense_m(6, 6), eigenvalue, quotient, residual
    real(dp) :: resolution, diagnostic_residual, diagnostic_resolution
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :), vector(:)
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
    call check_long_chain_resolution()

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
        real(dp) :: quotient, residual, resolution
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
        call require(resolution >= 4.0_dp * real(n, dp) &
            * epsilon(1.0_dp) * abs(quotient), &
            "long-chain norm reduction is missing from resolution")
    end subroutine check_long_chain_resolution

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
