program test_variable_block_factorization
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use variable_block_tridiagonal, only: factorize_variable_shifted, &
        pack_variable_blocks, solve_variable_factored, variable_block_factor_t, &
        variable_block_invalid, variable_block_ok, variable_block_singular, &
        variable_block_tridiagonal_t
    implicit none

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    integer, parameter :: widths(3) = [2, 3, 1]
    type(variable_block_tridiagonal_t) :: blocks
    type(variable_block_factor_t) :: factor, corrupt_factor
    real(dp) :: dense(6, 6), eigenvalues(6), eigenvectors(6, 6)
    real(dp) :: residual(6), rhs(6), solution(6), original_rhs(6), shift
    real(dp) :: pivot_dense(2, 2), pivot_rhs(2), short_rhs(1), singular(1, 1)
    integer :: info, i

    call build_dense_fixture(dense)
    call dense_eigenvalues(dense, eigenvalues, eigenvectors)
    call pack_variable_blocks(dense, widths, blocks, info)
    call require(info == variable_block_ok, "factor fixture packing failed")
    call factorize_variable_shifted(blocks, eigenvalues(1) - 1.0_dp, &
        factor, info)
    call require(info == variable_block_ok .and. factor%negative_count == 0, &
        "variable inertia below the spectrum is wrong")
    do i = 1, size(eigenvalues) - 1
        shift = 0.5_dp * (eigenvalues(i) + eigenvalues(i + 1))
        call factorize_variable_shifted(blocks, shift, factor, info)
        call require(info == variable_block_ok, &
            "variable shifted factorization failed")
        call require(factor%negative_count == i, &
            "variable inertia between eigenvalues is wrong")
    end do
    call factorize_variable_shifted(blocks, eigenvalues(6) + 1.0_dp, &
        factor, info)
    call require(info == variable_block_ok .and. factor%negative_count == 6, &
        "variable inertia above the spectrum is wrong")

    rhs = [0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, -0.2_dp, 0.5_dp]
    original_rhs = rhs
    call factorize_variable_shifted(blocks, 0.0_dp, factor, info)
    call require(info == variable_block_ok, "variable factorization failed")
    blocks%lower(1)%values = 99.0_dp
    solution = rhs
    call solve_variable_factored(factor, solution, info)
    call require(info == variable_block_ok, "variable factored solve failed")
    call matrix_residual(dense, solution, original_rhs, residual)
    call require(norm2(residual) < 1.0e-12_dp, &
        "variable factored solve disagrees with dense residual")

    pivot_dense = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
    call pack_variable_blocks(pivot_dense, [2], blocks, info)
    call factorize_variable_shifted(blocks, 0.0_dp, factor, info)
    call require(info == variable_block_ok .and. factor%negative_count == 1, &
        "two-by-two pivot inertia is wrong")
    pivot_rhs = [1.0_dp, 2.0_dp]
    call solve_variable_factored(factor, pivot_rhs, info)
    call require(info == variable_block_ok, "two-by-two pivot solve failed")
    call require(maxval(abs(pivot_rhs - [2.0_dp, 1.0_dp])) < 1.0e-14_dp, &
        "two-by-two pivot solution is wrong")
    short_rhs = 1.0_dp
    call solve_variable_factored(factor, short_rhs, info)
    call require(info == variable_block_invalid, "wrong-sized RHS was accepted")
    pivot_rhs = ieee_value(0.0_dp, ieee_quiet_nan)
    call solve_variable_factored(factor, pivot_rhs, info)
    call require(info == variable_block_invalid, "nonfinite RHS was accepted")

    corrupt_factor = factor
    corrupt_factor%widths(1) = 3
    call solve_variable_factored(corrupt_factor, pivot_rhs, info)
    call require(info == variable_block_invalid, "malformed factor was accepted")
    corrupt_factor = factor
    corrupt_factor%pivots(1)%values = [-1, -2]
    call solve_variable_factored(corrupt_factor, pivot_rhs, info)
    call require(info == variable_block_invalid, &
        "mismatched negative pivots were accepted")
    corrupt_factor = factor
    corrupt_factor%pivots(1)%values = [-1, 2]
    call solve_variable_factored(corrupt_factor, pivot_rhs, info)
    call require(info == variable_block_invalid, &
        "unpaired negative pivot was accepted")
    corrupt_factor = factor
    corrupt_factor%pivots(1)%values = [2, 2]
    call solve_variable_factored(corrupt_factor, pivot_rhs, info)
    call require(info == variable_block_invalid, &
        "positionally invalid positive pivot was accepted")
    corrupt_factor = factor
    corrupt_factor%pivots(1)%values = 0
    call solve_variable_factored(corrupt_factor, pivot_rhs, info)
    call require(info == variable_block_invalid, "zero pivots were accepted")
    shift = ieee_value(0.0_dp, ieee_quiet_nan)
    call factorize_variable_shifted(blocks, shift, factor, info)
    call require(info == variable_block_invalid, "nonfinite shift was accepted")
    singular = 0.0_dp
    call pack_variable_blocks(singular, [1], blocks, info)
    call factorize_variable_shifted(blocks, 0.0_dp, factor, info)
    call require(info == variable_block_singular, &
        "singular variable block was not identified")
    short_rhs = 1.0_dp
    call solve_variable_factored(factor, short_rhs, info)
    call require(info == variable_block_invalid, &
        "partial singular factor was accepted for solve")

    write (*, "(a)") "PASS"

contains

    pure subroutine matrix_residual(matrix, vector, rhs, residual)
        real(dp), intent(in) :: matrix(:, :), vector(:), rhs(:)
        real(dp), intent(out) :: residual(:)
        integer :: column, row

        do row = 1, size(matrix, 1)
            residual(row) = -rhs(row)
            do column = 1, size(matrix, 2)
                residual(row) = residual(row) &
                    + matrix(row, column) * vector(column)
            end do
        end do
    end subroutine matrix_residual

    subroutine build_dense_fixture(matrix)
        real(dp), intent(out) :: matrix(:, :)
        integer :: column, row

        matrix = 0.0_dp
        matrix(1:2, 1:2) = reshape([3.0_dp, 0.2_dp, 0.2_dp, 2.5_dp], [2, 2])
        matrix(3:5, 3:5) = reshape([4.0_dp, 0.1_dp, -0.2_dp, 0.1_dp, &
            3.5_dp, 0.3_dp, -0.2_dp, 0.3_dp, 2.8_dp], [3, 3])
        matrix(6, 6) = 2.2_dp
        matrix(3:5, 1:2) = reshape([-0.4_dp, 0.1_dp, 0.0_dp, 0.2_dp, &
            -0.3_dp, 0.05_dp], [3, 2])
        do column = 3, 5
            do row = 1, 2
                matrix(row, column) = matrix(column, row)
            end do
        end do
        matrix(6, 3:5) = [0.1_dp, -0.2_dp, 0.3_dp]
        matrix(3:5, 6) = matrix(6, 3:5)
    end subroutine build_dense_fixture

    subroutine dense_eigenvalues(matrix, eigenvalues, eigenvectors)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), contiguous, intent(out) :: eigenvalues(:)
        real(dp), contiguous, intent(out) :: eigenvectors(:, :)
        real(dp) :: work(8 * size(matrix, 1))
        integer :: info

        eigenvectors = matrix
        call dsyev("V", "U", size(matrix, 1), eigenvectors, size(matrix, 1), &
            eigenvalues, work, size(work), info)
        call require(info == 0, "dense variable-factor oracle failed")
    end subroutine dense_eigenvalues

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_variable_block_factorization
