module block_tridiagonal
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use symmetric_pivot_inertia, only: pivot_negative_count
    implicit none
    private

    type, public :: block_tridiagonal_t
        real(dp), allocatable :: diag(:, :, :)
        real(dp), allocatable :: off(:, :, :)
    end type block_tridiagonal_t

    type, public :: block_factor_t
        real(dp), allocatable :: schur(:, :, :)
        integer, allocatable :: pivots(:, :)
        integer :: negative_count = 0
    end type block_factor_t

    public :: factorize_shifted
    public :: solve_factored
    public :: apply_block_tridiagonal
    public :: apply_block_tridiagonal_into

    interface
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta
            real(dp), intent(in) :: a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: dp
            character(len=1), intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(dp), intent(in) :: alpha, beta
            real(dp), intent(in) :: a(lda, *), x(*)
            real(dp), intent(inout) :: y(*)
        end subroutine dgemv
        subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine dsytrf
        subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dsytrs
    end interface

contains

    subroutine factorize_shifted(matrix, shift, factor, info)
        type(block_tridiagonal_t), intent(in) :: matrix
        real(dp), intent(in) :: shift
        type(block_factor_t), intent(out) :: factor
        integer, intent(out) :: info
        real(dp), allocatable :: coupled(:, :), work(:)
        integer :: column, k, nb, i, j, row

        k = size(matrix%diag, 1)
        nb = size(matrix%diag, 3)
        allocate (factor%schur(k, k, nb), factor%pivots(k, nb))
        allocate (coupled(k, k), work(64 * k))
        factor%negative_count = 0
        do i = 1, nb
            factor%schur(:, :, i) = matrix%diag(:, :, i)
            do j = 1, k
                factor%schur(j, j, i) = factor%schur(j, j, i) - shift
            end do
            if (i > 1) then
                do column = 1, k
                    do row = 1, k
                        coupled(row, column) = matrix%off(column, row, i - 1)
                    end do
                end do
                call dsytrs("U", k, k, factor%schur(:, :, i - 1), k, &
                    factor%pivots(:, i - 1), coupled, k, info)
                if (info /= 0) return
                call dgemm("N", "N", k, k, k, -1.0_dp, &
                    matrix%off(:, :, i - 1), k, coupled, k, 1.0_dp, &
                    factor%schur(:, :, i), k)
            end if
            call dsytrf("U", k, factor%schur(:, :, i), k, &
                factor%pivots(:, i), work, size(work), info)
            if (info /= 0) return
            factor%negative_count = factor%negative_count &
                + pivot_negative_count(factor%schur(:, :, i), &
                factor%pivots(:, i))
        end do
        info = 0
    end subroutine factorize_shifted

    subroutine solve_factored(matrix, factor, rhs, info)
        type(block_tridiagonal_t), intent(in) :: matrix
        type(block_factor_t), intent(in) :: factor
        real(dp), contiguous, intent(inout) :: rhs(:, :)
        integer, intent(out) :: info
        real(dp) :: partial(size(rhs, 1))
        integer :: k, nb, i

        k = size(matrix%diag, 1)
        nb = size(matrix%diag, 3)
        do i = 1, nb
            if (i > 1) then
                call dgemv("N", k, k, -1.0_dp, matrix%off(:, :, i - 1), &
                    k, rhs(:, i - 1), 1, 1.0_dp, rhs(:, i), 1)
            end if
            call dsytrs("U", k, 1, factor%schur(:, :, i), k, &
                factor%pivots(:, i), rhs(:, i), k, info)
            if (info /= 0) return
        end do
        do i = nb - 1, 1, -1
            call dgemv("T", k, k, 1.0_dp, matrix%off(:, :, i), k, &
                rhs(:, i + 1), 1, 0.0_dp, partial, 1)
            call dsytrs("U", k, 1, factor%schur(:, :, i), k, &
                factor%pivots(:, i), partial, k, info)
            if (info /= 0) return
            rhs(:, i) = rhs(:, i) - partial
        end do
        info = 0
    end subroutine solve_factored

    function apply_block_tridiagonal(matrix, vector) result(image)
        type(block_tridiagonal_t), intent(in) :: matrix
        real(dp), contiguous, intent(in) :: vector(:, :)
        real(dp) :: image(size(vector, 1), size(vector, 2))

        call apply_block_tridiagonal_into(matrix, vector, image)
    end function apply_block_tridiagonal

    subroutine apply_block_tridiagonal_into(matrix, vector, image)
        type(block_tridiagonal_t), intent(in) :: matrix
        real(dp), contiguous, intent(in) :: vector(:, :)
        real(dp), contiguous, intent(out) :: image(:, :)
        integer :: k, nb, i

        k = size(matrix%diag, 1)
        nb = size(matrix%diag, 3)
        do i = 1, nb
            call dgemv("N", k, k, 1.0_dp, matrix%diag(:, :, i), k, &
                vector(:, i), 1, 0.0_dp, image(:, i), 1)
            if (i > 1) then
                call dgemv("N", k, k, 1.0_dp, matrix%off(:, :, i - 1), &
                    k, vector(:, i - 1), 1, 1.0_dp, image(:, i), 1)
            end if
            if (i < nb) then
                call dgemv("T", k, k, 1.0_dp, matrix%off(:, :, i), k, &
                    vector(:, i + 1), 1, 1.0_dp, image(:, i), 1)
            end if
        end do
    end subroutine apply_block_tridiagonal_into

end module block_tridiagonal
