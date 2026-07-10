module block_tridiagonal
    use, intrinsic :: iso_fortran_env, only: dp => real64
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

    interface
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
        integer :: k, nb, i, j

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
                coupled = transpose(matrix%off(:, :, i - 1))
                call dsytrs("U", k, k, factor%schur(:, :, i - 1), k, &
                    factor%pivots(:, i - 1), coupled, k, info)
                if (info /= 0) return
                factor%schur(:, :, i) = factor%schur(:, :, i) &
                    - matmul(matrix%off(:, :, i - 1), coupled)
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

    pure function pivot_negative_count(factored, pivots) result(count)
        real(dp), intent(in) :: factored(:, :)
        integer, intent(in) :: pivots(:)
        integer :: count
        real(dp) :: diagonal, determinant, trace
        integer :: j, k

        k = size(pivots)
        count = 0
        j = 1
        do while (j <= k)
            if (pivots(j) > 0) then
                diagonal = factored(j, j)
                if (diagonal < 0.0_dp) count = count + 1
                j = j + 1
            else
                determinant = factored(j, j) * factored(j + 1, j + 1) &
                    - factored(j, j + 1)**2
                trace = factored(j, j) + factored(j + 1, j + 1)
                if (determinant < 0.0_dp) then
                    count = count + 1
                else if (trace < 0.0_dp) then
                    count = count + 2
                end if
                j = j + 2
            end if
        end do
    end function pivot_negative_count

    subroutine solve_factored(matrix, factor, rhs, info)
        type(block_tridiagonal_t), intent(in) :: matrix
        type(block_factor_t), intent(in) :: factor
        real(dp), intent(inout) :: rhs(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: partial(:)
        integer :: k, nb, i

        k = size(matrix%diag, 1)
        nb = size(matrix%diag, 3)
        allocate (partial(k))
        do i = 1, nb
            if (i > 1) then
                partial = rhs(:, i - 1)
                rhs(:, i) = rhs(:, i) &
                    - matmul(matrix%off(:, :, i - 1), partial)
            end if
            call dsytrs("U", k, 1, factor%schur(:, :, i), k, &
                factor%pivots(:, i), rhs(:, i), k, info)
            if (info /= 0) return
        end do
        do i = nb - 1, 1, -1
            partial = matmul(transpose(matrix%off(:, :, i)), &
                rhs(:, i + 1))
            call dsytrs("U", k, 1, factor%schur(:, :, i), k, &
                factor%pivots(:, i), partial, k, info)
            if (info /= 0) return
            rhs(:, i) = rhs(:, i) - partial
        end do
        info = 0
    end subroutine solve_factored

    pure function apply_block_tridiagonal(matrix, vector) result(image)
        type(block_tridiagonal_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:, :)
        real(dp) :: image(size(vector, 1), size(vector, 2))
        integer :: nb, i

        nb = size(matrix%diag, 3)
        do i = 1, nb
            image(:, i) = matmul(matrix%diag(:, :, i), vector(:, i))
            if (i > 1) then
                image(:, i) = image(:, i) &
                    + matmul(matrix%off(:, :, i - 1), vector(:, i - 1))
            end if
            if (i < nb) then
                image(:, i) = image(:, i) &
                    + matmul(transpose(matrix%off(:, :, i)), &
                    vector(:, i + 1))
            end if
        end do
    end function apply_block_tridiagonal

end module block_tridiagonal
