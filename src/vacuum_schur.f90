module vacuum_schur
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: vacuum_schur_ok = 0
    integer, parameter, public :: vacuum_schur_invalid_input = 1
    integer, parameter, public :: vacuum_schur_not_spd = 2

    public :: eliminate_vacuum

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf

        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
    end interface

contains

    subroutine eliminate_vacuum(plasma, vacuum, coupling, effective, &
            response, info)
        real(dp), intent(in) :: plasma(:, :), vacuum(:, :), coupling(:, :)
        real(dp), allocatable, intent(out) :: effective(:, :), response(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: factor(:, :)
        integer :: n_plasma, n_vacuum, lapack_info

        info = vacuum_schur_invalid_input
        n_plasma = size(plasma, 1)
        n_vacuum = size(vacuum, 1)
        if (n_plasma < 1 .or. n_vacuum < 1) return
        if (size(plasma, 2) /= n_plasma) return
        if (size(vacuum, 2) /= n_vacuum) return
        if (size(coupling, 1) /= n_plasma) return
        if (size(coupling, 2) /= n_vacuum) return
        if (.not. all(ieee_is_finite(plasma))) return
        if (.not. all(ieee_is_finite(vacuum))) return
        if (.not. all(ieee_is_finite(coupling))) return
        if (.not. symmetric_matrix(plasma)) return
        if (.not. symmetric_matrix(vacuum)) return

        allocate (factor(n_vacuum, n_vacuum))
        factor = vacuum
        call dpotrf("L", n_vacuum, factor, n_vacuum, lapack_info)
        if (lapack_info /= 0) then
            if (lapack_info > 0) info = vacuum_schur_not_spd
            return
        end if

        allocate (response(n_vacuum, n_plasma))
        response = transpose(coupling)
        call dpotrs("L", n_vacuum, n_plasma, factor, n_vacuum, response, &
            n_vacuum, lapack_info)
        if (lapack_info /= 0) return
        response = -response
        effective = plasma + matmul(coupling, response)
        info = vacuum_schur_ok
    end subroutine eliminate_vacuum

    pure function symmetric_matrix(matrix) result(symmetric)
        real(dp), intent(in) :: matrix(:, :)
        logical :: symmetric
        real(dp) :: scale

        scale = max(1.0_dp, maxval(abs(matrix)))
        symmetric = maxval(abs(matrix - transpose(matrix))) &
            <= 128.0_dp * epsilon(1.0_dp) * scale
    end function symmetric_matrix

end module vacuum_schur
