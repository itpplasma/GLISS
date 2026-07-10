module symmetric_eigensolver
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: solve_three_component_modes
    public :: solve_symmetric_generalized

    interface
        subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, &
                lwork, info)
            import :: dp
            integer, intent(in) :: itype, n, lda, ldb, lwork
            character(len=1), intent(in) :: jobz, uplo
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: w(*)
            real(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine dsygv
    end interface

contains

    subroutine solve_three_component_modes(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        real(dp), intent(in) :: stiffness(3, 3), mass(3, 3)
        real(dp), intent(out) :: eigenvalues(3), eigenvectors(3, 3)
        integer, intent(out) :: info
        real(dp), allocatable :: values(:), vectors(:, :)

        call solve_symmetric_generalized(stiffness, mass, values, vectors, &
            info)
        if (info /= 0) return
        eigenvalues = values
        eigenvectors = vectors
    end subroutine solve_three_component_modes

    subroutine solve_symmetric_generalized(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :)
        real(dp), allocatable, intent(out) :: eigenvalues(:)
        real(dp), allocatable, intent(out) :: eigenvectors(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: mass_copy(:, :), work(:)
        integer :: n

        info = -1
        n = size(stiffness, 1)
        if (n < 1 .or. size(stiffness, 2) /= n) return
        if (size(mass, 1) /= n .or. size(mass, 2) /= n) return
        if (.not. all(ieee_is_finite(stiffness))) return
        if (.not. all(ieee_is_finite(mass))) return
        if (.not. symmetric_matrix(stiffness)) return
        if (.not. symmetric_matrix(mass)) return
        allocate (eigenvalues(n), eigenvectors(n, n), mass_copy(n, n))
        allocate (work(max(1, 8 * n)))
        eigenvectors = stiffness
        mass_copy = mass
        call dsygv(1, "V", "U", n, eigenvectors, n, mass_copy, n, &
            eigenvalues, work, size(work), info)
    end subroutine solve_symmetric_generalized

    pure function symmetric_matrix(matrix) result(symmetric)
        real(dp), intent(in) :: matrix(:, :)
        logical :: symmetric
        real(dp) :: scale

        scale = max(1.0_dp, maxval(abs(matrix)))
        symmetric = maxval(abs(matrix - transpose(matrix))) &
            <= 128.0_dp * epsilon(1.0_dp) * scale
    end function symmetric_matrix

end module symmetric_eigensolver
