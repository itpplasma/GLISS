module symmetric_eigensolver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: solve_three_component_modes

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
        real(dp) :: mass_copy(3, 3), work(32)

        eigenvectors = stiffness
        mass_copy = mass
        call dsygv(1, "V", "U", 3, eigenvectors, 3, mass_copy, 3, &
            eigenvalues, work, size(work), info)
    end subroutine solve_three_component_modes

end module symmetric_eigensolver
