module symmetric_eigensolver
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: symmetric_eigensolver_ok = 0
    integer, parameter, public :: symmetric_eigensolver_invalid = -1
    integer, parameter, public :: symmetric_eigensolver_allocation = -2

    public :: solve_three_component_modes
    public :: solve_symmetric_generalized
    public :: solve_symmetric_generalized_allocated

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
        real(dp), allocatable :: stiffness_copy(:, :), mass_copy(:, :)
        integer :: allocation_status

        info = symmetric_eigensolver_invalid
        if (.not. valid_generalized_problem(stiffness, mass)) return
        info = symmetric_eigensolver_allocation
        allocate (stiffness_copy, source=stiffness, stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (mass_copy, source=mass, stat=allocation_status)
        if (allocation_status /= 0) return
        call solve_symmetric_generalized_allocated(stiffness_copy, mass_copy, &
            eigenvalues, eigenvectors, info)
    end subroutine solve_symmetric_generalized

    subroutine solve_symmetric_generalized_allocated(stiffness, mass, &
            eigenvalues, eigenvectors, info)
        real(dp), allocatable, intent(inout) :: stiffness(:, :), mass(:, :)
        real(dp), allocatable, intent(out) :: eigenvalues(:)
        real(dp), allocatable, intent(out) :: eigenvectors(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: work(:)
        integer :: allocation_status, n

        info = symmetric_eigensolver_invalid
        if (.not. allocated(stiffness) .or. .not. allocated(mass)) return
        if (.not. valid_generalized_problem(stiffness, mass)) return
        n = size(stiffness, 1)
        if (n > ishft(huge(n), -3)) return
        info = symmetric_eigensolver_allocation
        allocate (eigenvalues(n), work(max(1, 8 * n)), stat=allocation_status)
        if (allocation_status /= 0) return
        call move_alloc(stiffness, eigenvectors)
        call dsygv(1, "V", "U", n, eigenvectors, n, mass, n, eigenvalues, &
            work, size(work), info)
    end subroutine solve_symmetric_generalized_allocated

    function valid_generalized_problem(stiffness, mass) result(valid)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :)
        logical :: valid
        integer :: n

        valid = .false.
        n = size(stiffness, 1)
        if (n < 1 .or. size(stiffness, 2) /= n) return
        if (size(mass, 1) /= n .or. size(mass, 2) /= n) return
        if (.not. all(ieee_is_finite(stiffness))) return
        if (.not. all(ieee_is_finite(mass))) return
        if (.not. symmetric_matrix(stiffness)) return
        if (.not. symmetric_matrix(mass)) return
        valid = .true.
    end function valid_generalized_problem

    pure function symmetric_matrix(matrix) result(symmetric)
        real(dp), intent(in) :: matrix(:, :)
        logical :: symmetric
        real(dp) :: scale

        scale = max(1.0_dp, maxval(abs(matrix)))
        symmetric = maxval(abs(matrix - transpose(matrix))) &
            <= 128.0_dp * epsilon(1.0_dp) * scale
    end function symmetric_matrix

end module symmetric_eigensolver
