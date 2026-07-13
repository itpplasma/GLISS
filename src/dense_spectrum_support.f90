module dense_spectrum_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    use variable_generalized_solver, only: variable_generalized_diagnostics, &
        variable_generalized_ok
    implicit none
    private

    integer, parameter, public :: dense_spectrum_ok = 0
    integer, parameter, public :: dense_spectrum_invalid = -1
    integer, parameter, public :: dense_spectrum_allocation = -2

    public :: diagnose_dense_spectrum
    public :: unpermute_dense_vectors

contains

    subroutine diagnose_dense_spectrum(stiffness, mass, eigenvalues, &
            eigenvectors, rayleigh_quotients, residuals, resolutions, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: eigenvalues(:), eigenvectors(:, :)
        real(dp), allocatable, intent(out) :: rayleigh_quotients(:)
        real(dp), allocatable, intent(out) :: residuals(:), resolutions(:)
        integer, intent(out) :: info
        integer :: allocation_status, index, count

        info = dense_spectrum_invalid
        count = size(eigenvalues)
        if (count < 1) return
        if (size(eigenvectors, 1) /= count &
            .or. size(eigenvectors, 2) /= count) return
        info = dense_spectrum_allocation
        allocate (rayleigh_quotients(count), residuals(count), &
            resolutions(count), stat=allocation_status)
        if (allocation_status /= 0) return
        do index = 1, count
            call variable_generalized_diagnostics(stiffness, mass, &
                eigenvectors(:, index), eigenvalues(index), &
                rayleigh_quotients(index), residuals(index), &
                resolutions(index), info)
            if (info /= variable_generalized_ok) then
                info = dense_spectrum_invalid
                return
            end if
        end do
        info = dense_spectrum_ok
    end subroutine diagnose_dense_spectrum

    subroutine unpermute_dense_vectors(vectors, permutation, info)
        real(dp), intent(inout) :: vectors(:, :)
        integer, intent(in) :: permutation(:)
        integer, intent(out) :: info
        logical, allocatable :: visited(:)
        real(dp), allocatable :: temporary(:)
        real(dp) :: value
        integer :: allocation_status, column, current, destination, start

        info = dense_spectrum_invalid
        if (size(vectors, 1) /= size(permutation)) return
        info = dense_spectrum_allocation
        allocate (visited(size(permutation)), source=.false., &
            stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (temporary(size(vectors, 2)), stat=allocation_status)
        if (allocation_status /= 0) return
        do start = 1, size(permutation)
            destination = permutation(start)
            if (destination < 1 .or. destination > size(permutation)) then
                info = dense_spectrum_invalid
                return
            end if
            if (visited(destination)) then
                info = dense_spectrum_invalid
                return
            end if
            visited(destination) = .true.
        end do
        visited = .false.
        do start = 1, size(permutation)
            if (visited(start)) cycle
            temporary = vectors(start, :)
            current = start
            do
                destination = permutation(current)
                do column = 1, size(vectors, 2)
                    value = vectors(destination, column)
                    vectors(destination, column) = temporary(column)
                    temporary(column) = value
                end do
                visited(current) = .true.
                current = destination
                if (current == start) exit
            end do
        end do
        info = dense_spectrum_ok
    end subroutine unpermute_dense_vectors

end module dense_spectrum_support
