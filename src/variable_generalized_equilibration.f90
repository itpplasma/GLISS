module variable_generalized_equilibration
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use variable_block_tridiagonal, only: validate_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: variable_equilibration_ok = 0
    integer, parameter, public :: variable_equilibration_invalid = -1
    integer, parameter :: equilibration_passes = 4

    public :: equilibrate_variable_generalized
    public :: undo_variable_congruence

contains

    subroutine equilibrate_variable_generalized(stiffness, mass, &
            balanced_stiffness, balanced_mass, scales, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        type(variable_block_tridiagonal_t), intent(out) :: balanced_stiffness
        type(variable_block_tridiagonal_t), intent(out) :: balanced_mass
        real(dp), allocatable, intent(out) :: scales(:)
        integer, intent(out) :: info
        type(variable_block_tridiagonal_t) :: next_mass, next_stiffness
        real(dp), allocatable :: row_scale(:), step(:)
        integer :: index, pass

        info = variable_equilibration_invalid
        call validate_variable_blocks(stiffness, info)
        if (info /= variable_block_ok) return
        call validate_variable_blocks(mass, info)
        if (info /= variable_block_ok) return
        if (size(stiffness%widths) /= size(mass%widths)) return
        if (any(stiffness%widths /= mass%widths)) return
        allocate (row_scale(sum(stiffness%widths)))
        allocate (scales(sum(stiffness%widths)), source=1.0_dp)
        allocate (step(sum(stiffness%widths)))
        balanced_stiffness = stiffness
        balanced_mass = mass
        do pass = 1, equilibration_passes
            row_scale = 0.0_dp
            call accumulate_row_scale(balanced_stiffness, row_scale)
            call accumulate_row_scale(balanced_mass, row_scale)
            if (.not. all(ieee_is_finite(row_scale))) return
            if (any(row_scale <= 0.0_dp)) return
            do index = 1, size(row_scale)
                step(index) = 1.0_dp / sqrt(row_scale(index))
            end do
            if (.not. all(ieee_is_finite(step))) return
            call apply_variable_congruence(balanced_stiffness, step, &
                next_stiffness)
            call apply_variable_congruence(balanced_mass, step, next_mass)
            balanced_stiffness = next_stiffness
            balanced_mass = next_mass
            do index = 1, size(scales)
                scales(index) = scales(index) * step(index)
            end do
            if (.not. all(ieee_is_finite(scales))) return
        end do
        call validate_variable_blocks(balanced_stiffness, info)
        if (info /= variable_block_ok) return
        call validate_variable_blocks(balanced_mass, info)
        if (info /= variable_block_ok) return
        info = variable_equilibration_ok
    end subroutine equilibrate_variable_generalized

    subroutine accumulate_row_scale(blocks, row_scale)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        real(dp), intent(inout) :: row_scale(:)
        integer :: block, first, next, row

        first = 1
        do block = 1, size(blocks%widths)
            do row = 1, blocks%widths(block)
                row_scale(first + row - 1) = max( &
                    row_scale(first + row - 1), &
                    maxval(abs(blocks%diagonal(block)%values(row, :))))
            end do
            if (block < size(blocks%widths)) then
                next = first + blocks%widths(block)
                do row = 1, blocks%widths(block + 1)
                    row_scale(next + row - 1) = max( &
                        row_scale(next + row - 1), &
                        maxval(abs(blocks%lower(block)%values(row, :))))
                end do
                do row = 1, blocks%widths(block)
                    row_scale(first + row - 1) = max( &
                        row_scale(first + row - 1), &
                        maxval(abs(blocks%lower(block)%values(:, row))))
                end do
            end if
            first = first + blocks%widths(block)
        end do
    end subroutine accumulate_row_scale

    subroutine apply_variable_congruence(blocks, scales, balanced)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        real(dp), intent(in) :: scales(:)
        type(variable_block_tridiagonal_t), intent(out) :: balanced
        integer :: block, column, first, next, row

        allocate (balanced%widths, source=blocks%widths)
        allocate (balanced%diagonal(size(blocks%diagonal)))
        allocate (balanced%lower(size(blocks%lower)))
        first = 1
        do block = 1, size(blocks%widths)
            allocate (balanced%diagonal(block)%values, &
                source=blocks%diagonal(block)%values)
            do column = 1, blocks%widths(block)
                do row = 1, blocks%widths(block)
                    balanced%diagonal(block)%values(row, column) = &
                        scales(first + row - 1) &
                        * blocks%diagonal(block)%values(row, column) &
                        * scales(first + column - 1)
                end do
            end do
            if (block < size(blocks%widths)) then
                next = first + blocks%widths(block)
                allocate (balanced%lower(block)%values, &
                    source=blocks%lower(block)%values)
                do column = 1, blocks%widths(block)
                    do row = 1, blocks%widths(block + 1)
                        balanced%lower(block)%values(row, column) = &
                            scales(next + row - 1) &
                            * blocks%lower(block)%values(row, column) &
                            * scales(first + column - 1)
                    end do
                end do
            end if
            first = first + blocks%widths(block)
        end do
    end subroutine apply_variable_congruence

    subroutine undo_variable_congruence(scales, balanced_vector, vector, info)
        real(dp), intent(in) :: scales(:), balanced_vector(:)
        real(dp), allocatable, intent(out) :: vector(:)
        integer, intent(out) :: info
        integer :: index

        info = variable_equilibration_invalid
        if (size(scales) /= size(balanced_vector)) return
        if (.not. all(ieee_is_finite(scales))) return
        if (.not. all(scales > 0.0_dp)) return
        if (.not. all(ieee_is_finite(balanced_vector))) return
        allocate (vector(size(scales)))
        do index = 1, size(scales)
            vector(index) = scales(index) * balanced_vector(index)
        end do
        if (.not. all(ieee_is_finite(vector))) return
        info = variable_equilibration_ok
    end subroutine undo_variable_congruence

end module variable_generalized_equilibration
