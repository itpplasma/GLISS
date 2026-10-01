! Interval elements of a dynamic family layout scattered directly into
! block-tridiagonal storage, without the dense matrix of the whole problem.
!
! The blocks are those of build_dynamic_block_permutation. Entries are
! summed in element order into the same slots that the dense assembly
! would fill, and the storage is completed from the upper triangle as
! pack_permuted_variable_blocks does, so both routes agree bit for bit.
module dynamic_block_scatter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        dynamic_family_layout_t, dynamic_layout_ok
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: dynamic_block_scatter_ok = 0
    integer, parameter, public :: dynamic_block_scatter_invalid = -1
    integer, parameter, public :: dynamic_block_scatter_allocation = -2

    ! Block and position within it of every global unknown.
    type, public :: dynamic_block_map_t
        integer, allocatable :: widths(:), permutation(:)
        integer, allocatable :: block(:), local(:)
    end type dynamic_block_map_t

    public :: build_dynamic_block_map
    public :: allocate_dynamic_blocks
    public :: add_mapped_block_element
    public :: complete_dynamic_blocks

contains

    subroutine build_dynamic_block_map(layout, map, info)
        type(dynamic_family_layout_t), intent(in) :: layout
        type(dynamic_block_map_t), intent(out) :: map
        integer, intent(out) :: info
        integer :: allocation_status, block, local, position, status

        info = dynamic_block_scatter_invalid
        call build_dynamic_block_permutation(layout, map%widths, &
            map%permutation, status)
        if (status /= dynamic_layout_ok) return
        allocate (map%block(size(map%permutation)), &
            map%local(size(map%permutation)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = dynamic_block_scatter_allocation
            return
        end if
        position = 0
        do block = 1, size(map%widths)
            do local = 1, map%widths(block)
                position = position + 1
                map%block(map%permutation(position)) = block
                map%local(map%permutation(position)) = local
            end do
        end do
        info = dynamic_block_scatter_ok
    end subroutine build_dynamic_block_map

    subroutine allocate_dynamic_blocks(map, blocks, info)
        type(dynamic_block_map_t), intent(in) :: map
        type(variable_block_tridiagonal_t), intent(out) :: blocks
        integer, intent(out) :: info
        integer :: allocation_status, block, count

        info = dynamic_block_scatter_allocation
        count = size(map%widths)
        allocate (blocks%widths, source=map%widths, stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (blocks%diagonal(count), blocks%lower(max(0, count - 1)), &
            stat=allocation_status)
        if (allocation_status /= 0) return
        do block = 1, count
            allocate (blocks%diagonal(block)%values(map%widths(block), &
                map%widths(block)), source=0.0_dp, stat=allocation_status)
            if (allocation_status /= 0) return
            if (block == count) cycle
            allocate (blocks%lower(block)%values(map%widths(block + 1), &
                map%widths(block)), source=0.0_dp, stat=allocation_status)
            if (allocation_status /= 0) return
        end do
        info = dynamic_block_scatter_ok
    end subroutine allocate_dynamic_blocks

    ! Adds the upper triangle of each diagonal block and the upper
    ! off-diagonal blocks, stored transposed as lower(k); the symmetric lower
    ! entries are completed by complete_dynamic_blocks. An entry coupling
    ! blocks that are not neighbours is rejected.
    subroutine add_mapped_block_element(local_to_global, element, map, &
            blocks, info)
        integer, intent(in) :: local_to_global(:)
        real(dp), intent(in) :: element(:, :)
        type(dynamic_block_map_t), intent(in) :: map
        type(variable_block_tridiagonal_t), intent(inout) :: blocks
        integer, intent(out) :: info
        integer :: block_a, block_b, global_a, global_b, local_a, local_b
        integer :: row, column

        info = dynamic_block_scatter_invalid
        if (any(shape(element) /= size(local_to_global))) return
        if (any(local_to_global < 0)) return
        if (any(local_to_global > size(map%block))) return
        do local_b = 1, size(element, 2)
            global_b = local_to_global(local_b)
            if (global_b == 0) cycle
            block_b = map%block(global_b)
            column = map%local(global_b)
            do local_a = 1, size(element, 1)
                global_a = local_to_global(local_a)
                if (global_a == 0) cycle
                block_a = map%block(global_a)
                row = map%local(global_a)
                if (block_a == block_b) then
                    if (row > column) cycle
                    blocks%diagonal(block_a)%values(row, column) = &
                        blocks%diagonal(block_a)%values(row, column) &
                        + element(local_a, local_b)
                else if (block_b == block_a + 1) then
                    blocks%lower(block_a)%values(column, row) = &
                        blocks%lower(block_a)%values(column, row) &
                        + element(local_a, local_b)
                else if (block_a /= block_b + 1) then
                    return
                end if
            end do
        end do
        info = dynamic_block_scatter_ok
    end subroutine add_mapped_block_element

    subroutine complete_dynamic_blocks(blocks)
        type(variable_block_tridiagonal_t), intent(inout) :: blocks
        integer :: block, i, j

        do block = 1, size(blocks%diagonal)
            associate (values => blocks%diagonal(block)%values)
                do j = 1, size(values, 2)
                    do i = j + 1, size(values, 1)
                        values(i, j) = values(j, i)
                    end do
                end do
            end associate
        end do
    end subroutine complete_dynamic_blocks

end module dynamic_block_scatter
