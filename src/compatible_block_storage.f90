module compatible_block_storage
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: compatible_block_ok = 0
    integer, parameter, public :: compatible_block_invalid = -1
    integer, parameter, public :: compatible_block_allocation = -2

    public :: initialize_compatible_block_pencil
    public :: scatter_symmetric_compatible_block
    public :: symmetrize_compatible_blocks

contains

    subroutine initialize_compatible_block_pencil(h1_dofs, l2_dofs, &
            normal_modes, eta_modes, degree, stiffness, mass, block_index, &
            local_index, info)
        integer, intent(in) :: h1_dofs, l2_dofs, normal_modes, eta_modes
        integer, intent(in) :: degree
        type(variable_block_tridiagonal_t), intent(out) :: stiffness, mass
        integer, allocatable, intent(out) :: block_index(:), local_index(:)
        integer, intent(out) :: info
        integer, allocatable :: widths(:)
        integer(int64) :: unknowns64
        integer :: allocation_status, basis, block, groups, local, mode
        integer :: normal_width, unknowns

        info = compatible_block_invalid
        if (h1_dofs < 0 .or. l2_dofs < 0) return
        if (normal_modes < 0 .or. eta_modes < 0) return
        if (degree < 1) return
        unknowns64 = int(h1_dofs, int64) * int(normal_modes, int64) &
            + int(l2_dofs, int64) * int(eta_modes, int64)
        if (unknowns64 < 1_int64) return
        if (unknowns64 > int(huge(unknowns), int64)) return
        unknowns = int(unknowns64)
        groups = max(ceiling_ratio(h1_dofs, degree), &
            ceiling_ratio(l2_dofs, degree))
        if (groups < 1) return
        allocate (widths(groups), source=0, stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_block_allocation
            return
        end if
        do block = 1, groups
            widths(block) = group_basis_count(h1_dofs, degree, block) &
                * normal_modes &
                + group_basis_count(l2_dofs, degree, block) * eta_modes
            if (widths(block) < 1) return
        end do
        allocate (block_index(unknowns), local_index(unknowns), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_block_allocation
            return
        end if
        do basis = 1, h1_dofs
            block = (basis - 1) / degree + 1
                do mode = 1, normal_modes
                    local = modulo(basis - 1, degree) * normal_modes + mode
                    block_index((basis - 1) * normal_modes + mode) = block
                    local_index((basis - 1) * normal_modes + mode) = local
                end do
            end do
            do basis = 1, l2_dofs
                block = (basis - 1) / degree + 1
                    normal_width = group_basis_count(h1_dofs, degree, block) &
                        * normal_modes
                    do mode = 1, eta_modes
                        local = normal_width &
                            + modulo(basis - 1, degree) * eta_modes + mode
                        block_index(h1_dofs * normal_modes &
                            + (basis - 1) * eta_modes + mode) = block
                        local_index(h1_dofs * normal_modes &
                            + (basis - 1) * eta_modes + mode) = local
                    end do
                end do
                call allocate_blocks(widths, stiffness, info)
                if (info /= compatible_block_ok) return
                call allocate_blocks(widths, mass, info)
            end subroutine initialize_compatible_block_pencil

            subroutine allocate_blocks(widths, blocks, info)
                integer, intent(in) :: widths(:)
                type(variable_block_tridiagonal_t), intent(out) :: blocks
                integer, intent(out) :: info
                integer :: allocation_status, block

                info = compatible_block_allocation
                allocate (blocks%widths, source=widths, stat=allocation_status)
                if (allocation_status /= 0) return
                allocate (blocks%diagonal(size(widths)), stat=allocation_status)
                if (allocation_status /= 0) return
                allocate (blocks%lower(max(0, size(widths) - 1)), &
                    stat=allocation_status)
                if (allocation_status /= 0) return
                do block = 1, size(widths)
                    allocate (blocks%diagonal(block)%values(widths(block), &
                        widths(block)), source=0.0_dp, stat=allocation_status)
                    if (allocation_status /= 0) return
                    if (block < size(widths)) then
                        allocate (blocks%lower(block)%values(widths(block + 1), &
                            widths(block)), source=0.0_dp, stat=allocation_status)
                        if (allocation_status /= 0) return
                    end if
                end do
                info = compatible_block_ok
            end subroutine allocate_blocks

            subroutine scatter_symmetric_compatible_block(map, source, scale, &
                    block_index, local_index, target, info)
                integer, intent(in) :: map(:), block_index(:), local_index(:)
                real(dp), intent(in) :: source(:, :), scale
                type(variable_block_tridiagonal_t), intent(inout) :: target
                integer, intent(out) :: info
                real(dp) :: symmetric_value
                integer :: a, b, block_a, block_b, global_a, global_b
                integer :: local_a, local_b

                info = compatible_block_invalid
                if (size(source, 1) /= size(map)) return
                if (size(source, 2) /= size(map)) return
                if (size(block_index) /= size(local_index)) return
                if (.not. ieee_is_finite(scale)) return
                do b = 1, size(map)
                    global_b = map(b)
                    if (global_b == 0) cycle
                    if (global_b < 1 .or. global_b > size(block_index)) return
                    block_b = block_index(global_b)
                    local_b = local_index(global_b)
                    do a = 1, size(map)
                        global_a = map(a)
                        if (global_a == 0) cycle
                        if (global_a < 1 .or. global_a > size(block_index)) return
                        block_a = block_index(global_a)
                        local_a = local_index(global_a)
                        if (abs(block_a - block_b) > 1) return
                        if (block_a == block_b) then
                            target%diagonal(block_a)%values(local_a, local_b) = &
                                target%diagonal(block_a)%values(local_a, local_b) &
                                + scale * source(a, b)
                        else if (block_a == block_b + 1) then
                            symmetric_value = 0.5_dp * (source(a, b) + source(b, a))
                            target%lower(block_b)%values(local_a, local_b) = &
                                target%lower(block_b)%values(local_a, local_b) &
                                + scale * symmetric_value
                        end if
                    end do
                end do
                info = compatible_block_ok
            end subroutine scatter_symmetric_compatible_block

            subroutine symmetrize_compatible_blocks(blocks)
                type(variable_block_tridiagonal_t), intent(inout) :: blocks
                real(dp) :: average
                integer :: block, column, row

                do block = 1, size(blocks%diagonal)
                    do column = 1, size(blocks%diagonal(block)%values, 2)
                        do row = column + 1, &
                                size(blocks%diagonal(block)%values, 1)
                            average = 0.5_dp &
                                * (blocks%diagonal(block)%values(row, column) &
                                + blocks%diagonal(block)%values(column, row))
                            blocks%diagonal(block)%values(row, column) = average
                            blocks%diagonal(block)%values(column, row) = average
                        end do
                    end do
                end do
            end subroutine symmetrize_compatible_blocks

            pure function ceiling_ratio(numerator, denominator) result(quotient)
                integer, intent(in) :: numerator, denominator
                integer :: quotient

                quotient = 0
                if (numerator > 0) quotient = 1 + (numerator - 1) / denominator
            end function ceiling_ratio

            pure function group_basis_count(dofs, degree, block) result(count)
                integer, intent(in) :: dofs, degree, block
                integer :: count, first

                first = (block - 1) * degree + 1
                count = max(0, min(degree, dofs - first + 1))
            end function group_basis_count

        end module compatible_block_storage
