! Blockwise CAS3D2MN congruence of a sparse compatible pencil (#11).
!
! Oracle: the retained dense congruence of the same pencil. A random
! symmetric degree-one block pencil (three H1 and four L2 functions, three
! physical modes) is transformed by a colliding labeled table (five labels
! on three modes) blockwise and densely; mapped to the component order, the
! labeled stiffness agrees entrywise. The quotient certificate counts the
! physical rank exactly and the labeled null space of the coincident labels.
program test_cas3d_block_congruence
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cas3d_phase_envelope_transform, only: &
        apply_cas3d_phase_envelope_block_congruence, &
        apply_cas3d_phase_envelope_congruence, &
        build_cas3d_phase_envelope_map, cas3d_phase_envelope_map_t, &
        cas3d_phase_transform_ok, cas3d_quotient_t
    use compatible_block_storage, only: build_compatible_block_indices, &
        compatible_block_ok
    use fourier_phase_kind, only: phase_cosine
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    implicit none

    integer, parameter :: h1 = 3, l2 = 4, physical = 3, labels = 5
    real(dp), parameter :: scale = 0.25_dp
    type(cas3d_phase_envelope_map_t) :: map
    type(variable_block_tridiagonal_t) :: stiffness, mass
    type(cas3d_quotient_t) :: quotient
    real(dp), allocatable :: dense(:, :), terms(:, :, :), dense_mass(:, :)
    real(dp), allocatable :: blockwise(:, :)
    integer :: info

    call build_cas3d_phase_envelope_map([3, 4, 2, 2, 4], [2, 2, 2, 2, 2], &
        [1, 1, 1, 1, 1], [3, 4, 2], [2, 2, 2], phase_cosine, map, info)
    call require(info == cas3d_phase_transform_ok, 'colliding map failed')
    call random_pencil(stiffness, mass)
    call to_component_order(stiffness, physical, dense)
    allocate (terms(size(dense, 1), size(dense, 2), 1))
    terms(:, :, 1) = dense
    allocate (dense_mass(size(dense, 1), size(dense, 2)), source=0.0_dp)
    call apply_cas3d_phase_envelope_congruence(map, h1, l2, scale, dense, &
        terms, dense_mass, info)
    call require(info == cas3d_phase_transform_ok, 'dense congruence failed')

    call apply_cas3d_phase_envelope_block_congruence(map, h1, l2, scale, &
        stiffness, mass, quotient, info)
    call require(info == cas3d_phase_transform_ok, 'block congruence failed')
    call require(all(stiffness%widths == [2 * labels, 2 * labels, &
        2 * labels, labels]), 'labeled block widths differ')
    call to_component_order(stiffness, labels, blockwise)
    call require(all(shape(blockwise) == shape(dense)), &
        'labeled dimension differs from the dense congruence')
    call require(maxval(abs(blockwise - dense)) <= 1.0e-14_dp &
        * maxval(abs(dense)), 'blockwise congruence differs from the dense one')
    call to_component_order(mass, labels, blockwise)
    call require(all(blockwise == dense_mass), &
        'labeled coefficient mass differs')
    call require(quotient%physical_unknowns == (h1 + l2) * physical &
        .and. quotient%envelope_unknowns == (h1 + l2) * labels, &
        'quotient dimensions differ')
    call require(quotient%quotient_rank == (h1 + l2) * physical, &
        'quotient rank differs')
    call require(quotient%nullity == (h1 + l2) * (labels - physical), &
        'labeled nullity differs')
    call require(quotient%peak_block_width == 2 * labels, &
        'peak block width differs')
    write (*, '(a)') 'blockwise CAS3D2MN congruence matches the dense oracle'

contains

    subroutine random_pencil(stiffness, mass)
        type(variable_block_tridiagonal_t), intent(out) :: stiffness, mass
        integer, allocatable :: widths(:), block_index(:), local_index(:)
        integer :: block, status, column

        call build_compatible_block_indices(h1, l2, physical, physical, 1, &
            widths, block_index, local_index, status)
        call require(status == compatible_block_ok, 'physical layout failed')
        allocate (stiffness%widths, source=widths)
        allocate (mass%widths, source=widths)
        allocate (stiffness%diagonal(size(widths)), &
            stiffness%lower(size(widths) - 1), mass%diagonal(size(widths)), &
            mass%lower(size(widths) - 1))
        call random_seed(put=[(17 + column, column = 1, 64)])
        do block = 1, size(widths)
            allocate (stiffness%diagonal(block)%values(widths(block), &
                widths(block)), mass%diagonal(block)%values(widths(block), &
                widths(block)))
            call random_number(stiffness%diagonal(block)%values)
            call symmetrize(stiffness%diagonal(block)%values)
            mass%diagonal(block)%values = 0.0_dp
            do column = 1, widths(block)
                mass%diagonal(block)%values(column, column) = 1.0_dp
            end do
            if (block == size(widths)) cycle
            allocate (stiffness%lower(block)%values(widths(block + 1), &
                widths(block)), mass%lower(block)%values(widths(block + 1), &
                widths(block)), source=0.0_dp)
            call random_number(stiffness%lower(block)%values)
        end do
    end subroutine random_pencil

    pure subroutine symmetrize(matrix)
        real(dp), intent(inout) :: matrix(:, :)
        integer :: row, column

        do column = 1, size(matrix, 2)
            do row = column + 1, size(matrix, 1)
                matrix(row, column) = matrix(row, column) + matrix(column, row)
                matrix(column, row) = matrix(row, column)
            end do
            matrix(column, column) = 2.0_dp * matrix(column, column)
        end do
    end subroutine symmetrize

    subroutine to_component_order(blocks, modes, matrix)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        integer, intent(in) :: modes
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, allocatable :: widths(:), block_index(:), local_index(:)
        integer :: first, second, status, block_a, block_b

        ! Components first (normal, then eta), basis-major, mode-minor.
        call build_compatible_block_indices(h1, l2, modes, modes, 1, widths, &
            block_index, local_index, status)
        call require(status == compatible_block_ok, 'layout failed')
        allocate (matrix(size(block_index), size(block_index)), source=0.0_dp)
        do second = 1, size(block_index)
            block_b = block_index(second)
            do first = 1, size(block_index)
                block_a = block_index(first)
                if (block_a == block_b) then
                    matrix(first, second) = blocks%diagonal(block_a)%values( &
                        local_index(first), local_index(second))
                else if (block_a == block_b + 1) then
                    matrix(first, second) = blocks%lower(block_b)%values( &
                        local_index(first), local_index(second))
                else if (block_b == block_a + 1) then
                    matrix(first, second) = blocks%lower(block_a)%values( &
                        local_index(second), local_index(first))
                end if
            end do
        end do
    end subroutine to_component_order

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') 'FAIL: ' // message
            error stop 1
        end if
    end subroutine require

end program test_cas3d_block_congruence
