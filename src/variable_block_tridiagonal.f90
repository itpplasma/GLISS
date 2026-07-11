module variable_block_tridiagonal
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: variable_block_ok = 0
    integer, parameter, public :: variable_block_invalid = -1

    type, public :: variable_matrix_block_t
        real(dp), allocatable :: values(:, :)
    end type variable_matrix_block_t

    type, public :: variable_block_tridiagonal_t
        integer, allocatable :: widths(:)
        type(variable_matrix_block_t), allocatable :: diagonal(:)
        type(variable_matrix_block_t), allocatable :: lower(:)
    end type variable_block_tridiagonal_t

    public :: apply_variable_block_tridiagonal
    public :: pack_permuted_variable_blocks
    public :: pack_variable_blocks
    public :: validate_variable_blocks
    public :: variable_block_to_dense

contains

    subroutine pack_variable_blocks(dense, widths, blocks, info)
        real(dp), intent(in) :: dense(:, :)
        integer, intent(in) :: widths(:)
        type(variable_block_tridiagonal_t), intent(out) :: blocks
        integer, intent(out) :: info
        real(dp), allocatable :: reconstructed(:, :)
        real(dp) :: scale
        integer :: block, first, next, status

        info = variable_block_invalid
        if (.not. valid_dense_input(dense, widths)) return
        allocate (blocks%widths, source=widths)
        allocate (blocks%diagonal(size(widths)))
        allocate (blocks%lower(max(0, size(widths) - 1)))
        first = 1
        do block = 1, size(widths)
            allocate (blocks%diagonal(block)%values(widths(block), &
                widths(block)))
            blocks%diagonal(block)%values = dense(first:first + widths(block) &
                - 1, first:first + widths(block) - 1)
            if (block < size(widths)) then
                next = first + widths(block)
                allocate (blocks%lower(block)%values(widths(block + 1), &
                    widths(block)))
                blocks%lower(block)%values = dense(next:next &
                    + widths(block + 1) - 1, first:first + widths(block) - 1)
            end if
            first = first + widths(block)
        end do
        call variable_block_to_dense(blocks, reconstructed, status)
        if (status /= variable_block_ok) return
        scale = max(1.0_dp, maxval(abs(dense)))
        if (maxval(abs(reconstructed - dense)) &
            > 128.0_dp * epsilon(1.0_dp) * scale) return
        info = variable_block_ok
    end subroutine pack_variable_blocks

    subroutine pack_permuted_variable_blocks(dense, permutation, widths, &
            blocks, info)
        real(dp), intent(in) :: dense(:, :)
        integer, intent(in) :: permutation(:), widths(:)
        type(variable_block_tridiagonal_t), intent(out) :: blocks
        integer, intent(out) :: info
        logical, allocatable :: seen(:)
        real(dp), allocatable :: permuted(:, :)
        integer :: i, j, n

        info = variable_block_invalid
        n = size(dense, 1)
        if (n < 1 .or. size(dense, 2) /= n) return
        if (size(permutation) /= n) return
        allocate (seen(n), source=.false.)
        do i = 1, n
            if (permutation(i) < 1 .or. permutation(i) > n) return
            if (seen(permutation(i))) return
            seen(permutation(i)) = .true.
        end do
        allocate (permuted(n, n))
        do j = 1, n
            do i = 1, n
                permuted(i, j) = dense(permutation(i), permutation(j))
            end do
        end do
        call pack_variable_blocks(permuted, widths, blocks, info)
    end subroutine pack_permuted_variable_blocks

    subroutine variable_block_to_dense(blocks, dense, info)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        real(dp), allocatable, intent(out) :: dense(:, :)
        integer, intent(out) :: info
        integer :: block, first, n, next

        call validate_variable_blocks(blocks, info)
        if (info /= variable_block_ok) return
        n = sum(blocks%widths)
        allocate (dense(n, n), source=0.0_dp)
        first = 1
        do block = 1, size(blocks%widths)
            dense(first:first + blocks%widths(block) - 1, &
                first:first + blocks%widths(block) - 1) = &
                blocks%diagonal(block)%values
            if (block < size(blocks%widths)) then
                next = first + blocks%widths(block)
                dense(next:next + blocks%widths(block + 1) - 1, &
                    first:first + blocks%widths(block) - 1) = &
                    blocks%lower(block)%values
                dense(first:first + blocks%widths(block) - 1, &
                    next:next + blocks%widths(block + 1) - 1) = &
                    transpose(blocks%lower(block)%values)
            end if
            first = first + blocks%widths(block)
        end do
        info = variable_block_ok
    end subroutine variable_block_to_dense

    subroutine apply_variable_block_tridiagonal(blocks, vector, image, info)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        real(dp), intent(in) :: vector(:)
        real(dp), intent(out) :: image(:)
        integer, intent(out) :: info
        integer :: block, first, last, next_first, next_last

        call validate_variable_blocks(blocks, info)
        if (info /= variable_block_ok) return
        info = variable_block_invalid
        if (size(vector) /= sum(blocks%widths)) return
        if (size(image) /= size(vector)) return
        if (.not. all(ieee_is_finite(vector))) return
        first = 1
        do block = 1, size(blocks%widths)
            last = first + blocks%widths(block) - 1
            image(first:last) = matmul(blocks%diagonal(block)%values, &
                vector(first:last))
            if (block > 1) then
                next_first = first - blocks%widths(block - 1)
                next_last = first - 1
                image(first:last) = image(first:last) &
                    + matmul(blocks%lower(block - 1)%values, &
                    vector(next_first:next_last))
            end if
            if (block < size(blocks%widths)) then
                next_first = last + 1
                next_last = last + blocks%widths(block + 1)
                image(first:last) = image(first:last) &
                    + matmul(transpose(blocks%lower(block)%values), &
                    vector(next_first:next_last))
            end if
            first = last + 1
        end do
        info = variable_block_ok
    end subroutine apply_variable_block_tridiagonal

    subroutine validate_variable_blocks(blocks, info)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        integer, intent(out) :: info
        real(dp) :: scale
        integer :: block

        info = variable_block_invalid
        if (.not. allocated(blocks%widths)) return
        if (size(blocks%widths) < 1) return
        if (any(blocks%widths < 1)) return
        if (.not. allocated(blocks%diagonal)) return
        if (size(blocks%diagonal) /= size(blocks%widths)) return
        if (.not. allocated(blocks%lower)) return
        if (size(blocks%lower) /= size(blocks%widths) - 1) return
        do block = 1, size(blocks%widths)
            if (.not. allocated(blocks%diagonal(block)%values)) return
            if (any(shape(blocks%diagonal(block)%values) &
                /= blocks%widths(block))) return
            if (.not. all(ieee_is_finite( &
                blocks%diagonal(block)%values))) return
            scale = max(1.0_dp, maxval(abs( &
                blocks%diagonal(block)%values)))
            if (maxval(abs(blocks%diagonal(block)%values &
                - transpose(blocks%diagonal(block)%values))) &
                > 128.0_dp * epsilon(1.0_dp) * scale) return
            if (block == size(blocks%widths)) cycle
            if (.not. allocated(blocks%lower(block)%values)) return
            if (size(blocks%lower(block)%values, 1) &
                /= blocks%widths(block + 1)) return
            if (size(blocks%lower(block)%values, 2) &
                /= blocks%widths(block)) return
            if (.not. all(ieee_is_finite(blocks%lower(block)%values))) return
        end do
        info = variable_block_ok
    end subroutine validate_variable_blocks

    pure function valid_dense_input(dense, widths) result(valid)
        real(dp), intent(in) :: dense(:, :)
        integer, intent(in) :: widths(:)
        logical :: valid
        real(dp) :: scale

        valid = size(widths) > 0
        if (.not. valid) return
        valid = all(widths > 0)
        if (.not. valid) return
        valid = size(dense, 1) == size(dense, 2)
        if (.not. valid) return
        valid = size(dense, 1) == sum(widths)
        if (.not. valid) return
        valid = all(ieee_is_finite(dense))
        if (.not. valid) return
        scale = max(1.0_dp, maxval(abs(dense)))
        valid = maxval(abs(dense - transpose(dense))) &
            <= 128.0_dp * epsilon(1.0_dp) * scale
    end function valid_dense_input

end module variable_block_tridiagonal
