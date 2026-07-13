module variable_block_tridiagonal
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use symmetric_pivot_inertia, only: pivot_negative_count
    implicit none
    private

    integer, parameter, public :: variable_block_ok = 0
    integer, parameter, public :: variable_block_invalid = -1
    integer, parameter, public :: variable_block_singular = -2
    integer, parameter, public :: variable_block_allocation = -3

    type, public :: variable_matrix_block_t
        real(dp), allocatable :: values(:, :)
    end type variable_matrix_block_t

    type, public :: variable_block_tridiagonal_t
        integer, allocatable :: widths(:)
        type(variable_matrix_block_t), allocatable :: diagonal(:)
        type(variable_matrix_block_t), allocatable :: lower(:)
    end type variable_block_tridiagonal_t

    type, public :: variable_integer_block_t
        integer, allocatable :: values(:)
    end type variable_integer_block_t

    type, public :: variable_block_factor_t
        integer, allocatable :: widths(:)
        type(variable_matrix_block_t), allocatable :: schur(:)
        type(variable_matrix_block_t), allocatable :: lower(:)
        type(variable_integer_block_t), allocatable :: pivots(:)
        integer :: negative_count = 0
        logical :: complete = .false.
    end type variable_block_factor_t

    public :: apply_variable_block_tridiagonal
    public :: factorize_variable_shifted
    public :: pack_permuted_variable_blocks
    public :: pack_variable_blocks
    public :: solve_variable_factored
    public :: validate_variable_blocks
    public :: variable_block_to_dense

    interface
        subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine dsytrf
        subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dsytrs
    end interface

contains

    subroutine factorize_variable_shifted(blocks, shift, factor, info)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        real(dp), intent(in) :: shift
        type(variable_block_factor_t), intent(out) :: factor
        integer, intent(out) :: info
        real(dp), allocatable :: coupled(:, :), work(:)
        integer :: block, current_width, j, maximum_width, previous_width

        call validate_variable_blocks(blocks, info)
        if (info /= variable_block_ok) return
        info = variable_block_invalid
        if (.not. ieee_is_finite(shift)) return
        maximum_width = maxval(blocks%widths)
        call initialize_variable_factor(blocks, factor)
        allocate (coupled(maximum_width, maximum_width))
        allocate (work(64 * maximum_width))
        factor%negative_count = 0
        do block = 1, size(blocks%widths)
            current_width = blocks%widths(block)
            factor%schur(block)%values = blocks%diagonal(block)%values
            do j = 1, current_width
                factor%schur(block)%values(j, j) = &
                    factor%schur(block)%values(j, j) - shift
            end do
            if (block > 1) then
                previous_width = blocks%widths(block - 1)
                call update_variable_schur(block, previous_width, &
                    current_width, factor, blocks%lower(block - 1)%values, &
                    coupled, info)
                if (info /= variable_block_ok) return
            end if
            call dsytrf("U", current_width, factor%schur(block)%values, &
                current_width, factor%pivots(block)%values, work, size(work), &
                info)
            if (info > 0) then
                info = variable_block_singular
                return
            end if
            if (info < 0) then
                info = variable_block_invalid
                return
            end if
            factor%negative_count = factor%negative_count &
                + pivot_negative_count(factor%schur(block)%values, &
                factor%pivots(block)%values)
        end do
        factor%complete = .true.
        info = variable_block_ok
    end subroutine factorize_variable_shifted

    subroutine initialize_variable_factor(blocks, factor)
        type(variable_block_tridiagonal_t), intent(in) :: blocks
        type(variable_block_factor_t), intent(out) :: factor
        integer :: block, width

        allocate (factor%widths, source=blocks%widths)
        allocate (factor%schur(size(blocks%widths)))
        allocate (factor%lower(max(0, size(blocks%widths) - 1)))
        allocate (factor%pivots(size(blocks%widths)))
        do block = 1, size(blocks%widths)
            width = blocks%widths(block)
            allocate (factor%schur(block)%values(width, width))
            allocate (factor%pivots(block)%values(width))
            factor%pivots(block)%values = 0
            if (block < size(blocks%widths)) then
                allocate (factor%lower(block)%values, &
                    source=blocks%lower(block)%values)
            end if
        end do
    end subroutine initialize_variable_factor

    subroutine update_variable_schur(block, previous_width, current_width, &
            factor, lower, coupled, info)
        integer, intent(in) :: block, previous_width, current_width
        type(variable_block_factor_t), intent(inout) :: factor
        real(dp), intent(in) :: lower(:, :)
        real(dp), intent(inout) :: coupled(:, :)
        integer, intent(out) :: info

        coupled(1:previous_width, 1:current_width) = transpose(lower)
        call dsytrs("U", previous_width, current_width, &
            factor%schur(block - 1)%values, previous_width, &
            factor%pivots(block - 1)%values, coupled, size(coupled, 1), info)
        if (info /= 0) then
            info = variable_block_invalid
            return
        end if
        factor%schur(block)%values = factor%schur(block)%values &
            - matmul(lower, coupled(1:previous_width, 1:current_width))
        info = variable_block_ok
    end subroutine update_variable_schur

    subroutine solve_variable_factored(factor, rhs, info)
        type(variable_block_factor_t), intent(in) :: factor
        real(dp), intent(inout) :: rhs(:)
        integer, intent(out) :: info
        integer, allocatable :: offsets(:)
        real(dp), allocatable :: partial(:)
        integer :: block, first, last, next_first, next_last, width

        call validate_variable_factor(factor, info)
        if (info /= variable_block_ok) return
        info = variable_block_invalid
        if (size(rhs) /= sum(factor%widths)) return
        if (.not. all(ieee_is_finite(rhs))) return
        call build_block_offsets(factor%widths, offsets)
        allocate (partial(maxval(factor%widths)))
        do block = 1, size(factor%widths)
            first = offsets(block)
            last = offsets(block + 1) - 1
            width = factor%widths(block)
            if (block > 1) then
                rhs(first:last) = rhs(first:last) &
                    - matmul(factor%lower(block - 1)%values, &
                    rhs(offsets(block - 1):first - 1))
            end if
            call dsytrs("U", width, 1, factor%schur(block)%values, width, &
                factor%pivots(block)%values, rhs(first:last), width, info)
            if (info /= 0) then
                info = variable_block_invalid
                return
            end if
        end do
        do block = size(factor%widths) - 1, 1, -1
            first = offsets(block)
            last = offsets(block + 1) - 1
            next_first = offsets(block + 1)
            next_last = offsets(block + 2) - 1
            width = factor%widths(block)
            partial(1:width) = matmul(transpose( &
                factor%lower(block)%values), rhs(next_first:next_last))
            call dsytrs("U", width, 1, factor%schur(block)%values, width, &
                factor%pivots(block)%values, partial, width, info)
            if (info /= 0) then
                info = variable_block_invalid
                return
            end if
            rhs(first:last) = rhs(first:last) - partial(1:width)
        end do
        info = variable_block_ok
    end subroutine solve_variable_factored

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
            call copy_upper_to_lower(blocks%diagonal(block)%values)
            if (block < size(widths)) then
                next = first + widths(block)
                allocate (blocks%lower(block)%values(widths(block + 1), &
                    widths(block)))
                blocks%lower(block)%values = transpose(dense(first:first &
                    + widths(block) - 1, next:next + widths(block + 1) - 1))
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

    pure subroutine copy_upper_to_lower(matrix)
        real(dp), intent(inout) :: matrix(:, :)
        integer :: i, j

        do j = 1, size(matrix, 2)
            do i = j + 1, size(matrix, 1)
                matrix(i, j) = matrix(j, i)
            end do
        end do
    end subroutine copy_upper_to_lower

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
        integer :: allocation_status, block, first, n, next

        call validate_variable_blocks(blocks, info)
        if (info /= variable_block_ok) return
        n = sum(blocks%widths)
        allocate (dense(n, n), source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) then
            info = variable_block_allocation
            return
        end if
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

    subroutine validate_variable_factor(factor, info)
        type(variable_block_factor_t), intent(in) :: factor
        integer, intent(out) :: info
        integer :: block

        info = variable_block_invalid
        if (.not. factor%complete) return
        if (.not. allocated(factor%widths)) return
        if (size(factor%widths) < 1 .or. any(factor%widths < 1)) return
        if (.not. allocated(factor%schur)) return
        if (size(factor%schur) /= size(factor%widths)) return
        if (.not. allocated(factor%pivots)) return
        if (size(factor%pivots) /= size(factor%widths)) return
        if (.not. allocated(factor%lower)) return
        if (size(factor%lower) /= size(factor%widths) - 1) return
        do block = 1, size(factor%widths)
            if (.not. allocated(factor%schur(block)%values)) return
            if (any(shape(factor%schur(block)%values) &
                /= factor%widths(block))) return
            if (.not. all(ieee_is_finite(factor%schur(block)%values))) return
            if (.not. allocated(factor%pivots(block)%values)) return
            if (size(factor%pivots(block)%values) /= factor%widths(block)) &
                return
            if (.not. valid_upper_pivots(factor%pivots(block)%values)) return
            if (block == size(factor%widths)) cycle
            if (.not. allocated(factor%lower(block)%values)) return
            if (size(factor%lower(block)%values, 1) &
                /= factor%widths(block + 1)) return
            if (size(factor%lower(block)%values, 2) &
                /= factor%widths(block)) return
            if (.not. all(ieee_is_finite(factor%lower(block)%values))) return
        end do
        info = variable_block_ok
    end subroutine validate_variable_factor

    pure function valid_upper_pivots(pivots) result(valid)
        integer, intent(in) :: pivots(:)
        logical :: valid
        integer :: position

        valid = .false.
        position = size(pivots)
        do while (position >= 1)
            if (pivots(position) == 0) return
            if (pivots(position) > 0) then
                if (pivots(position) > position) return
                position = position - 1
            else
                if (position == 1) return
                if (pivots(position) /= pivots(position - 1)) return
                if (pivots(position) < -(position - 1)) return
                position = position - 2
            end if
        end do
        valid = .true.
    end function valid_upper_pivots

    pure subroutine build_block_offsets(widths, offsets)
        integer, intent(in) :: widths(:)
        integer, allocatable, intent(out) :: offsets(:)
        integer :: block

        allocate (offsets(size(widths) + 1))
        offsets(1) = 1
        do block = 1, size(widths)
            offsets(block + 1) = offsets(block) + widths(block)
        end do
    end subroutine build_block_offsets

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
