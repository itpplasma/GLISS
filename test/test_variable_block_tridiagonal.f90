program test_variable_block_tridiagonal
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        build_dynamic_family_layout, dynamic_family_layout_t, dynamic_layout_ok
    use variable_block_tridiagonal, only: &
        apply_variable_block_tridiagonal, pack_permuted_variable_blocks, &
        pack_variable_blocks, variable_block_invalid, variable_block_ok, &
        variable_block_to_dense, variable_block_tridiagonal_t
    !$ use omp_lib, only: omp_set_num_threads
    implicit none

    integer, parameter :: widths(3) = [2, 3, 1]
    type(variable_block_tridiagonal_t) :: blocks
    type(dynamic_family_layout_t) :: layout
    real(dp) :: dense(6, 6), corrupt(6, 6), vector(6), image(6), reference(6)
    real(dp), allocatable :: reconstructed(:, :), canonical(:, :)
    integer, allocatable :: dynamic_widths(:), permutation(:)
    integer :: expected_permutation(22), info

    vector = 0.0_dp
    call apply_variable_block_tridiagonal(blocks, vector, image, info)
    call require(info == variable_block_invalid, "uninitialized blocks accepted")
    call build_dense_fixture(dense)
    call pack_variable_blocks(dense, widths, blocks, info)
    call require(info == variable_block_ok, "valid variable blocks rejected")
    call variable_block_to_dense(blocks, reconstructed, info)
    call require(info == variable_block_ok, "variable block reconstruction failed")
    call require(maxval(abs(reconstructed - dense)) < 1.0e-14_dp, &
        "variable block roundtrip changed the matrix")
    vector = [0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, -0.2_dp, 0.5_dp]
    call apply_variable_block_tridiagonal(blocks, vector, image, info)
    call require(info == variable_block_ok, "variable block apply failed")
    call matrix_vector_product(dense, vector, reference)
    call require(maxval(abs(image - reference)) < 1.0e-14_dp, &
        "variable block apply disagrees with dense matmul")

    call build_dynamic_family_layout(2, 4, layout, info)
    call require(info == dynamic_layout_ok, "dynamic layout setup failed")
    call build_dynamic_block_permutation(layout, dynamic_widths, permutation, &
        info)
    call require(info == dynamic_layout_ok, "dynamic block permutation failed")
    call require(all(dynamic_widths == [6, 6, 6, 4]), &
        "dynamic block widths are wrong")
    expected_permutation = [1, 2, 7, 8, 15, 16, 3, 4, 9, 10, 17, 18, &
        5, 6, 11, 12, 19, 20, 13, 14, 21, 22]
    call require(all(permutation == expected_permutation), &
        "dynamic block permutation is wrong")

    allocate (canonical(6, 6), source=0.0_dp)
    call place_in_canonical_order(dense, [3, 1, 6, 2, 5, 4], canonical)
    call pack_permuted_variable_blocks(canonical, [3, 1, 6, 2, 5, 4], &
        widths, blocks, info)
    call require(info == variable_block_ok, "valid permuted packing failed")
    call variable_block_to_dense(blocks, reconstructed, info)
    call require(maxval(abs(reconstructed - dense)) < 1.0e-14_dp, &
        "permuted variable block packing is wrong")

    corrupt = dense
    corrupt(1, 6) = 0.1_dp
    corrupt(6, 1) = 0.1_dp
    call pack_variable_blocks(corrupt, widths, blocks, info)
    call require(info == variable_block_invalid, "off-band fill was accepted")
    corrupt = dense
    corrupt(1, 2) = corrupt(1, 2) + 0.1_dp
    call pack_variable_blocks(corrupt, widths, blocks, info)
    call require(info == variable_block_invalid, &
        "nonsymmetric variable matrix was accepted")
    call pack_variable_blocks(dense, [2, 2, 1], blocks, info)
    call require(info == variable_block_invalid, &
        "widths with the wrong total were accepted")
    call pack_permuted_variable_blocks(dense, [1, 2, 2, 4, 5, 6], widths, &
        blocks, info)
    call require(info == variable_block_invalid, &
        "duplicate permutation index was accepted")

    call check_threaded_apply()
    write (*, "(a)") "PASS"

contains

    subroutine check_threaded_apply()
        integer, parameter :: widths(7) = [1, 50, 2, 250, 5, 200, 3]
        integer, parameter :: n = sum(widths)
        type(variable_block_tridiagonal_t) :: chain
        real(dp) :: dense(n, n), vector(n), image(n), expected(n)
        integer :: i, threads, status

        dense = 0.0_dp
        do i = 1, n
            dense(i, i) = 2.0_dp
            vector(i) = real(i, dp) / 128.0_dp
        end do
        do i = 1, n - 1
            dense(i + 1, i) = -0.25_dp
            dense(i, i + 1) = -0.25_dp
        end do
        expected = 2.0_dp * vector
        expected(2:) = expected(2:) - 0.25_dp * vector(:n - 1)
        expected(:n - 1) = expected(:n - 1) - 0.25_dp * vector(2:)
        call pack_variable_blocks(dense, widths, chain, status)
        call require(status == variable_block_ok, "threaded chain packing")
        ! Unequal chunks, with more threads than blocks at the last count.
        do threads = 1, 8
            !$ call omp_set_num_threads(threads)
            call apply_variable_block_tridiagonal(chain, vector, image, status)
            call require(status == variable_block_ok, "threaded chain apply")
            call require(maxval(abs(image - expected)) &
                <= 4.0_dp * epsilon(1.0_dp) * maxval(abs(expected)), &
                "threaded block apply differs from exact dyadic stencil")
        end do
        !$ call omp_set_num_threads(1)
    end subroutine check_threaded_apply


    pure subroutine matrix_vector_product(matrix, vector, image)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp), intent(out) :: image(:)
        integer :: column, row

        do row = 1, size(matrix, 1)
            image(row) = 0.0_dp
            do column = 1, size(matrix, 2)
                image(row) = image(row) + matrix(row, column) * vector(column)
            end do
        end do
    end subroutine matrix_vector_product

    subroutine build_dense_fixture(matrix)
        real(dp), intent(out) :: matrix(:, :)
        integer :: column, row

        matrix = 0.0_dp
        matrix(1:2, 1:2) = reshape([3.0_dp, 0.2_dp, 0.2_dp, 2.5_dp], [2, 2])
        matrix(3:5, 3:5) = reshape([4.0_dp, 0.1_dp, -0.2_dp, 0.1_dp, &
            3.5_dp, 0.3_dp, -0.2_dp, 0.3_dp, 2.8_dp], [3, 3])
        matrix(6, 6) = 2.2_dp
        matrix(3:5, 1:2) = reshape([-0.4_dp, 0.1_dp, 0.0_dp, 0.2_dp, &
            -0.3_dp, 0.05_dp], [3, 2])
        do column = 3, 5
            do row = 1, 2
                matrix(row, column) = matrix(column, row)
            end do
        end do
        matrix(6, 3:5) = [0.1_dp, -0.2_dp, 0.3_dp]
        matrix(3:5, 6) = matrix(6, 3:5)
    end subroutine build_dense_fixture

    subroutine place_in_canonical_order(blocked, order, canonical)
        real(dp), intent(in) :: blocked(:, :)
        integer, intent(in) :: order(:)
        real(dp), intent(out) :: canonical(:, :)
        integer :: i, j

        canonical = 0.0_dp
        do j = 1, size(order)
            do i = 1, size(order)
                canonical(order(i), order(j)) = blocked(i, j)
            end do
        end do
    end subroutine place_in_canonical_order

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_variable_block_tridiagonal
