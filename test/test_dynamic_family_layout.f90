program test_dynamic_family_layout
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: add_dynamic_element, &
        add_mapped_dynamic_element, build_dynamic_block_permutation, &
        build_dynamic_element_map, build_dynamic_family_layout, &
        build_resolved_dynamic_family_layout, dynamic_family_layout_t, &
        dynamic_element_map_is_valid, dynamic_layout_ok, eta_global_index, &
        mu_global_index, normal_global_index
    use trial_space_topology, only: trial_space_topology_t
    implicit none

    type(dynamic_family_layout_t) :: layout
    real(dp) :: element(12, 12), malformed_matrix(32, 32)
    real(dp) :: ragged_element(8, 8), ragged_matrix(11, 11)
    real(dp) :: ragged_expected(11, 11)
    integer, allocatable :: element_map(:, :), permutation(:), widths(:)
    type(trial_space_topology_t) :: topology
    integer, parameter :: local_to_global(8) = [1, 0, 2, 0, 0, 5, 0, 9]
    integer :: i, info, j

    call build_dynamic_family_layout(3, 4, layout, info)
    call require(info == dynamic_layout_ok, "valid dynamic layout failed")
    call require(layout%normal_unknowns == 9, "normal count is wrong")
    call require(layout%eta_unknowns == 12, "eta count is wrong")
    call require(layout%mu_unknowns == 12, "mu count is wrong")
    call require(layout%total_unknowns == 33, "total count is wrong")
    call require(normal_global_index(layout, 0, 1) == 0, &
        "axis normal coefficient was retained")
    call require(normal_global_index(layout, 4, 1) == 0, &
        "edge normal coefficient was retained")
    call require(normal_global_index(layout, 1, 1) == 1, &
        "first interior normal index is wrong")
    call require(normal_global_index(layout, 3, 3) == 9, &
        "last interior normal index is wrong")
    call require(eta_global_index(layout, 1, 1) == 10, &
        "first eta index is wrong")
    call require(eta_global_index(layout, 4, 3) == 21, &
        "last eta index is wrong")
    call require(mu_global_index(layout, 1, 1) == 22, &
        "first mu index is wrong")
    call require(mu_global_index(layout, 4, 3) == 33, &
        "last mu index is wrong")
    call check_retained_outer_normal()
    call check_malformed_indices(layout)
    element = 1.0_dp
    malformed_matrix = 0.0_dp
    layout%total_unknowns = 32
    call add_dynamic_element(layout, 4, element, malformed_matrix, info)
    call require(info /= dynamic_layout_ok, &
        "inconsistent public dynamic layout was accepted")
    call require(all(malformed_matrix == 0.0_dp), &
        "rejected dynamic layout modified the matrix")
    call build_dynamic_family_layout(0, 4, layout, info)
    call require(info /= dynamic_layout_ok, "zero-trial layout was accepted")
    allocate (topology%active(3, 2))
    topology%active(:, 1) = [.true., .false., .false.]
    topology%active(:, 2) = [.false., .true., .true.]
    call build_resolved_dynamic_family_layout(topology, 4, layout, info)
    call require(info == dynamic_layout_ok, "ragged dynamic layout failed")
    call require(layout%normal_unknowns == 3, &
        "ragged normal count is wrong")
    call require(layout%eta_unknowns == 4, "ragged eta count is wrong")
    call require(layout%mu_unknowns == 4, "ragged mu count is wrong")
    call require(layout%total_unknowns == 11, "ragged total count is wrong")
    call require(normal_global_index(layout, 1, 1) == 1 .and. &
        normal_global_index(layout, 1, 2) == 0, &
        "ragged normal indices are wrong")
    call require(eta_global_index(layout, 1, 1) == 0 .and. &
        eta_global_index(layout, 1, 2) == 4, &
        "ragged eta indices are wrong")
    call require(mu_global_index(layout, 4, 1) == 0 .and. &
        mu_global_index(layout, 4, 2) == 11, &
        "ragged mu indices are wrong")
    call build_dynamic_element_map(layout, element_map, info)
    call require(info == dynamic_layout_ok, "ragged element map failed")
    call require(all(element_map(:, 2) == local_to_global), &
        "ragged element map changed component coordinates")
    call require(dynamic_element_map_is_valid(element_map, 2, 4, 11), &
        "valid ragged element map was rejected")
    call build_dynamic_block_permutation(layout, widths, permutation, info)
    call require(info == dynamic_layout_ok, &
        "ragged block permutation failed")
    call require(all(widths == [3, 3, 3, 2]), &
        "ragged block widths are wrong")
    call require(all(permutation == [1, 4, 8, 2, 5, 9, 3, 6, 10, 7, 11]), &
        "ragged block permutation is wrong")
    do j = 1, 8
        do i = 1, 8
            ragged_element(i, j) = real(100 * i + j, dp)
        end do
    end do
    ragged_matrix = 0.0_dp
    ragged_expected = 0.0_dp
    call add_dynamic_element(layout, 2, ragged_element, ragged_matrix, info)
    call require(info == dynamic_layout_ok, "ragged element gather failed")
    do j = 1, 8
        if (local_to_global(j) == 0) cycle
        do i = 1, 8
            if (local_to_global(i) == 0) cycle
            ragged_expected(local_to_global(i), local_to_global(j)) = &
                ragged_element(i, j)
        end do
    end do
    call require(all(ragged_matrix == ragged_expected), &
        "ragged element gather changed component coordinates")
    ragged_matrix = 0.0_dp
    element_map(1, 2) = -1
    call require(.not. dynamic_element_map_is_valid(element_map, 2, 4, 11), &
        "negative element map was accepted")
    call add_mapped_dynamic_element(element_map(:, 2), ragged_element, &
        ragged_matrix, info)
    call require(info /= dynamic_layout_ok, &
        "negative mapped element index was accepted")
    call require(all(ragged_matrix == 0.0_dp), &
        "rejected mapped element modified the matrix")
    topology%active(:, 2) = [.false., .false., .false.]
    call build_resolved_dynamic_family_layout(topology, 4, layout, info)
    call require(info == dynamic_layout_ok, "normal-only layout failed")
    call build_dynamic_block_permutation(layout, widths, permutation, info)
    call require(info == dynamic_layout_ok, &
        "normal-only block permutation failed")
    call require(all(widths == [1, 1, 1]), &
        "empty edge block was not removed")
    call require(all(permutation == [1, 2, 3]), &
        "normal-only block permutation is wrong")

    write (*, "(a)") "PASS"

contains

    subroutine check_retained_outer_normal()
        type(dynamic_family_layout_t) :: free_layout
        integer, allocatable :: free_permutation(:), free_widths(:)
        integer, parameter :: expected(36) = [ &
            1, 2, 3, 13, 14, 15, 25, 26, 27, &
            4, 5, 6, 16, 17, 18, 28, 29, 30, &
            7, 8, 9, 19, 20, 21, 31, 32, 33, &
            10, 11, 12, 22, 23, 24, 34, 35, 36]
        integer :: local_info

        call build_dynamic_family_layout(3, 4, free_layout, local_info, &
            retain_outer_normal=.true.)
        call require(local_info == dynamic_layout_ok, &
            "retained-edge dynamic layout failed")
        call require(free_layout%normal_unknowns == 12 &
            .and. free_layout%total_unknowns == 36, &
            "retained-edge dynamic counts are wrong")
        call require(normal_global_index(free_layout, 4, 1) == 10 &
            .and. normal_global_index(free_layout, 4, 3) == 12, &
            "retained edge normal indices are wrong")
        call build_dynamic_block_permutation(free_layout, free_widths, &
            free_permutation, local_info)
        call require(local_info == dynamic_layout_ok, &
            "retained-edge block permutation failed")
        call require(all(free_widths == [9, 9, 9, 9]), &
            "retained-edge block widths are wrong")
        call require(all(free_permutation == expected), &
            "retained-edge block permutation is wrong")
    end subroutine check_retained_outer_normal

    subroutine check_malformed_indices(valid_layout)
        type(dynamic_family_layout_t), intent(in) :: valid_layout
        type(dynamic_family_layout_t) :: corrupt

        corrupt = valid_layout
        corrupt%total_unknowns = corrupt%total_unknowns - 1
        call require_indices_rejected(corrupt, "total count")
        corrupt = valid_layout
        corrupt%normal_unknowns = corrupt%normal_unknowns - 1
        call require_indices_rejected(corrupt, "normal count")
        corrupt = valid_layout
        corrupt%eta_unknowns = corrupt%eta_unknowns - 1
        call require_indices_rejected(corrupt, "eta count")
        corrupt = valid_layout
        corrupt%trials = corrupt%trials - 1
        call require_indices_rejected(corrupt, "trial count")
        corrupt = valid_layout
        corrupt%intervals = corrupt%intervals - 1
        call require_indices_rejected(corrupt, "interval count")
        corrupt = valid_layout
        corrupt%active_count(1) = corrupt%active_count(1) - 1
        call require_indices_rejected(corrupt, "active count")
        corrupt = valid_layout
        corrupt%active_rank(1, 1) = corrupt%active_rank(1, 1) + 1
        call require_indices_rejected(corrupt, "active rank")
        corrupt = valid_layout
        deallocate (corrupt%active)
        allocate (corrupt%active(2, valid_layout%trials), source=.true.)
        call require_indices_rejected(corrupt, "activity shape")
        corrupt = valid_layout
        deallocate (corrupt%active_rank)
        allocate (corrupt%active_rank(2, valid_layout%trials), source=1)
        call require_indices_rejected(corrupt, "rank shape")
    end subroutine check_malformed_indices

    subroutine require_indices_rejected(layout, label)
        type(dynamic_family_layout_t), intent(in) :: layout
        character(len=*), intent(in) :: label

        call require(normal_global_index(layout, 1, 1) == 0, &
            "malformed " // label // " produced a normal index")
        call require(eta_global_index(layout, 1, 1) == 0, &
            "malformed " // label // " produced an eta index")
        call require(mu_global_index(layout, 1, 1) == 0, &
            "malformed " // label // " produced a mu index")
    end subroutine require_indices_rejected

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_dynamic_family_layout
