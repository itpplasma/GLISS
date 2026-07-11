module dynamic_family_layout
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: dynamic_layout_ok = 0
    integer, parameter, public :: dynamic_layout_invalid = -1

    type, public :: dynamic_family_layout_t
        integer :: trials = 0
        integer :: intervals = 0
        integer :: normal_unknowns = 0
        integer :: eta_unknowns = 0
        integer :: mu_unknowns = 0
        integer :: total_unknowns = 0
    end type dynamic_family_layout_t

    public :: build_dynamic_family_layout
    public :: build_dynamic_block_permutation
    public :: add_dynamic_element
    public :: eta_global_index
    public :: mu_global_index
    public :: normal_global_index

contains

    subroutine add_dynamic_element(layout, interval, element, matrix, info)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: interval
        real(dp), intent(in) :: element(:, :)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(out) :: info
        integer :: local_a, local_b, global_a, global_b

        info = dynamic_layout_invalid
        if (.not. dynamic_layout_is_consistent(layout)) return
        if (interval < 1 .or. interval > layout%intervals) return
        if (any(shape(element) /= 4 * layout%trials)) return
        if (any(shape(matrix) /= layout%total_unknowns)) return
        do local_b = 1, size(element, 2)
            global_b = element_global_index(layout, interval, local_b)
            if (global_b == 0) cycle
            do local_a = 1, size(element, 1)
                global_a = element_global_index(layout, interval, local_a)
                if (global_a == 0) cycle
                matrix(global_a, global_b) = matrix(global_a, global_b) &
                    + element(local_a, local_b)
            end do
        end do
        info = dynamic_layout_ok
    end subroutine add_dynamic_element

    subroutine build_dynamic_block_permutation(layout, widths, permutation, &
            info)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, allocatable, intent(out) :: widths(:), permutation(:)
        integer, intent(out) :: info
        integer :: base, cell, trial

        info = dynamic_layout_invalid
        if (.not. dynamic_layout_is_consistent(layout)) return
        allocate (widths(layout%intervals), &
            permutation(layout%total_unknowns))
        widths(1:layout%intervals - 1) = 3 * layout%trials
        widths(layout%intervals) = 2 * layout%trials
        do cell = 1, layout%intervals - 1
            base = 3 * layout%trials * (cell - 1)
            do trial = 1, layout%trials
                permutation(base + trial) = &
                    normal_global_index(layout, cell, trial)
                permutation(base + layout%trials + trial) = &
                    eta_global_index(layout, cell, trial)
                permutation(base + 2 * layout%trials + trial) = &
                    mu_global_index(layout, cell, trial)
            end do
        end do
        base = 3 * layout%trials * (layout%intervals - 1)
        do trial = 1, layout%trials
            permutation(base + trial) = &
                eta_global_index(layout, layout%intervals, trial)
            permutation(base + layout%trials + trial) = &
                mu_global_index(layout, layout%intervals, trial)
        end do
        info = dynamic_layout_ok
    end subroutine build_dynamic_block_permutation

    pure function element_global_index(layout, interval, local) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: interval, local
        integer :: index, component_block, trial

        component_block = (local - 1) / layout%trials + 1
        trial = modulo(local - 1, layout%trials) + 1
        select case (component_block)
        case (1)
            index = normal_global_index(layout, interval - 1, trial)
        case (2)
            index = normal_global_index(layout, interval, trial)
        case (3)
            index = eta_global_index(layout, interval, trial)
        case (4)
            index = mu_global_index(layout, interval, trial)
        case default
            index = 0
        end select
    end function element_global_index

    pure subroutine build_dynamic_family_layout(trials, intervals, layout, info)
        integer, intent(in) :: trials, intervals
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info

        layout = dynamic_family_layout_t()
        info = dynamic_layout_invalid
        if (trials < 1 .or. intervals < 2) return
        layout%trials = trials
        layout%intervals = intervals
        layout%normal_unknowns = trials * (intervals - 1)
        layout%eta_unknowns = trials * intervals
        layout%mu_unknowns = trials * intervals
        layout%total_unknowns = layout%normal_unknowns &
            + layout%eta_unknowns + layout%mu_unknowns
        info = dynamic_layout_ok
    end subroutine build_dynamic_family_layout

    pure function dynamic_layout_is_consistent(layout) result(consistent)
        type(dynamic_family_layout_t), intent(in) :: layout
        type(dynamic_family_layout_t) :: expected
        logical :: consistent
        integer :: info

        call build_dynamic_family_layout(layout%trials, layout%intervals, &
            expected, info)
        consistent = info == dynamic_layout_ok
        if (.not. consistent) return
        consistent = layout%normal_unknowns == expected%normal_unknowns
        if (.not. consistent) return
        consistent = layout%eta_unknowns == expected%eta_unknowns
        if (.not. consistent) return
        consistent = layout%mu_unknowns == expected%mu_unknowns
        if (.not. consistent) return
        consistent = layout%total_unknowns == expected%total_unknowns
    end function dynamic_layout_is_consistent

    pure function normal_global_index(layout, node, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: node, trial
        integer :: index

        index = 0
        if (node <= 0 .or. node >= layout%intervals) return
        if (trial < 1 .or. trial > layout%trials) return
        index = (node - 1) * layout%trials + trial
    end function normal_global_index

    pure function eta_global_index(layout, cell, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: cell, trial
        integer :: index

        index = 0
        if (cell < 1 .or. cell > layout%intervals) return
        if (trial < 1 .or. trial > layout%trials) return
        index = layout%normal_unknowns + (cell - 1) * layout%trials + trial
    end function eta_global_index

    pure function mu_global_index(layout, cell, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: cell, trial
        integer :: index

        index = 0
        if (cell < 1 .or. cell > layout%intervals) return
        if (trial < 1 .or. trial > layout%trials) return
        index = layout%normal_unknowns + layout%eta_unknowns &
            + (cell - 1) * layout%trials + trial
    end function mu_global_index

end module dynamic_family_layout
