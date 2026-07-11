module dynamic_family_layout
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
    public :: eta_global_index
    public :: mu_global_index
    public :: normal_global_index

contains

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
