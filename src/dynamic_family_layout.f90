module dynamic_family_layout
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use trial_space_topology, only: trial_space_topology_t
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
        integer :: active_count(3) = 0
        logical :: outer_normal_retained = .false.
        logical, allocatable :: active(:, :)
        integer, allocatable :: active_rank(:, :)
    end type dynamic_family_layout_t

    public :: build_dynamic_family_layout
    public :: build_resolved_dynamic_family_layout
    public :: build_dynamic_block_permutation
    public :: build_dynamic_element_map
    public :: add_dynamic_element
    public :: add_mapped_dynamic_element
    public :: dynamic_element_map_is_valid
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
        integer :: local_a
        integer :: local_to_global(size(element, 1))

        info = dynamic_layout_invalid
        if (.not. dynamic_layout_is_consistent(layout)) return
        if (interval < 1 .or. interval > layout%intervals) return
        if (any(shape(element) /= 4 * layout%trials)) return
        if (any(shape(matrix) /= layout%total_unknowns)) return
        do local_a = 1, size(element, 1)
            local_to_global(local_a) = &
                element_global_index(layout, interval, local_a)
        end do
        call add_mapped_dynamic_element(local_to_global, element, matrix, info)
    end subroutine add_dynamic_element

    subroutine add_mapped_dynamic_element(local_to_global, element, matrix, &
            info)
        integer, intent(in) :: local_to_global(:)
        real(dp), intent(in) :: element(:, :)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(out) :: info
        integer :: local_a, local_b, global_a, global_b

        info = dynamic_layout_invalid
        if (any(shape(element) /= size(local_to_global))) return
        if (size(matrix, 1) /= size(matrix, 2)) return
        if (any(local_to_global < 0)) return
        if (any(local_to_global > size(matrix, 1))) return
        do local_b = 1, size(element, 2)
            global_b = local_to_global(local_b)
            if (global_b == 0) cycle
            do local_a = 1, size(element, 1)
                global_a = local_to_global(local_a)
                if (global_a == 0) cycle
                matrix(global_a, global_b) = matrix(global_a, global_b) &
                    + element(local_a, local_b)
            end do
        end do
        info = dynamic_layout_ok
    end subroutine add_mapped_dynamic_element

    subroutine build_dynamic_element_map(layout, element_to_global, info)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, allocatable, intent(out) :: element_to_global(:, :)
        integer, intent(out) :: info
        integer :: interval, local

        info = dynamic_layout_invalid
        if (.not. dynamic_layout_is_consistent(layout)) return
        allocate (element_to_global(4 * layout%trials, layout%intervals))
        do interval = 1, layout%intervals
            do local = 1, size(element_to_global, 1)
                element_to_global(local, interval) = &
                    element_global_index(layout, interval, local)
            end do
        end do
        info = dynamic_layout_ok
    end subroutine build_dynamic_element_map

    pure function dynamic_element_map_is_valid(element_to_global, trials, &
            intervals, total_unknowns) result(valid)
        integer, intent(in) :: element_to_global(:, :)
        integer, intent(in) :: trials, intervals, total_unknowns
        logical :: valid
        integer :: global

        valid = .false.
        if (trials < 1 .or. intervals < 2 .or. total_unknowns < 1) return
        if (size(element_to_global, 1) /= 4 * trials) return
        if (size(element_to_global, 2) /= intervals) return
        if (any(element_to_global < 0)) return
        if (any(element_to_global > total_unknowns)) return
        if (maxval(element_to_global) /= total_unknowns) return
        do global = 1, total_unknowns
            if (.not. any(element_to_global == global)) return
        end do
        valid = .true.
    end function dynamic_element_map_is_valid

    subroutine build_dynamic_block_permutation(layout, widths, permutation, &
            info)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, allocatable, intent(out) :: widths(:), permutation(:)
        integer, intent(out) :: info
        integer :: block_count, cell, edge_width, position, regular_width

        info = dynamic_layout_invalid
        if (.not. dynamic_layout_is_consistent(layout)) return
        regular_width = sum(layout%active_count)
        edge_width = sum(layout%active_count(2:3))
        block_count = layout%intervals - 1
        if (edge_width > 0 .or. layout%outer_normal_retained) &
            block_count = layout%intervals
        allocate (widths(block_count), permutation(layout%total_unknowns))
        widths(1:layout%intervals - 1) = regular_width
        if (block_count == layout%intervals) then
            widths(block_count) = edge_width
            if (layout%outer_normal_retained) widths(block_count) = regular_width
        end if
        position = 0
        do cell = 1, layout%intervals - 1
            call append_cell_indices(layout, cell, .true., permutation, &
                position)
        end do
        if (block_count == layout%intervals) call append_cell_indices(layout, &
            layout%intervals, layout%outer_normal_retained, permutation, &
            position)
        if (position /= layout%total_unknowns) return
        info = dynamic_layout_ok
    end subroutine build_dynamic_block_permutation

    pure subroutine append_cell_indices(layout, cell, include_normal, &
            permutation, position)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: cell
        logical, intent(in) :: include_normal
        integer, intent(inout) :: permutation(:), position
        integer :: component, index, trial

        do component = 1, 3
            if (component == 1 .and. .not. include_normal) cycle
            do trial = 1, layout%trials
                select case (component)
                case (1)
                    index = normal_global_index_unchecked(layout, cell, trial)
                case (2)
                    index = eta_global_index_unchecked(layout, cell, trial)
                case (3)
                    index = mu_global_index_unchecked(layout, cell, trial)
                end select
                if (index == 0) cycle
                position = position + 1
                permutation(position) = index
            end do
        end do
    end subroutine append_cell_indices

    pure function element_global_index(layout, interval, local) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: interval, local
        integer :: index, component_block, trial

        component_block = (local - 1) / layout%trials + 1
        trial = modulo(local - 1, layout%trials) + 1
        select case (component_block)
        case (1)
            index = normal_global_index_unchecked(layout, interval - 1, trial)
        case (2)
            index = normal_global_index_unchecked(layout, interval, trial)
        case (3)
            index = eta_global_index_unchecked(layout, interval, trial)
        case (4)
            index = mu_global_index_unchecked(layout, interval, trial)
        case default
            index = 0
        end select
    end function element_global_index

    pure subroutine build_dynamic_family_layout(trials, intervals, layout, &
            info, retain_outer_normal)
        integer, intent(in) :: trials, intervals
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        logical, intent(in), optional :: retain_outer_normal

        type(trial_space_topology_t) :: topology

        info = dynamic_layout_invalid
        if (trials < 1) return
        allocate (topology%active(3, trials), source=.true.)
        call build_resolved_dynamic_family_layout(topology, intervals, &
            layout, info, retain_outer_normal)
    end subroutine build_dynamic_family_layout

    pure subroutine build_resolved_dynamic_family_layout(topology, intervals, &
            layout, info, retain_outer_normal)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: intervals
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        logical, intent(in), optional :: retain_outer_normal
        integer :: normal_nodes

        layout = dynamic_family_layout_t()
        info = dynamic_layout_invalid
        if (.not. allocated(topology%active)) return
        if (size(topology%active, 1) /= 3 &
            .or. size(topology%active, 2) < 1) return
        if (intervals < 2 .or. .not. any(topology%active)) return
        layout%trials = size(topology%active, 2)
        layout%intervals = intervals
        if (present(retain_outer_normal)) &
            layout%outer_normal_retained = retain_outer_normal
        allocate (layout%active, source=topology%active)
        call build_activity_ranks(layout)
        normal_nodes = intervals - 1
        if (layout%outer_normal_retained) normal_nodes = intervals
        layout%normal_unknowns = layout%active_count(1) * normal_nodes
        layout%eta_unknowns = layout%active_count(2) * intervals
        layout%mu_unknowns = layout%active_count(3) * intervals
        layout%total_unknowns = layout%normal_unknowns &
            + layout%eta_unknowns + layout%mu_unknowns
        info = dynamic_layout_ok
    end subroutine build_resolved_dynamic_family_layout

    pure subroutine build_activity_ranks(layout)
        type(dynamic_family_layout_t), intent(inout) :: layout
        integer :: component, rank, trial

        allocate (layout%active_rank(3, layout%trials), source=0)
        do component = 1, 3
            rank = 0
            do trial = 1, layout%trials
                if (.not. layout%active(component, trial)) cycle
                rank = rank + 1
                layout%active_rank(component, trial) = rank
            end do
            layout%active_count(component) = rank
        end do
    end subroutine build_activity_ranks

    pure function dynamic_layout_is_consistent(layout) result(consistent)
        type(dynamic_family_layout_t), intent(in) :: layout
        logical :: consistent

        consistent = layout_metadata_is_consistent(layout)
    end function dynamic_layout_is_consistent

    pure function normal_global_index(layout, node, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: node, trial
        integer :: index

        index = 0
        if (.not. layout_metadata_is_consistent(layout)) return
        index = normal_global_index_unchecked(layout, node, trial)
    end function normal_global_index

    pure function normal_global_index_unchecked(layout, node, trial) &
            result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: node, trial
        integer :: index

        index = 0
        if (node <= 0 .or. node > layout%intervals) return
        if (node == layout%intervals .and. &
            .not. layout%outer_normal_retained) return
        if (trial < 1 .or. trial > layout%trials) return
        if (.not. layout%active(1, trial)) return
        index = (node - 1) * layout%active_count(1) &
            + layout%active_rank(1, trial)
    end function normal_global_index_unchecked

    pure function eta_global_index(layout, cell, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: cell, trial
        integer :: index

        index = 0
        if (.not. layout_metadata_is_consistent(layout)) return
        index = eta_global_index_unchecked(layout, cell, trial)
    end function eta_global_index

    pure function eta_global_index_unchecked(layout, cell, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: cell, trial
        integer :: index

        index = 0
        if (cell < 1 .or. cell > layout%intervals) return
        if (trial < 1 .or. trial > layout%trials) return
        if (.not. layout%active(2, trial)) return
        index = layout%normal_unknowns &
            + (cell - 1) * layout%active_count(2) &
            + layout%active_rank(2, trial)
    end function eta_global_index_unchecked

    pure function mu_global_index(layout, cell, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: cell, trial
        integer :: index

        index = 0
        if (.not. layout_metadata_is_consistent(layout)) return
        index = mu_global_index_unchecked(layout, cell, trial)
    end function mu_global_index

    pure function mu_global_index_unchecked(layout, cell, trial) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: cell, trial
        integer :: index

        index = 0
        if (cell < 1 .or. cell > layout%intervals) return
        if (trial < 1 .or. trial > layout%trials) return
        if (.not. layout%active(3, trial)) return
        index = layout%normal_unknowns + layout%eta_unknowns &
            + (cell - 1) * layout%active_count(3) &
            + layout%active_rank(3, trial)
    end function mu_global_index_unchecked

    pure function layout_activity_is_shaped(layout) result(valid)
        type(dynamic_family_layout_t), intent(in) :: layout
        logical :: valid

        valid = allocated(layout%active)
        if (.not. valid) return
        valid = size(layout%active, 1) == 3
        if (.not. valid) return
        valid = size(layout%active, 2) == layout%trials
        if (.not. valid) return
        valid = allocated(layout%active_rank)
        if (.not. valid) return
        valid = size(layout%active_rank, 1) == 3 &
            .and. size(layout%active_rank, 2) == layout%trials
    end function layout_activity_is_shaped

    pure function layout_metadata_is_consistent(layout) result(valid)
        type(dynamic_family_layout_t), intent(in) :: layout
        logical :: valid

        valid = layout_activity_is_shaped(layout)
        if (.not. valid) return
        valid = layout%trials >= 1 .and. layout%intervals >= 2
        if (.not. valid) return
        valid = any(layout%active)
        if (.not. valid) return
        valid = activity_ranks_are_consistent(layout)
        if (.not. valid) return
        if (layout%outer_normal_retained) then
            valid = layout%normal_unknowns == layout%active_count(1) &
                * layout%intervals
        else
            valid = layout%normal_unknowns == layout%active_count(1) &
                * (layout%intervals - 1)
        end if
        if (.not. valid) return
        valid = layout%eta_unknowns == layout%active_count(2) &
            * layout%intervals
        if (.not. valid) return
        valid = layout%mu_unknowns == layout%active_count(3) &
            * layout%intervals
        if (.not. valid) return
        valid = layout%total_unknowns == layout%normal_unknowns &
            + layout%eta_unknowns + layout%mu_unknowns
    end function layout_metadata_is_consistent

    pure function activity_ranks_are_consistent(layout) result(valid)
        type(dynamic_family_layout_t), intent(in) :: layout
        logical :: valid
        integer :: component, rank, trial

        valid = .false.
        do component = 1, 3
            rank = 0
            do trial = 1, layout%trials
                if (layout%active(component, trial)) rank = rank + 1
                if (layout%active_rank(component, trial) /= merge(rank, 0, &
                    layout%active(component, trial))) return
            end do
            if (layout%active_count(component) /= rank) return
        end do
        valid = .true.
    end function activity_ranks_are_consistent

end module dynamic_family_layout
