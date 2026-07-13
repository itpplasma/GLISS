module terpsichore_reduced_layout
    use dynamic_family_layout, only: build_dynamic_element_map, &
        build_resolved_dynamic_family_layout, dynamic_family_layout_t, &
        dynamic_layout_ok
    use trial_space_topology, only: build_trial_space_topology, &
        trial_space_topology_t, trial_topology_ok
    implicit none
    private

    integer, parameter, public :: terpsichore_reduced_layout_ok = 0
    integer, parameter, public :: terpsichore_reduced_layout_invalid = -1

    public :: build_terpsichore_reduced_fixed_boundary_layout
    public :: build_terpsichore_reduced_free_boundary_layout

contains

    subroutine build_terpsichore_reduced_fixed_boundary_layout(mode_m, mode_n, &
            parity, intervals, layout, element_to_global, info)
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:), intervals
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, allocatable, intent(out) :: element_to_global(:, :)
        integer, intent(out) :: info
        call build_terpsichore_reduced_layout(mode_m, mode_n, parity, &
            intervals, .false., layout, element_to_global, info)
    end subroutine build_terpsichore_reduced_fixed_boundary_layout

    subroutine build_terpsichore_reduced_free_boundary_layout(mode_m, mode_n, &
            parity, intervals, layout, element_to_global, info)
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:), intervals
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, allocatable, intent(out) :: element_to_global(:, :)
        integer, intent(out) :: info

        call build_terpsichore_reduced_layout(mode_m, mode_n, parity, &
            intervals, .true., layout, element_to_global, info)
    end subroutine build_terpsichore_reduced_free_boundary_layout

    subroutine build_terpsichore_reduced_layout(mode_m, mode_n, parity, &
            intervals, retain_outer_normal, layout, element_to_global, info)
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:), intervals
        logical, intent(in) :: retain_outer_normal
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, allocatable, intent(out) :: element_to_global(:, :)
        integer, intent(out) :: info
        integer, allocatable :: full_map(:, :)
        type(trial_space_topology_t) :: topology

        info = terpsichore_reduced_layout_invalid
        call build_trial_space_topology(mode_m, mode_n, parity, topology, info)
        if (info /= trial_topology_ok) return
        topology%active(3, :) = .false.
        call build_resolved_dynamic_family_layout(topology, intervals, layout, &
            info, retain_outer_normal=retain_outer_normal)
        if (info /= dynamic_layout_ok) return
        call build_dynamic_element_map(layout, full_map, info)
        if (info /= dynamic_layout_ok) return
        allocate (element_to_global(3 * size(mode_m), intervals), &
            source=full_map(:3 * size(mode_m), :))
        info = terpsichore_reduced_layout_ok
    end subroutine build_terpsichore_reduced_layout

end module terpsichore_reduced_layout
