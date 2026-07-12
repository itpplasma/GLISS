module physical_mass_family_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: add_mapped_dynamic_element, &
        build_dynamic_element_map, build_resolved_dynamic_family_layout, &
        dynamic_element_map_is_valid, dynamic_family_layout_t, &
        dynamic_layout_ok
    use mass_density_policy, only: evaluate_mass_density, mass_density_ok, &
        mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use physical_mass_assembly, only: &
        assemble_physical_mass_surface_resolved
    use radial_slice_interpolation, only: gauss_two_lower, &
        gauss_two_upper, interpolate_slice_pair
    use radial_space_policy, only: radial_space_config_t
    use trial_space_topology, only: build_trial_space_topology, &
        trial_space_topology_t, trial_topology_ok
    implicit none
    private

    real(dp), parameter :: radial_tolerance = 1.0e-12_dp

    public :: assemble_physical_family_mass
    public :: assemble_physical_family_mass_resolved
    public :: assemble_physical_family_mass_fixed_layout

contains

    subroutine assemble_physical_family_mass(fields, density_profile, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_space, radial_step, phase_assembly, mass, layout, info)
        real(dp), intent(in) :: fields(:, :, :, :)
        type(mass_density_profile_t), intent(in) :: density_profile
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_step
        real(dp), allocatable, intent(out) :: mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        type(trial_space_topology_t) :: topology
        real(dp), allocatable :: density_kg_m3(:)
        real(dp) :: radial_coordinate
        integer :: interval

        call build_trial_space_topology(trial_m, trial_n, trial_parity, &
            topology, info)
        if (info /= trial_topology_ok) return
        call build_resolved_dynamic_family_layout(topology, size(fields, 4), &
            layout, info)
        if (info /= dynamic_layout_ok) return
        allocate (density_kg_m3(layout%intervals))
        do interval = 1, layout%intervals
            radial_coordinate = (real(interval, dp) - 0.5_dp) * radial_step
            call evaluate_mass_density(density_profile, radial_coordinate, &
                density_kg_m3(interval), info)
            if (info /= mass_density_ok) then
                info = -1
                return
            end if
        end do
        allocate (mass(layout%total_unknowns, layout%total_unknowns))
        call assemble_physical_family_mass_resolved(fields, density_kg_m3, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_space, radial_step, phase_assembly, mass, info)
    end subroutine assemble_physical_family_mass

    subroutine assemble_physical_family_mass_resolved(fields, density_kg_m3, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_space, radial_step, phase_assembly, mass, info)
        real(dp), intent(in) :: fields(:, :, :, :), density_kg_m3(:)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: info
        type(dynamic_family_layout_t) :: layout
        type(trial_space_topology_t) :: topology
        integer, allocatable :: element_to_global(:, :)

        call build_trial_space_topology(trial_m, trial_n, trial_parity, &
            topology, info)
        if (info /= trial_topology_ok) return
        call build_resolved_dynamic_family_layout(topology, size(fields, 4), &
            layout, info)
        if (info /= dynamic_layout_ok) return
        if (any(shape(mass) /= layout%total_unknowns)) then
            info = -1
            return
        end if
        call build_dynamic_element_map(layout, element_to_global, info)
        if (info /= dynamic_layout_ok) return
        call assemble_physical_family_mass_fixed_layout(fields, density_kg_m3, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_space, radial_step, phase_assembly, element_to_global, &
            mass, info)
    end subroutine assemble_physical_family_mass_resolved

    subroutine assemble_physical_family_mass_fixed_layout(fields, &
            density_kg_m3, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_step, phase_assembly, &
            element_to_global, mass, info)
        real(dp), intent(in) :: fields(:, :, :, :), density_kg_m3(:)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_step
        integer, intent(in) :: element_to_global(:, :)
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: info
        type(radial_space_config_t) :: node_space
        real(dp) :: element(4 * size(trial_m), 4 * size(trial_m))
        real(dp) :: element_total(4 * size(trial_m), 4 * size(trial_m))
        real(dp) :: fields_node(size(fields, 1), size(fields, 2), &
            size(fields, 3))
        real(dp) :: nodes(2), radial_coordinate, blend, density_node
        integer :: interval, node, node_count, left, right

        call validate_fixed_inputs(fields, density_kg_m3, trial_m, &
            trial_n, trial_parity, stored_power, field_periods, radial_step, &
            phase_assembly, element_to_global, mass, info)
        if (info /= 0) return
        if (radial_space%quadrature_points == 2) then
            node_count = 2
            nodes(1) = gauss_two_lower
            nodes(2) = gauss_two_upper
        else
            node_count = 1
            nodes(1) = radial_space%evaluation_coordinate
        end if
        mass = 0.0_dp
        do interval = 1, size(fields, 4)
            element_total = 0.0_dp
            do node = 1, node_count
                node_space = radial_space
                node_space%evaluation_coordinate = nodes(node)
                node_space%weight_fraction = radial_space%weight_fraction &
                    / real(node_count, dp)
                radial_coordinate = (real(interval - 1, dp) + nodes(node)) &
                    * radial_step
                call interpolate_slice_pair(interval, size(fields, 4), &
                    nodes(node), left, right, blend)
                fields_node = (1.0_dp - blend) * fields(:, :, :, left) &
                    + blend * fields(:, :, :, right)
                density_node = (1.0_dp - blend) * density_kg_m3(left) &
                    + blend * density_kg_m3(right)
                call assemble_physical_mass_surface_resolved( &
                    fields_node, density_node, trial_m, &
                    trial_n, trial_parity, stored_power, field_periods, &
                    node_space, radial_coordinate, radial_step, &
                    phase_assembly, element, info)
                if (info /= 0) return
                element_total = element_total + element
            end do
            call add_mapped_dynamic_element(element_to_global(:, interval), &
                element_total, mass, info)
            if (info /= dynamic_layout_ok) return
        end do
        info = 0
    end subroutine assemble_physical_family_mass_fixed_layout

    subroutine validate_fixed_inputs(fields, density_kg_m3, trial_m, &
            trial_n, trial_parity, stored_power, field_periods, radial_step, &
            phase_assembly, element_to_global, mass, info)
        real(dp), intent(in) :: fields(:, :, :, :), density_kg_m3(:)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        integer, intent(in) :: field_periods, phase_assembly
        integer, intent(in) :: element_to_global(:, :)
        real(dp), intent(in) :: mass(:, :)
        integer, intent(out) :: info
        integer :: trials, intervals

        info = -1
        trials = size(trial_m)
        intervals = size(fields, 4)
        if (trials < 1 .or. intervals < 1) return
        if (size(trial_n) /= trials .or. size(trial_parity) /= trials) return
        if (size(stored_power) /= trials) return
        if (size(density_kg_m3) /= intervals) return
        if (.not. all(ieee_is_finite(density_kg_m3))) return
        if (any(density_kg_m3 <= 0.0_dp)) return
        if (.not. ieee_is_finite(radial_step)) return
        if (abs(radial_step * real(intervals, dp) - 1.0_dp) &
            > radial_tolerance) return
        if (field_periods < 1) return
        if (phase_assembly /= phase_assembly_transformed .and. &
            phase_assembly /= phase_assembly_direct) return
        if (size(mass, 1) /= size(mass, 2)) return
        if (.not. dynamic_element_map_is_valid(element_to_global, trials, &
            intervals, size(mass, 1))) return
        info = 0
    end subroutine validate_fixed_inputs

end module physical_mass_family_assembly
