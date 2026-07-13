module compressible_stiffness_family_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compressible_stiffness_assembly, only: &
        assemble_compressible_stiffness_surface_resolved, &
        compressible_stiffness_term_count
    use dynamic_family_layout, only: add_mapped_dynamic_element, &
        build_dynamic_element_map, build_resolved_dynamic_family_layout, &
        dynamic_element_map_is_valid, dynamic_family_layout_t, &
        dynamic_layout_ok
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use radial_slice_interpolation, only: gauss_two_lower, &
        gauss_two_upper, interpolate_slice_pair
    use radial_space_policy, only: radial_space_config_t
    use trial_space_topology, only: build_trial_space_topology, &
        trial_space_topology_t, trial_topology_ok
    implicit none
    private

    real(dp), parameter :: radial_tolerance = 1.0e-12_dp
    integer, parameter, public :: compressible_family_ok = 0
    integer, parameter, public :: compressible_family_invalid = -1
    integer, parameter, public :: compressible_family_allocation_error = -2

    public :: assemble_compressible_family_stiffness
    public :: assemble_compressible_family_stiffness_with_terms
    public :: assemble_compressible_family_stiffness_resolved
    public :: assemble_compressible_family_stiffness_fixed_layout

contains

    subroutine assemble_compressible_family_stiffness(fields, drive, &
            signed_sqrtg_radial, signed_sqrtg_theta, signed_sqrtg_zeta, &
            gamma_pressure_pa, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_step, phase_assembly, &
            stiffness, layout, info)
        real(dp), intent(in) :: fields(:, :, :, :), drive(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_radial(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_theta(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_zeta(:, :, :)
        real(dp), intent(in) :: gamma_pressure_pa(:, :, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_step
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        type(trial_space_topology_t) :: topology

        call build_trial_space_topology(trial_m, trial_n, trial_parity, &
            topology, info)
        if (info /= trial_topology_ok) return
        call build_resolved_dynamic_family_layout(topology, size(fields, 4), &
            layout, info)
        if (info /= dynamic_layout_ok) return
        allocate (stiffness(layout%total_unknowns, layout%total_unknowns))
        call assemble_compressible_family_stiffness_resolved(fields, drive, &
            signed_sqrtg_radial, signed_sqrtg_theta, signed_sqrtg_zeta, &
            gamma_pressure_pa, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_step, phase_assembly, &
            stiffness, info)
    end subroutine assemble_compressible_family_stiffness

    subroutine assemble_compressible_family_stiffness_with_terms(fields, &
            drive, signed_sqrtg_radial, signed_sqrtg_theta, &
            signed_sqrtg_zeta, gamma_pressure_pa, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_step, phase_assembly, stiffness, stiffness_terms, layout, &
            info)
        real(dp), intent(in) :: fields(:, :, :, :), drive(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_radial(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_theta(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_zeta(:, :, :)
        real(dp), intent(in) :: gamma_pressure_pa(:, :, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        real(dp), allocatable, intent(out) :: stiffness_terms(:, :, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        type(trial_space_topology_t) :: topology
        integer :: allocation_status

        call build_trial_space_topology(trial_m, trial_n, trial_parity, &
            topology, info)
        if (info /= trial_topology_ok) return
        call build_resolved_dynamic_family_layout(topology, size(fields, 4), &
            layout, info)
        if (info /= dynamic_layout_ok) return
        allocate (stiffness(layout%total_unknowns, layout%total_unknowns), &
            stiffness_terms(layout%total_unknowns, layout%total_unknowns, &
            compressible_stiffness_term_count), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compressible_family_allocation_error
            return
        end if
        call assemble_compressible_family_stiffness_resolved(fields, drive, &
            signed_sqrtg_radial, signed_sqrtg_theta, signed_sqrtg_zeta, &
            gamma_pressure_pa, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_step, phase_assembly, &
            stiffness, info, stiffness_terms)
    end subroutine assemble_compressible_family_stiffness_with_terms

    subroutine assemble_compressible_family_stiffness_resolved(fields, &
            drive, signed_sqrtg_radial, signed_sqrtg_theta, &
            signed_sqrtg_zeta, gamma_pressure_pa, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_step, phase_assembly, stiffness, info, stiffness_terms)
        real(dp), intent(in) :: fields(:, :, :, :), drive(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_radial(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_theta(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_zeta(:, :, :)
        real(dp), intent(in) :: gamma_pressure_pa(:, :, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(out) :: stiffness_terms(:, :, :)
        type(dynamic_family_layout_t) :: layout
        type(trial_space_topology_t) :: topology
        integer, allocatable :: element_to_global(:, :)

        call build_trial_space_topology(trial_m, trial_n, trial_parity, &
            topology, info)
        if (info /= trial_topology_ok) return
        call build_resolved_dynamic_family_layout(topology, size(fields, 4), &
            layout, info)
        if (info /= dynamic_layout_ok) return
        if (any(shape(stiffness) /= layout%total_unknowns)) then
            info = compressible_family_invalid
            return
        end if
        if (present(stiffness_terms)) then
            if (any(shape(stiffness_terms) /= [layout%total_unknowns, &
                layout%total_unknowns, compressible_stiffness_term_count])) then
                info = compressible_family_invalid
                return
            end if
        end if
        call build_dynamic_element_map(layout, element_to_global, info)
        if (info /= dynamic_layout_ok) return
        if (present(stiffness_terms)) then
            call assemble_compressible_family_stiffness_fixed_layout(fields, &
                drive, signed_sqrtg_radial, signed_sqrtg_theta, &
                signed_sqrtg_zeta, gamma_pressure_pa, trial_m, trial_n, &
                trial_parity, stored_power, field_periods, radial_space, &
                radial_step, phase_assembly, element_to_global, stiffness, &
                info, stiffness_terms)
        else
            call assemble_compressible_family_stiffness_fixed_layout(fields, &
                drive, signed_sqrtg_radial, signed_sqrtg_theta, &
                signed_sqrtg_zeta, gamma_pressure_pa, trial_m, trial_n, &
                trial_parity, stored_power, field_periods, radial_space, &
                radial_step, phase_assembly, element_to_global, stiffness, info)
        end if
    end subroutine assemble_compressible_family_stiffness_resolved

    subroutine assemble_compressible_family_stiffness_fixed_layout(fields, &
            drive, signed_sqrtg_radial, signed_sqrtg_theta, &
            signed_sqrtg_zeta, gamma_pressure_pa, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_step, phase_assembly, element_to_global, stiffness, info, &
            stiffness_terms)
        real(dp), intent(in) :: fields(:, :, :, :), drive(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_radial(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_theta(:, :, :)
        real(dp), intent(in) :: signed_sqrtg_zeta(:, :, :)
        real(dp), intent(in) :: gamma_pressure_pa(:, :, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_step
        integer, intent(in) :: element_to_global(:, :)
        real(dp), intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(out) :: stiffness_terms(:, :, :)
        type(radial_space_config_t) :: node_space
        real(dp) :: element(4 * size(trial_m), 4 * size(trial_m))
        real(dp) :: element_total(4 * size(trial_m), 4 * size(trial_m))
        real(dp) :: element_terms(4 * size(trial_m), 4 * size(trial_m), &
            compressible_stiffness_term_count)
        real(dp) :: element_terms_total(4 * size(trial_m), &
            4 * size(trial_m), compressible_stiffness_term_count)
        real(dp) :: fields_node(size(fields, 1), size(fields, 2), &
            size(fields, 3))
        real(dp) :: drive_node(size(fields, 1), size(fields, 2))
        real(dp) :: sqrtg_radial_node(size(fields, 1), size(fields, 2))
        real(dp) :: sqrtg_theta_node(size(fields, 1), size(fields, 2))
        real(dp) :: sqrtg_zeta_node(size(fields, 1), size(fields, 2))
        real(dp) :: gamma_node(size(fields, 1), size(fields, 2))
        real(dp) :: nodes(2), radial_coordinate, blend
        integer :: interval, node, node_count, left, right, term

        call validate_inputs(fields, drive, signed_sqrtg_radial, &
            signed_sqrtg_theta, signed_sqrtg_zeta, gamma_pressure_pa, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_step, phase_assembly, element_to_global, stiffness, info)
        if (info /= 0) return
        if (present(stiffness_terms)) then
            if (any(shape(stiffness_terms) /= [size(stiffness, 1), &
                size(stiffness, 2), compressible_stiffness_term_count])) then
                info = compressible_family_invalid
                return
            end if
            stiffness_terms = 0.0_dp
        end if
        if (radial_space%quadrature_points == 2) then
            node_count = 2
            nodes(1) = gauss_two_lower
            nodes(2) = gauss_two_upper
        else
            node_count = 1
            nodes(1) = radial_space%evaluation_coordinate
        end if
        stiffness = 0.0_dp
        do interval = 1, size(fields, 4)
            element_total = 0.0_dp
            if (present(stiffness_terms)) element_terms_total = 0.0_dp
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
                drive_node = (1.0_dp - blend) * drive(:, :, left) &
                    + blend * drive(:, :, right)
                sqrtg_radial_node = (1.0_dp - blend) &
                    * signed_sqrtg_radial(:, :, left) &
                    + blend * signed_sqrtg_radial(:, :, right)
                sqrtg_theta_node = (1.0_dp - blend) &
                    * signed_sqrtg_theta(:, :, left) &
                    + blend * signed_sqrtg_theta(:, :, right)
                sqrtg_zeta_node = (1.0_dp - blend) &
                    * signed_sqrtg_zeta(:, :, left) &
                    + blend * signed_sqrtg_zeta(:, :, right)
                gamma_node = (1.0_dp - blend) &
                    * gamma_pressure_pa(:, :, left) &
                    + blend * gamma_pressure_pa(:, :, right)
                if (present(stiffness_terms)) then
                    call assemble_compressible_stiffness_surface_resolved( &
                        fields_node, drive_node, sqrtg_radial_node, &
                        sqrtg_theta_node, sqrtg_zeta_node, gamma_node, &
                        trial_m, trial_n, trial_parity, stored_power, &
                        field_periods, node_space, radial_coordinate, &
                        radial_step, phase_assembly, element, info, &
                        element_terms)
                else
                    call assemble_compressible_stiffness_surface_resolved( &
                        fields_node, drive_node, sqrtg_radial_node, &
                        sqrtg_theta_node, sqrtg_zeta_node, gamma_node, &
                        trial_m, trial_n, trial_parity, stored_power, &
                        field_periods, node_space, radial_coordinate, &
                        radial_step, phase_assembly, element, info)
                end if
                if (info /= 0) return
                element_total = element_total + element
                if (present(stiffness_terms)) &
                    element_terms_total = element_terms_total + element_terms
            end do
            call add_mapped_dynamic_element(element_to_global(:, interval), &
                element_total, stiffness, info)
            if (info /= dynamic_layout_ok) return
            if (present(stiffness_terms)) then
                do term = 1, compressible_stiffness_term_count
                    call add_mapped_dynamic_element( &
                        element_to_global(:, interval), &
                        element_terms_total(:, :, term), &
                        stiffness_terms(:, :, term), info)
                    if (info /= dynamic_layout_ok) return
                end do
            end if
        end do
        info = 0
    end subroutine assemble_compressible_family_stiffness_fixed_layout

    subroutine validate_inputs(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_step, &
            phase_assembly, element_to_global, stiffness, info)
        real(dp), intent(in) :: fields(:, :, :, :), drive(:, :, :)
        real(dp), intent(in) :: jacobian_radial(:, :, :)
        real(dp), intent(in) :: jacobian_theta(:, :, :)
        real(dp), intent(in) :: jacobian_zeta(:, :, :), gamma_pressure(:, :, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        integer, intent(in) :: field_periods, phase_assembly
        integer, intent(in) :: element_to_global(:, :)
        real(dp), intent(in) :: stiffness(:, :)
        integer, intent(out) :: info
        integer :: angular_radial_shape(3), intervals, trials

        info = compressible_family_invalid
        trials = size(trial_m)
        intervals = size(fields, 4)
        angular_radial_shape = [size(fields, 1), size(fields, 2), intervals]
        if (size(trial_n) /= trials .or. size(trial_parity) /= trials) return
        if (size(stored_power) /= trials) return
        if (size(fields, 3) < 13) return
        if (any(shape(drive) /= angular_radial_shape)) return
        if (any(shape(jacobian_radial) /= angular_radial_shape)) return
        if (any(shape(jacobian_theta) /= angular_radial_shape)) return
        if (any(shape(jacobian_zeta) /= angular_radial_shape)) return
        if (any(shape(gamma_pressure) /= angular_radial_shape)) return
        if (.not. ieee_is_finite(radial_step)) return
        if (abs(radial_step * real(intervals, dp) - 1.0_dp) &
            > radial_tolerance) return
        if (field_periods < 1) return
        if (phase_assembly /= phase_assembly_transformed .and. &
            phase_assembly /= phase_assembly_direct) return
        if (size(stiffness, 1) /= size(stiffness, 2)) return
        if (.not. dynamic_element_map_is_valid(element_to_global, trials, &
            intervals, size(stiffness, 1))) return
        info = 0
    end subroutine validate_inputs

end module compressible_stiffness_family_assembly
