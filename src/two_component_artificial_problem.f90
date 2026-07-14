module two_component_artificial_problem
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: add_mapped_dynamic_element, &
        build_dynamic_block_permutation, build_dynamic_element_map, &
        build_resolved_dynamic_family_layout, dynamic_family_layout_t, &
        dynamic_layout_ok
    use family_assembly, only: family_assembly_options_t, surface_geometry_t
    use family_point_assembly, only: assemble_direct_surface_resolved, &
        assemble_transformed_surface_resolved
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use trial_space_topology, only: build_trial_space_topology, &
        trial_space_topology_t, trial_topology_ok
    use variable_block_tridiagonal, only: pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: artificial_problem_ok = 0
    integer, parameter, public :: artificial_problem_invalid = -1
    integer, parameter, public :: artificial_problem_assembly_error = -2

    type, public :: two_component_artificial_problem_t
        type(variable_block_tridiagonal_t) :: stiffness
        type(variable_block_tridiagonal_t) :: mass
        integer :: unknowns = 0
        integer :: normal_unknowns = 0
        integer :: eta_unknowns = 0
    end type two_component_artificial_problem_t

    public :: build_two_component_artificial_problem

contains

    subroutine build_two_component_artificial_problem(geometry, mode_m, &
            mode_n, stored_power, radial_step, options, problem, info)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        type(family_assembly_options_t), intent(in) :: options
        type(two_component_artificial_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        type(dynamic_family_layout_t) :: layout
        type(trial_space_topology_t) :: topology
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        integer, allocatable :: mode_parity(:), element_map(:, :)
        integer, allocatable :: permutation(:), widths(:)

        info = artificial_problem_invalid
        if (.not. valid_inputs(geometry, mode_m, mode_n, stored_power, &
            radial_step, options)) return
        allocate (mode_parity(size(mode_m)), source=options%parity_class)
        call build_trial_space_topology(mode_m, mode_n, mode_parity, &
            topology, info)
        if (info /= trial_topology_ok) return
        topology%active(3, :) = .false.
        call build_resolved_dynamic_family_layout(topology, size(geometry), &
            layout, info)
        if (info /= dynamic_layout_ok) return
        call build_dynamic_element_map(layout, element_map, info)
        if (info /= dynamic_layout_ok) return
        allocate (stiffness(layout%total_unknowns, layout%total_unknowns), &
            source=0.0_dp)
        call assemble_stiffness(geometry, mode_m, mode_n, mode_parity, &
            stored_power, radial_step, options, element_map, stiffness, info)
        if (info /= artificial_problem_ok) return
        allocate (mass(layout%total_unknowns, layout%total_unknowns), &
            source=0.0_dp)
        call assemble_artificial_mass(radial_step, mass)
        call build_dynamic_block_permutation(layout, widths, permutation, info)
        if (info /= dynamic_layout_ok) then
            info = artificial_problem_assembly_error
            return
        end if
        call pack_permuted_variable_blocks(stiffness, permutation, widths, &
            problem%stiffness, info)
        if (info /= variable_block_ok) then
            info = artificial_problem_assembly_error
            return
        end if
        call pack_permuted_variable_blocks(mass, permutation, widths, &
            problem%mass, info)
        if (info /= variable_block_ok) then
            info = artificial_problem_assembly_error
            return
        end if
        problem%unknowns = layout%total_unknowns
        problem%normal_unknowns = layout%normal_unknowns
        problem%eta_unknowns = layout%eta_unknowns
        info = artificial_problem_ok
    end subroutine build_two_component_artificial_problem

    subroutine assemble_stiffness(geometry, mode_m, mode_n, mode_parity, &
            stored_power, radial_step, options, element_map, stiffness, info)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:), mode_parity(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        type(family_assembly_options_t), intent(in) :: options
        integer, intent(in) :: element_map(:, :)
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp) :: full(3 * size(mode_m), 3 * size(mode_m))
        real(dp) :: element(4 * size(mode_m), 4 * size(mode_m))
        real(dp) :: radial_coordinate
        integer :: interval

        do interval = 1, size(geometry)
            full = 0.0_dp
            radial_coordinate = (real(interval, dp) - 0.5_dp) * radial_step
            if (options%phase_assembly == phase_assembly_direct) then
                call assemble_direct_surface_resolved( &
                    geometry(interval)%fields, geometry(interval)%drive, &
                    mode_m, mode_n, mode_parity, stored_power, &
                    options%field_periods, options%radial_space, &
                    radial_coordinate, radial_step, full, info)
            else
                call assemble_transformed_surface_resolved( &
                    geometry(interval)%fields, geometry(interval)%drive, &
                    mode_m, mode_n, mode_parity, stored_power, &
                    options%field_periods, options%radial_space, &
                    radial_coordinate, radial_step, full, info)
            end if
            if (info /= 0) then
                info = artificial_problem_assembly_error
                return
            end if
            element = 0.0_dp
            element(1:size(full, 1), 1:size(full, 2)) = radial_step * full
            call add_mapped_dynamic_element(element_map(:, interval), &
                element, stiffness, info)
            if (info /= dynamic_layout_ok) then
                info = artificial_problem_assembly_error
                return
            end if
        end do
        info = artificial_problem_ok
    end subroutine assemble_stiffness

    pure subroutine assemble_artificial_mass(radial_step, mass)
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: mass(:, :)
        integer :: coefficient

        mass = 0.0_dp
        do coefficient = 1, size(mass, 1)
            mass(coefficient, coefficient) = radial_step
        end do
    end subroutine assemble_artificial_mass

    pure function valid_inputs(geometry, mode_m, mode_n, stored_power, &
            radial_step, options) result(valid)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        type(family_assembly_options_t), intent(in) :: options
        logical :: valid

        valid = size(geometry) >= 2 .and. size(mode_m) >= 1
        if (.not. valid) return
        valid = size(mode_n) == size(mode_m) &
            .and. size(stored_power) == size(mode_m)
        if (.not. valid) return
        valid = radial_step > 0.0_dp
        if (.not. valid) return
        valid = options%field_periods >= 1
        if (.not. valid) return
        valid = options%parity_class >= 1 .and. options%parity_class <= 2
        if (.not. valid) return
        valid = options%phase_assembly == phase_assembly_direct &
            .or. options%phase_assembly == phase_assembly_transformed
    end function valid_inputs

end module two_component_artificial_problem
