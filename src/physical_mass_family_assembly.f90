module physical_mass_family_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: build_dynamic_family_layout, &
        dynamic_family_layout_t, dynamic_layout_ok, eta_global_index, &
        mu_global_index, normal_global_index
    use mass_density_policy, only: evaluate_mass_density, mass_density_ok, &
        mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use physical_mass_assembly, only: &
        assemble_physical_mass_surface_resolved
    use radial_space_policy, only: radial_space_config_t
    implicit none
    private

    real(dp), parameter :: radial_tolerance = 1.0e-12_dp

    public :: assemble_physical_family_mass
    public :: assemble_physical_family_mass_resolved

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
        real(dp), allocatable :: density_kg_m3(:)
        real(dp) :: radial_coordinate
        integer :: interval

        call build_dynamic_family_layout(size(trial_m), size(fields, 4), &
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
        real(dp) :: element(4 * size(trial_m), 4 * size(trial_m))
        real(dp) :: radial_coordinate
        integer :: interval

        call validate_resolved_inputs(fields, density_kg_m3, trial_m, &
            trial_n, trial_parity, stored_power, field_periods, radial_step, &
            phase_assembly, mass, layout, info)
        if (info /= 0) return
        mass = 0.0_dp
        do interval = 1, layout%intervals
            radial_coordinate = (real(interval, dp) - 0.5_dp) * radial_step
            call assemble_physical_mass_surface_resolved( &
                fields(:, :, :, interval), density_kg_m3(interval), trial_m, &
                trial_n, trial_parity, stored_power, field_periods, &
                radial_space, radial_coordinate, radial_step, phase_assembly, &
                element, info)
            if (info /= 0) return
            call add_element(layout, interval, element, mass)
        end do
        info = 0
    end subroutine assemble_physical_family_mass_resolved

    subroutine add_element(layout, interval, element, mass)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: interval
        real(dp), intent(in) :: element(:, :)
        real(dp), intent(inout) :: mass(:, :)
        integer :: local_a, local_b, global_a, global_b

        do local_b = 1, size(element, 2)
            global_b = element_global_index(layout, interval, local_b)
            if (global_b == 0) cycle
            do local_a = 1, size(element, 1)
                global_a = element_global_index(layout, interval, local_a)
                if (global_a == 0) cycle
                mass(global_a, global_b) = mass(global_a, global_b) &
                    + element(local_a, local_b)
            end do
        end do
    end subroutine add_element

    pure function element_global_index(layout, interval, local) result(index)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: interval, local
        integer :: index, block, trial

        block = (local - 1) / layout%trials + 1
            trial = modulo(local - 1, layout%trials) + 1
            select case (block)
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

        subroutine validate_resolved_inputs(fields, density_kg_m3, trial_m, &
                trial_n, trial_parity, stored_power, field_periods, radial_step, &
                phase_assembly, mass, layout, info)
            real(dp), intent(in) :: fields(:, :, :, :), density_kg_m3(:)
            integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
            real(dp), intent(in) :: stored_power(:), radial_step
            integer, intent(in) :: field_periods, phase_assembly
            real(dp), intent(in) :: mass(:, :)
            type(dynamic_family_layout_t), intent(out) :: layout
            integer, intent(out) :: info
            integer :: trials, intervals

            info = -1
            trials = size(trial_m)
            intervals = size(fields, 4)
            call build_dynamic_family_layout(trials, intervals, layout, info)
            if (info /= dynamic_layout_ok) then
                info = -1
                return
            end if
            info = -1
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
            if (size(mass, 1) /= layout%total_unknowns) return
            if (size(mass, 2) /= layout%total_unknowns) return
            info = 0
        end subroutine validate_resolved_inputs

    end module physical_mass_family_assembly
