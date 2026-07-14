module physical_mass_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use mass_density_policy, only: evaluate_mass_density, mass_density_ok, &
        mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    implicit none
    private

    public :: assemble_physical_mass_surface
    public :: assemble_physical_mass_surface_resolved

contains

    subroutine assemble_physical_mass_surface(fields, density_profile, &
            trial_m, trial_n, trial_parity, normal_stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            phase_assembly, mass, info)
        real(dp), intent(in) :: fields(:, :, :)
        type(mass_density_profile_t), intent(in) :: density_profile
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: density_kg_m3

        call evaluate_mass_density(density_profile, radial_coordinate, &
            density_kg_m3, info)
        if (info /= mass_density_ok) then
            info = -1
            return
        end if
        allocate (mass(4 * size(trial_m), 4 * size(trial_m)))
        call assemble_physical_mass_surface_resolved(fields, density_kg_m3, &
            trial_m, trial_n, trial_parity, normal_stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            phase_assembly, mass, info)
    end subroutine assemble_physical_mass_surface

    subroutine assemble_physical_mass_surface_resolved(fields, &
            density_kg_m3, trial_m, trial_n, trial_parity, &
            normal_stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, phase_assembly, mass, info)
        real(dp), intent(in) :: fields(:, :, :), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: h1_values(2, size(trial_m))
        real(dp) :: l2_values(1, size(trial_m)), derivatives(2)
        integer :: trial

        call validate_inputs(fields, density_kg_m3, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, radial_step, &
            phase_assembly, mass, info)
        if (info /= 0) return
        do trial = 1, size(trial_m)
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, &
                radial_space%evaluation_coordinate, h1_values(:, trial), &
                derivatives, info, normal_stored_power(trial))
            if (info /= radial_space_ok) return
        end do
        l2_values = 1.0_dp
        mass = 0.0_dp
        call assemble_compatible_physical_mass_surface(fields, density_kg_m3, &
            trial_m, trial_n, trial_parity, field_periods, h1_values, &
            l2_values, radial_step * radial_space%weight_fraction, &
            phase_assembly, mass, info)
    end subroutine assemble_physical_mass_surface_resolved

    subroutine validate_inputs(fields, density, trial_m, trial_n, parity, &
            stored_power, field_periods, radial_step, phase_assembly, mass, info)
        real(dp), intent(in) :: fields(:, :, :), density, stored_power(:)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:)
        integer, intent(in) :: field_periods, phase_assembly
        real(dp), intent(in) :: radial_step, mass(:, :)
        integer, intent(out) :: info
        integer :: trials

        info = -1
        trials = size(trial_m)
        if (trials < 1 .or. size(trial_n) /= trials) return
        if (size(parity) /= trials .or. size(stored_power) /= trials) return
        if (any(trial_m < 0) .or. any(parity < 1) .or. any(parity > 2)) return
        if (.not. all(ieee_is_finite(stored_power))) return
        if (.not. ieee_is_finite(density) .or. density <= 0.0_dp) return
        if (field_periods < 1) return
        if (phase_assembly /= phase_assembly_transformed .and. &
            phase_assembly /= phase_assembly_direct) return
        if (any(shape(mass) /= 4 * trials)) return
        if (size(fields, 1) < 1 .or. size(fields, 2) < 1 &
            .or. size(fields, 3) < 13) return
        if (.not. all(ieee_is_finite(fields(:, :, 1:13)))) return
        if (.not. ieee_is_finite(radial_step) .or. radial_step <= 0.0_dp) &
            return
        if (any(fields(:, :, 7) == 0.0_dp) &
            .or. any(fields(:, :, 8) <= 0.0_dp) &
            .or. any(fields(:, :, 9) <= 0.0_dp)) return
        if (any(fields(:, :, 1)**2 + fields(:, :, 2)**2 <= 0.0_dp)) return
        info = 0
    end subroutine validate_inputs

end module physical_mass_assembly
