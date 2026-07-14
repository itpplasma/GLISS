module compressible_stiffness_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface, &
        compatible_stiffness_term_count
    use compressible_stiffness_validation, only: &
        validate_compressible_stiffness_inputs
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    implicit none
    private

    integer, parameter, public :: compressible_stiffness_term_count = &
        compatible_stiffness_term_count

    public :: assemble_compressible_stiffness_surface
    public :: assemble_compressible_stiffness_surface_resolved

contains

    subroutine assemble_compressible_stiffness_surface(fields, drive, &
            signed_sqrtg_radial, signed_sqrtg_theta, signed_sqrtg_zeta, &
            gamma_pressure_pa, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            phase_assembly, stiffness, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: signed_sqrtg_radial(:, :)
        real(dp), intent(in) :: signed_sqrtg_theta(:, :)
        real(dp), intent(in) :: signed_sqrtg_zeta(:, :)
        real(dp), intent(in) :: gamma_pressure_pa(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        integer :: dimension

        call validate_compressible_stiffness_inputs(fields, drive, &
            signed_sqrtg_radial, signed_sqrtg_theta, signed_sqrtg_zeta, &
            gamma_pressure_pa, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_step, phase_assembly, info)
        if (info /= 0) return
        dimension = 4 * size(trial_m)
        allocate (stiffness(dimension, dimension))
        call assemble_compressible_stiffness_surface_resolved(fields, drive, &
            signed_sqrtg_radial, signed_sqrtg_theta, signed_sqrtg_zeta, &
            gamma_pressure_pa, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            phase_assembly, stiffness, info)
    end subroutine assemble_compressible_stiffness_surface

    subroutine assemble_compressible_stiffness_surface_resolved(fields, &
            drive, signed_sqrtg_radial, signed_sqrtg_theta, &
            signed_sqrtg_zeta, gamma_pressure_pa, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, phase_assembly, stiffness, info, &
            stiffness_terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: signed_sqrtg_radial(:, :)
        real(dp), intent(in) :: signed_sqrtg_theta(:, :)
        real(dp), intent(in) :: signed_sqrtg_zeta(:, :)
        real(dp), intent(in) :: gamma_pressure_pa(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(out) :: stiffness_terms(:, :, :)
        real(dp) :: h1_derivatives(2, size(trial_m))
        real(dp) :: h1_values(2, size(trial_m))
        real(dp) :: l2_values(1, size(trial_m))

        call validate_compressible_stiffness_inputs(fields, drive, &
            signed_sqrtg_radial, signed_sqrtg_theta, signed_sqrtg_zeta, &
            gamma_pressure_pa, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_step, phase_assembly, info)
        if (info /= 0) return
        if (any(shape(stiffness) /= 4 * size(trial_m))) then
            info = -1
            return
        end if
        if (present(stiffness_terms)) then
            if (any(shape(stiffness_terms) /= [size(stiffness, 1), &
                size(stiffness, 2), compressible_stiffness_term_count])) then
                info = -1
                return
            end if
            stiffness_terms = 0.0_dp
        end if
        call build_legacy_factors(radial_space, trial_m, stored_power, &
            radial_coordinate, radial_step, h1_values, h1_derivatives, &
            l2_values, info)
        if (info /= 0) return
        stiffness = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, &
            drive, signed_sqrtg_radial, signed_sqrtg_theta, &
            signed_sqrtg_zeta, gamma_pressure_pa, trial_m, trial_n, &
            trial_parity, field_periods, h1_values, h1_derivatives, &
            l2_values, radial_step * radial_space%weight_fraction, &
            phase_assembly, stiffness, info, stiffness_terms)
    end subroutine assemble_compressible_stiffness_surface_resolved

    subroutine build_legacy_factors(radial_space, trial_m, stored_power, &
            radial_coordinate, radial_step, h1_values, h1_derivatives, &
            l2_values, info)
        type(radial_space_config_t), intent(in) :: radial_space
        integer, intent(in) :: trial_m(:)
        real(dp), intent(in) :: stored_power(:), radial_coordinate, radial_step
        real(dp), intent(out) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(out) :: l2_values(:, :)
        integer, intent(out) :: info
        integer :: trial

        info = -1
        do trial = 1, size(trial_m)
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, &
                radial_space%evaluation_coordinate, h1_values(:, trial), &
                h1_derivatives(:, trial), info, stored_power(trial))
            if (info /= radial_space_ok) return
        end do
        l2_values = 1.0_dp
        info = 0
    end subroutine build_legacy_factors

end module compressible_stiffness_assembly
