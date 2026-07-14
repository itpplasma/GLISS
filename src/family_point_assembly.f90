module family_point_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_family_point_assembly, only: &
        assemble_compatible_direct_surface, &
        assemble_compatible_transformed_surface
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    implicit none
    private

    public :: assemble_direct_surface
    public :: assemble_direct_surface_resolved
    public :: assemble_transformed_surface
    public :: assemble_transformed_surface_resolved
    public :: resolve_normal_stored_power

contains

    subroutine assemble_direct_surface(fields, drive, trial_m, trial_n, &
            trial_parity, field_periods, radial_space, radial_coordinate, &
            radial_step, full, info, normal_stored_power)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in), optional :: normal_stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: stored_power(:)
        call resolve_normal_stored_power(normal_stored_power, size(trial_m), &
            stored_power, info)
        if (info /= 0) return
        call assemble_direct_surface_resolved(fields, drive, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, full, info)
    end subroutine assemble_direct_surface

    subroutine assemble_direct_surface_resolved(fields, drive, trial_m, &
            trial_n, trial_parity, normal_stored_power, field_periods, &
            radial_space, radial_coordinate, radial_step, full, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp) :: h1_derivatives(2, size(trial_m))
        real(dp) :: h1_values(2, size(trial_m))
        real(dp) :: l2_values(1, size(trial_m))

        call validate_surface_inputs(fields, drive, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, &
            radial_coordinate, radial_step, full, info)
        if (info /= 0) return
        call build_compatibility_factors(radial_space, trial_m, &
            normal_stored_power, radial_coordinate, radial_step, h1_values, &
            h1_derivatives, l2_values, info)
        if (info /= 0) return
        call assemble_compatible_direct_surface(fields, drive, trial_m, &
            trial_n, trial_parity, field_periods, h1_values, h1_derivatives, &
            l2_values, full, info)
    end subroutine assemble_direct_surface_resolved

    subroutine assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            trial_parity, field_periods, radial_space, radial_coordinate, &
            radial_step, full, info, normal_stored_power)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in), optional :: normal_stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: stored_power(:)
        call resolve_normal_stored_power(normal_stored_power, size(trial_m), &
            stored_power, info)
        if (info /= 0) return
        call assemble_transformed_surface_resolved(fields, drive, trial_m, &
            trial_n, trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, full, info)
    end subroutine assemble_transformed_surface

    subroutine assemble_transformed_surface_resolved(fields, drive, trial_m, &
            trial_n, trial_parity, normal_stored_power, field_periods, &
            radial_space, radial_coordinate, radial_step, full, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp) :: h1_derivatives(2, size(trial_m))
        real(dp) :: h1_values(2, size(trial_m))
        real(dp) :: l2_values(1, size(trial_m))

        call validate_surface_inputs(fields, drive, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, &
            radial_coordinate, radial_step, full, info)
        if (info /= 0) return
        call build_compatibility_factors(radial_space, trial_m, &
            normal_stored_power, radial_coordinate, radial_step, h1_values, &
            h1_derivatives, l2_values, info)
        if (info /= 0) return
        call assemble_compatible_transformed_surface(fields, drive, trial_m, &
            trial_n, trial_parity, field_periods, h1_values, h1_derivatives, &
            l2_values, full, info)
    end subroutine assemble_transformed_surface_resolved

    subroutine resolve_normal_stored_power(input, trials, power, info)
        real(dp), intent(in), optional :: input(:)
        integer, intent(in) :: trials
        real(dp), allocatable, intent(out) :: power(:)
        integer, intent(out) :: info

        info = -1
        allocate (power(trials), source=0.0_dp)
        if (present(input)) then
            if (size(input) /= trials) return
            if (.not. all(ieee_is_finite(input))) return
            power = input
        end if
        info = 0
    end subroutine resolve_normal_stored_power

    subroutine validate_surface_inputs(fields, drive, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, &
            radial_coordinate, radial_step, full, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(in) :: full(:, :)
        integer, intent(out) :: info
        integer :: trials

        info = -1
        trials = size(trial_m)
        if (field_periods < 1 .or. trials < 1) return
        if (size(trial_n) /= trials .or. size(trial_parity) /= trials) return
        if (size(normal_stored_power) /= trials) return
        if (size(fields, 3) < 13) return
        if (size(drive, 1) /= size(fields, 1)) return
        if (size(drive, 2) /= size(fields, 2)) return
        if (size(full, 1) /= 3 * trials .or. &
            size(full, 2) /= 3 * trials) return
        if (any(trial_m < 0)) return
        if (any(trial_parity < 1) .or. any(trial_parity > 2)) return
        if (.not. all(ieee_is_finite(normal_stored_power))) return
        if (.not. ieee_is_finite(radial_coordinate)) return
        if (.not. ieee_is_finite(radial_step)) return
        if (radial_step <= 0.0_dp) return
        if (.not. all(ieee_is_finite(fields(:, :, 1:13)))) return
        if (.not. all(ieee_is_finite(drive))) return
        info = 0
    end subroutine validate_surface_inputs

    subroutine build_compatibility_factors(radial_space, trial_m, &
            normal_stored_power, radial_coordinate, radial_step, h1_values, &
            h1_derivatives, l2_values, info)
        type(radial_space_config_t), intent(in) :: radial_space
        integer, intent(in) :: trial_m(:)
        real(dp), intent(in) :: normal_stored_power(:)
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(out) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(out) :: l2_values(:, :)
        integer, intent(out) :: info
        integer :: trial

        info = -1
        if (size(h1_values, 1) /= 2 .or. size(l2_values, 1) /= 1) return
        if (size(h1_values, 2) /= size(trial_m)) return
        if (any(shape(h1_derivatives) /= shape(h1_values))) return
        if (size(l2_values, 2) /= size(trial_m)) return
        do trial = 1, size(trial_m)
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, 0.5_dp, &
                h1_values(:, trial), h1_derivatives(:, trial), info, &
                normal_stored_power(trial))
            if (info /= radial_space_ok) return
        end do
        l2_values = 1.0_dp
        info = 0
    end subroutine build_compatibility_factors

end module family_point_assembly
