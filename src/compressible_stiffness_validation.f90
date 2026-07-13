module compressible_stiffness_validation
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    implicit none
    private

    public :: validate_compressible_stiffness_inputs

contains

    subroutine validate_compressible_stiffness_inputs(fields, drive, &
            jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_step, phase_assembly, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), gamma_pressure(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        integer, intent(in) :: field_periods, phase_assembly
        integer, intent(out) :: info
        integer :: angular_shape(2), trials

        info = -1
        trials = size(trial_m)
        angular_shape = [size(fields, 1), size(fields, 2)]
        if (trials < 1 .or. size(trial_n) /= trials) return
        if (size(trial_parity) /= trials .or. size(stored_power) /= trials) &
            return
        if (any(trial_m < 0)) return
        if (any(trial_parity < 1) .or. any(trial_parity > 2)) return
        if (.not. all(ieee_is_finite(stored_power))) return
        if (field_periods < 1) return
        if (phase_assembly /= phase_assembly_transformed .and. &
            phase_assembly /= phase_assembly_direct) return
        if (size(fields, 1) < 1 .or. size(fields, 2) < 1) return
        if (size(fields, 3) < 13) return
        if (any(shape(drive) /= angular_shape)) return
        if (any(shape(jacobian_radial) /= angular_shape)) return
        if (any(shape(jacobian_theta) /= angular_shape)) return
        if (any(shape(jacobian_zeta) /= angular_shape)) return
        if (any(shape(gamma_pressure) /= angular_shape)) return
        if (.not. all(ieee_is_finite(fields(:, :, 1:13)))) return
        if (.not. all(ieee_is_finite(drive))) return
        if (.not. all(ieee_is_finite(jacobian_radial))) return
        if (.not. all(ieee_is_finite(jacobian_theta))) return
        if (.not. all(ieee_is_finite(jacobian_zeta))) return
        if (.not. all(ieee_is_finite(gamma_pressure))) return
        if (any(gamma_pressure < 0.0_dp)) return
        if (.not. ieee_is_finite(radial_step) .or. radial_step <= 0.0_dp) &
            return
        if (any(fields(:, :, 7) == 0.0_dp)) return
        if (any(fields(:, :, 8) <= 0.0_dp)) return
        if (any(fields(:, :, 9) <= 0.0_dp)) return
        if (any(fields(:, :, 1)**2 + fields(:, :, 2)**2 <= 0.0_dp)) return
        info = 0
    end subroutine validate_compressible_stiffness_inputs

end module compressible_stiffness_validation
