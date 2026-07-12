module compressible_stiffness_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use phase_factor_topology, only: phase_cosine, &
        phase_product_coefficients, phase_sine
    use physical_constants, only: vacuum_permeability
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    use three_component_kernel, only: compressible_divergence_value
    use two_component_kernel, only: bending_component_value, &
        compression_component_value, shear_component_value
    implicit none
    private

    integer, parameter :: response_count = 5
    integer, parameter :: xi_value = 1, xi_radial = 2
    integer, parameter :: xi_theta = 3, xi_zeta = 4, eta_value = 5
    integer, parameter :: eta_theta = 6, eta_zeta = 7
    integer, parameter :: mu_theta = 8, mu_zeta = 9
    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

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

        call validate_inputs(fields, drive, signed_sqrtg_radial, &
            signed_sqrtg_theta, signed_sqrtg_zeta, gamma_pressure_pa, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_step, phase_assembly, info)
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
            radial_coordinate, radial_step, phase_assembly, stiffness, info)
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

        call validate_inputs(fields, drive, signed_sqrtg_radial, &
            signed_sqrtg_theta, signed_sqrtg_zeta, gamma_pressure_pa, &
            trial_m, trial_n, trial_parity, stored_power, field_periods, &
            radial_step, phase_assembly, info)
        if (info /= 0) return
        if (any(shape(stiffness) /= 4 * size(trial_m))) then
            info = -1
            return
        end if
        stiffness = 0.0_dp
        if (phase_assembly == phase_assembly_direct) then
            call assemble_direct(fields, drive, signed_sqrtg_radial, &
                signed_sqrtg_theta, signed_sqrtg_zeta, gamma_pressure_pa, &
                trial_m, trial_n, trial_parity, stored_power, field_periods, &
                radial_space, radial_coordinate, radial_step, stiffness, info)
        else
            call assemble_transformed(fields, drive, signed_sqrtg_radial, &
                signed_sqrtg_theta, signed_sqrtg_zeta, gamma_pressure_pa, &
                trial_m, trial_n, trial_parity, stored_power, field_periods, &
                radial_space, radial_coordinate, radial_step, stiffness, info)
        end if
    end subroutine assemble_compressible_stiffness_surface_resolved

    subroutine assemble_direct(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, stiffness, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), gamma_pressure(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp) :: theta, zeta, weight
        integer :: j, k, period

        weight = radial_step * radial_space%weight_fraction &
            / real(size(fields, 1) * size(fields, 2) &
            * field_periods, dp)
        do period = 0, field_periods - 1
            do k = 1, size(fields, 2)
                zeta = real(k - 1, dp) / real(size(fields, 2), dp) &
                    + real(period, dp)
                do j = 1, size(fields, 1)
                    theta = real(j - 1, dp) / real(size(fields, 1), dp)
                    call accumulate_direct_point(fields(j, k, :), &
                        drive(j, k), jacobian_radial(j, k), &
                        jacobian_theta(j, k), jacobian_zeta(j, k), &
                        gamma_pressure(j, k), trial_m, trial_n, trial_parity, &
                        stored_power, field_periods, radial_space, &
                        radial_coordinate, radial_step, theta, zeta, weight, &
                        stiffness, info)
                    if (info /= 0) return
                end do
            end do
        end do
        info = 0
    end subroutine assemble_direct

    subroutine accumulate_direct_point(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, theta, zeta, weight, stiffness, &
            info)
        real(dp), intent(in) :: fields(:), drive, jacobian_radial
        real(dp), intent(in) :: jacobian_theta, jacobian_zeta, gamma_pressure
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step, theta, zeta
        real(dp), intent(in) :: weight
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp) :: responses(response_count, 2, size(stiffness, 1))
        real(dp) :: values(response_count, size(stiffness, 1)), phase
        integer :: degree, trial, trials

        call build_point_responses(fields, jacobian_radial, jacobian_theta, &
            jacobian_zeta, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            responses, info)
        if (info /= 0) return
        trials = size(trial_m)
        do degree = 1, size(stiffness, 1)
            trial = modulo(degree - 1, trials) + 1
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            values(:, degree) = responses(:, phase_cosine, degree) * cos(phase) &
                + responses(:, phase_sine, degree) * sin(phase)
        end do
        call accumulate_response_values(values, drive, gamma_pressure, &
            fields(7), weight, stiffness)
        info = 0
    end subroutine accumulate_direct_point

    subroutine assemble_transformed(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, stiffness, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), gamma_pressure(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp) :: theta, zeta, weight
        integer :: j, k

        weight = radial_step * radial_space%weight_fraction &
            / real(size(fields, 1) * size(fields, 2), dp)
        do k = 1, size(fields, 2)
            zeta = real(k - 1, dp) / real(size(fields, 2), dp)
            do j = 1, size(fields, 1)
                theta = real(j - 1, dp) / real(size(fields, 1), dp)
                call accumulate_transformed_point(fields(j, k, :), &
                    drive(j, k), jacobian_radial(j, k), &
                    jacobian_theta(j, k), jacobian_zeta(j, k), &
                    gamma_pressure(j, k), trial_m, trial_n, trial_parity, &
                    stored_power, field_periods, radial_space, &
                    radial_coordinate, radial_step, theta, zeta, weight, &
                    stiffness, info)
                if (info /= 0) return
            end do
        end do
        info = 0
    end subroutine assemble_transformed

    subroutine accumulate_transformed_point(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, theta, zeta, weight, stiffness, &
            info)
        real(dp), intent(in) :: fields(:), drive, jacobian_radial
        real(dp), intent(in) :: jacobian_theta, jacobian_zeta, gamma_pressure
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step, theta, zeta
        real(dp), intent(in) :: weight
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp) :: responses(response_count, 2, size(stiffness, 1))
        real(dp) :: cosine(size(trial_m)), sine(size(trial_m)), phase
        integer :: trial

        call build_point_responses(fields, jacobian_radial, jacobian_theta, &
            jacobian_zeta, trial_m, trial_n, trial_parity, stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            responses, info)
        if (info /= 0) return
        do trial = 1, size(trial_m)
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            cosine(trial) = cos(phase)
            sine(trial) = sin(phase)
        end do
        call accumulate_response_coefficients(responses, cosine, sine, &
            trial_n, field_periods, drive, gamma_pressure, fields(7), weight, &
            stiffness)
        info = 0
    end subroutine accumulate_transformed_point

    subroutine build_point_responses(fields, jacobian_radial, &
            jacobian_theta, jacobian_zeta, trial_m, trial_n, trial_parity, &
            stored_power, field_periods, radial_space, radial_coordinate, &
            radial_step, responses, info)
        real(dp), intent(in) :: fields(:), jacobian_radial, jacobian_theta
        real(dp), intent(in) :: jacobian_zeta, stored_power(:)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(out) :: responses(:, :, :)
        integer, intent(out) :: info
        real(dp) :: basis(9, 2), normal_values(2), normal_derivatives(2)
        real(dp) :: point_responses(response_count, 2)
        integer :: block, degree, trial, trials

        trials = size(trial_m)
        responses = 0.0_dp
        do trial = 1, trials
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, &
                radial_space%evaluation_coordinate, normal_values, &
                normal_derivatives, info, stored_power(trial))
            if (info /= radial_space_ok) return
            do block = 1, 4
                degree = (block - 1) * trials + trial
                call build_basis_coefficients(block, trial_parity(trial), &
                    trial_m(trial), trial_n(trial), field_periods, &
                    normal_values, normal_derivatives, basis)
                call build_energy_responses(fields, jacobian_radial, &
                    jacobian_theta, jacobian_zeta, basis, point_responses)
                responses(:, :, degree) = point_responses
            end do
        end do
        info = 0
    end subroutine build_point_responses

    pure subroutine build_basis_coefficients(block, parity, mode_m, mode_n, &
            field_periods, normal_values, normal_derivatives, basis)
        integer, intent(in) :: block, parity, mode_m, mode_n, field_periods
        real(dp), intent(in) :: normal_values(2), normal_derivatives(2)
        real(dp), intent(out) :: basis(9, 2)
        real(dp) :: xi_phase(2), tangential_phase(2)
        integer :: normal_index

        xi_phase = 0.0_dp
        tangential_phase = 0.0_dp
        xi_phase(parity) = 1.0_dp
        if (parity == phase_cosine) then
            tangential_phase(phase_sine) = 1.0_dp
        else
            tangential_phase(phase_cosine) = 1.0_dp
        end if
        basis = 0.0_dp
        if (block <= 2) then
            normal_index = block
            basis(xi_value, :) = normal_values(normal_index) * xi_phase
            basis(xi_radial, :) = normal_derivatives(normal_index) * xi_phase
            basis(xi_theta, :) = angular_derivative(xi_phase, &
                two_pi * real(mode_m, dp)) * normal_values(normal_index)
            basis(xi_zeta, :) = angular_derivative(xi_phase, &
                -two_pi * real(mode_n, dp) / real(field_periods, dp)) &
                * normal_values(normal_index)
        else if (block == 3) then
            basis(eta_value, :) = tangential_phase
            basis(eta_theta, :) = angular_derivative(tangential_phase, &
                two_pi * real(mode_m, dp))
            basis(eta_zeta, :) = angular_derivative(tangential_phase, &
                -two_pi * real(mode_n, dp) / real(field_periods, dp))
        else
            basis(mu_theta, :) = angular_derivative(tangential_phase, &
                two_pi * real(mode_m, dp))
            basis(mu_zeta, :) = angular_derivative(tangential_phase, &
                -two_pi * real(mode_n, dp) / real(field_periods, dp))
        end if
    end subroutine build_basis_coefficients

    pure function angular_derivative(coefficients, wavenumber) result(values)
        real(dp), intent(in) :: coefficients(2), wavenumber
        real(dp) :: values(2)

        values(phase_cosine) = wavenumber * coefficients(phase_sine)
        values(phase_sine) = -wavenumber * coefficients(phase_cosine)
    end function angular_derivative

    pure subroutine build_energy_responses(fields, jacobian_radial, &
            jacobian_theta, jacobian_zeta, basis, responses)
        real(dp), intent(in) :: fields(:), jacobian_radial, jacobian_theta
        real(dp), intent(in) :: jacobian_zeta, basis(9, 2)
        real(dp), intent(out) :: responses(response_count, 2)
        real(dp) :: sqrtg_xi_radial, sqrtg_eta_theta, sqrtg_eta_zeta
        integer :: kind

        do kind = phase_cosine, phase_sine
            sqrtg_xi_radial = jacobian_radial * basis(xi_value, kind) &
                + fields(7) * basis(xi_radial, kind)
            sqrtg_eta_theta = jacobian_theta * basis(eta_value, kind) &
                + fields(7) * basis(eta_theta, kind)
            sqrtg_eta_zeta = jacobian_zeta * basis(eta_value, kind) &
                + fields(7) * basis(eta_zeta, kind)
            responses(1, kind) = bending_component_value(fields(1), &
                fields(2), fields(7), fields(9), basis(xi_theta, kind), &
                basis(xi_zeta, kind))
            responses(2, kind) = shear_component_value(fields(1), fields(2), &
                fields(3), fields(4), fields(7), fields(8), fields(9), &
                fields(10), fields(12), basis(xi_value, kind), &
                basis(xi_theta, kind), basis(xi_zeta, kind), &
                basis(eta_theta, kind), basis(eta_zeta, kind))
            responses(3, kind) = compression_component_value(fields(1), &
                fields(2), fields(3), fields(4), fields(5), fields(6), &
                fields(7), fields(8), fields(11), fields(13), &
                basis(xi_value, kind), basis(xi_radial, kind), &
                basis(xi_theta, kind), basis(xi_zeta, kind), &
                basis(eta_theta, kind), basis(eta_zeta, kind))
            responses(4, kind) = basis(xi_value, kind)
            responses(5, kind) = compressible_divergence_value(fields(1), &
                fields(2), fields(7), sqrtg_xi_radial, sqrtg_eta_theta, &
                sqrtg_eta_zeta, basis(mu_theta, kind), basis(mu_zeta, kind))
        end do
    end subroutine build_energy_responses

    subroutine accumulate_response_values(responses, drive, gamma_pressure, &
            signed_sqrtg, weight, stiffness)
        real(dp), intent(in) :: responses(:, :), drive, gamma_pressure
        real(dp), intent(in) :: signed_sqrtg, weight
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp) :: factors(response_count)
        integer :: a, b

        factors = response_factors(drive, gamma_pressure, signed_sqrtg)
        do b = 1, size(stiffness, 2)
            do a = 1, size(stiffness, 1)
                stiffness(a, b) = stiffness(a, b) + weight &
                    * sum(factors * responses(:, a) * responses(:, b))
            end do
        end do
    end subroutine accumulate_response_values

    subroutine accumulate_response_coefficients(responses, cosine, sine, &
            trial_n, field_periods, drive, gamma_pressure, signed_sqrtg, &
            weight, stiffness)
        real(dp), intent(in) :: responses(:, :, :), cosine(:), sine(:)
        integer, intent(in) :: trial_n(:), field_periods
        real(dp), intent(in) :: drive, gamma_pressure, signed_sqrtg, weight
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp) :: factors(response_count), products(2, 2), contribution
        integer :: a, b, kind_a, kind_b, trial_a, trial_b, trials

        trials = size(trial_n)
        factors = response_factors(drive, gamma_pressure, signed_sqrtg)
        do b = 1, size(stiffness, 2)
            trial_b = modulo(b - 1, trials) + 1
            do a = 1, size(stiffness, 1)
                trial_a = modulo(a - 1, trials) + 1
                call phase_product_coefficients(cosine(trial_a), &
                    sine(trial_a), cosine(trial_b), sine(trial_b), &
                    trial_n(trial_a), trial_n(trial_b), field_periods, products)
                contribution = 0.0_dp
                do kind_b = phase_cosine, phase_sine
                    do kind_a = phase_cosine, phase_sine
                        contribution = contribution + products(kind_a, kind_b) &
                            * sum(factors * responses(:, kind_a, a) &
                            * responses(:, kind_b, b))
                    end do
                end do
                stiffness(a, b) = stiffness(a, b) + weight * contribution
            end do
        end do
    end subroutine accumulate_response_coefficients

    pure function response_factors(drive, gamma_pressure, signed_sqrtg) &
            result(factors)
        real(dp), intent(in) :: drive, gamma_pressure, signed_sqrtg
        real(dp) :: factors(response_count)

        factors(1:3) = abs(signed_sqrtg) / vacuum_permeability
        factors(4) = -drive * abs(signed_sqrtg) / vacuum_permeability
        factors(5) = gamma_pressure * abs(signed_sqrtg)
    end function response_factors

    subroutine validate_inputs(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_step, &
            phase_assembly, info)
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
    end subroutine validate_inputs

end module compressible_stiffness_assembly
