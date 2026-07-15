module compatible_compressible_stiffness_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use phase_factor_topology, only: phase_cosine, &
        phase_product_coefficients, phase_sine
    use physical_constants, only: vacuum_permeability
    use three_component_kernel, only: compressible_divergence_value
    use two_component_kernel, only: bending_component_value, &
        compression_component_value, shear_component_value
    implicit none
    private

    integer, parameter, public :: compatible_stiffness_term_count = 5
    integer, parameter :: xi_value = 1, xi_radial = 2
    integer, parameter :: xi_theta = 3, xi_zeta = 4, eta_value = 5
    integer, parameter :: eta_theta = 6, eta_zeta = 7
    integer, parameter :: mu_theta = 8, mu_zeta = 9
    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: assemble_compatible_compressible_stiffness_surface

contains

    subroutine assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_radial, jacobian_theta, jacobian_zeta, &
            gamma_pressure, trial_m, trial_n, trial_parity, field_periods, &
            h1_values, h1_derivatives, l2_values, radial_weight, &
            phase_assembly, stiffness, info, stiffness_terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), gamma_pressure(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods, phase_assembly
        real(dp), intent(in) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(in) :: l2_values(:, :), radial_weight
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(inout) :: stiffness_terms(:, :, :)
        real(dp) :: angular_weight
        integer :: j, k, period

        call validate_inputs(fields, drive, jacobian_radial, jacobian_theta, &
            jacobian_zeta, gamma_pressure, trial_m, trial_n, trial_parity, &
            field_periods, h1_values, h1_derivatives, l2_values, &
            radial_weight, phase_assembly, stiffness, info, stiffness_terms)
        if (info /= 0) return
        if (phase_assembly == phase_assembly_direct) then
            angular_weight = radial_weight / real(size(fields, 1) &
                * size(fields, 2) * field_periods, dp)
            do period = 0, field_periods - 1
                do k = 1, size(fields, 2)
                    do j = 1, size(fields, 1)
                        call accumulate_direct(fields(j, k, :), drive(j, k), &
                            jacobian_radial(j, k), jacobian_theta(j, k), &
                            jacobian_zeta(j, k), gamma_pressure(j, k), &
                            trial_m, trial_n, trial_parity, field_periods, &
                            h1_values, h1_derivatives, l2_values, &
                            real(j - 1, dp) / real(size(fields, 1), dp), &
                            real(k - 1, dp) / real(size(fields, 2), dp) &
                            + real(period, dp), angular_weight, stiffness, &
                            stiffness_terms)
                    end do
                end do
            end do
        else
            angular_weight = radial_weight / real(size(fields, 1) &
                * size(fields, 2), dp)
            do k = 1, size(fields, 2)
                do j = 1, size(fields, 1)
                    call accumulate_transformed(fields(j, k, :), drive(j, k), &
                        jacobian_radial(j, k), jacobian_theta(j, k), &
                        jacobian_zeta(j, k), gamma_pressure(j, k), trial_m, &
                        trial_n, trial_parity, field_periods, h1_values, &
                        h1_derivatives, l2_values, &
                        real(j - 1, dp) / real(size(fields, 1), dp), &
                        real(k - 1, dp) / real(size(fields, 2), dp), &
                        angular_weight, stiffness, stiffness_terms)
                end do
            end do
        end if
        info = 0
    end subroutine assemble_compatible_compressible_stiffness_surface

    subroutine accumulate_direct(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            parity, field_periods, h1, dh1, l2, theta, zeta, weight, &
            stiffness, stiffness_terms)
        real(dp), intent(in) :: fields(:), drive, jacobian_radial
        real(dp), intent(in) :: jacobian_theta, jacobian_zeta, gamma_pressure
        real(dp), intent(in) :: h1(:, :), dh1(:, :), l2(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp), optional, intent(inout) :: stiffness_terms(:, :, :)
        real(dp) :: coefficients(5, 2, size(stiffness, 1))
        real(dp) :: responses(5, size(stiffness, 1)), factors(5)
        real(dp) :: phase, cosine, sine
        integer :: column, trial, trials

        call build_response_coefficients(fields, jacobian_radial, &
            jacobian_theta, jacobian_zeta, trial_m, trial_n, parity, &
            field_periods, h1, dh1, l2, coefficients)
        trials = size(trial_m)
        do column = 1, size(stiffness, 1)
            trial = modulo(column - 1, trials) + 1
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            cosine = cos(phase)
            sine = sin(phase)
            responses(:, column) = coefficients(:, phase_cosine, column) &
                * cosine + coefficients(:, phase_sine, column) * sine
        end do
        call build_response_factors(drive, gamma_pressure, fields(7), factors)
        call rank_update(responses, factors, weight, stiffness, stiffness_terms)
    end subroutine accumulate_direct

    subroutine accumulate_transformed(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            parity, field_periods, h1, dh1, l2, theta, zeta, weight, &
            stiffness, stiffness_terms)
        real(dp), intent(in) :: fields(:), drive, jacobian_radial
        real(dp), intent(in) :: jacobian_theta, jacobian_zeta, gamma_pressure
        real(dp), intent(in) :: h1(:, :), dh1(:, :), l2(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp), optional, intent(inout) :: stiffness_terms(:, :, :)
        real(dp) :: coefficients(5, 2, size(stiffness, 1))
        real(dp) :: cosine(size(trial_m)), sine(size(trial_m)), factors(5), phase
        integer :: trial

        call build_response_coefficients(fields, jacobian_radial, &
            jacobian_theta, jacobian_zeta, trial_m, trial_n, parity, &
            field_periods, h1, dh1, l2, coefficients)
        do trial = 1, size(trial_m)
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            cosine(trial) = cos(phase)
            sine(trial) = sin(phase)
        end do
        call build_response_factors(drive, gamma_pressure, fields(7), factors)
        call transformed_rank_update(coefficients, cosine, sine, trial_n, &
            field_periods, factors, weight, stiffness, stiffness_terms)
    end subroutine accumulate_transformed

    pure subroutine build_response_coefficients(fields, jacobian_radial, &
            jacobian_theta, jacobian_zeta, trial_m, trial_n, parity, &
            field_periods, h1, dh1, l2, responses)
        real(dp), intent(in) :: fields(:), jacobian_radial, jacobian_theta
        real(dp), intent(in) :: jacobian_zeta, h1(:, :), dh1(:, :), l2(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), contiguous, intent(out) :: responses(:, :, :)
        real(dp) :: basis(9, 2)
        real(dp) :: phase_coefficients(2)
        integer :: basis_index, column, trial, trials

        responses = 0.0_dp
        trials = size(trial_m)
        do basis_index = 1, size(h1, 1)
            do trial = 1, trials
                basis = 0.0_dp
                basis(xi_value, parity(trial)) = h1(basis_index, trial)
                basis(xi_radial, parity(trial)) = dh1(basis_index, trial)
                phase_coefficients = basis(xi_value, :)
                call angular_derivative(phase_coefficients(phase_cosine), &
                    phase_coefficients(phase_sine), &
                    two_pi * real(trial_m(trial), dp), &
                    basis(xi_theta, phase_cosine), &
                    basis(xi_theta, phase_sine))
                call angular_derivative(phase_coefficients(phase_cosine), &
                    phase_coefficients(phase_sine), &
                    -two_pi * real(trial_n(trial), dp) &
                    / real(field_periods, dp), basis(xi_zeta, phase_cosine), &
                    basis(xi_zeta, phase_sine))
                column = (basis_index - 1) * trials + trial
                call build_energy_responses(fields, jacobian_radial, &
                    jacobian_theta, jacobian_zeta, basis, &
                    responses(:, :, column))
            end do
        end do
        do basis_index = 1, size(l2, 1)
            do trial = 1, trials
                call build_tangential_basis(trial_m(trial), trial_n(trial), &
                    parity(trial), field_periods, l2(basis_index, trial), &
                    .true., basis)
                column = size(h1, 1) * trials &
                    + (basis_index - 1) * trials + trial
                call build_energy_responses(fields, jacobian_radial, &
                    jacobian_theta, jacobian_zeta, basis, &
                    responses(:, :, column))
                call build_tangential_basis(trial_m(trial), trial_n(trial), &
                    parity(trial), field_periods, l2(basis_index, trial), &
                    .false., basis)
                column = (size(h1, 1) + size(l2, 1)) * trials &
                    + (basis_index - 1) * trials + trial
                call build_energy_responses(fields, jacobian_radial, &
                    jacobian_theta, jacobian_zeta, basis, &
                    responses(:, :, column))
            end do
        end do
    end subroutine build_response_coefficients

    pure subroutine build_tangential_basis(mode_m, mode_n, parity, &
            field_periods, value, is_eta, basis)
        integer, intent(in) :: mode_m, mode_n, parity, field_periods
        real(dp), intent(in) :: value
        logical, intent(in) :: is_eta
        real(dp), intent(out) :: basis(9, 2)
        real(dp) :: phase_coefficients(2)
        integer :: kind

        if (parity == phase_cosine) then
            kind = phase_sine
        else
            kind = phase_cosine
        end if
        basis = 0.0_dp
        phase_coefficients = 0.0_dp
        phase_coefficients(kind) = value
        if (is_eta) then
            basis(eta_value, :) = phase_coefficients
            call angular_derivative(phase_coefficients(phase_cosine), &
                phase_coefficients(phase_sine), two_pi * real(mode_m, dp), &
                basis(eta_theta, phase_cosine), &
                basis(eta_theta, phase_sine))
            call angular_derivative(phase_coefficients(phase_cosine), &
                phase_coefficients(phase_sine), -two_pi * real(mode_n, dp) &
                / real(field_periods, dp), basis(eta_zeta, phase_cosine), &
                basis(eta_zeta, phase_sine))
        else
            call angular_derivative(phase_coefficients(phase_cosine), &
                phase_coefficients(phase_sine), two_pi * real(mode_m, dp), &
                basis(mu_theta, phase_cosine), basis(mu_theta, phase_sine))
            call angular_derivative(phase_coefficients(phase_cosine), &
                phase_coefficients(phase_sine), -two_pi * real(mode_n, dp) &
                / real(field_periods, dp), basis(mu_zeta, phase_cosine), &
                basis(mu_zeta, phase_sine))
        end if
    end subroutine build_tangential_basis

    pure subroutine angular_derivative(cosine, sine, wavenumber, &
            derivative_cosine, derivative_sine)
        real(dp), intent(in) :: cosine, sine, wavenumber
        real(dp), intent(out) :: derivative_cosine, derivative_sine

        derivative_cosine = wavenumber * sine
        derivative_sine = -wavenumber * cosine
    end subroutine angular_derivative

    pure subroutine build_energy_responses(fields, jacobian_radial, &
            jacobian_theta, jacobian_zeta, basis, responses)
        real(dp), intent(in) :: fields(:), jacobian_radial, jacobian_theta
        real(dp), intent(in) :: jacobian_zeta, basis(9, 2)
        real(dp), intent(out) :: responses(5, 2)
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

    subroutine rank_update(responses, factors, weight, stiffness, terms)
        real(dp), intent(in) :: responses(:, :), factors(:), weight
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp), optional, intent(inout) :: terms(:, :, :)
        real(dp) :: contributions(5)
        integer :: a, b

        do b = 1, size(stiffness, 2)
            do a = 1, size(stiffness, 1)
                contributions = factors * responses(:, a) * responses(:, b)
                stiffness(a, b) = stiffness(a, b) &
                    + weight * sum(contributions)
                if (present(terms)) &
                    terms(a, b, :) = terms(a, b, :) + weight * contributions
            end do
        end do
    end subroutine rank_update

    subroutine transformed_rank_update(responses, cosine, sine, trial_n, &
            field_periods, factors, weight, stiffness, terms)
        real(dp), intent(in) :: responses(:, :, :), cosine(:), sine(:)
        integer, intent(in) :: trial_n(:), field_periods
        real(dp), intent(in) :: factors(:), weight
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp), optional, intent(inout) :: terms(:, :, :)
        real(dp) :: products(2, 2), contributions(5)
        integer :: a, b, kind_a, kind_b, trial_a, trial_b, trials

        trials = size(trial_n)
        do b = 1, size(stiffness, 2)
            trial_b = modulo(b - 1, trials) + 1
            do a = 1, size(stiffness, 1)
                trial_a = modulo(a - 1, trials) + 1
                call phase_product_coefficients(cosine(trial_a), &
                    sine(trial_a), cosine(trial_b), sine(trial_b), &
                    trial_n(trial_a), trial_n(trial_b), field_periods, products)
                contributions = 0.0_dp
                do kind_b = phase_cosine, phase_sine
                    do kind_a = phase_cosine, phase_sine
                        contributions = contributions &
                            + products(kind_a, kind_b) * factors &
                            * responses(:, kind_a, a) * responses(:, kind_b, b)
                    end do
                end do
                stiffness(a, b) = stiffness(a, b) &
                    + weight * sum(contributions)
                if (present(terms)) &
                    terms(a, b, :) = terms(a, b, :) + weight * contributions
            end do
        end do
    end subroutine transformed_rank_update

    pure subroutine build_response_factors(drive, gamma_pressure, signed_sqrtg, &
            factors)
        real(dp), intent(in) :: drive, gamma_pressure, signed_sqrtg
        real(dp), intent(out) :: factors(5)

        factors(1:3) = abs(signed_sqrtg) / vacuum_permeability
        factors(4) = -drive * abs(signed_sqrtg) / vacuum_permeability
        factors(5) = gamma_pressure * abs(signed_sqrtg)
    end subroutine build_response_factors

    subroutine validate_inputs(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            parity, field_periods, h1, dh1, l2, radial_weight, phase_assembly, &
            stiffness, info, terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), gamma_pressure(:, :)
        real(dp), intent(in) :: h1(:, :), dh1(:, :), l2(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:)
        integer, intent(in) :: field_periods, phase_assembly
        real(dp), intent(in) :: radial_weight, stiffness(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(in) :: terms(:, :, :)
        integer :: expected, trials

        info = -1
        trials = size(trial_m)
        if (trials < 1 .or. field_periods < 1) return
        if (size(trial_n) /= trials .or. size(parity) /= trials) return
        if (any(trial_m < 0) .or. any(parity < 1) .or. any(parity > 2)) return
        if (size(h1, 1) < 1 .or. size(l2, 1) < 1) return
        if (size(h1, 2) /= trials .or. any(shape(dh1) /= shape(h1)) &
            .or. size(l2, 2) /= trials) return
        expected = trials * (size(h1, 1) + 2 * size(l2, 1))
        if (any(shape(stiffness) /= expected)) return
        if (present(terms)) then
            if (size(terms, 1) /= expected .or. size(terms, 2) /= expected &
                .or. size(terms, 3) /= compatible_stiffness_term_count) return
        end if
        if (size(fields, 1) < 1 .or. size(fields, 2) < 1 &
            .or. size(fields, 3) < 13) return
        if (.not. has_angular_shape(drive, fields) &
            .or. .not. has_angular_shape(jacobian_radial, fields) &
            .or. .not. has_angular_shape(jacobian_theta, fields) &
            .or. .not. has_angular_shape(jacobian_zeta, fields) &
            .or. .not. has_angular_shape(gamma_pressure, fields)) return
        if (.not. all(ieee_is_finite(fields(:, :, 1:13))) &
            .or. .not. all(ieee_is_finite(drive)) &
            .or. .not. all(ieee_is_finite(jacobian_radial)) &
            .or. .not. all(ieee_is_finite(jacobian_theta)) &
            .or. .not. all(ieee_is_finite(jacobian_zeta)) &
            .or. .not. all(ieee_is_finite(gamma_pressure)) &
            .or. .not. all(ieee_is_finite(h1)) &
            .or. .not. all(ieee_is_finite(dh1)) &
            .or. .not. all(ieee_is_finite(l2))) return
        if (any(gamma_pressure < 0.0_dp)) return
        if (.not. ieee_is_finite(radial_weight) .or. radial_weight <= 0.0_dp) &
            return
        if (phase_assembly /= phase_assembly_direct .and. &
            phase_assembly /= phase_assembly_transformed) return
        if (any(fields(:, :, 7) == 0.0_dp) &
            .or. any(fields(:, :, 8) <= 0.0_dp) &
            .or. any(fields(:, :, 9) <= 0.0_dp)) return
        if (any(fields(:, :, 1)**2 + fields(:, :, 2)**2 <= 0.0_dp)) return
        info = 0
    end subroutine validate_inputs

    pure logical function has_angular_shape(values, fields) result(valid)
        real(dp), intent(in) :: values(:, :), fields(:, :, :)

        valid = size(values, 1) == size(fields, 1) &
            .and. size(values, 2) == size(fields, 2)
    end function has_angular_shape

end module compatible_compressible_stiffness_assembly
