module compatible_compressible_stiffness_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use phase_factor_topology, only: phase_cosine, phase_sine
    use period_averaged_assembly, only: accumulate_period_averaged, &
        accumulate_period_averaged_tangent, period_masks
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
    public :: assemble_compatible_pressure_stiffness_tangent

contains

    subroutine assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_radial, jacobian_theta, jacobian_zeta, &
            gamma_pressure, trial_m, trial_n, trial_parity, field_periods, &
            h1_values, h1_derivatives, eta_values, l2_values, radial_weight, &
            phase_assembly, stiffness, info, stiffness_terms, sine_stiffness, &
            sine_terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), gamma_pressure(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods, phase_assembly
        real(dp), intent(in) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(in) :: eta_values(:, :), l2_values(:, :)
        real(dp), intent(in) :: radial_weight
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(inout) :: stiffness_terms(:, :, :)
        ! Cosine-parity trials only: also the stiffness and terms of the same
        ! modes with sine parity, from the same angular products.
        real(dp), optional, intent(inout) :: sine_stiffness(:, :)
        real(dp), optional, intent(inout) :: sine_terms(:, :, :)
        real(dp) :: angular_weight
        integer :: j, k, period

        call validate_inputs(fields, drive, jacobian_radial, jacobian_theta, &
            jacobian_zeta, gamma_pressure, trial_m, trial_n, trial_parity, &
            field_periods, h1_values, h1_derivatives, l2_values, &
            radial_weight, phase_assembly, stiffness, info, stiffness_terms)
        if (info /= 0) return
        info = -1
        if (any(shape(eta_values) /= shape(l2_values))) return
        if (.not. all(ieee_is_finite(eta_values))) return
        if (present(sine_stiffness) .neqv. present(sine_terms)) return
        if (present(sine_stiffness)) then
            if (phase_assembly /= phase_assembly_transformed) return
            if (.not. present(stiffness_terms)) return
            if (any(trial_parity /= phase_cosine)) return
            if (any(shape(sine_stiffness) /= shape(stiffness))) return
            if (any(shape(sine_terms) /= shape(stiffness_terms))) return
        end if
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
                            h1_values, h1_derivatives, eta_values, &
                            l2_values, &
                            real(j - 1, dp) / real(size(fields, 1), dp), &
                            real(k - 1, dp) / real(size(fields, 2), dp) &
                            + real(period, dp), angular_weight, stiffness, &
                            stiffness_terms)
                    end do
                end do
            end do
        else
            call assemble_transformed(fields, drive, jacobian_radial, &
                jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, &
                trial_n, trial_parity, field_periods, h1_values, &
                h1_derivatives, eta_values, l2_values, radial_weight, &
                stiffness, stiffness_terms, sine_stiffness, sine_terms)
        end if
        info = 0
    end subroutine assemble_compatible_compressible_stiffness_surface

    ! Contract v1: only pressure fields 10, 11 and 13, drive and gamma*p
    ! vary. Geometry, phases, basis and quadrature stay fixed. Differentiate
    ! each response square with the full bilinear product rule.
    subroutine assemble_compatible_pressure_stiffness_tangent(fields, drive, &
            jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure, &
            fields_tangent, drive_tangent, gamma_tangent, trial_m, trial_n, &
            parity, field_periods, h1, dh1, eta, l2, radial_weight, &
            stiffness, terms, info)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :), gamma_pressure(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), fields_tangent(:, :, :)
        real(dp), intent(in) :: drive_tangent(:, :), gamma_tangent
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: h1(:, :), dh1(:, :), eta(:, :), l2(:, :)
        real(dp), intent(in) :: radial_weight
        real(dp), intent(inout) :: stiffness(:, :), terms(:, :, :)
        integer, intent(out) :: info
        integer, parameter :: chunk_limit = 256
        real(dp), allocatable :: cosine(:, :, :), sine(:, :, :)
        real(dp), allocatable :: dc(:, :, :), ds(:, :, :), pc(:, :), ps(:, :)
        real(dp), allocatable :: weight(:, :), dweight(:, :), term(:, :)
        real(dp), allocatable :: plus(:, :), minus(:, :)
        real(dp) :: response(5, 2, size(stiffness, 1))
        real(dp) :: tangent(5, 2, size(stiffness, 1)), factors(5), angle
        real(dp) :: angular_weight
        logical :: mixed
        integer :: columns, trials, points, chunk, first, count, point
        integer :: j, k, trial, column, component

        call validate_inputs(fields, drive, jacobian_radial, jacobian_theta, &
            jacobian_zeta, gamma_pressure, trial_m, trial_n, parity, field_periods, &
            h1, dh1, l2, radial_weight, phase_assembly_transformed, stiffness, &
            info, terms)
        if (info /= 0) return
        info = -1
        if (any(shape(fields_tangent) /= shape(fields))) return
        if (any(shape(drive_tangent) /= shape(drive))) return
        if (any(shape(eta) /= shape(l2))) return
        if (.not. all(ieee_is_finite(eta))) return
        if (.not. all(ieee_is_finite(fields_tangent))) return
        if (.not. all(ieee_is_finite(drive_tangent))) return
        if (.not. ieee_is_finite(gamma_tangent)) return
        if (any(fields_tangent(:, :, :9) /= 0.0_dp)) return
        if (any(fields_tangent(:, :, 12) /= 0.0_dp)) return
        columns = size(stiffness, 1)
        trials = size(trial_m)
        points = size(fields, 1) * size(fields, 2)
        chunk = min(chunk_limit, points)
        angular_weight = radial_weight / real(points, dp)
        allocate (cosine(chunk, columns, 5), sine(chunk, columns, 5), &
            dc(chunk, columns, 5), ds(chunk, columns, 5), pc(chunk, columns), &
            ps(chunk, columns), weight(chunk, 5), dweight(chunk, 5), &
            term(columns, columns), plus(trials, trials), minus(trials, trials))
        call period_masks(trial_n, field_periods, plus, minus, mixed)
        do first = 1, points, chunk
            count = min(chunk, points - first + 1)
            cosine = 0.0_dp
            sine = 0.0_dp
            dc = 0.0_dp
            ds = 0.0_dp
            pc = 0.0_dp
            ps = 0.0_dp
            weight = 0.0_dp
            dweight = 0.0_dp
            do point = 1, count
                j = modulo(first + point - 2, size(fields, 1)) + 1
                k = (first + point - 2) / size(fields, 1) + 1
                call build_response_coefficients(fields(j, k, :), &
                    jacobian_radial(j, k), jacobian_theta(j, k), jacobian_zeta(j, k), &
                    trial_m, trial_n, parity, field_periods, h1, dh1, eta, l2, response)
                call build_response_coefficients(fields(j, k, :), &
                    jacobian_radial(j, k), jacobian_theta(j, k), jacobian_zeta(j, k), &
                    trial_m, trial_n, parity, field_periods, h1, dh1, eta, l2, &
                    tangent, &
                    fields_tangent(j, k, :))
                call build_response_factors(drive(j, k), gamma_pressure(j, k), &
                    fields(j, k, 7), factors)
                weight(point, :) = angular_weight * factors
                call build_response_factors(drive_tangent(j, k), gamma_tangent, &
                    fields(j, k, 7), factors)
                dweight(point, 4:5) = angular_weight * factors(4:5)
                do column = 1, columns
                    trial = modulo(column - 1, trials) + 1
                    angle = two_pi * (real(trial_m(trial), dp) &
                        * real(j - 1, dp) / real(size(fields, 1), dp) &
                        - real(trial_n(trial), dp) * real(k - 1, dp) &
                        / real(size(fields, 2) * field_periods, dp))
                    pc(point, column) = cos(angle)
                    ps(point, column) = sin(angle)
                    do component = 1, 5
                        cosine(point, column, component) = &
                            response(component, 1, column)
                        sine(point, column, component) = response(component, 2, column)
                        dc(point, column, component) = tangent(component, 1, column)
                        ds(point, column, component) = tangent(component, 2, column)
                    end do
                end do
            end do
            do component = 2, 5
                term = 0.0_dp
                if (component <= 3) then
                    call accumulate_period_averaged_tangent(cosine(:, :, component), &
                        sine(:, :, component), dc(:, :, component), &
                        ds(:, :, component), &
                        pc, ps, weight(:, component), plus, minus, mixed, term)
                else
                    call accumulate_period_averaged(cosine(:, :, component), &
                        sine(:, :, component), pc, ps, dweight(:, component), &
                        plus, minus, mixed, term)
                end if
                stiffness = stiffness + term
                terms(:, :, component) = terms(:, :, component) + term
            end do
        end do
        info = 0
    end subroutine assemble_compatible_pressure_stiffness_tangent

    ! Period-averaged angular quadrature as matrix products over chunks of
    ! angular points; each energy term is one weighted channel.
    subroutine assemble_transformed(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            parity, field_periods, h1, dh1, eta, l2, radial_weight, &
            stiffness, stiffness_terms, sine_stiffness, sine_terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        real(dp), intent(in) :: jacobian_radial(:, :), jacobian_theta(:, :)
        real(dp), intent(in) :: jacobian_zeta(:, :), gamma_pressure(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: h1(:, :), dh1(:, :), eta(:, :), l2(:, :)
        real(dp), intent(in) :: radial_weight
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp), optional, intent(inout) :: stiffness_terms(:, :, :)
        real(dp), optional, intent(inout) :: sine_stiffness(:, :)
        real(dp), optional, intent(inout) :: sine_terms(:, :, :)
        integer, parameter :: chunk_limit = 256
        real(dp), allocatable :: cosine_part(:, :, :), sine_part(:, :, :)
        real(dp), allocatable :: cosine_phase(:, :), sine_phase(:, :)
        real(dp), allocatable :: weight(:, :), term(:, :), sine_term(:, :)
        real(dp) :: coefficients(5, 2, size(stiffness, 1)), factors(5)
        real(dp) :: angular_weight, phase, theta, zeta
        real(dp) :: trial_cosine(size(trial_m)), trial_sine(size(trial_m))
        real(dp) :: rotation(size(stiffness, 1))
        real(dp), allocatable :: plus(:, :), minus(:, :)
        logical :: mixed
        integer :: chunk, column, columns, component, count, first, j, k
        integer :: point, points, trial, trials

        columns = size(stiffness, 1)
        trials = size(trial_m)
        points = size(fields, 1) * size(fields, 2)
        angular_weight = radial_weight / real(points, dp)
        chunk = min(chunk_limit, points)
        allocate (cosine_part(chunk, columns, 5), sine_part(chunk, columns, 5), &
            cosine_phase(chunk, columns), sine_phase(chunk, columns), &
            weight(chunk, 5), term(columns, columns), &
            sine_term(columns, columns), plus(trials, trials), &
            minus(trials, trials))
        call period_masks(trial_n, field_periods, plus, minus, mixed)
        ! The sine-parity basis turns the normal component a quarter period
        ! forward and the tangential ones a quarter period back.
        rotation = -1.0_dp
        rotation(:size(h1, 1) * trials) = 1.0_dp
        do first = 1, points, chunk
            count = min(chunk, points - first + 1)
            ! Rows past the last point carry zero weight and zero response.
            cosine_part(count + 1:, :, :) = 0.0_dp
            sine_part(count + 1:, :, :) = 0.0_dp
            cosine_phase(count + 1:, :) = 0.0_dp
            sine_phase(count + 1:, :) = 0.0_dp
            weight(count + 1:, :) = 0.0_dp
            do point = 1, count
                j = modulo(first + point - 2, size(fields, 1)) + 1
                k = (first + point - 2) / size(fields, 1) + 1
                call build_response_coefficients(fields(j, k, :), &
                    jacobian_radial(j, k), jacobian_theta(j, k), &
                    jacobian_zeta(j, k), trial_m, trial_n, parity, &
                    field_periods, h1, dh1, eta, l2, coefficients)
                call build_response_factors(drive(j, k), gamma_pressure(j, k), &
                    fields(j, k, 7), factors)
                weight(point, :) = angular_weight * factors
                theta = real(j - 1, dp) / real(size(fields, 1), dp)
                zeta = real(k - 1, dp) / real(size(fields, 2), dp)
                do trial = 1, trials
                    phase = two_pi * (real(trial_m(trial), dp) * theta &
                        - real(trial_n(trial), dp) * zeta &
                        / real(field_periods, dp))
                    trial_cosine(trial) = cos(phase)
                    trial_sine(trial) = sin(phase)
                end do
                do column = 1, columns
                    trial = modulo(column - 1, trials) + 1
                    cosine_phase(point, column) = trial_cosine(trial)
                    sine_phase(point, column) = trial_sine(trial)
                    do component = 1, 5
                        cosine_part(point, column, component) = &
                            coefficients(component, phase_cosine, column)
                        sine_part(point, column, component) = &
                            coefficients(component, phase_sine, column)
                    end do
                end do
            end do
            do component = 1, 5
                term = 0.0_dp
                if (present(sine_stiffness)) then
                    sine_term = 0.0_dp
                    call accumulate_period_averaged( &
                        cosine_part(:, :, component), &
                        sine_part(:, :, component), cosine_phase, &
                        sine_phase, weight(:, component), plus, minus, &
                        mixed, term, sine_term, rotation)
                    sine_stiffness = sine_stiffness + sine_term
                    sine_terms(:, :, component) = sine_terms(:, :, component) &
                        + sine_term
                else
                    call accumulate_period_averaged( &
                        cosine_part(:, :, component), &
                        sine_part(:, :, component), cosine_phase, &
                        sine_phase, weight(:, component), plus, minus, &
                        mixed, term)
                end if
                stiffness = stiffness + term
                if (present(stiffness_terms)) stiffness_terms(:, :, component) &
                    = stiffness_terms(:, :, component) + term
            end do
        end do
    end subroutine assemble_transformed

    subroutine accumulate_direct(fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, trial_m, trial_n, &
            parity, field_periods, h1, dh1, eta, l2, theta, zeta, weight, &
            stiffness, stiffness_terms)
        real(dp), intent(in) :: fields(:), drive, jacobian_radial
        real(dp), intent(in) :: jacobian_theta, jacobian_zeta, gamma_pressure
        real(dp), intent(in) :: h1(:, :), dh1(:, :), eta(:, :), l2(:, :)
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
            field_periods, h1, dh1, eta, l2, coefficients)
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

    ! The responses are linear in the radial basis values, so each trial
    ! evaluates the kernel once per unit channel (xi^s value, its radial
    ! derivative, eta and mu) and every basis function scales them.
    pure subroutine build_response_coefficients(fields, jacobian_radial, &
            jacobian_theta, jacobian_zeta, trial_m, trial_n, parity, &
            field_periods, h1, dh1, eta, l2, responses, pressure_tangent)
        real(dp), intent(in) :: fields(:), jacobian_radial, jacobian_theta
        real(dp), intent(in) :: jacobian_zeta, h1(:, :), dh1(:, :)
        real(dp), intent(in) :: eta(:, :), l2(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), contiguous, intent(out) :: responses(:, :, :)
        real(dp), optional, intent(in) :: pressure_tangent(:)
        real(dp) :: basis(9, 2), phase_coefficients(2)
        real(dp) :: unit_value(5, 2), unit_radial(5, 2), unit_eta(5, 2)
        real(dp) :: unit_mu(5, 2)
        integer :: basis_index, column, trial, trials

        trials = size(trial_m)
        do trial = 1, trials
            basis = 0.0_dp
            basis(xi_value, parity(trial)) = 1.0_dp
            phase_coefficients = basis(xi_value, :)
            call angular_derivative(phase_coefficients(phase_cosine), &
                phase_coefficients(phase_sine), &
                two_pi * real(trial_m(trial), dp), &
                basis(xi_theta, phase_cosine), basis(xi_theta, phase_sine))
            call angular_derivative(phase_coefficients(phase_cosine), &
                phase_coefficients(phase_sine), &
                -two_pi * real(trial_n(trial), dp) &
                / real(field_periods, dp), basis(xi_zeta, phase_cosine), &
                basis(xi_zeta, phase_sine))
            call build_energy_responses(fields, jacobian_radial, &
                jacobian_theta, jacobian_zeta, basis, unit_value, pressure_tangent)
            basis = 0.0_dp
            basis(xi_radial, parity(trial)) = 1.0_dp
            call build_energy_responses(fields, jacobian_radial, &
                jacobian_theta, jacobian_zeta, basis, unit_radial, pressure_tangent)
            call build_tangential_basis(trial_m(trial), trial_n(trial), &
                parity(trial), field_periods, 1.0_dp, .true., basis)
            call build_energy_responses(fields, jacobian_radial, &
                jacobian_theta, jacobian_zeta, basis, unit_eta, pressure_tangent)
            call build_tangential_basis(trial_m(trial), trial_n(trial), &
                parity(trial), field_periods, 1.0_dp, .false., basis)
            call build_energy_responses(fields, jacobian_radial, &
                jacobian_theta, jacobian_zeta, basis, unit_mu, pressure_tangent)
            do basis_index = 1, size(h1, 1)
                column = (basis_index - 1) * trials + trial
                responses(:, :, column) = h1(basis_index, trial) * unit_value &
                    + dh1(basis_index, trial) * unit_radial
            end do
            do basis_index = 1, size(l2, 1)
                column = size(h1, 1) * trials + (basis_index - 1) * trials &
                    + trial
                responses(:, :, column) = eta(basis_index, trial) * unit_eta
                column = (size(h1, 1) + size(l2, 1)) * trials &
                    + (basis_index - 1) * trials + trial
                responses(:, :, column) = l2(basis_index, trial) * unit_mu
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
            jacobian_theta, jacobian_zeta, basis, responses, pressure_tangent)
        real(dp), intent(in) :: fields(:), jacobian_radial, jacobian_theta
        real(dp), intent(in) :: jacobian_zeta, basis(9, 2)
        real(dp), intent(out) :: responses(5, 2)
        real(dp), optional, intent(in) :: pressure_tangent(:)
        real(dp) :: sqrtg_xi_radial, sqrtg_eta_theta, sqrtg_eta_zeta
        integer :: kind

        if (present(pressure_tangent)) then
            responses = 0.0_dp
            do kind = phase_cosine, phase_sine
                responses(2, kind) = -pressure_tangent(10) * basis(xi_value, kind) &
                    / (fields(8) * sqrt(fields(9)))
                responses(3, kind) = (-pressure_tangent(11) * basis(xi_value, kind) &
                    + pressure_tangent(13) * (fields(1) * basis(xi_theta, kind) &
                    + fields(2) * basis(xi_zeta, kind)) / fields(7)) / fields(8)
            end do
            return
        end if

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
            ! The third unknown is nu=mu-(FP'/FT') sqrtg eta, proportional
            ! to sqrtg xi^zeta and regular at the axis; mu itself inherits
            ! the singular s^(-1/2) part of eta for |m|=1.
            responses(5, kind) = compressible_divergence_value(fields(1), &
                fields(2), fields(7), sqrtg_xi_radial, sqrtg_eta_theta, &
                sqrtg_eta_zeta, basis(mu_theta, kind) &
                + fields(2) / fields(1) * sqrtg_eta_theta, &
                basis(mu_zeta, kind) + fields(2) / fields(1) * sqrtg_eta_zeta)
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
