module compatible_family_point_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use period_averaged_assembly, only: accumulate_period_averaged, &
        period_masks
    use two_component_kernel, only: two_component_components
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
    integer, parameter, public :: compatible_two_component_term_count = 4

    public :: assemble_compatible_direct_surface
    public :: assemble_compatible_transformed_surface

contains

    subroutine assemble_compatible_direct_surface(fields, drive, trial_m, &
            trial_n, trial_parity, field_periods, h1_values, h1_derivatives, &
            l2_values, full, info, terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(in) :: l2_values(:, :)
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(inout) :: terms(:, :, :)
        real(dp) :: weight
        integer :: j, k, period

        call validate_inputs(fields, drive, trial_m, trial_n, trial_parity, &
            field_periods, h1_values, h1_derivatives, l2_values, full, info, &
            terms)
        if (info /= 0) return
        weight = 1.0_dp / real(size(fields, 1) * size(fields, 2) &
            * field_periods, dp)
        do period = 0, field_periods - 1
            do k = 1, size(fields, 2)
                do j = 1, size(fields, 1)
                    call accumulate_direct(fields(j, k, :), drive(j, k), &
                        trial_m, trial_n, trial_parity, field_periods, &
                        h1_values, h1_derivatives, l2_values, &
                        real(j - 1, dp) / real(size(fields, 1), dp), &
                        real(k - 1, dp) / real(size(fields, 2), dp) &
                        + real(period, dp), weight, full, terms)
                end do
            end do
        end do
        info = 0
    end subroutine assemble_compatible_direct_surface

    ! Period-averaged angular quadrature as matrix products over chunks of
    ! angular points (period_averaged_assembly): the three kernel
    ! components and the drive are four weighted channels.
    subroutine assemble_compatible_transformed_surface(fields, drive, &
            trial_m, trial_n, trial_parity, field_periods, h1_values, &
            h1_derivatives, l2_values, full, info, terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(in) :: l2_values(:, :)
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(inout) :: terms(:, :, :)
        integer, parameter :: chunk_limit = 256
        real(dp), allocatable :: cosine_part(:, :, :), sine_part(:, :, :)
        real(dp), allocatable :: cosine_phase(:, :), sine_phase(:, :)
        real(dp), allocatable :: weight(:, :), term(:, :)
        real(dp), allocatable :: plus(:, :), minus(:, :)
        real(dp), allocatable :: cosine_rows(:, :), sine_rows(:, :)
        real(dp) :: angular_weight, phase, theta, zeta
        logical :: mixed
        integer :: channel, chunk, column, columns, count, first, j, k
        integer :: point, points, trial, trials

        call validate_inputs(fields, drive, trial_m, trial_n, trial_parity, &
            field_periods, h1_values, h1_derivatives, l2_values, full, info, &
            terms)
        if (info /= 0) return
        columns = size(full, 1)
        trials = size(trial_m)
        points = size(fields, 1) * size(fields, 2)
        angular_weight = 1.0_dp / real(points, dp)
        chunk = min(chunk_limit, points)
        allocate (cosine_part(chunk, columns, 4), sine_part(chunk, columns, 4), &
            cosine_phase(chunk, columns), sine_phase(chunk, columns), &
            weight(chunk, 4), term(columns, columns), &
            plus(columns, columns), minus(columns, columns), &
            cosine_rows(4, columns), sine_rows(4, columns))
        call period_masks(trial_n, field_periods, columns, plus, minus, mixed)
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
                call phase_rows(fields(j, k, :), trial_m, trial_n, &
                    trial_parity, field_periods, h1_values, h1_derivatives, &
                    l2_values, cosine_rows, sine_rows)
                weight(point, 1:3) = angular_weight * abs(fields(j, k, 7))
                weight(point, 4) = -drive(j, k) * angular_weight &
                    * abs(fields(j, k, 7))
                theta = real(j - 1, dp) / real(size(fields, 1), dp)
                zeta = real(k - 1, dp) / real(size(fields, 2), dp)
                do column = 1, columns
                    trial = modulo(column - 1, trials) + 1
                    phase = two_pi * (real(trial_m(trial), dp) * theta &
                        - real(trial_n(trial), dp) * zeta &
                        / real(field_periods, dp))
                    cosine_phase(point, column) = cos(phase)
                    sine_phase(point, column) = sin(phase)
                    cosine_part(point, column, :) = cosine_rows(:, column)
                    sine_part(point, column, :) = sine_rows(:, column)
                end do
            end do
            do channel = 1, 4
                term = 0.0_dp
                call accumulate_period_averaged(cosine_part(:, :, channel), &
                    sine_part(:, :, channel), cosine_phase, sine_phase, &
                    weight(:, channel), plus, minus, mixed, term)
                full = full + term
                if (present(terms)) terms(:, :, channel) = &
                    terms(:, :, channel) + term
            end do
        end do
        info = 0
    end subroutine assemble_compatible_transformed_surface

    ! Kernel rows multiplying cos and sin of each trial phase.
    subroutine phase_rows(fields, trial_m, trial_n, trial_parity, &
            field_periods, h1_values, h1_derivatives, l2_values, cosine_rows, &
            sine_rows)
        real(dp), intent(in) :: fields(:), h1_values(:, :)
        real(dp), intent(in) :: h1_derivatives(:, :), l2_values(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods
        real(dp), intent(out) :: cosine_rows(:, :), sine_rows(:, :)
        real(dp) :: coefficients(3, 6), toroidal_wave
        integer :: trial

        call kernel_coefficients(fields, coefficients)
        cosine_rows = 0.0_dp
        sine_rows = 0.0_dp
        do trial = 1, size(trial_m)
            toroidal_wave = real(trial_n(trial), dp) &
                / real(field_periods, dp)
            if (trial_parity(trial) == 1) then
                call add_trial_columns(cosine_rows, trial, trial_m(trial), &
                    toroidal_wave, h1_values(:, trial), &
                    h1_derivatives(:, trial), l2_values(:, trial), &
                    1.0_dp, 0.0_dp, 1.0_dp, coefficients)
                call add_trial_columns(sine_rows, trial, trial_m(trial), &
                    toroidal_wave, h1_values(:, trial), &
                    h1_derivatives(:, trial), l2_values(:, trial), &
                    0.0_dp, -1.0_dp, 0.0_dp, coefficients)
            else
                call add_trial_columns(cosine_rows, trial, trial_m(trial), &
                    toroidal_wave, h1_values(:, trial), &
                    h1_derivatives(:, trial), l2_values(:, trial), &
                    0.0_dp, 1.0_dp, 0.0_dp, coefficients)
                call add_trial_columns(sine_rows, trial, trial_m(trial), &
                    toroidal_wave, h1_values(:, trial), &
                    h1_derivatives(:, trial), l2_values(:, trial), &
                    1.0_dp, 0.0_dp, -1.0_dp, coefficients)
            end if
        end do
    end subroutine phase_rows

    subroutine accumulate_direct(fields, drive, trial_m, trial_n, &
            trial_parity, field_periods, h1_values, h1_derivatives, &
            l2_values, theta, zeta, weight, full, terms)
        real(dp), intent(in) :: fields(:), drive, h1_values(:, :)
        real(dp), intent(in) :: h1_derivatives(:, :), l2_values(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: full(:, :)
        real(dp), optional, intent(inout) :: terms(:, :, :)
        real(dp) :: rows(4, size(full, 1)), coefficients(3, 6)
        real(dp) :: phase, toroidal_wave, value, dvalue, dother
        integer :: trial

        call kernel_coefficients(fields, coefficients)
        rows = 0.0_dp
        do trial = 1, size(trial_m)
            toroidal_wave = real(trial_n(trial), dp) &
                / real(field_periods, dp)
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - toroidal_wave * zeta)
            call phase_factors(phase, trial_parity(trial), value, dvalue, &
                dother)
            call add_trial_columns(rows, trial, trial_m(trial), &
                toroidal_wave, h1_values(:, trial), &
                h1_derivatives(:, trial), l2_values(:, trial), value, &
                dvalue, dother, coefficients)
        end do
        call rank_update(rows, drive, weight * abs(fields(7)), full, terms)
    end subroutine accumulate_direct

    subroutine add_trial_columns(rows, trial, m, toroidal_wave, h1_values, &
            h1_derivatives, l2_values, value, dvalue, dother, coefficients)
        real(dp), intent(inout) :: rows(:, :)
        integer, intent(in) :: trial, m
        real(dp), intent(in) :: toroidal_wave, h1_values(:)
        real(dp), intent(in) :: h1_derivatives(:), l2_values(:)
        real(dp), intent(in) :: value, dvalue, dother, coefficients(:, :)
        real(dp) :: inputs(6)
        integer :: basis, column, trials

        trials = size(rows, 2) / (size(h1_values) + size(l2_values))
        do basis = 1, size(h1_values)
            column = (basis - 1) * trials + trial
            inputs(1) = value * h1_values(basis)
            inputs(2) = value * h1_derivatives(basis)
            inputs(3) = two_pi * real(m, dp) * dvalue * h1_values(basis)
            inputs(4) = -two_pi * toroidal_wave * dvalue * h1_values(basis)
            inputs(5:6) = 0.0_dp
            rows(1:3, column) = matmul(coefficients, inputs)
            rows(4, column) = value * h1_values(basis)
        end do
        do basis = 1, size(l2_values)
            column = size(h1_values) * trials + (basis - 1) * trials + trial
            inputs(1:4) = 0.0_dp
            inputs(5) = two_pi * real(m, dp) * dother * l2_values(basis)
            inputs(6) = -two_pi * toroidal_wave * dother * l2_values(basis)
            rows(1:3, column) = matmul(coefficients, inputs)
        end do
    end subroutine add_trial_columns

    subroutine kernel_coefficients(fields, coefficients)
        real(dp), intent(in) :: fields(:)
        real(dp), intent(out) :: coefficients(3, 6)
        real(dp) :: inputs(6)
        integer :: entry

        do entry = 1, 6
            inputs = 0.0_dp
            inputs(entry) = 1.0_dp
            call two_component_components(fields(1), fields(2), fields(3), &
                fields(4), fields(5), fields(6), fields(7), fields(8), &
                fields(9), fields(10), fields(11), fields(12), fields(13), &
                inputs(1), inputs(2), inputs(3), inputs(4), inputs(5), &
                inputs(6), coefficients(1, entry), coefficients(2, entry), &
                coefficients(3, entry))
        end do
    end subroutine kernel_coefficients

    pure subroutine phase_factors(phase, parity, value, dvalue, dother)
        real(dp), intent(in) :: phase
        integer, intent(in) :: parity
        real(dp), intent(out) :: value, dvalue, dother

        if (parity == 1) then
            value = cos(phase)
            dvalue = -sin(phase)
            dother = cos(phase)
        else
            value = sin(phase)
            dvalue = cos(phase)
            dother = -sin(phase)
        end if
    end subroutine phase_factors

    subroutine rank_update(rows, drive, weight, full, terms)
        real(dp), intent(in) :: rows(:, :), drive, weight
        real(dp), intent(inout) :: full(:, :)
        real(dp), optional, intent(inout) :: terms(:, :, :)
        real(dp) :: contributions(4)
        integer :: a, b

        do b = 1, size(full, 2)
            do a = 1, size(full, 1)
                contributions(1) = rows(1, a) * rows(1, b)
                contributions(2) = rows(2, a) * rows(2, b)
                contributions(3) = rows(3, a) * rows(3, b)
                contributions(4) = -drive * rows(4, a) * rows(4, b)
                full(a, b) = full(a, b) + weight * sum(contributions)
                if (present(terms)) terms(a, b, :) = terms(a, b, :) &
                    + weight * contributions
            end do
        end do
    end subroutine rank_update

    subroutine validate_inputs(fields, drive, trial_m, trial_n, trial_parity, &
            field_periods, h1_values, h1_derivatives, l2_values, full, info, &
            terms)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(in) :: l2_values(:, :), full(:, :)
        integer, intent(out) :: info
        real(dp), optional, intent(in) :: terms(:, :, :)
        integer :: expected, trials

        info = -1
        trials = size(trial_m)
        if (trials < 1 .or. field_periods < 1) return
        if (size(trial_n) /= trials .or. size(trial_parity) /= trials) return
        if (any(trial_m < 0) .or. any(trial_parity < 1) &
            .or. any(trial_parity > 2)) return
        if (size(h1_values, 1) < 1 .or. size(l2_values, 1) < 1) return
        if (size(h1_values, 2) /= trials &
            .or. any(shape(h1_derivatives) /= shape(h1_values)) &
            .or. size(l2_values, 2) /= trials) return
        expected = trials * (size(h1_values, 1) + size(l2_values, 1))
        if (any(shape(full) /= expected)) return
        if (present(terms)) then
            if (size(terms, 1) /= expected .or. size(terms, 2) /= expected &
                .or. size(terms, 3) /= compatible_two_component_term_count) &
                return
        end if
        if (size(fields, 1) < 1 .or. size(fields, 2) < 1 &
            .or. size(fields, 3) < 13) return
        if (any(shape(drive) /= shape(fields(:, :, 1)))) return
        if (.not. all(ieee_is_finite(fields(:, :, 1:13))) &
            .or. .not. all(ieee_is_finite(drive)) &
            .or. .not. all(ieee_is_finite(h1_values)) &
            .or. .not. all(ieee_is_finite(h1_derivatives)) &
            .or. .not. all(ieee_is_finite(l2_values))) return
        info = 0
    end subroutine validate_inputs

end module compatible_family_point_assembly
