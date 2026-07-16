module compatible_family_point_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_factor_topology, only: phase_product_coefficients
    use two_component_kernel, only: two_component_components
    !$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
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
        real(dp) :: weight
        integer :: first_column, j, k, last_column, thread, threads

        call validate_inputs(fields, drive, trial_m, trial_n, trial_parity, &
            field_periods, h1_values, h1_derivatives, l2_values, full, info, &
            terms)
        if (info /= 0) return
        weight = 1.0_dp / real(size(fields, 1) * size(fields, 2), dp)
        !$omp parallel default(none) &
        !$omp shared(fields, drive, trial_m, trial_n, trial_parity, field_periods) &
        !$omp shared(h1_values, h1_derivatives, l2_values, weight, full, terms) &
        !$omp private(first_column, last_column, thread, threads, j, k)
        thread = 0
        threads = 1
        !$      thread = omp_get_thread_num()
        !$      threads = omp_get_num_threads()
        first_column = 1 + thread * size(full, 2) / threads
        last_column = (thread + 1) * size(full, 2) / threads
        do k = 1, size(fields, 2)
            do j = 1, size(fields, 1)
                call accumulate_transformed(fields(j, k, :), drive(j, k), &
                    trial_m, trial_n, trial_parity, field_periods, h1_values, &
                    h1_derivatives, l2_values, &
                    real(j - 1, dp) / real(size(fields, 1), dp), &
                    real(k - 1, dp) / real(size(fields, 2), dp), weight, &
                    first_column, last_column, full, terms)
            end do
        end do
        !$omp end parallel
        info = 0
    end subroutine assemble_compatible_transformed_surface

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

    subroutine accumulate_transformed(fields, drive, trial_m, trial_n, &
            trial_parity, field_periods, h1_values, h1_derivatives, &
            l2_values, theta, zeta, weight, first_column, last_column, full, &
            terms)
        real(dp), intent(in) :: fields(:), drive, h1_values(:, :)
        real(dp), intent(in) :: h1_derivatives(:, :), l2_values(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: theta, zeta, weight
        integer, intent(in) :: first_column, last_column
        real(dp), intent(inout) :: full(:, :)
        real(dp), optional, intent(inout) :: terms(:, :, :)
        real(dp) :: cosine_rows(4, size(full, 1))
        real(dp) :: sine_rows(4, size(full, 1)), coefficients(3, 6)
        real(dp) :: phases(size(trial_m)), toroidal_wave
        integer :: trial

        call kernel_coefficients(fields, coefficients)
        cosine_rows = 0.0_dp
        sine_rows = 0.0_dp
        do trial = 1, size(trial_m)
            toroidal_wave = real(trial_n(trial), dp) &
                / real(field_periods, dp)
            phases(trial) = two_pi * (real(trial_m(trial), dp) * theta &
                - toroidal_wave * zeta)
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
        call transformed_rank_update(cosine_rows, sine_rows, phases, trial_n, &
            field_periods, drive, weight * abs(fields(7)), first_column, &
            last_column, full, terms)
    end subroutine accumulate_transformed

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

    subroutine transformed_rank_update(cosine_rows, sine_rows, phases, &
            trial_n, field_periods, drive, weight, first_column, last_column, &
            full, terms)
        real(dp), intent(in) :: cosine_rows(:, :), sine_rows(:, :), phases(:)
        integer, intent(in) :: trial_n(:), field_periods
        integer, intent(in) :: first_column, last_column
        real(dp), intent(in) :: drive, weight
        real(dp), intent(inout) :: full(:, :)
        real(dp), optional, intent(inout) :: terms(:, :, :)
        real(dp) :: products(2, 2), cosine(size(trial_n)), sine(size(trial_n))
        real(dp) :: contributions(4)
        integer :: a, b, first_trial, second_trial, trials

        trials = size(trial_n)
        cosine = cos(phases)
        sine = sin(phases)
        do b = first_column, last_column
            second_trial = modulo(b - 1, trials) + 1
            do a = 1, size(full, 1)
                first_trial = modulo(a - 1, trials) + 1
                call phase_product_coefficients(cosine(first_trial), &
                    sine(first_trial), cosine(second_trial), &
                    sine(second_trial), trial_n(first_trial), &
                    trial_n(second_trial), field_periods, products)
                call component_products(cosine_rows, sine_rows, products, &
                    a, b, contributions(1:3))
                contributions(4) = -drive * (products(1, 1) &
                    * cosine_rows(4, a) * cosine_rows(4, b) &
                    + products(1, 2) * cosine_rows(4, a) * sine_rows(4, b) &
                    + products(2, 1) * sine_rows(4, a) * cosine_rows(4, b) &
                    + products(2, 2) * sine_rows(4, a) * sine_rows(4, b))
                full(a, b) = full(a, b) + weight * sum(contributions)
                if (present(terms)) terms(a, b, :) = terms(a, b, :) &
                    + weight * contributions
            end do
        end do
    end subroutine transformed_rank_update

    pure subroutine component_products(cosine_rows, sine_rows, products, &
            a, b, contributions)
        real(dp), intent(in) :: cosine_rows(:, :), sine_rows(:, :)
        real(dp), intent(in) :: products(2, 2)
        integer, intent(in) :: a, b
        real(dp), intent(out) :: contributions(3)
        integer :: component

        do component = 1, 3
            contributions(component) = products(1, 1) &
                * cosine_rows(component, a) * cosine_rows(component, b) &
                + products(1, 2) * cosine_rows(component, a) &
                * sine_rows(component, b) + products(2, 1) &
                * sine_rows(component, a) * cosine_rows(component, b) &
                + products(2, 2) * sine_rows(component, a) &
                * sine_rows(component, b)
        end do
    end subroutine component_products

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
