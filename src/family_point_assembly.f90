module family_point_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_factor_topology, only: phase_product_coefficients
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    use two_component_kernel, only: two_component_components
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

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
        real(dp) :: weight
        integer :: j, l, period

        call validate_surface_inputs(fields, drive, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, &
            radial_coordinate, radial_step, full, info)
        if (info /= 0) return
        weight = 1.0_dp / real(size(fields, 1) * size(fields, 2) &
            * field_periods, dp)
        do period = 0, field_periods - 1
            do l = 1, size(fields, 2)
                do j = 1, size(fields, 1)
                    call accumulate_direct_point(fields(j, l, :), drive(j, l), &
                        trial_m, trial_n, trial_parity, normal_stored_power, &
                        field_periods, radial_space, radial_coordinate, &
                        (real(j, dp) - 1.0_dp) / real(size(fields, 1), dp), &
                        (real(l, dp) - 1.0_dp) / real(size(fields, 2), dp) &
                        + real(period, dp), radial_step, weight, full, info)
                    if (info /= 0) return
                end do
            end do
        end do
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
        real(dp) :: weight
        integer :: j, l

        call validate_surface_inputs(fields, drive, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, &
            radial_coordinate, radial_step, full, info)
        if (info /= 0) return
        weight = 1.0_dp / real(size(fields, 1) * size(fields, 2), dp)
        do l = 1, size(fields, 2)
            do j = 1, size(fields, 1)
                call accumulate_transformed_point(fields(j, l, :), &
                    drive(j, l), trial_m, trial_n, trial_parity, &
                    normal_stored_power, field_periods, radial_space, &
                    radial_coordinate, &
                    (real(j, dp) - 1.0_dp) / real(size(fields, 1), dp), &
                    (real(l, dp) - 1.0_dp) / real(size(fields, 2), dp), &
                    radial_step, weight, full, info)
                if (info /= 0) return
            end do
        end do
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

    subroutine accumulate_direct_point(fields, drive, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, radial_space, &
            radial_coordinate, theta, zeta, radial_step, weight, full, info)
        real(dp), intent(in) :: fields(:), drive, theta, zeta
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate
        real(dp), intent(in) :: radial_step, weight
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp) :: rows(4, 3 * size(trial_m))
        real(dp) :: phase, cosine, sine, toroidal_wave
        real(dp) :: value, dvalue, dother
        real(dp) :: normal_values(2), normal_derivatives(2)
        real(dp) :: c1_of(6), c2_of(6), c3_of(6)
        integer :: trials, trial

        call kernel_linear_coefficients(fields, c1_of, c2_of, c3_of)
        trials = size(trial_m)
        rows = 0.0_dp
        do trial = 1, trials
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, 0.5_dp, normal_values, &
                normal_derivatives, info, normal_stored_power(trial))
            if (info /= radial_space_ok) return
            toroidal_wave = real(trial_n(trial), dp) &
                / real(field_periods, dp)
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - toroidal_wave * zeta)
            cosine = cos(phase)
            sine = sin(phase)
            if (trial_parity(trial) == 1) then
                value = cosine
                dvalue = -sine
                dother = cosine
            else
                value = sine
                dvalue = cosine
                dother = -sine
            end if
            call add_trial_rows(rows, trial, trials, trial_m(trial), &
                toroidal_wave, normal_values, normal_derivatives, value, &
                dvalue, dother, c1_of, c2_of, c3_of)
        end do
        call rank_updates(rows, drive, weight * abs(fields(7)), full)
        info = 0
    end subroutine accumulate_direct_point

    subroutine accumulate_transformed_point(fields, drive, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, radial_space, &
            radial_coordinate, theta, zeta, radial_step, weight, full, info)
        real(dp), intent(in) :: fields(:), drive, theta, zeta
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step, weight
        real(dp), intent(inout) :: full(:, :)
        integer, intent(out) :: info
        real(dp) :: cosine_rows(4, 3 * size(trial_m))
        real(dp) :: sine_rows(4, 3 * size(trial_m))
        real(dp) :: phases(size(trial_m)), normal_values(2)
        real(dp) :: normal_derivatives(2), c1_of(6), c2_of(6), c3_of(6)
        real(dp) :: toroidal_wave
        integer :: trials, trial

        call kernel_linear_coefficients(fields, c1_of, c2_of, c3_of)
        trials = size(trial_m)
        cosine_rows = 0.0_dp
        sine_rows = 0.0_dp
        do trial = 1, trials
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, 0.5_dp, normal_values, &
                normal_derivatives, info, normal_stored_power(trial))
            if (info /= radial_space_ok) return
            toroidal_wave = real(trial_n(trial), dp) &
                / real(field_periods, dp)
            phases(trial) = two_pi * (real(trial_m(trial), dp) * theta &
                - toroidal_wave * zeta)
            if (trial_parity(trial) == 1) then
                call add_trial_rows(cosine_rows, trial, trials, &
                    trial_m(trial), toroidal_wave, normal_values, &
                    normal_derivatives, 1.0_dp, 0.0_dp, 1.0_dp, c1_of, &
                    c2_of, c3_of)
                call add_trial_rows(sine_rows, trial, trials, &
                    trial_m(trial), toroidal_wave, normal_values, &
                    normal_derivatives, 0.0_dp, -1.0_dp, 0.0_dp, c1_of, &
                    c2_of, c3_of)
            else
                call add_trial_rows(cosine_rows, trial, trials, &
                    trial_m(trial), toroidal_wave, normal_values, &
                    normal_derivatives, 0.0_dp, 1.0_dp, 0.0_dp, c1_of, &
                    c2_of, c3_of)
                call add_trial_rows(sine_rows, trial, trials, &
                    trial_m(trial), toroidal_wave, normal_values, &
                    normal_derivatives, 1.0_dp, 0.0_dp, -1.0_dp, c1_of, &
                    c2_of, c3_of)
            end if
        end do
        call rank_updates_transformed(cosine_rows, sine_rows, phases, &
            trial_n, field_periods, drive, weight * abs(fields(7)), full)
        info = 0
    end subroutine accumulate_transformed_point

    subroutine kernel_linear_coefficients(fields, c1_of, c2_of, c3_of)
        real(dp), intent(in) :: fields(:)
        real(dp), intent(out) :: c1_of(6), c2_of(6), c3_of(6)
        real(dp) :: unit_inputs(6)
        integer :: entry_index

        do entry_index = 1, 6
            unit_inputs = 0.0_dp
            unit_inputs(entry_index) = 1.0_dp
            call two_component_components(fields(1), fields(2), &
                fields(3), fields(4), fields(5), fields(6), &
                fields(7), fields(8), fields(9), fields(10), &
                fields(11), fields(12), fields(13), &
                unit_inputs(1), unit_inputs(2), unit_inputs(3), &
                unit_inputs(4), unit_inputs(5), unit_inputs(6), &
                c1_of(entry_index), c2_of(entry_index), &
                c3_of(entry_index))
        end do
    end subroutine kernel_linear_coefficients

    subroutine add_trial_rows(rows, trial, trials, m, toroidal_wave, &
            normal_values, normal_derivatives, value, dvalue, dother, &
            c1_of, c2_of, c3_of)
        real(dp), intent(inout) :: rows(:, :)
        integer, intent(in) :: trial, trials, m
        real(dp), intent(in) :: toroidal_wave, normal_values(2)
        real(dp), intent(in) :: normal_derivatives(2)
        real(dp), intent(in) :: value, dvalue, dother
        real(dp), intent(in) :: c1_of(6), c2_of(6), c3_of(6)
        integer :: entry_index

        do entry_index = 1, 6
            call add_linear(rows, trial, trials, entry_index, value, &
                dvalue, dother, m, toroidal_wave, normal_values, &
                normal_derivatives, c1_of(entry_index), c2_of(entry_index), &
                c3_of(entry_index))
        end do
        rows(4, trial) = rows(4, trial) + normal_values(1) * value
        rows(4, trials + trial) = rows(4, trials + trial) &
            + normal_values(2) * value
    end subroutine add_trial_rows

    subroutine add_linear(rows, trial, trials, entry_index, value, &
            dvalue, dother, m, toroidal_wave, normal_values, &
            normal_derivatives, c1, c2, c3)
        real(dp), intent(inout) :: rows(:, :)
        integer, intent(in) :: trial, trials, entry_index, m
        real(dp), intent(in) :: value, dvalue, dother, toroidal_wave
        real(dp), intent(in) :: normal_values(2), normal_derivatives(2)
        real(dp), intent(in) :: c1, c2, c3

        select case (entry_index)
        case (1)
            call apply(rows, trial, trials, value * normal_values(1), &
                value * normal_values(2), 0.0_dp, c1, c2, c3)
        case (2)
            call apply(rows, trial, trials, value * normal_derivatives(1), &
                value * normal_derivatives(2), 0.0_dp, c1, c2, c3)
        case (3)
            call apply(rows, trial, trials, &
                two_pi * real(m, dp) * dvalue * normal_values(1), &
                two_pi * real(m, dp) * dvalue * normal_values(2), 0.0_dp, &
                c1, c2, c3)
        case (4)
            call apply(rows, trial, trials, &
                -two_pi * toroidal_wave * dvalue * normal_values(1), &
                -two_pi * toroidal_wave * dvalue * normal_values(2), &
                0.0_dp, c1, c2, c3)
        case (5)
            call apply(rows, trial, trials, 0.0_dp, 0.0_dp, &
                two_pi * real(m, dp) * dother, c1, c2, c3)
        case (6)
            call apply(rows, trial, trials, 0.0_dp, 0.0_dp, &
                -two_pi * toroidal_wave * dother, c1, c2, c3)
        end select
    end subroutine add_linear

    subroutine apply(rows, trial, trials, left_factor, right_factor, &
            tangential_factor, c1, c2, c3)
        real(dp), intent(inout) :: rows(:, :)
        integer, intent(in) :: trial, trials
        real(dp), intent(in) :: left_factor, right_factor
        real(dp), intent(in) :: tangential_factor, c1, c2, c3

        rows(1, trial) = rows(1, trial) + c1 * left_factor
        rows(1, trials + trial) = rows(1, trials + trial) &
            + c1 * right_factor
        rows(1, 2 * trials + trial) = rows(1, 2 * trials + trial) &
            + c1 * tangential_factor
        rows(2, trial) = rows(2, trial) + c2 * left_factor
        rows(2, trials + trial) = rows(2, trials + trial) &
            + c2 * right_factor
        rows(2, 2 * trials + trial) = rows(2, 2 * trials + trial) &
            + c2 * tangential_factor
        rows(3, trial) = rows(3, trial) + c3 * left_factor
        rows(3, trials + trial) = rows(3, trials + trial) &
            + c3 * right_factor
        rows(3, 2 * trials + trial) = rows(3, 2 * trials + trial) &
            + c3 * tangential_factor
    end subroutine apply

    subroutine rank_updates(rows, drive, weight, full)
        real(dp), intent(in) :: rows(:, :), drive, weight
        real(dp), intent(inout) :: full(:, :)
        integer :: a, b, component

        do b = 1, size(full, 2)
            do a = 1, size(full, 1)
                do component = 1, 3
                    full(a, b) = full(a, b) + weight &
                        * rows(component, a) * rows(component, b)
                end do
                full(a, b) = full(a, b) - weight * drive &
                    * rows(4, a) * rows(4, b)
            end do
        end do
    end subroutine rank_updates

    subroutine rank_updates_transformed(cosine_rows, sine_rows, phases, &
            trial_n, field_periods, drive, weight, full)
        real(dp), intent(in) :: cosine_rows(:, :), sine_rows(:, :), phases(:)
        integer, intent(in) :: trial_n(:), field_periods
        real(dp), intent(in) :: drive, weight
        real(dp), intent(inout) :: full(:, :)
        real(dp) :: cc, cs, sc, ss, product
        real(dp) :: products(2, 2)
        real(dp) :: cosine(size(trial_n)), sine(size(trial_n))
        integer :: a, b, component, first_trial, second_trial, trials
        integer :: first_block, second_block

        trials = size(trial_n)
        cosine = cos(phases)
        sine = sin(phases)
        do second_trial = 1, trials
            do first_trial = 1, trials
                call phase_product_coefficients(cosine(first_trial), &
                    sine(first_trial), cosine(second_trial), &
                    sine(second_trial), trial_n(first_trial), &
                    trial_n(second_trial), field_periods, products)
                cc = products(1, 1)
                cs = products(1, 2)
                sc = products(2, 1)
                ss = products(2, 2)
                do second_block = 0, 2
                    b = second_block * trials + second_trial
                    do first_block = 0, 2
                        a = first_block * trials + first_trial
                        product = 0.0_dp
                        do component = 1, 3
                            product = product &
                                + cc * cosine_rows(component, a) &
                                * cosine_rows(component, b) &
                                + cs * cosine_rows(component, a) &
                                * sine_rows(component, b) &
                                + sc * sine_rows(component, a) &
                                * cosine_rows(component, b) &
                                + ss * sine_rows(component, a) &
                                * sine_rows(component, b)
                        end do
                        product = product - drive &
                            * (cc * cosine_rows(4, a) * cosine_rows(4, b) &
                            + cs * cosine_rows(4, a) * sine_rows(4, b) &
                            + sc * sine_rows(4, a) * cosine_rows(4, b) &
                            + ss * sine_rows(4, a) * sine_rows(4, b))
                        full(a, b) = full(a, b) + weight * product
                    end do
                end do
            end do
        end do
    end subroutine rank_updates_transformed

end module family_point_assembly
