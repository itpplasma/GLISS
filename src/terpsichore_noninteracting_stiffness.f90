module terpsichore_noninteracting_stiffness
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: add_mapped_dynamic_element, &
        dynamic_family_layout_t, dynamic_layout_ok
    use fourier_phase_kind, only: phase_sine
    use terpsichore_matrix_fixture, only: terpsichore_matrix_fixture_t
    use terpsichore_noninteracting_coefficients, only: &
        build_terpsichore_noninteracting_coefficients, &
        terpsichore_coefficients_ok
    use terpsichore_reduced_layout, only: &
        build_terpsichore_reduced_fixed_boundary_layout, &
        terpsichore_reduced_layout_ok
    implicit none
    private

    integer, parameter, public :: terpsichore_noninteracting_ok = 0
    integer, parameter, public :: terpsichore_noninteracting_invalid = -1
    real(dp), parameter :: negative_jacobian_to_physical_sign = -1.0_dp

    type :: interval_factors_t
        real(dp) :: derivative_scale
        real(dp) :: radial_weight
        real(dp) :: current_flux
        real(dp) :: current_curve
        real(dp) :: flux_cross
        real(dp) :: radial_metric
    end type interval_factors_t

    public :: assemble_terpsichore_noninteracting_fixed_boundary_stiffness

contains

    subroutine assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
            fixture, stiffness, layout, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        real(dp), allocatable :: coefficients(:, :, :, :), element(:, :)
        integer, allocatable :: element_to_global(:, :), parity(:)
        integer :: allocation_status, interval

        info = terpsichore_noninteracting_invalid
        call build_terpsichore_noninteracting_coefficients(fixture, &
            coefficients, info)
        if (info /= terpsichore_coefficients_ok) return
        allocate (parity(fixture%modes), source=phase_sine, &
            stat=allocation_status)
        if (allocation_status /= 0) return
        call build_terpsichore_reduced_fixed_boundary_layout(fixture%mode_m, &
            fixture%mode_n, parity, fixture%intervals, layout, &
            element_to_global, info)
        if (info /= terpsichore_reduced_layout_ok) return
        allocate (stiffness(layout%total_unknowns, layout%total_unknowns), &
            source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (element(3 * fixture%modes, 3 * fixture%modes), &
            stat=allocation_status)
        if (allocation_status /= 0) return
        do interval = 1, fixture%intervals
            call build_local_element(fixture, coefficients(:, :, :, interval), &
                interval, element)
            if (.not. all(ieee_is_finite(element))) return
            call add_mapped_dynamic_element(element_to_global(:, interval), &
                element, stiffness, info)
            if (info /= dynamic_layout_ok) then
                info = terpsichore_noninteracting_invalid
                return
            end if
        end do
        if (.not. all(ieee_is_finite(stiffness))) return
        info = terpsichore_noninteracting_ok
    end subroutine assemble_terpsichore_noninteracting_fixed_boundary_stiffness

    pure subroutine build_local_element(fixture, coefficient, interval, element)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval
        real(dp), intent(out) :: element(:, :)
        type(interval_factors_t) :: factors

        call build_interval_factors(fixture, interval, factors)
        element = 0.0_dp
        call build_diagonal_blocks(fixture, coefficient, interval, factors, &
            element)
        call build_coupling_blocks(fixture, coefficient, interval, factors, &
            element)
        element = negative_jacobian_to_physical_sign * element
    end subroutine build_local_element

    pure subroutine build_interval_factors(fixture, interval, factors)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        type(interval_factors_t), intent(out) :: factors
        real(dp) :: cell_width, midpoint

        cell_width = fixture%s(interval) - fixture%s(interval - 1)
        midpoint = 0.5_dp * (fixture%s(interval) + fixture%s(interval - 1))
        factors%radial_weight = real(fixture%intervals, dp) * cell_width
        factors%derivative_scale = 2.0_dp / cell_width
        factors%current_flux = fixture%current_j(interval) &
            * fixture%flux_p_slope(interval) &
            - fixture%current_i(interval) * fixture%flux_t_slope(interval)
        factors%flux_cross = fixture%flux_p_slope(interval) &
            * fixture%flux_t_curve(interval) &
            - fixture%flux_t_slope(interval) * fixture%flux_p_curve(interval)
        factors%current_curve = fixture%current_j(interval) &
            * fixture%flux_p_curve(interval) &
            - fixture%current_i(interval) * fixture%flux_t_curve(interval)
        factors%radial_metric = factors%current_flux / midpoint
    end subroutine build_interval_factors

    pure subroutine build_diagonal_blocks(fixture, coefficient, interval, &
            factors, element)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval
        type(interval_factors_t), intent(in) :: factors
        real(dp), intent(inout) :: element(:, :)
        real(dp) :: left_diagonal, right_diagonal
        integer :: a, b, modes

        modes = fixture%modes
        do b = 1, fixture%modes
            do a = 1, b
                left_diagonal = normal_diagonal(fixture, coefficient, interval, &
                    a, b, factors%derivative_scale, factors%radial_weight, &
                    factors%current_flux, factors%current_curve, &
                    factors%radial_metric, 1.0_dp)
                right_diagonal = normal_diagonal(fixture, coefficient, interval, &
                    a, b, factors%derivative_scale, factors%radial_weight, &
                    factors%current_flux, factors%current_curve, &
                    factors%radial_metric, -1.0_dp)
                call set_symmetric(element, a, b, left_diagonal)
                call set_symmetric(element, modes + a, modes + b, &
                    right_diagonal)
                call set_symmetric(element, 2 * modes + a, 2 * modes + b, &
                    tangential_diagonal(fixture, coefficient, interval, a, b, &
                    factors%radial_weight, factors%current_flux))
            end do
        end do
    end subroutine build_diagonal_blocks

    pure subroutine build_coupling_blocks(fixture, coefficient, interval, &
            factors, element)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval
        type(interval_factors_t), intent(in) :: factors
        real(dp), intent(inout) :: element(:, :)
        integer :: a, b, modes

        modes = fixture%modes
        do b = 1, fixture%modes
            do a = 1, modes
                call set_transpose_pair(element, a, modes + b, &
                    normal_cross(fixture, coefficient, interval, a, b, &
                    factors%derivative_scale, factors%radial_weight, &
                    factors%current_flux, factors%current_curve, &
                    factors%radial_metric))
                call set_transpose_pair(element, a, 2 * modes + b, &
                    left_tangential(fixture, coefficient, interval, a, b, &
                    factors%derivative_scale, factors%radial_weight, &
                    factors%current_flux, factors%current_curve, &
                    factors%radial_metric, factors%flux_cross))
                call set_transpose_pair(element, 2 * modes + a, modes + b, &
                    tangential_right(fixture, coefficient, interval, a, b, &
                    factors%derivative_scale, factors%radial_weight, &
                    factors%current_flux, factors%current_curve, &
                    factors%radial_metric, factors%flux_cross))
            end do
        end do
    end subroutine build_coupling_blocks

    pure function normal_diagonal(fixture, coefficient, interval, a, b, &
            derivative_scale, radial_weight, current_flux, current_curve, &
            radial_metric, orientation) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative_scale, radial_weight, current_flux
        real(dp), intent(in) :: current_curve, radial_metric, orientation
        real(dp) :: left_parallel, right_parallel, scale, value
        real(dp) :: derivative_term, bending_term, midpoint

        midpoint = 0.5_dp * (fixture%s(interval) + fixture%s(interval - 1))
        left_parallel = parallel_mode(fixture, interval, a)
        right_parallel = parallel_mode(fixture, interval, b)
        scale = 0.25_dp * radial_weight &
            / midpoint**(fixture%radial_power(a) + fixture%radial_power(b))
        derivative_term = derivative_scale * (right_parallel &
            * coefficient(1, a, b) + left_parallel * coefficient(1, b, a) &
            - 2.0_dp * coefficient(2, a, b) &
            - (2.0_dp * current_curve - (fixture%radial_power(a) &
            + fixture%radial_power(b)) * radial_metric) &
            * coefficient(0, a, b))
        bending_term = derivative_scale**2 * current_flux &
            * coefficient(0, a, b)
        value = scale * (bending_term + orientation * derivative_term &
            + coefficient(9, a, b))
    end function normal_diagonal

    pure function normal_cross(fixture, coefficient, interval, a, b, &
            derivative_scale, radial_weight, current_flux, current_curve, &
            radial_metric) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative_scale, radial_weight, current_flux
        real(dp), intent(in) :: current_curve, radial_metric
        real(dp) :: left_parallel, right_parallel, scale, value
        real(dp) :: coupling, midpoint

        midpoint = 0.5_dp * (fixture%s(interval) + fixture%s(interval - 1))
        left_parallel = parallel_mode(fixture, interval, a)
        right_parallel = parallel_mode(fixture, interval, b)
        scale = 0.25_dp * radial_weight &
            / midpoint**(fixture%radial_power(a) + fixture%radial_power(b))
        coupling = right_parallel * coefficient(1, a, b) &
            - left_parallel * coefficient(1, b, a) &
            + (fixture%radial_power(b) - fixture%radial_power(a)) &
            * radial_metric * coefficient(0, a, b)
        value = scale * (-derivative_scale**2 * current_flux &
            * coefficient(0, a, b) + coefficient(9, a, b) &
            + derivative_scale * coupling)
    end function normal_cross

    pure function tangential_diagonal(fixture, coefficient, interval, a, b, &
            radial_weight, current_flux) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: radial_weight, current_flux
        real(dp) :: value

        value = radial_weight / current_flux &
            * (parallel_mode(fixture, interval, a) &
            * parallel_mode(fixture, interval, b) * coefficient(6, a, b) &
            + current_mode(fixture, interval, a) &
            * current_mode(fixture, interval, b) * coefficient(0, a, b))
    end function tangential_diagonal

    pure function left_tangential(fixture, coefficient, interval, a, b, &
            derivative_scale, radial_weight, current_flux, current_curve, &
            radial_metric, flux_cross) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative_scale, radial_weight, current_flux
        real(dp), intent(in) :: current_curve, radial_metric, flux_cross
        real(dp) :: right_parallel, right_current, scale, value, midpoint

        midpoint = 0.5_dp * (fixture%s(interval) + fixture%s(interval - 1))
        right_parallel = parallel_mode(fixture, interval, b)
        right_current = current_mode(fixture, interval, b)
        scale = 0.5_dp * radial_weight / midpoint**fixture%radial_power(a)
        value = scale * (-derivative_scale * right_current &
            * coefficient(0, a, b) + (right_parallel &
            * (parallel_mode(fixture, interval, a) * coefficient(4, b, a) &
            - current_flux * coefficient(5, a, b) &
            - flux_cross * coefficient(6, a, b)) &
            + right_current * (coefficient(2, a, b) &
            - parallel_mode(fixture, interval, a) * coefficient(1, b, a) &
            + (current_curve - fixture%radial_power(a) * radial_metric) &
            * coefficient(0, a, b))) / current_flux)
    end function left_tangential

    pure function tangential_right(fixture, coefficient, interval, a, b, &
            derivative_scale, radial_weight, current_flux, current_curve, &
            radial_metric, flux_cross) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(in) :: coefficient(0:, :, :)
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative_scale, radial_weight, current_flux
        real(dp), intent(in) :: current_curve, radial_metric, flux_cross
        real(dp) :: left_parallel, left_current, scale, value, midpoint

        midpoint = 0.5_dp * (fixture%s(interval) + fixture%s(interval - 1))
        left_parallel = parallel_mode(fixture, interval, a)
        left_current = current_mode(fixture, interval, a)
        scale = 0.5_dp * radial_weight / midpoint**fixture%radial_power(b)
        value = scale * (derivative_scale * left_current &
            * coefficient(0, a, b) + (left_parallel &
            * (parallel_mode(fixture, interval, b) * coefficient(4, a, b) &
            - current_flux * coefficient(5, a, b) &
            - flux_cross * coefficient(6, a, b)) &
            + left_current * (coefficient(2, a, b) &
            - parallel_mode(fixture, interval, b) * coefficient(1, a, b) &
            + (current_curve - fixture%radial_power(b) * radial_metric) &
            * coefficient(0, a, b))) / current_flux)
    end function tangential_right

    pure function parallel_mode(fixture, interval, mode) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, mode
        real(dp) :: value

        value = real(fixture%mode_m(mode), dp) &
            * fixture%flux_p_slope(interval) &
            - real(fixture%mode_n(mode), dp) &
            * fixture%flux_t_slope(interval)
    end function parallel_mode

    pure function current_mode(fixture, interval, mode) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, mode
        real(dp) :: value

        value = real(fixture%mode_m(mode), dp) * fixture%current_i(interval) &
            - real(fixture%mode_n(mode), dp) * fixture%current_j(interval)
    end function current_mode

    pure subroutine set_symmetric(matrix, row, column, value)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(in) :: row, column
        real(dp), intent(in) :: value

        matrix(row, column) = value
        matrix(column, row) = value
    end subroutine set_symmetric

    pure subroutine set_transpose_pair(matrix, row, column, value)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(in) :: row, column
        real(dp), intent(in) :: value

        matrix(row, column) = value
        matrix(column, row) = value
    end subroutine set_transpose_pair

end module terpsichore_noninteracting_stiffness
