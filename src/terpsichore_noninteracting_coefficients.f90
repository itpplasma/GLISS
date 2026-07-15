module terpsichore_noninteracting_coefficients
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use terpsichore_matrix_fixture, only: &
        terpsichore_matrix_fixture_t, terpsichore_potential_fixture_is_valid
    use terpsichore_model_policy, only: decode_terpsichore_model, &
        potential_model_non_interacting, terpsichore_model_config_t, &
        terpsichore_model_ok
    use terpsichore_pair_average, only: terpsichore_pair_averages, &
        terpsichore_pair_ok
    implicit none
    private

    integer, parameter, public :: terpsichore_coefficients_ok = 0
    integer, parameter, public :: terpsichore_coefficients_invalid = -1
    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: build_terpsichore_noninteracting_coefficients
    public :: build_terpsichore_noninteracting_coefficients_direct

contains

    subroutine build_terpsichore_noninteracting_coefficients(fixture, &
            coefficients, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: coefficients(:, :, :, :)
        integer, intent(out) :: info
        integer :: allocation_status

        info = terpsichore_coefficients_invalid
        if (.not. fixture_selects_noninteracting(fixture)) return
        allocate (coefficients(0:9, fixture%modes, fixture%modes, &
            fixture%intervals), stat=allocation_status)
        if (allocation_status /= 0) return
        call build_coefficients_resolved(fixture, coefficients, info)
    end subroutine build_terpsichore_noninteracting_coefficients

    subroutine build_coefficients_resolved(fixture, coefficients, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), intent(out) :: coefficients(0:, :, :, :)
        integer, intent(out) :: info
        integer :: interval

        info = terpsichore_coefficients_invalid
        if (.not. fixture_selects_noninteracting(fixture)) return
        if (size(coefficients, 1) /= 10) return
        if (size(coefficients, 2) /= fixture%modes) return
        if (size(coefficients, 3) /= fixture%modes) return
        if (size(coefficients, 4) /= fixture%intervals) return
        coefficients = 0.0_dp
        do interval = 1, fixture%intervals
            call build_interval_coefficients_fast(fixture, interval, &
                coefficients(:, :, :, interval), info)
            if (info /= terpsichore_coefficients_ok) return
        end do
        info = terpsichore_coefficients_invalid
        if (.not. all(ieee_is_finite(coefficients))) return
        info = terpsichore_coefficients_ok
    end subroutine build_coefficients_resolved

    ! naive point-pair oracle, O(points modes^2) per interval; retained
    ! for the equivalence test against the transform path.
    subroutine build_terpsichore_noninteracting_coefficients_direct( &
            fixture, coefficients, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: coefficients(:, :, :, :)
        integer, intent(out) :: info
        integer :: allocation_status, interval

        info = terpsichore_coefficients_invalid
        if (.not. fixture_selects_noninteracting(fixture)) return
        allocate (coefficients(0:9, fixture%modes, fixture%modes, &
            fixture%intervals), stat=allocation_status)
        if (allocation_status /= 0) return
        coefficients = 0.0_dp
        do interval = 1, fixture%intervals
            call build_interval_coefficients(fixture, interval, &
                coefficients(:, :, :, interval))
        end do
        if (.not. all(ieee_is_finite(coefficients))) return
        info = terpsichore_coefficients_ok
    end subroutine build_terpsichore_noninteracting_coefficients_direct

    subroutine build_interval_coefficients_fast(fixture, interval, &
            coefficient, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        real(dp), intent(out) :: coefficient(0:, :, :)
        integer, intent(out) :: info
        real(dp) :: field_grid(fixture%poloidal_points &
            * fixture%toroidal_points, 0:7)
        real(dp) :: normal_normal(fixture%modes, fixture%modes)
        real(dp) :: normal_tangent(fixture%modes, fixture%modes)
        real(dp) :: tangent_tangent(fixture%modes, fixture%modes)
        real(dp) :: point_field(0:7)
        integer :: kind, point, pair_status

        do point = 1, fixture%poloidal_points * fixture%toroidal_points
            call build_point_fields(fixture, interval, point, point_field)
            field_grid(point, :) = point_field
        end do
        coefficient = 0.0_dp
        do kind = 0, 7
            call terpsichore_pair_averages(field_grid(:, kind), &
                fixture%poloidal_points, fixture%toroidal_points, &
                fixture%stability_periods, fixture%field_periods, &
                fixture%mode_m, fixture%mode_n, &
                normal_normal, normal_tangent, tangent_tangent, &
                pair_status)
            if (pair_status /= terpsichore_pair_ok) then
                info = terpsichore_coefficients_invalid
                return
            end if
            select case (kind)
            case (1, 4)
                coefficient(kind, :, :) = normal_tangent
            case (3)
                coefficient(kind, :, :) = tangent_tangent
            case default
                coefficient(kind, :, :) = normal_normal
            end select
        end do
        call build_composite_coefficient(fixture, interval, coefficient)
        info = terpsichore_coefficients_ok
    end subroutine build_interval_coefficients_fast

    pure function fixture_selects_noninteracting(fixture) result(valid)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        type(terpsichore_model_config_t) :: model
        integer :: model_status
        logical :: valid

        valid = terpsichore_potential_fixture_is_valid(fixture)
        if (.not. valid) return
        valid = fixture%parity == 0.0_dp
        if (.not. valid) return
        valid = denominators_are_resolved(fixture)
        if (.not. valid) return
        call decode_terpsichore_model(fixture%legacy_modelk, model, model_status)
        valid = model_status == terpsichore_model_ok &
            .and. model%potential_model == potential_model_non_interacting
    end function fixture_selects_noninteracting

    pure function denominators_are_resolved(fixture) result(valid)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp) :: current_flux(fixture%intervals)
        real(dp) :: current_scale, sigma_scale
        real(dp), parameter :: reciprocal_product_floor = 1.0_dp &
            / sqrt(huge(1.0_dp))
        logical :: valid

        current_flux = fixture%current_j * &
            fixture%flux_p_slope(:fixture%intervals) &
            - fixture%current_i * fixture%flux_t_slope(:fixture%intervals)
        sigma_scale = maxval(abs(fixture%sigma_b))
        current_scale = maxval(abs(current_flux))
        valid = minval(abs(fixture%sigma_b)) > max( &
            epsilon(1.0_dp) * sigma_scale, reciprocal_product_floor) &
            .and. minval(abs(current_flux)) > max( &
            epsilon(1.0_dp) * current_scale, reciprocal_product_floor)
    end function denominators_are_resolved

    subroutine build_interval_coefficients(fixture, interval, coefficient)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        real(dp), intent(out) :: coefficient(0:, :, :)
        real(dp) :: field(0:7), normal_a, normal_b, tangent_a, tangent_b
        real(dp) :: scale
        integer :: a, b, point, replica

        coefficient = 0.0_dp
        scale = 2.0_dp / real(fixture%poloidal_points &
            * fixture%toroidal_points * fixture%stability_periods, dp)
        do replica = 0, fixture%stability_periods - 1
            do point = 1, fixture%poloidal_points * fixture%toroidal_points
                call build_point_fields(fixture, interval, point, field)
                do b = 1, fixture%modes
                    call phase_values(fixture, point, replica, b, normal_b, &
                        tangent_b)
                    do a = 1, fixture%modes
                        call phase_values(fixture, point, replica, a, &
                            normal_a, tangent_a)
                        coefficient(0, a, b) = coefficient(0, a, b) &
                            + scale * field(0) * normal_a * normal_b
                        coefficient(1, a, b) = coefficient(1, a, b) &
                            + scale * field(1) * normal_a * tangent_b
                        coefficient(2, a, b) = coefficient(2, a, b) &
                            + scale * field(2) * normal_a * normal_b
                        coefficient(3, a, b) = coefficient(3, a, b) &
                            + scale * field(3) * tangent_a * tangent_b
                        coefficient(4, a, b) = coefficient(4, a, b) &
                            + scale * field(4) * normal_a * tangent_b
                        coefficient(5, a, b) = coefficient(5, a, b) &
                            + scale * field(5) * normal_a * normal_b
                        coefficient(6, a, b) = coefficient(6, a, b) &
                            + scale * field(6) * normal_a * normal_b
                        coefficient(7, a, b) = coefficient(7, a, b) &
                            + scale * field(7) * normal_a * normal_b
                    end do
                end do
            end do
        end do
        call build_composite_coefficient(fixture, interval, coefficient)
    end subroutine build_interval_coefficients

    pure subroutine build_point_fields(fixture, interval, point, field)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, point
        real(dp), intent(out) :: field(0:7)
        real(dp) :: g_sp, gssu, sigma_b, t7, t8

        sigma_b = fixture%sigma_b(point, interval)
        g_sp = (fixture%sigma_b_s(point, interval) / sigma_b &
            - fixture%flux_p_slope(interval) &
            * fixture%metric_st_over_jacobian(point, interval)) &
            / fixture%flux_t_slope(interval)
        t7 = (fixture%current_j(interval) * fixture%flux_p_slope(interval) &
            - fixture%current_i(interval) * fixture%flux_t_slope(interval)) &
            / fixture%flux_t_slope(interval)**2
        t8 = (fixture%current_j(interval) &
            / fixture%flux_t_slope(interval))**2
        gssu = t7 * sigma_b &
            * fixture%metric_tt_over_jacobian(point, interval) - t8
        field(0) = 1.0_dp / sigma_b
        field(1) = fixture%sigma_b_s(point, interval) / sigma_b
        field(2) = fixture%pressure_slope(interval) &
            * fixture%signed_bjac(point, interval)
        field(3) = fixture%metric_ss_over_jacobian(point, interval)
        field(4) = fixture%current_i(interval) &
            * fixture%metric_st_over_jacobian(point, interval) &
            + fixture%current_j(interval) * g_sp
        field(5) = fixture%parallel_current(point, interval)
        field(6) = gssu / sigma_b
        field(7) = fixture%pressure_slope(interval) &
            * fixture%signed_bjac(point, interval) &
            * fixture%signed_bjac_radial(point, interval)
    end subroutine build_point_fields

    ! phases at absolute toroidal mode numbers; the replica argument
    ! shifts the point into later stability periods so the direct
    ! oracle can quadrature the full torus (the one-period sum aliases
    ! pair channels whose toroidal number is not stability-divisible).
    pure subroutine phase_values(fixture, point, replica, mode, normal, &
            tangent)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: point, replica, mode
        real(dp), intent(out) :: normal, tangent
        real(dp) :: angle, theta, zeta
        integer :: j, k

        j = modulo(point - 1, fixture%poloidal_points)
        k = (point - 1) / fixture%poloidal_points
        theta = two_pi * real(j, dp) / real(fixture%poloidal_points, dp)
        zeta = two_pi * (real(k, dp) &
            / real(fixture%toroidal_points, dp) + real(replica, dp)) &
            / real(fixture%field_periods, dp)
        angle = real(fixture%mode_m(mode), dp) * theta &
            - real(fixture%mode_n(mode), dp) * zeta
        normal = sin(angle)
        tangent = cos(angle)
    end subroutine phase_values

    pure subroutine build_composite_coefficient(fixture, interval, coefficient)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        real(dp), intent(inout) :: coefficient(0:, :, :)
        real(dp) :: flux_cross, flux_ratio, flux_square, midpoint
        real(dp) :: current_flux, left_parallel, right_parallel
        integer :: a, b

        midpoint = 0.5_dp * (fixture%s(interval - 1) + fixture%s(interval))
        flux_cross = fixture%flux_p_slope(interval) &
            * fixture%flux_t_curve(interval) &
            - fixture%flux_t_slope(interval) * fixture%flux_p_curve(interval)
        current_flux = fixture%current_j(interval) &
            * fixture%flux_p_slope(interval) &
            - fixture%current_i(interval) * fixture%flux_t_slope(interval)
        flux_ratio = flux_cross / current_flux
        flux_square = flux_cross**2 / current_flux
        do b = 1, fixture%modes
            right_parallel = real(fixture%mode_m(b), dp) &
                * fixture%flux_p_slope(interval) &
                - real(fixture%mode_n(b), dp) &
                * fixture%flux_t_slope(interval)
            do a = 1, fixture%modes
                left_parallel = real(fixture%mode_m(a), dp) &
                    * fixture%flux_p_slope(interval) &
                    - real(fixture%mode_n(a), dp) &
                    * fixture%flux_t_slope(interval)
                coefficient(9, a, b) = composite_value(fixture, interval, &
                    midpoint, flux_ratio, flux_square, current_flux, &
                    left_parallel, right_parallel, a, b, coefficient)
            end do
        end do
    end subroutine build_composite_coefficient

    pure function composite_value(fixture, interval, midpoint, flux_ratio, &
            flux_square, current_flux, left_parallel, right_parallel, a, b, &
            coefficient) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: midpoint, flux_ratio, flux_square, current_flux
        real(dp), intent(in) :: left_parallel, right_parallel
        real(dp), intent(in) :: coefficient(0:, :, :)
        real(dp) :: left_shift, right_shift, value
        real(dp) :: current_curve

        current_curve = (fixture%current_j(interval) &
            * fixture%flux_p_curve(interval) &
            - fixture%current_i(interval) * fixture%flux_t_curve(interval)) &
            / current_flux
        left_shift = current_curve - fixture%radial_power(a) / midpoint
        right_shift = current_curve - fixture%radial_power(b) / midpoint
        value = coefficient(7, a, b) + flux_cross_value(fixture, interval) &
            * coefficient(5, a, b) + flux_square * coefficient(6, a, b) &
            + left_shift * right_shift * current_flux * coefficient(0, a, b) &
            + (current_curve - (fixture%radial_power(a) &
            + fixture%radial_power(b)) / midpoint) * coefficient(2, a, b) &
            + left_parallel * right_parallel * coefficient(3, a, b) &
            - right_parallel * (flux_ratio * coefficient(4, a, b) &
            + left_shift * coefficient(1, a, b)) &
            - left_parallel * (flux_ratio * coefficient(4, b, a) &
            + right_shift * coefficient(1, b, a))
    end function composite_value

    pure function flux_cross_value(fixture, interval) result(value)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        real(dp) :: value

        value = fixture%flux_p_slope(interval) &
            * fixture%flux_t_curve(interval) &
            - fixture%flux_t_slope(interval) * fixture%flux_p_curve(interval)
    end function flux_cross_value

end module terpsichore_noninteracting_coefficients
