module terpsichore_pseudoplasma_stiffness
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use terpsichore_pseudoplasma_fixture, only: &
        terpsichore_pseudoplasma_fixture_is_valid, &
        terpsichore_pseudoplasma_fixture_t
    implicit none
    private

    integer, parameter, public :: pseudoplasma_stiffness_ok = 0
    integer, parameter, public :: pseudoplasma_stiffness_invalid = -1
    real(dp), parameter :: negative_jacobian_to_physical_sign = -1.0_dp

    public :: assemble_terpsichore_pseudoplasma_stiffness

contains

    subroutine assemble_terpsichore_pseudoplasma_stiffness(fixture, &
            stiffness, info)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: element(:, :)
        integer, allocatable :: element_to_global(:)
        integer :: allocation_status, interval, order

        info = pseudoplasma_stiffness_invalid
        if (.not. terpsichore_pseudoplasma_fixture_is_valid(fixture)) return
        order = 2 * fixture%modes * fixture%vacuum_intervals
        allocate (stiffness(order, order), source=0.0_dp, &
            stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (element(3 * fixture%modes, 3 * fixture%modes), &
            element_to_global(3 * fixture%modes), stat=allocation_status)
        if (allocation_status /= 0) return
        do interval = 1, fixture%vacuum_intervals
            call build_local_element(fixture, interval, element)
            if (.not. all(ieee_is_finite(element))) return
            call build_element_map(fixture, interval, element_to_global)
            call add_element(element, element_to_global, stiffness)
        end do
        if (.not. all(ieee_is_finite(stiffness))) return
        info = pseudoplasma_stiffness_ok
    end subroutine assemble_terpsichore_pseudoplasma_stiffness

    pure subroutine build_local_element(fixture, interval, element)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        real(dp), intent(out) :: element(:, :)
        real(dp) :: derivative_scale, radial_weight

        derivative_scale = 2.0_dp / &
            (fixture%s(interval) - fixture%s(interval - 1))
        radial_weight = real(fixture%plasma_intervals, dp) &
            * (fixture%s(interval) - fixture%s(interval - 1))
        element = 0.0_dp
        call build_normal_blocks(fixture, interval, derivative_scale, &
            radial_weight, element)
        call build_tangential_blocks(fixture, interval, derivative_scale, &
            radial_weight, element)
        element = negative_jacobian_to_physical_sign * element
    end subroutine build_local_element

    pure subroutine build_normal_blocks(fixture, interval, derivative_scale, &
            radial_weight, element)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        real(dp), intent(in) :: derivative_scale, radial_weight
        real(dp), intent(inout) :: element(:, :)
        real(dp) :: radial, forward, backward, bending
        integer :: a, b, modes

        modes = fixture%modes
        do b = 1, modes
            do a = 1, b
                radial = radial_normal_term(fixture, interval, a, b, &
                    derivative_scale)
                forward = forward_normal_term(fixture, interval, a, b, &
                    derivative_scale)
                backward = forward_normal_term(fixture, interval, b, a, &
                    derivative_scale)
                bending = parallel_mode(fixture, a) &
                    * parallel_mode(fixture, b) &
                    * fixture%coefficient(3, a, b, interval)
                call set_symmetric(element, a, b, 0.25_dp * radial_weight &
                    * (radial + forward + backward + bending))
                call set_symmetric(element, modes + a, modes + b, &
                    0.25_dp * radial_weight &
                    * (radial - forward - backward + bending))
                element(a, modes + b) = 0.25_dp * radial_weight &
                    * (-radial + forward - backward + bending)
                element(modes + b, a) = element(a, modes + b)
                if (a == b) cycle
                element(b, modes + a) = 0.25_dp * radial_weight &
                    * (-radial - forward + backward + bending)
                element(modes + a, b) = element(b, modes + a)
            end do
        end do
    end subroutine build_normal_blocks

    pure subroutine build_tangential_blocks(fixture, interval, &
            derivative_scale, radial_weight, element)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        real(dp), intent(in) :: derivative_scale, radial_weight
        real(dp), intent(inout) :: element(:, :)
        real(dp) :: left_coupling, right_coupling
        integer :: a, b, modes

        modes = fixture%modes
        do b = 1, modes
            do a = 1, b
                call set_symmetric(element, 2 * modes + a, 2 * modes + b, &
                    radial_weight * tangential_term(fixture, interval, a, b))
            end do
            do a = 1, modes
                left_coupling = 0.5_dp * radial_weight &
                    * left_tangential_term(fixture, interval, a, b, &
                    derivative_scale)
                right_coupling = 0.5_dp * radial_weight &
                    * right_tangential_term(fixture, interval, a, b, &
                    derivative_scale)
                element(a, 2 * modes + b) = left_coupling
                element(2 * modes + b, a) = left_coupling
                element(2 * modes + a, modes + b) = right_coupling
                element(modes + b, 2 * modes + a) = right_coupling
            end do
        end do
    end subroutine build_tangential_blocks

    pure function radial_normal_term(fixture, interval, a, b, derivative) &
            result(value)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative
        real(dp) :: value, fp, ft

        fp = fixture%flux_p_slope
        ft = fixture%flux_t_slope
        value = derivative**2 * (fp**2 &
            * fixture%coefficient(6, a, b, interval) + 2.0_dp * fp * ft &
            * fixture%coefficient(2, a, b, interval) + ft**2 &
            * fixture%coefficient(1, a, b, interval))
    end function radial_normal_term

    pure function forward_normal_term(fixture, interval, a, b, derivative) &
            result(value)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative
        real(dp) :: value

        value = derivative * parallel_mode(fixture, b) &
            * (fixture%flux_p_slope &
            * fixture%coefficient(4, a, b, interval) &
            + fixture%flux_t_slope &
            * fixture%coefficient(5, a, b, interval))
    end function forward_normal_term

    pure function tangential_term(fixture, interval, a, b) result(value)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, a, b
        real(dp) :: ma, mb, na, nb, value

        ma = real(fixture%mode_m(a), dp)
        mb = real(fixture%mode_m(b), dp)
        na = real(fixture%mode_n(a), dp)
        nb = real(fixture%mode_n(b), dp)
        value = ma * mb &
            * fixture%coefficient(1, a, b, interval) &
            + na * nb &
            * fixture%coefficient(6, a, b, interval) &
            + (ma * nb + mb * na) &
            * fixture%coefficient(2, a, b, interval)
    end function tangential_term

    pure function left_tangential_term(fixture, interval, a, b, derivative) &
            result(value)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative
        real(dp) :: value, fp, ft

        fp = fixture%flux_p_slope
        ft = fixture%flux_t_slope
        value = derivative * (real(fixture%mode_m(b), dp) &
            * (fp * fixture%coefficient(2, a, b, interval) &
            + ft * fixture%coefficient(1, a, b, interval)) &
            + real(fixture%mode_n(b), dp) &
            * (fp * fixture%coefficient(6, a, b, interval) &
            + ft * fixture%coefficient(2, a, b, interval))) &
            + parallel_mode(fixture, a) &
            * (real(fixture%mode_m(b), dp) &
            * fixture%coefficient(5, b, a, interval) &
            + real(fixture%mode_n(b), dp) &
            * fixture%coefficient(4, b, a, interval))
    end function left_tangential_term

    pure function right_tangential_term(fixture, interval, a, b, derivative) &
            result(value)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval, a, b
        real(dp), intent(in) :: derivative
        real(dp) :: value, fp, ft

        fp = fixture%flux_p_slope
        ft = fixture%flux_t_slope
        value = -derivative * (real(fixture%mode_m(a), dp) &
            * (fp * fixture%coefficient(2, a, b, interval) &
            + ft * fixture%coefficient(1, a, b, interval)) &
            + real(fixture%mode_n(a), dp) &
            * (fp * fixture%coefficient(6, a, b, interval) &
            + ft * fixture%coefficient(2, a, b, interval))) &
            + parallel_mode(fixture, b) &
            * (real(fixture%mode_m(a), dp) &
            * fixture%coefficient(5, a, b, interval) &
            + real(fixture%mode_n(a), dp) &
            * fixture%coefficient(4, a, b, interval))
    end function right_tangential_term

    pure function parallel_mode(fixture, mode) result(value)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: mode
        real(dp) :: value

        value = real(fixture%mode_m(mode), dp) * fixture%flux_p_slope &
            - real(fixture%mode_n(mode), dp) * fixture%flux_t_slope
    end function parallel_mode

    pure subroutine build_element_map(fixture, interval, mapping)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        integer, intent(in) :: interval
        integer, intent(out) :: mapping(:)
        integer :: mode, modes

        modes = fixture%modes
        do mode = 1, modes
            mapping(mode) = (interval - 1) * modes + mode
            mapping(modes + mode) = interval * modes + mode
            if (interval == fixture%vacuum_intervals) &
                mapping(modes + mode) = 0
            mapping(2 * modes + mode) = modes * fixture%vacuum_intervals &
                + (interval - 1) * modes + mode
        end do
    end subroutine build_element_map

    pure subroutine add_element(element, mapping, stiffness)
        real(dp), intent(in) :: element(:, :)
        integer, intent(in) :: mapping(:)
        real(dp), intent(inout) :: stiffness(:, :)
        integer :: a, b

        do b = 1, size(mapping)
            if (mapping(b) == 0) cycle
            do a = 1, size(mapping)
                if (mapping(a) == 0) cycle
                stiffness(mapping(a), mapping(b)) = &
                    stiffness(mapping(a), mapping(b)) + element(a, b)
            end do
        end do
    end subroutine add_element

    pure subroutine set_symmetric(matrix, a, b, value)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(in) :: a, b
        real(dp), intent(in) :: value

        matrix(a, b) = value
        matrix(b, a) = value
    end subroutine set_symmetric

end module terpsichore_pseudoplasma_stiffness
