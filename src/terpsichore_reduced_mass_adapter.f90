module terpsichore_reduced_mass_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: add_mapped_dynamic_element, &
        dynamic_family_layout_t, dynamic_layout_ok
    use fourier_phase_kind, only: phase_sine
    use terpsichore_matrix_fixture, only: &
        terpsichore_fixed_fixture_is_valid, terpsichore_matrix_fixture_t
    use terpsichore_pair_average, only: terpsichore_pair_averages, &
        terpsichore_pair_ok
    use terpsichore_reduced_layout, only: &
        build_terpsichore_reduced_fixed_boundary_layout, &
        terpsichore_reduced_layout_ok
    use terpsichore_reduced_mass, only: add_reduced_values
    implicit none
    private

    integer, parameter, public :: terpsichore_reduced_adapter_ok = 0
    integer, parameter, public :: terpsichore_reduced_adapter_invalid = -1

    public :: assemble_terpsichore_fixture_reduced_mass

contains

    ! Fixture-driven reduced mass through the pair-average transform:
    ! per interval one |BJAC| transform on the difference/sum mode
    ! table replaces the O(points modes^2) point-pair loop.  The
    ! phase-table family assembly remains the small-size oracle
    ! (test_terpsichore_reduced_mass_adapter compares both routes).
    subroutine assemble_terpsichore_fixture_reduced_mass(fixture, mass, &
            layout, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        real(dp), allocatable :: radial_factor(:, :), radial_weight(:)
        real(dp), allocatable :: normal_normal(:, :), normal_tangent(:, :)
        real(dp), allocatable :: tangent_tangent(:, :), element(:, :)
        integer, allocatable :: element_to_global(:, :), parity(:)
        real(dp) :: normal_value, tangential_value
        integer :: allocation_status, interval, first, second, modes

        info = terpsichore_reduced_adapter_invalid
        if (.not. terpsichore_fixed_fixture_is_valid(fixture)) return
        if (fixture%parity /= 0.0_dp) return
        call build_radial_values(fixture, radial_factor, radial_weight, &
            allocation_status)
        if (allocation_status /= 0) return
        modes = fixture%modes
        allocate (parity(modes), source=phase_sine, stat=allocation_status)
        if (allocation_status /= 0) return
        call build_terpsichore_reduced_fixed_boundary_layout( &
            fixture%mode_m, fixture%mode_n, parity, fixture%intervals, &
            layout, element_to_global, info)
        if (info /= terpsichore_reduced_layout_ok) then
            info = terpsichore_reduced_adapter_invalid
            return
        end if
        allocate (mass(layout%total_unknowns, layout%total_unknowns), &
            normal_normal(modes, modes), normal_tangent(modes, modes), &
            tangent_tangent(modes, modes), element(3 * modes, 3 * modes), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = terpsichore_reduced_adapter_invalid
            return
        end if
        mass = 0.0_dp
        do interval = 1, fixture%intervals
            call terpsichore_pair_averages( &
                abs(fixture%signed_bjac(:, interval)), &
                fixture%poloidal_points, fixture%toroidal_points, &
                fixture%field_periods, fixture%mode_m, fixture%mode_n, &
                normal_normal, normal_tangent, tangent_tangent, info)
            if (info /= terpsichore_pair_ok) then
                info = terpsichore_reduced_adapter_invalid
                return
            end if
            element = 0.0_dp
            do second = 1, modes
                do first = 1, modes
                    normal_value = 0.25_dp * radial_weight(interval) &
                        * radial_factor(first, interval) &
                        * radial_factor(second, interval) &
                        * normal_normal(first, second)
                    tangential_value = 0.25_dp * radial_weight(interval) &
                        * tangent_tangent(first, second) &
                        / fixture%flux_t_slope(interval)**2
                    call add_reduced_values(element, modes, first, second, &
                        normal_value, tangential_value)
                end do
            end do
            call add_mapped_dynamic_element(element_to_global(:, interval), &
                element, mass, info)
            if (info /= dynamic_layout_ok) then
                info = terpsichore_reduced_adapter_invalid
                return
            end if
        end do
        info = terpsichore_reduced_adapter_ok
    end subroutine assemble_terpsichore_fixture_reduced_mass

    subroutine build_radial_values(fixture, factor, weight, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: factor(:, :), weight(:)
        integer, intent(out) :: info
        real(dp) :: midpoint
        integer :: interval

        allocate (factor(fixture%modes, fixture%intervals), stat=info)
        if (info /= 0) return
        allocate (weight(fixture%intervals), stat=info)
        if (info /= 0) return
        do interval = 1, fixture%intervals
            midpoint = 0.5_dp * (fixture%s(interval - 1) + fixture%s(interval))
            factor(:, interval) = midpoint**(-fixture%radial_power)
            weight(interval) = real(fixture%intervals, dp) &
                * (fixture%s(interval) - fixture%s(interval - 1))
        end do
    end subroutine build_radial_values

end module terpsichore_reduced_mass_adapter
