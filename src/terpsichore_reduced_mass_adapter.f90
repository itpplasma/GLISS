module terpsichore_reduced_mass_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: dynamic_family_layout_t
    use fourier_phase_kind, only: phase_sine
    use terpsichore_matrix_fixture, only: &
        terpsichore_fixed_fixture_is_valid, terpsichore_matrix_fixture_t
    use terpsichore_reduced_mass_family_assembly, only: &
        assemble_terpsichore_reduced_fixed_boundary_mass, &
        terpsichore_reduced_family_ok
    implicit none
    private

    integer, parameter, public :: terpsichore_reduced_adapter_ok = 0
    integer, parameter, public :: terpsichore_reduced_adapter_invalid = -1

    public :: assemble_terpsichore_fixture_reduced_mass

contains

    subroutine assemble_terpsichore_fixture_reduced_mass(fixture, mass, &
            layout, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        real(dp), allocatable :: normal_phase(:, :, :)
        real(dp), allocatable :: radial_factor(:, :), radial_weight(:)
        real(dp), allocatable :: tangential_phase(:, :, :)
        integer, allocatable :: parity(:)
        integer :: allocation_status

        info = terpsichore_reduced_adapter_invalid
        if (.not. terpsichore_fixed_fixture_is_valid(fixture)) return
        if (fixture%parity /= 0.0_dp) return
        call build_phase_values(fixture, normal_phase, tangential_phase, &
            allocation_status)
        if (allocation_status /= 0) return
        call build_radial_values(fixture, radial_factor, radial_weight, &
            allocation_status)
        if (allocation_status /= 0) return
        allocate (parity(fixture%modes), source=phase_sine, &
            stat=allocation_status)
        if (allocation_status /= 0) return
        call assemble_terpsichore_reduced_fixed_boundary_mass( &
            fixture%signed_bjac(:, 1:fixture%intervals), &
            fixture%flux_t_slope(:fixture%intervals), normal_phase, &
            tangential_phase, radial_factor, radial_weight, fixture%mode_m, &
            fixture%mode_n, parity, mass, layout, info)
        if (info /= terpsichore_reduced_family_ok) then
            info = terpsichore_reduced_adapter_invalid
        else
            info = terpsichore_reduced_adapter_ok
        end if
    end subroutine assemble_terpsichore_fixture_reduced_mass

    subroutine build_phase_values(fixture, normal, tangential, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: normal(:, :, :)
        real(dp), allocatable, intent(out) :: tangential(:, :, :)
        integer, intent(out) :: info
        real(dp) :: angle, theta, zeta
        integer :: interval, j, k, mode, point

        allocate (normal(fixture%modes, &
            fixture%poloidal_points * fixture%toroidal_points, &
            fixture%intervals), stat=info)
        if (info /= 0) return
        allocate (tangential(fixture%modes, &
            fixture%poloidal_points * fixture%toroidal_points, &
            fixture%intervals), stat=info)
        if (info /= 0) return
        do k = 1, fixture%toroidal_points
            zeta = 2.0_dp * acos(-1.0_dp) * real(k - 1, dp) &
                / real(fixture%toroidal_points * fixture%field_periods, dp)
            do j = 1, fixture%poloidal_points
                theta = 2.0_dp * acos(-1.0_dp) * real(j - 1, dp) &
                    / real(fixture%poloidal_points, dp)
                point = j + fixture%poloidal_points * (k - 1)
                do mode = 1, fixture%modes
                    angle = real(fixture%mode_m(mode), dp) * theta &
                        - real(fixture%mode_n(mode), dp) * zeta
                    do interval = 1, fixture%intervals
                        normal(mode, point, interval) = sin(angle)
                        tangential(mode, point, interval) = cos(angle)
                    end do
                end do
            end do
        end do
    end subroutine build_phase_values

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
