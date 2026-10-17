module primitive_cartesian_field
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cartesian_harmonic_spline, only: cartesian_harmonic_ok, &
        cartesian_jet_grid_t, evaluate_cartesian_harmonic_spline
    use field_periodic_cartesian, only: convert_field_periodic_jet, &
        field_periodic_cartesian_ok
    use gliss_cartesian_field_generated, only: generated_cartesian_field, &
        generated_cartesian_field_frame
    use primitive_equilibrium_spline, only: primitive_equilibrium_spline_t
    use radial_cubic_spline, only: evaluate_radial_cubic_spline_field, &
        radial_cubic_spline_ok
    implicit none
    private

    integer, parameter, public :: primitive_cartesian_field_ok = 0
    integer, parameter, public :: primitive_cartesian_field_invalid = -1
    integer, parameter, public :: primitive_cartesian_field_allocation_error = -2

    ! Cartesian SI components: B in T; gradient, curl and divergence in T/m.
    ! gradient(i,j) is d B_i / d x_j, independent of chart orientation.
    type, public :: primitive_cartesian_field_point_t
        real(dp) :: position(3) = 0.0_dp
        real(dp) :: signed_jacobian = 0.0_dp
        real(dp) :: b(3) = 0.0_dp
        real(dp) :: gradient(3, 3) = 0.0_dp
        real(dp) :: curl(3) = 0.0_dp
        real(dp) :: divergence = 0.0_dp
    end type primitive_cartesian_field_point_t

    type, public :: primitive_cartesian_field_grid_t
        type(primitive_cartesian_field_point_t), allocatable :: point(:, :)
    end type primitive_cartesian_field_grid_t

    public :: build_primitive_cartesian_field_point
    public :: evaluate_primitive_cartesian_field

contains

    pure subroutine build_primitive_cartesian_field_point(position, position_s, &
            position_theta, position_zeta, position_ss, position_s_theta, &
            position_s_zeta, position_theta_theta, position_theta_zeta, &
            position_zeta_zeta, field_periods, toroidal_flux_slope, &
            poloidal_flux_slope, toroidal_flux_curvature, &
            poloidal_flux_curvature, field, info)
        real(dp), intent(in) :: position(3), position_s(3), position_theta(3)
        real(dp), intent(in) :: position_zeta(3), position_ss(3)
        real(dp), intent(in) :: position_s_theta(3), position_s_zeta(3)
        real(dp), intent(in) :: position_theta_theta(3)
        real(dp), intent(in) :: position_theta_zeta(3), position_zeta_zeta(3)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: toroidal_flux_slope, poloidal_flux_slope
        real(dp), intent(in) :: toroidal_flux_curvature, poloidal_flux_curvature
        type(primitive_cartesian_field_point_t), intent(out) :: field
        integer, intent(out) :: info
        real(dp) :: packed(32), values(16), frame(2)

        field = primitive_cartesian_field_point_t()
        info = primitive_cartesian_field_invalid
        if (field_periods < 1) return
        if (.not. all(ieee_is_finite(position))) return
        packed(1:3) = position_s
        packed(4:6) = position_theta
        packed(7:9) = position_zeta
        packed(10:12) = position_ss
        packed(13:15) = position_s_theta
        packed(16:18) = position_s_zeta
        packed(19:21) = position_theta_theta
        packed(22:24) = position_theta_zeta
        packed(25:27) = position_zeta_zeta
        packed(28) = toroidal_flux_slope
        packed(29) = poloidal_flux_slope
        packed(30) = toroidal_flux_curvature
        packed(31) = poloidal_flux_curvature
        packed(32) = 1.0_dp / real(field_periods, dp)
        if (.not. all(ieee_is_finite(packed))) return
        call generated_cartesian_field_frame(packed(1:9), frame)
        if (.not. all(ieee_is_finite(frame))) return
        if (frame(2) <= tiny(1.0_dp)) return
        if (abs(frame(1)) <= 128.0_dp * epsilon(1.0_dp) * frame(2)) return
        call generated_cartesian_field(packed, values)
        if (.not. all(ieee_is_finite(values))) return
        field%position = position
        field%signed_jacobian = frame(1)
        field%b = values(1:3)
        field%gradient(:, 1) = values(4:6)
        field%gradient(:, 2) = values(7:9)
        field%gradient(:, 3) = values(10:12)
        field%curl = values(13:15)
        field%divergence = values(16)
        info = primitive_cartesian_field_ok
    end subroutine build_primitive_cartesian_field_point

    subroutine evaluate_primitive_cartesian_field(spline, coordinate, theta, &
            zeta_period, field, pressure, pressure_slope, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        real(dp), intent(in) :: coordinate, theta(:), zeta_period(:)
        type(primitive_cartesian_field_grid_t), intent(out) :: field
        real(dp), intent(out) :: pressure, pressure_slope
        integer, intent(out) :: info
        type(cartesian_jet_grid_t) :: jet
        real(dp) :: profiles(3), slopes(3), seconds(3)
        real(dp) :: position(3), point_jet(3, 9)
        integer :: local_info, allocation_status, point_info, point_status, j, k

        field = primitive_cartesian_field_grid_t()
        pressure = 0.0_dp
        pressure_slope = 0.0_dp
        info = primitive_cartesian_field_invalid
        if (spline%field_periods < 1) return
        call evaluate_cartesian_harmonic_spline(spline%radial_grid, &
            spline%position, coordinate, theta, zeta_period, jet, local_info)
        if (local_info /= cartesian_harmonic_ok) return
        call convert_field_periodic_jet(zeta_period, spline%field_periods, &
            spline%winding, jet, local_info)
        if (local_info /= field_periodic_cartesian_ok) return
        call evaluate_radial_cubic_spline_field(spline%radial_grid, &
            spline%profiles, coordinate, profiles, slopes, seconds, local_info)
        if (local_info /= radial_cubic_spline_ok) return
        if (.not. all(ieee_is_finite(profiles))) return
        if (.not. all(ieee_is_finite(slopes))) return
        if (.not. all(ieee_is_finite(seconds))) return
        allocate (field%point(size(theta), size(zeta_period)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = primitive_cartesian_field_allocation_error
            return
        end if
        point_status = primitive_cartesian_field_ok
        !$omp parallel do collapse(2) default(none) &
        !$omp shared(jet, field, spline, slopes, seconds) &
        !$omp private(j, k, position, point_jet, point_info) &
        !$omp reduction(min:point_status)
        do k = 1, size(jet%value, 2)
            do j = 1, size(jet%value, 1)
                position = jet%value(j, k, :)
                point_jet(:, 1) = jet%radial(j, k, :)
                point_jet(:, 2) = jet%poloidal(j, k, :)
                point_jet(:, 3) = jet%toroidal(j, k, :)
                point_jet(:, 4) = jet%radial_radial(j, k, :)
                point_jet(:, 5) = jet%radial_poloidal(j, k, :)
                point_jet(:, 6) = jet%radial_toroidal(j, k, :)
                point_jet(:, 7) = jet%poloidal_poloidal(j, k, :)
                point_jet(:, 8) = jet%poloidal_toroidal(j, k, :)
                point_jet(:, 9) = jet%toroidal_toroidal(j, k, :)
                call build_primitive_cartesian_field_point(position, &
                    point_jet(:, 1), point_jet(:, 2), point_jet(:, 3), &
                    point_jet(:, 4), point_jet(:, 5), point_jet(:, 6), &
                    point_jet(:, 7), point_jet(:, 8), point_jet(:, 9), &
                    spline%field_periods, slopes(1), slopes(2), &
                    seconds(1), seconds(2), field%point(j, k), point_info)
                point_status = min(point_status, point_info)
            end do
        end do
        !$omp end parallel do
        if (point_status /= primitive_cartesian_field_ok) then
            field = primitive_cartesian_field_grid_t()
            return
        end if
        pressure = profiles(3)
        pressure_slope = slopes(3)
        info = primitive_cartesian_field_ok
    end subroutine evaluate_primitive_cartesian_field

end module primitive_cartesian_field
