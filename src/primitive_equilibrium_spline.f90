module primitive_equilibrium_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cartesian_harmonic_spline, only: cartesian_harmonic_allocation_error, &
        cartesian_harmonic_ok, cartesian_harmonic_spline_t, &
        cartesian_jet_grid_t, evaluate_cartesian_harmonic_spline, &
        fit_cartesian_harmonic_spline
    use field_periodic_cartesian, only: convert_field_periodic_jet, &
        field_periodic_cartesian_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, radial_grid_half
    use primitive_geometry_grid, only: build_primitive_geometry_grid, &
        primitive_geometry_grid_allocation_error, primitive_geometry_grid_ok, &
        primitive_geometry_grid_t
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        evaluate_radial_cubic_spline_field, fit_radial_cubic_spline_field, &
        radial_cubic_spline_allocation_error, radial_cubic_spline_field_t, &
        radial_cubic_spline_grid_t, radial_cubic_spline_ok
    implicit none
    private

    integer, parameter, public :: primitive_equilibrium_ok = 0
    integer, parameter, public :: primitive_equilibrium_invalid = -1
    integer, parameter, public :: primitive_equilibrium_allocation_error = -2

    type, public :: primitive_equilibrium_spline_t
        integer :: field_periods = 0
        integer :: winding = 0
        type(radial_cubic_spline_grid_t) :: radial_grid
        type(cartesian_harmonic_spline_t) :: position
        type(radial_cubic_spline_field_t) :: profiles
    end type primitive_equilibrium_spline_t

    public :: evaluate_primitive_equilibrium
    public :: fit_primitive_equilibrium

contains

    subroutine fit_primitive_equilibrium(equilibrium, spline, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(primitive_equilibrium_spline_t), intent(out) :: spline
        integer, intent(out) :: info
        real(dp), allocatable :: profiles(:, :)
        integer :: allocation_status, local_info

        spline = primitive_equilibrium_spline_t()
        info = primitive_equilibrium_invalid
        if (.not. equilibrium_primitives_are_valid(equilibrium)) return
        call build_radial_cubic_spline_grid(equilibrium%s, 0.0_dp, 1.0_dp, &
            spline%radial_grid, local_info)
        if (local_info == radial_cubic_spline_allocation_error) then
            info = primitive_equilibrium_allocation_error
            return
        end if
        if (local_info /= radial_cubic_spline_ok) return
        call fit_cartesian_harmonic_spline(spline%radial_grid, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, &
            equilibrium%xhat, equilibrium%yhat, equilibrium%zhat, &
            spline%position, local_info)
        if (local_info == cartesian_harmonic_allocation_error) then
            info = primitive_equilibrium_allocation_error
            return
        end if
        if (local_info /= cartesian_harmonic_ok) return
        allocate (profiles(size(equilibrium%s), 3), stat=allocation_status)
        if (allocation_status /= 0) then
            info = primitive_equilibrium_allocation_error
            return
        end if
        profiles(:, 1) = equilibrium%toroidal_flux
        profiles(:, 2) = equilibrium%poloidal_flux
        profiles(:, 3) = equilibrium%pressure
        call fit_radial_cubic_spline_field(spline%radial_grid, profiles, &
            spline%profiles, local_info)
        if (local_info == radial_cubic_spline_allocation_error) then
            info = primitive_equilibrium_allocation_error
            return
        end if
        if (local_info /= radial_cubic_spline_ok) return
        spline%field_periods = equilibrium%field_periods
        spline%winding = equilibrium%winding
        info = primitive_equilibrium_ok
    end subroutine fit_primitive_equilibrium

    subroutine evaluate_primitive_equilibrium(spline, coordinate, theta, &
            zeta_period, geometry, pressure, pressure_slope, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        real(dp), intent(in) :: coordinate, theta(:), zeta_period(:)
        type(primitive_geometry_grid_t), intent(out) :: geometry
        real(dp), intent(out) :: pressure, pressure_slope
        integer, intent(out) :: info
        type(cartesian_jet_grid_t) :: jet
        real(dp) :: profiles(3), slopes(3), seconds(3)
        integer :: local_info

        geometry = primitive_geometry_grid_t()
        pressure = 0.0_dp
        pressure_slope = 0.0_dp
        info = primitive_equilibrium_invalid
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
        call build_primitive_geometry_grid(jet, spline%field_periods, &
            slopes(1), slopes(2), geometry, local_info)
        if (local_info == primitive_geometry_grid_allocation_error) then
            info = primitive_equilibrium_allocation_error
            geometry = primitive_geometry_grid_t()
            return
        end if
        if (local_info /= primitive_geometry_grid_ok) then
            geometry = primitive_geometry_grid_t()
            return
        end if
        pressure = profiles(3)
        pressure_slope = slopes(3)
        if (.not. ieee_is_finite(pressure) &
            .or. .not. ieee_is_finite(pressure_slope)) then
            geometry = primitive_geometry_grid_t()
            pressure = 0.0_dp
            pressure_slope = 0.0_dp
            return
        end if
        info = primitive_equilibrium_ok
    end subroutine evaluate_primitive_equilibrium

    function equilibrium_primitives_are_valid(equilibrium) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        logical :: valid
        integer :: count

        valid = equilibrium%radial_grid == radial_grid_half &
            .and. equilibrium%field_periods >= 1 &
            .and. (equilibrium%winding == 0 &
            .or. equilibrium%has_boozer_position_frame)
        if (.not. valid) return
        valid = allocated(equilibrium%s) &
            .and. allocated(equilibrium%poloidal_modes) &
            .and. allocated(equilibrium%toroidal_modes) &
            .and. allocated(equilibrium%toroidal_flux) &
            .and. allocated(equilibrium%poloidal_flux) &
            .and. allocated(equilibrium%pressure)
        if (.not. valid) return
        count = size(equilibrium%s)
        valid = count >= 4 &
            .and. size(equilibrium%toroidal_flux) == count &
            .and. size(equilibrium%poloidal_flux) == count &
            .and. size(equilibrium%pressure) == count
        if (.not. valid) return
        valid = all(ieee_is_finite(equilibrium%toroidal_flux)) &
            .and. all(ieee_is_finite(equilibrium%poloidal_flux)) &
            .and. all(ieee_is_finite(equilibrium%pressure))
    end function equilibrium_primitives_are_valid

end module primitive_equilibrium_spline
