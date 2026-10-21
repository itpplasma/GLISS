module primitive_boundary_geometry
    ! Resolve the interface from independently exported s=1 data when present.
    ! Schema-0 positions retain the legacy half-grid extrapolation; flux slopes
    ! require edge data. Neither path certifies exterior equilibrium matching.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cartesian_harmonic_spline, only: cartesian_harmonic_ok, &
        cartesian_jet_grid_t, evaluate_cartesian_harmonic_spline
    use field_periodic_cartesian, only: convert_field_periodic_jet, &
        field_periodic_cartesian_ok
    use gvec_cas3d_reconstruction, only: reconstruct_harmonic_grid, &
        reconstruction_ok
    use gvec_cas3d_types, only: harmonic_pair_t
    use primitive_equilibrium_spline, only: primitive_equilibrium_spline_t
    use radial_cubic_spline, only: evaluate_radial_cubic_spline_field, &
        radial_cubic_spline_ok
    implicit none
    private

    integer, parameter, public :: primitive_boundary_ok = 0
    integer, parameter, public :: primitive_boundary_invalid = -1
    integer, parameter, public :: primitive_boundary_missing_edge_slopes = -2

    public :: evaluate_primitive_boundary

contains

    subroutine evaluate_primitive_boundary(spline, theta, zeta, jet, fp, ft, &
            info, edge_flux_slopes)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        real(dp), intent(in) :: theta(:), zeta(:)
        ! Only values and angular derivatives are defined by the edge schema.
        ! Radial entries are zero on that path and must not be used as a jet.
        type(cartesian_jet_grid_t), intent(out) :: jet
        real(dp), intent(out) :: fp, ft
        integer, intent(out) :: info
        ! Independently supplied [Phi'(1), chi'(1)] for schema-0 inputs.
        real(dp), optional, intent(in) :: edge_flux_slopes(:)
        real(dp) :: profiles(3), slopes(3), seconds(3)
        integer :: local_info

        info = primitive_boundary_invalid
        fp = 0.0_dp
        ft = 0.0_dp
        if (spline%field_periods < 1) return
        if (size(theta) < 1 .or. size(zeta) < 1) return
        if (.not. all(ieee_is_finite(theta))) return
        if (.not. all(ieee_is_finite(zeta))) return
        select case (spline%boundary%schema_version)
        case (0)
            call evaluate_cartesian_harmonic_spline(spline%radial_grid, &
                spline%position, 1.0_dp, theta, zeta, jet, local_info)
            if (local_info /= cartesian_harmonic_ok) return
            if (present(edge_flux_slopes)) then
                if (size(edge_flux_slopes) /= 2) return
                if (.not. all(ieee_is_finite(edge_flux_slopes))) return
                if (edge_flux_slopes(1) == 0.0_dp) return
                ft = -edge_flux_slopes(1)
                fp = -edge_flux_slopes(2)
            else
                if (.not. allocated(spline%radial_grid%nodes)) return
                if (maxval(spline%radial_grid%nodes) < 1.0_dp) then
                    info = primitive_boundary_missing_edge_slopes
                    return
                end if
                call evaluate_radial_cubic_spline_field(spline%radial_grid, &
                    spline%profiles, 1.0_dp, profiles, slopes, seconds, local_info)
                if (local_info /= radial_cubic_spline_ok) return
                fp = -slopes(2)
                ft = -slopes(1)
            end if
        case (1)
            if (.not. allocated(spline%position%poloidal_modes)) return
            if (.not. allocated(spline%position%toroidal_modes)) return
            if (.not. ieee_is_finite(spline%boundary%toroidal_flux_slope)) return
            if (.not. ieee_is_finite(spline%boundary%poloidal_flux_slope)) return
            if (spline%boundary%toroidal_flux_slope == 0.0_dp) return
            call allocate_surface_jet(jet, size(theta), size(zeta))
            call reconstruct_component(spline%boundary%xhat, spline, &
                theta, zeta, 1, jet, local_info)
            if (local_info /= reconstruction_ok) return
            call reconstruct_component(spline%boundary%yhat, spline, &
                theta, zeta, 2, jet, local_info)
            if (local_info /= reconstruction_ok) return
            call reconstruct_component(spline%boundary%zhat, spline, &
                theta, zeta, 3, jet, local_info)
            if (local_info /= reconstruction_ok) return
            fp = -spline%boundary%poloidal_flux_slope
            ft = -spline%boundary%toroidal_flux_slope
        case default
            return
        end select
        call convert_field_periodic_jet(zeta, spline%field_periods, &
            spline%winding, jet, local_info)
        if (local_info /= field_periodic_cartesian_ok) return
        if (.not. ieee_is_finite(fp) .or. .not. ieee_is_finite(ft)) return
        if (fp == 0.0_dp .and. ft == 0.0_dp) return
        info = primitive_boundary_ok
    end subroutine evaluate_primitive_boundary

    subroutine allocate_surface_jet(jet, nu, nv)
        type(cartesian_jet_grid_t), intent(out) :: jet
        integer, intent(in) :: nu, nv

        allocate (jet%value(nu, nv, 3), jet%poloidal(nu, nv, 3), &
            jet%toroidal(nu, nv, 3), jet%poloidal_poloidal(nu, nv, 3), &
            jet%poloidal_toroidal(nu, nv, 3), jet%toroidal_toroidal(nu, nv, 3))
        allocate (jet%radial(nu, nv, 3), jet%radial_radial(nu, nv, 3), &
            jet%radial_poloidal(nu, nv, 3), jet%radial_toroidal(nu, nv, 3))
        jet%radial = 0.0_dp
        jet%radial_radial = 0.0_dp
        jet%radial_poloidal = 0.0_dp
        jet%radial_toroidal = 0.0_dp
    end subroutine allocate_surface_jet

    subroutine reconstruct_component(pair, spline, theta, zeta, component, jet, info)
        type(harmonic_pair_t), intent(in) :: pair
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        real(dp), intent(in) :: theta(:), zeta(:)
        integer, intent(in) :: component
        type(cartesian_jet_grid_t), intent(inout) :: jet
        integer, intent(out) :: info
        real(dp), allocatable :: value(:, :), dt(:, :), dz(:, :)
        real(dp), allocatable :: dtt(:, :), dtz(:, :), dzz(:, :)

        call reconstruct_harmonic_grid(pair, 1, spline%position%poloidal_modes, &
            spline%position%toroidal_modes, theta, zeta, value, dt, dz, &
            info, dtt, dtz, dzz)
        if (info /= reconstruction_ok) return
        jet%value(:, :, component) = value
        jet%poloidal(:, :, component) = dt
        jet%toroidal(:, :, component) = dz
        jet%poloidal_poloidal(:, :, component) = dtt
        jet%poloidal_toroidal(:, :, component) = dtz
        jet%toroidal_toroidal(:, :, component) = dzz
    end subroutine reconstruct_component

end module primitive_boundary_geometry
