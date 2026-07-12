module compressible_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_reconstruction, only: reconstruct_harmonic_grid, &
        reconstruction_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_half
    use mercier_diagnostic, only: build_angular_grids, differentiate_pair
    implicit none
    private

    integer, parameter, public :: compressible_geometry_ok = 0
    integer, parameter, public :: compressible_geometry_invalid_input = 1
    integer, parameter, public :: compressible_geometry_reconstruction_error = 2

    public :: build_compressible_geometry

contains

    subroutine build_compressible_geometry(equilibrium, n_theta, n_zeta, &
            adiabatic_index, jacobian_radial, jacobian_theta, &
            jacobian_zeta, gamma_pressure_pa, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        real(dp), intent(in) :: adiabatic_index
        real(dp), allocatable, intent(out) :: jacobian_radial(:, :, :)
        real(dp), allocatable, intent(out) :: jacobian_theta(:, :, :)
        real(dp), allocatable, intent(out) :: jacobian_zeta(:, :, :)
        real(dp), allocatable, intent(out) :: gamma_pressure_pa(:, :, :)
        integer, intent(out) :: info
        type(harmonic_pair_t) :: jacobian_slope
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: slope_grid(:, :), value_grid(:, :)
        real(dp), allocatable :: grid_theta(:, :), grid_zeta(:, :)
        real(dp), allocatable :: discard_a(:, :), discard_b(:, :)
        real(dp) :: pressure_scale
        integer :: ns, i, rec_info

        info = compressible_geometry_invalid_input
        if (n_theta < 8 .or. n_zeta < 8) return
        if (equilibrium%radial_grid /= radial_grid_half) return
        if (.not. ieee_is_finite(adiabatic_index)) return
        if (adiabatic_index < 0.0_dp) return
        ns = size(equilibrium%s)
        if (ns < 5) return
        if (.not. all(ieee_is_finite(equilibrium%pressure))) return
        pressure_scale = maxval(equilibrium%pressure)
        if (pressure_scale <= 0.0_dp) return
        ! spline exports undershoot slightly below zero where the edge
        ! pressure vanishes; accept undershoot within a relative noise
        ! floor and clamp it out of gamma*p so the compression energy
        ! stays positive semidefinite.
        if (any(equilibrium%pressure &
            < -1.0e-3_dp * pressure_scale)) return

        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        call differentiate_pair(equilibrium%s, equilibrium%jacobian, &
            jacobian_slope)
        allocate (jacobian_radial(n_theta, n_zeta, ns))
        allocate (jacobian_theta(n_theta, n_zeta, ns))
        allocate (jacobian_zeta(n_theta, n_zeta, ns))
        allocate (gamma_pressure_pa(n_theta, n_zeta, ns))
        do i = 1, ns
            call reconstruct_harmonic_grid(jacobian_slope, i, &
                equilibrium%poloidal_modes, equilibrium%toroidal_modes, &
                theta, zeta, slope_grid, discard_a, discard_b, rec_info)
            if (rec_info /= reconstruction_ok) then
                info = compressible_geometry_reconstruction_error
                return
            end if
            call reconstruct_harmonic_grid(equilibrium%jacobian, i, &
                equilibrium%poloidal_modes, equilibrium%toroidal_modes, &
                theta, zeta, value_grid, grid_theta, grid_zeta, rec_info)
            if (rec_info /= reconstruction_ok) then
                info = compressible_geometry_reconstruction_error
                return
            end if
            jacobian_radial(:, :, i) = slope_grid
            jacobian_theta(:, :, i) = grid_theta
            jacobian_zeta(:, :, i) = grid_zeta
            gamma_pressure_pa(:, :, i) = adiabatic_index &
                * max(equilibrium%pressure(i), 0.0_dp)
        end do
        info = compressible_geometry_ok
    end subroutine build_compressible_geometry

end module compressible_geometry
