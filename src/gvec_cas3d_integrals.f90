module gvec_cas3d_integrals
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, radial_grid_half
    implicit none
    private

    integer, parameter, public :: integration_ok = 0
    integer, parameter, public :: integration_invalid_grid = 1
    integer, parameter, public :: integration_shape_mismatch = 2
    integer, parameter, public :: integration_nonfinite_input = 3
    integer, parameter, public :: integration_invalid_orientation = 4

    public :: integrate_half_mesh_volume

contains

    pure subroutine integrate_half_mesh_volume(equilibrium, &
            signed_period_volume, full_device_volume, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(out) :: signed_period_volume, full_device_volume
        integer, intent(out) :: info
        integer :: poloidal_zero, toroidal_zero

        signed_period_volume = 0.0_dp
        full_device_volume = 0.0_dp
        call validate_equilibrium(equilibrium, info)
        if (info /= integration_ok) return
        poloidal_zero = findloc(equilibrium%poloidal_modes, 0, dim=1)
        toroidal_zero = findloc(equilibrium%toroidal_modes, 0, dim=1)
        signed_period_volume = sum(equilibrium%jacobian%cosine(:, &
            poloidal_zero, toroidal_zero)) / real(size(equilibrium%s), dp)
        full_device_volume = -real(equilibrium%field_periods, dp) * &
            signed_period_volume
        info = integration_invalid_orientation
        if (signed_period_volume >= 0.0_dp) return
        if (full_device_volume <= 0.0_dp) return
        info = integration_ok
    end subroutine integrate_half_mesh_volume

    pure subroutine validate_equilibrium(equilibrium, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(out) :: info
        integer :: radial

        info = integration_invalid_grid
        if (equilibrium%radial_grid /= radial_grid_half) return
        if (equilibrium%field_periods < 1) return
        info = integration_shape_mismatch
        if (.not. allocated(equilibrium%s)) return
        if (.not. allocated(equilibrium%poloidal_modes)) return
        if (.not. allocated(equilibrium%toroidal_modes)) return
        if (.not. allocated(equilibrium%jacobian%cosine)) return
        if (count(equilibrium%poloidal_modes == 0) /= 1) return
        if (count(equilibrium%toroidal_modes == 0) /= 1) return
        if (size(equilibrium%jacobian%cosine, 1) /= size(equilibrium%s)) return
        if (size(equilibrium%jacobian%cosine, 2) /= &
            size(equilibrium%poloidal_modes)) return
        if (size(equilibrium%jacobian%cosine, 3) /= &
            size(equilibrium%toroidal_modes)) return
        info = integration_nonfinite_input
        if (.not. all(ieee_is_finite(equilibrium%s))) return
        if (.not. all(ieee_is_finite(equilibrium%jacobian%cosine))) return
        info = integration_invalid_grid
        do radial = 1, size(equilibrium%s)
            if (abs(equilibrium%s(radial) - &
                (real(radial, dp) - 0.5_dp) / &
                real(size(equilibrium%s), dp)) > 1.0e-10_dp) return
        end do
        info = integration_ok
    end subroutine validate_equilibrium

end module gvec_cas3d_integrals
