module primitive_kernel_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use export_surface_geometry, only: build_surface_kernel_fields, &
        mercier_ok, solve_beta_derivatives_modes, surface_data_t, &
        surface_profiles_t
    use primitive_equilibrium_spline, only: evaluate_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_geometry_grid, only: primitive_geometry_grid_t
    implicit none
    private

    integer, parameter, public :: primitive_kernel_ok = 0
    integer, parameter, public :: primitive_kernel_invalid = -1
    integer, parameter, public :: primitive_kernel_allocation_error = -2

    public :: evaluate_primitive_kernel_surface

contains

    subroutine evaluate_primitive_kernel_surface(spline, coordinate, theta, &
            zeta_period, fields, drive, info, jacobian_radial, &
            jacobian_theta, jacobian_zeta, pressure_pa, geometric_drive)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        real(dp), intent(in) :: coordinate, theta(:), zeta_period(:)
        real(dp), allocatable, intent(out) :: fields(:, :, :), drive(:, :)
        integer, intent(out) :: info
        real(dp), allocatable, optional, intent(out) :: jacobian_radial(:, :)
        real(dp), allocatable, optional, intent(out) :: jacobian_theta(:, :)
        real(dp), allocatable, optional, intent(out) :: jacobian_zeta(:, :)
        real(dp), optional, intent(out) :: pressure_pa
        real(dp), allocatable, optional, intent(out) :: geometric_drive(:, :)
        type(primitive_geometry_grid_t) :: geometry
        type(surface_data_t) :: surface
        type(surface_profiles_t) :: profiles
        real(dp) :: pressure, pressure_slope
        integer :: allocation_status, local_info

        info = primitive_kernel_invalid
        if (present(pressure_pa)) pressure_pa = 0.0_dp
        call evaluate_primitive_equilibrium(spline, coordinate, theta, &
            zeta_period, geometry, pressure, pressure_slope, local_info)
        if (local_info /= primitive_equilibrium_ok) return
        if (.not. geometry%has_radial_field_derivatives) return
        call copy_surface(geometry, surface, allocation_status)
        if (allocation_status /= 0) then
            info = primitive_kernel_allocation_error
            return
        end if
        call compute_profiles(geometry, pressure_slope, profiles, local_info)
        if (local_info /= primitive_kernel_ok) return
        allocate (fields(size(theta), size(zeta_period), 13), &
            drive(size(theta), size(zeta_period)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = primitive_kernel_allocation_error
            return
        end if
        call build_surface_kernel_fields(spline%position%poloidal_modes, &
            spline%position%toroidal_modes, .true., surface, profiles, &
            geometry%jacobian_s, theta, zeta_period, fields, drive, local_info)
        if (local_info /= mercier_ok) then
            deallocate (fields, drive)
            return
        end if
        if (.not. all(ieee_is_finite(fields)) &
            .or. .not. all(ieee_is_finite(drive))) then
            deallocate (fields, drive)
            return
        end if
        if (present(jacobian_radial)) then
            allocate (jacobian_radial, source=geometry%jacobian_s, &
                stat=allocation_status)
            if (allocation_status /= 0) then
                deallocate (fields, drive)
                call clear_optional_outputs(jacobian_radial, jacobian_theta, &
                    jacobian_zeta, geometric_drive)
                info = primitive_kernel_allocation_error
                return
            end if
        end if
        if (present(jacobian_theta)) then
            allocate (jacobian_theta, source=geometry%jacobian_theta, &
                stat=allocation_status)
            if (allocation_status /= 0) then
                deallocate (fields, drive)
                call clear_optional_outputs(jacobian_radial, jacobian_theta, &
                    jacobian_zeta, geometric_drive)
                info = primitive_kernel_allocation_error
                return
            end if
        end if
        if (present(jacobian_zeta)) then
            allocate (jacobian_zeta, source=geometry%jacobian_zeta, &
                stat=allocation_status)
            if (allocation_status /= 0) then
                deallocate (fields, drive)
                call clear_optional_outputs(jacobian_radial, jacobian_theta, &
                    jacobian_zeta, geometric_drive)
                info = primitive_kernel_allocation_error
                return
            end if
        end if
        if (present(pressure_pa)) pressure_pa = pressure
        if (present(geometric_drive)) then
            call compute_geometric_drive(spline, geometry, surface, profiles, &
                theta, zeta_period, geometric_drive, local_info)
            if (local_info /= primitive_kernel_ok) then
                deallocate (fields, drive)
                call clear_optional_outputs(jacobian_radial, jacobian_theta, &
                    jacobian_zeta, geometric_drive)
                info = local_info
                return
            end if
        end if
        info = primitive_kernel_ok
    end subroutine evaluate_primitive_kernel_surface

    subroutine clear_optional_outputs(jacobian_radial, jacobian_theta, &
            jacobian_zeta, geometric_drive)
        real(dp), allocatable, optional, intent(inout) :: jacobian_radial(:, :)
        real(dp), allocatable, optional, intent(inout) :: jacobian_theta(:, :)
        real(dp), allocatable, optional, intent(inout) :: jacobian_zeta(:, :)
        real(dp), allocatable, optional, intent(inout) :: geometric_drive(:, :)

        if (present(jacobian_radial)) then
            if (allocated(jacobian_radial)) deallocate (jacobian_radial)
        end if
        if (present(jacobian_theta)) then
            if (allocated(jacobian_theta)) deallocate (jacobian_theta)
        end if
        if (present(jacobian_zeta)) then
            if (allocated(jacobian_zeta)) deallocate (jacobian_zeta)
        end if
        if (present(geometric_drive)) then
            if (allocated(geometric_drive)) deallocate (geometric_drive)
        end if
    end subroutine clear_optional_outputs

    subroutine compute_geometric_drive(spline, geometry, surface, profiles, &
            theta, zeta, drive, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(primitive_geometry_grid_t), intent(in) :: geometry
        type(surface_data_t), intent(in) :: surface
        type(surface_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: drive(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: beta(:, :), beta_theta(:, :)
        real(dp), allocatable :: beta_zeta(:, :), current_theta(:, :)
        real(dp), allocatable :: current_zeta(:, :), grad_s2(:, :)
        real(dp), allocatable :: curvature_theta(:, :)
        real(dp), allocatable :: curvature_zeta(:, :)

        info = primitive_kernel_invalid
        call solve_beta_derivatives_modes(spline%position%poloidal_modes, &
            spline%position%toroidal_modes, surface, theta, zeta, &
            profiles%covariant_theta_slope, &
            profiles%covariant_zeta_slope, profiles%pressure_slope, &
            profiles%poloidal_slope, profiles%flux_slope, beta, beta_theta, &
            beta_zeta, info=info)
        if (info /= mercier_ok) return
        grad_s2 = (surface%g_tt * surface%g_zz - surface%g_tz**2) &
            / surface%jacobian**2
        current_theta = (beta_zeta - profiles%covariant_zeta_slope) &
            / surface%jacobian
        current_zeta = (profiles%covariant_theta_slope - beta_theta) &
            / surface%jacobian
        beta = surface%g_tt * current_theta + surface%g_tz * current_zeta
        beta_theta = surface%g_tz * current_theta &
            + surface%g_zz * current_zeta
        curvature_theta = -sqrt(grad_s2) &
            * (geometry%b_contravariant(:, :, 1) &
            * geometry%second_form(:, :, 1, 1) &
            + geometry%b_contravariant(:, :, 2) &
            * geometry%second_form(:, :, 1, 2))
        curvature_zeta = -sqrt(grad_s2) &
            * (geometry%b_contravariant(:, :, 1) &
            * geometry%second_form(:, :, 1, 2) &
            + geometry%b_contravariant(:, :, 2) &
            * geometry%second_form(:, :, 2, 2))
        allocate (drive(size(theta), size(zeta)))
        drive = 2.0_dp * (beta_theta * curvature_theta &
            - beta * curvature_zeta) / (surface%jacobian * grad_s2**2)
        if (.not. all(ieee_is_finite(drive))) then
            deallocate (drive)
            info = primitive_kernel_invalid
            return
        end if
        info = primitive_kernel_ok
    end subroutine compute_geometric_drive

    subroutine copy_surface(geometry, surface, allocation_status)
        type(primitive_geometry_grid_t), intent(in) :: geometry
        type(surface_data_t), intent(out) :: surface
        integer, intent(out) :: allocation_status
        integer :: n_theta, n_zeta

        n_theta = size(geometry%signed_jacobian, 1)
        n_zeta = size(geometry%signed_jacobian, 2)
        allocate (surface%jacobian(n_theta, n_zeta), &
            surface%g_tt(n_theta, n_zeta), &
            surface%g_tz(n_theta, n_zeta), &
            surface%g_zz(n_theta, n_zeta), &
            surface%b_theta(n_theta, n_zeta), &
            surface%b_zeta(n_theta, n_zeta), &
            surface%g_st(n_theta, n_zeta), &
            surface%g_sz(n_theta, n_zeta), &
            surface%mod_b(n_theta, n_zeta), &
            surface%area_element(n_theta, n_zeta), stat=allocation_status)
        if (allocation_status /= 0) return
        surface%jacobian = geometry%signed_jacobian
        surface%g_tt = geometry%metric(:, :, 2, 2)
        surface%g_tz = geometry%metric(:, :, 2, 3)
        surface%g_zz = geometry%metric(:, :, 3, 3)
        surface%b_theta = geometry%b_contravariant(:, :, 1)
        surface%b_zeta = geometry%b_contravariant(:, :, 2)
        surface%g_st = geometry%metric(:, :, 1, 2)
        surface%g_sz = geometry%metric(:, :, 1, 3)
        surface%mod_b = geometry%mod_b
        surface%area_element = sqrt(surface%g_tt * surface%g_zz &
            - surface%g_tz**2)
    end subroutine copy_surface

    subroutine compute_profiles(geometry, pressure_slope, profiles, info)
        type(primitive_geometry_grid_t), intent(in) :: geometry
        real(dp), intent(in) :: pressure_slope
        type(surface_profiles_t), intent(out) :: profiles
        integer, intent(out) :: info
        real(dp) :: count

        info = primitive_kernel_invalid
        count = real(size(geometry%signed_jacobian), dp)
        if (count < 1.0_dp .or. .not. ieee_is_finite(pressure_slope)) return
        profiles%flux_slope = sum(geometry%signed_jacobian &
            * geometry%b_contravariant(:, :, 2)) / count
        profiles%poloidal_slope = sum(geometry%signed_jacobian &
            * geometry%b_contravariant(:, :, 1)) / count
        profiles%flux_curvature = sum(geometry%jacobian_s &
            * geometry%b_contravariant(:, :, 2) &
            + geometry%signed_jacobian &
            * geometry%b_contravariant_radial(:, :, 2)) / count
        profiles%poloidal_curvature = sum(geometry%jacobian_s &
            * geometry%b_contravariant(:, :, 1) &
            + geometry%signed_jacobian &
            * geometry%b_contravariant_radial(:, :, 1)) / count
        profiles%covariant_theta = &
            sum(geometry%b_covariant(:, :, 1)) / count
        profiles%covariant_zeta = &
            sum(geometry%b_covariant(:, :, 2)) / count
        profiles%covariant_theta_slope = &
            sum(geometry%b_covariant_radial(:, :, 1)) / count
        profiles%covariant_zeta_slope = &
            sum(geometry%b_covariant_radial(:, :, 2)) / count
        profiles%pressure_slope = pressure_slope
        if (.not. profiles_are_finite(profiles)) return
        info = primitive_kernel_ok
    end subroutine compute_profiles

    pure function profiles_are_finite(profiles) result(finite)
        type(surface_profiles_t), intent(in) :: profiles
        logical :: finite

        finite = ieee_is_finite(profiles%flux_slope) &
            .and. ieee_is_finite(profiles%poloidal_slope) &
            .and. ieee_is_finite(profiles%flux_curvature) &
            .and. ieee_is_finite(profiles%poloidal_curvature) &
            .and. ieee_is_finite(profiles%covariant_theta) &
            .and. ieee_is_finite(profiles%covariant_zeta) &
            .and. ieee_is_finite(profiles%covariant_theta_slope) &
            .and. ieee_is_finite(profiles%covariant_zeta_slope) &
            .and. ieee_is_finite(profiles%pressure_slope)
    end function profiles_are_finite

end module primitive_kernel_geometry
