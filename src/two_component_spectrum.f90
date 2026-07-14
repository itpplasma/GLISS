module two_component_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, &
        ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use eigenvalue_tracking, only: certified_lowest_eigenvalue
    use family_assembly, only: family_assembly_options_t, &
        family_negative_count, surface_geometry_t
    use field_profile_identities, only: compute_field_profile_identities, &
        field_profile_identities_ok, field_profile_identity_result_t
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none
    private

    integer, parameter, public :: two_component_spectrum_ok = 0
    integer, parameter, public :: two_component_spectrum_invalid = 1
    integer, parameter, public :: two_component_spectrum_compute_error = 2

    type, public :: two_component_spectrum_result_t
        logical :: has_eigenpair = .false.
        integer :: field_periods = 0
        integer :: mode_count = 0
        integer :: radial_surfaces = 0
        integer :: parity_class = 0
        integer :: radial_quadrature = 0
        integer :: negative_count = 0
        real(dp) :: lowest_eigenvalue = 0.0_dp
        real(dp) :: certificate = 0.0_dp
        real(dp) :: eigenpair_residual = 0.0_dp
        real(dp) :: force_balance_residual = 0.0_dp
    end type two_component_spectrum_result_t

    public :: compute_two_component_spectrum

contains

    subroutine compute_two_component_spectrum(equilibrium, mode_m, mode_n, &
            normal_stored_power, parity_class, radial_quadrature, n_theta, &
            n_zeta, solve_eigenpair, result, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: parity_class, radial_quadrature
        integer, intent(in) :: n_theta, n_zeta
        logical, intent(in) :: solve_eigenpair
        type(two_component_spectrum_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(two_component_spectrum_result_t) :: candidate
        type(family_assembly_options_t) :: options
        type(field_profile_identity_result_t) :: identities
        type(surface_geometry_t), allocatable :: geometry(:)
        real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
        real(dp) :: radial_step
        integer :: surface

        call validate_inputs(equilibrium, mode_m, mode_n, &
            normal_stored_power, parity_class, radial_quadrature, n_theta, &
            n_zeta, info, message)
        if (info /= two_component_spectrum_ok) return
        call build_kernel_geometry(equilibrium, n_theta, n_zeta, fields, &
            drive, info)
        if (info /= mercier_ok) then
            info = two_component_spectrum_compute_error
            message = "kernel geometry could not be built"
            return
        end if
        allocate (geometry(size(equilibrium%s)))
        do surface = 1, size(geometry)
            geometry(surface)%fields = fields(:, :, :, surface)
            geometry(surface)%drive = drive(:, :, surface)
        end do
        radial_step = 1.0_dp / real(size(geometry), dp)
        options%field_periods = equilibrium%field_periods
        options%parity_class = parity_class
        options%radial_space%quadrature_points = radial_quadrature
        call solve_assembled(geometry, mode_m, mode_n, normal_stored_power, &
            radial_step, options, solve_eigenpair, candidate, info, message)
        if (info /= two_component_spectrum_ok) return
        call compute_field_profile_identities(equilibrium, n_theta, n_zeta, &
            identities, info)
        if (info /= field_profile_identities_ok) then
            info = two_component_spectrum_compute_error
            message = "force-balance diagnostic failed"
            return
        end if
        candidate%has_eigenpair = solve_eigenpair
        candidate%field_periods = equilibrium%field_periods
        candidate%mode_count = size(mode_m)
        candidate%radial_surfaces = size(geometry)
        candidate%parity_class = parity_class
        candidate%radial_quadrature = radial_quadrature
        candidate%force_balance_residual = &
            maxval(identities%general_force_balance_deviation)
        result = candidate
        info = two_component_spectrum_ok
        message = ""
    end subroutine compute_two_component_spectrum

    subroutine solve_assembled(geometry, mode_m, mode_n, stored_power, &
            radial_step, options, solve_eigenpair, result, info, message)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:), radial_step
        type(family_assembly_options_t), intent(in) :: options
        logical, intent(in) :: solve_eigenpair
        type(two_component_spectrum_result_t), intent(inout) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        call family_negative_count(geometry, mode_m, mode_n, radial_step, &
            0.0_dp, result%negative_count, info, options, stored_power)
        if (info /= 0) then
            info = two_component_spectrum_compute_error
            message = "zero-shift inertia failed"
            return
        end if
        if (solve_eigenpair) then
            call certified_lowest_eigenvalue(geometry, mode_m, mode_n, &
                radial_step, result%lowest_eigenvalue, result%certificate, &
                info, options, stored_power, result%eigenpair_residual)
            if (info /= 0) then
                info = two_component_spectrum_compute_error
                message = "certified eigensolve failed"
                return
            end if
        else
            result%lowest_eigenvalue = &
                ieee_value(result%lowest_eigenvalue, ieee_quiet_nan)
            result%certificate = &
                ieee_value(result%certificate, ieee_quiet_nan)
            result%eigenpair_residual = &
                ieee_value(result%eigenpair_residual, ieee_quiet_nan)
        end if
        info = two_component_spectrum_ok
        message = ""
    end subroutine solve_assembled

    subroutine validate_inputs(equilibrium, mode_m, mode_n, stored_power, &
            parity_class, radial_quadrature, n_theta, n_zeta, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, radial_quadrature
        integer, intent(in) :: n_theta, n_zeta
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        info = two_component_spectrum_invalid
        if (.not. valid_mode_table(mode_m, mode_n, stored_power)) then
            message = "mode table is invalid"
        else if (parity_class < 1 .or. parity_class > 2) then
            message = "parity class must be 1 or 2"
        else if (radial_quadrature /= 1) then
            message = "radial quadrature must be midpoint (1)"
        else if (n_theta < 4 .or. n_zeta < 4) then
            message = "angular resolutions must be at least 4"
        else if (.not. equilibrium%has_chart_metric) then
            message = "equilibrium export lacks g_st/g_sz chart metrics"
        else if (equilibrium%field_periods < 1) then
            message = "field periods must be positive"
        else if (size(equilibrium%s) < 2) then
            message = "equilibrium requires at least two radial surfaces"
        else if (angular_grid_aliases(equilibrium, mode_m, mode_n, n_theta, &
                n_zeta)) then
            message = "mode table aliases the angular quadrature"
        else
            info = two_component_spectrum_ok
            message = ""
        end if
    end subroutine validate_inputs

    function valid_mode_table(mode_m, mode_n, stored_power) result(valid)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        logical :: valid
        integer :: first, second

        valid = .false.
        if (size(mode_m) < 1) return
        if (size(mode_n) /= size(mode_m)) return
        if (size(stored_power) /= size(mode_m)) return
        if (any(mode_m < 0)) return
        if (.not. all(ieee_is_finite(stored_power))) return
        do first = 1, size(mode_m)
            if (mode_m(first) == 0 .and. mode_n(first) < 0) return
            do second = 1, first - 1
                if (mode_m(first) == mode_m(second) &
                    .and. mode_n(first) == mode_n(second)) return
            end do
        end do
        valid = .true.
    end function valid_mode_table

    function angular_grid_aliases(equilibrium, mode_m, mode_n, n_theta, &
            n_zeta) result(aliases)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:), n_theta, n_zeta
        logical :: aliases
        integer(int64) :: poloidal_bandwidth, toroidal_bandwidth

        aliases = .true.
        if (.not. allocated(equilibrium%poloidal_modes)) return
        if (.not. allocated(equilibrium%toroidal_modes)) return
        if (size(equilibrium%poloidal_modes) < 1) return
        if (size(equilibrium%toroidal_modes) < 1) return
        poloidal_bandwidth = 2_int64 * &
            int(maxval(mode_m), int64) + int(maxval( &
            abs(equilibrium%poloidal_modes)), int64)
        toroidal_bandwidth = 2_int64 * &
            maxval(abs(int(mode_n, int64))) &
            + maxval(abs(int(equilibrium%toroidal_modes, int64)))
        aliases = poloidal_bandwidth >= int(n_theta, int64) &
            .or. toroidal_bandwidth >= int(n_zeta, int64)
    end function angular_grid_aliases

end module two_component_spectrum
