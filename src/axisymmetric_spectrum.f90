module axisymmetric_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use eigenvalue_tracking, only: certified_lowest_eigenvalue
    use family_assembly, only: family_assembly_options_t, &
        family_negative_count, surface_geometry_t
    use field_profile_identities, only: compute_field_profile_identities, &
        field_profile_identities_ok, field_profile_identity_result_t
    use gvec_cas3d_types, only: equilibrium_is_axisymmetric, &
        gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none
    private

    integer, parameter :: n_theta = 64, n_zeta = 8
    integer, parameter, public :: axisymmetric_spectrum_ok = 0
    integer, parameter, public :: axisymmetric_spectrum_invalid_input = 1
    integer, parameter, public :: axisymmetric_spectrum_compute_error = 2

    type, public :: axisymmetric_spectrum_result_t
        logical :: has_eigenpair = .false.
        integer :: field_periods = 0
        integer :: toroidal_mode = 0
        integer :: poloidal_max = 0
        integer :: mode_count = 0
        integer :: radial_surfaces = 0
        integer :: parity_class = 0
        integer :: radial_quadrature = 0
        integer :: negative_count = 0
        real(dp) :: lowest_eigenvalue = 0.0_dp
        real(dp) :: certificate = 0.0_dp
        real(dp) :: eigenpair_residual = 0.0_dp
        real(dp) :: force_balance_residual = 0.0_dp
    end type axisymmetric_spectrum_result_t

    public :: compute_axisymmetric_spectrum

contains

    subroutine compute_axisymmetric_spectrum(equilibrium, toroidal_mode, &
            poloidal_max, radial_quadrature, solve_eigenpair, result, info, &
            message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: toroidal_mode, poloidal_max, radial_quadrature
        logical, intent(in) :: solve_eigenpair
        type(axisymmetric_spectrum_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(axisymmetric_spectrum_result_t) :: candidate
        type(family_assembly_options_t) :: options
        type(field_profile_identity_result_t) :: identities
        type(surface_geometry_t), allocatable :: geometry(:)
        real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
        real(dp), allocatable :: normal_stored_power(:)
        integer, allocatable :: mode_m(:), mode_n(:)
        real(dp) :: step
        integer :: i, solver_info

        call validate_input(equilibrium, toroidal_mode, poloidal_max, &
            radial_quadrature, info, message)
        if (info /= axisymmetric_spectrum_ok) return
        call build_mode_table(toroidal_mode, poloidal_max, mode_m, mode_n, &
            normal_stored_power)
        call build_kernel_geometry(equilibrium, n_theta, n_zeta, fields, &
            drive, solver_info)
        if (solver_info /= mercier_ok) then
            info = axisymmetric_spectrum_compute_error
            message = "kernel geometry could not be built"
            return
        end if
        allocate (geometry(size(equilibrium%s)))
        do i = 1, size(geometry)
            geometry(i)%fields = fields(:, :, :, i)
            geometry(i)%drive = drive(:, :, i)
        end do
        step = 1.0_dp / real(size(geometry), dp)
        options%field_periods = 1
        options%parity_class = 1
        options%radial_space%quadrature_points = radial_quadrature
        call family_negative_count(geometry, mode_m, mode_n, step, 0.0_dp, &
            candidate%negative_count, solver_info, options, &
            normal_stored_power)
        if (solver_info /= 0) then
            info = axisymmetric_spectrum_compute_error
            message = "zero-shift inertia failed"
            return
        end if
        call compute_eigenpair(geometry, mode_m, mode_n, step, options, &
            normal_stored_power, solve_eigenpair, candidate, info, message)
        if (info /= axisymmetric_spectrum_ok) return
        call compute_field_profile_identities(equilibrium, n_theta, n_zeta, &
            identities, solver_info)
        if (solver_info /= field_profile_identities_ok) then
            info = axisymmetric_spectrum_compute_error
            message = "force-balance diagnostic failed"
            return
        end if
        candidate%has_eigenpair = solve_eigenpair
        candidate%field_periods = equilibrium%field_periods
        candidate%toroidal_mode = toroidal_mode
        candidate%poloidal_max = poloidal_max
        candidate%mode_count = size(mode_m)
        candidate%radial_surfaces = size(geometry)
        candidate%parity_class = options%parity_class
        candidate%radial_quadrature = radial_quadrature
        candidate%force_balance_residual = &
            maxval(identities%general_force_balance_deviation)
        result = candidate
        info = axisymmetric_spectrum_ok
        message = ""
    end subroutine compute_axisymmetric_spectrum

    subroutine validate_input(equilibrium, toroidal_mode, poloidal_max, &
            radial_quadrature, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: toroidal_mode, poloidal_max, radial_quadrature
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        info = axisymmetric_spectrum_invalid_input
        if (toroidal_mode <= 0) then
            message = "toroidal mode must be positive"
        else if (poloidal_max < 1) then
            message = "poloidal maximum must be positive"
        else if (radial_quadrature /= 1 .and. radial_quadrature /= 2) then
            message = "radial quadrature must be 1 or 2"
        else if (.not. equilibrium%has_chart_metric) then
            message = "equilibrium export lacks g_st/g_sz chart metrics"
        else if (equilibrium%field_periods /= 1) then
            message = "axisymmetric comparison requires N_FP=1"
        else if (.not. equilibrium_is_axisymmetric(equilibrium)) then
            message = "equilibrium contains nonaxisymmetric harmonics"
        else if (poloidal_max >= n_theta / 2) then
            message = "poloidal maximum aliases the fixed angular quadrature"
        else if (2 * poloidal_max + &
                maxval(abs(equilibrium%poloidal_modes)) >= n_theta) then
            message = "poloidal maximum aliases the fixed angular quadrature"
        else
            info = axisymmetric_spectrum_ok
            message = ""
        end if
    end subroutine validate_input

    subroutine build_mode_table(toroidal_mode, poloidal_max, mode_m, mode_n, &
            normal_stored_power)
        integer, intent(in) :: toroidal_mode, poloidal_max
        integer, allocatable, intent(out) :: mode_m(:), mode_n(:)
        real(dp), allocatable, intent(out) :: normal_stored_power(:)
        integer :: m, mode

        allocate (mode_m(2 * poloidal_max + 1))
        allocate (mode_n(2 * poloidal_max + 1))
        mode_m(1) = 0
        mode_n(1) = toroidal_mode
        mode = 1
        do m = 1, poloidal_max
            mode = mode + 1
            mode_m(mode) = m
            mode_n(mode) = -toroidal_mode
            mode = mode + 1
            mode_m(mode) = m
            mode_n(mode) = toroidal_mode
        end do
        allocate (normal_stored_power(size(mode_m)), source=0.0_dp)
        where (mode_m > 0)
            normal_stored_power = 1.0_dp - 0.5_dp * real(mode_m, dp)
        end where
    end subroutine build_mode_table

    subroutine compute_eigenpair(geometry, mode_m, mode_n, step, options, &
            normal_stored_power, solve_eigenpair, result, info, message)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: step, normal_stored_power(:)
        type(family_assembly_options_t), intent(in) :: options
        logical, intent(in) :: solve_eigenpair
        type(axisymmetric_spectrum_result_t), intent(inout) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer :: solver_info

        if (solve_eigenpair) then
            call certified_lowest_eigenvalue(geometry, mode_m, mode_n, step, &
                result%lowest_eigenvalue, result%certificate, solver_info, &
                options, normal_stored_power, result%eigenpair_residual)
            if (solver_info /= 0) then
                info = axisymmetric_spectrum_compute_error
                message = "certified eigensolve failed"
                return
            end if
        else
            result%lowest_eigenvalue = &
                ieee_value(result%lowest_eigenvalue, ieee_quiet_nan)
            result%certificate = ieee_value(result%certificate, ieee_quiet_nan)
            result%eigenpair_residual = &
                ieee_value(result%eigenpair_residual, ieee_quiet_nan)
        end if
        info = axisymmetric_spectrum_ok
        message = ""
    end subroutine compute_eigenpair

end module axisymmetric_spectrum
