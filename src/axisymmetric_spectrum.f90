module axisymmetric_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_types, only: equilibrium_is_axisymmetric, &
        gvec_cas3d_equilibrium_t
    use marginality_spectrum, only: compute_marginality_spectrum, &
        marginality_spectrum_invalid, marginality_spectrum_ok, &
        marginality_spectrum_result_t
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
        integer :: degree = 0
        integer :: negative_count = 0
        real(dp) :: lowest_eigenvalue = 0.0_dp
        real(dp) :: certificate = 0.0_dp
        real(dp) :: eigenpair_residual = 0.0_dp
        real(dp) :: force_balance_residual = 0.0_dp
    end type axisymmetric_spectrum_result_t

    public :: build_axisymmetric_mode_table
    public :: compute_axisymmetric_spectrum

contains

    subroutine compute_axisymmetric_spectrum(equilibrium, toroidal_mode, &
            poloidal_max, degree, solve_eigenpair, result, info, &
            message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: toroidal_mode, poloidal_max, degree
        logical, intent(in) :: solve_eigenpair
        type(axisymmetric_spectrum_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(marginality_spectrum_result_t) :: general
        integer, allocatable :: mode_m(:), mode_n(:)
        real(dp), allocatable :: normal_stored_power(:)
        integer :: general_info

        call validate_input(equilibrium, toroidal_mode, poloidal_max, &
            degree, info, message)
        if (info /= axisymmetric_spectrum_ok) return
        call build_axisymmetric_mode_table(toroidal_mode, poloidal_max, &
            mode_m, mode_n, normal_stored_power)
        call compute_marginality_spectrum(equilibrium, mode_m, mode_n, &
            normal_stored_power, 1, degree, n_theta, n_zeta, &
            solve_eigenpair, general, general_info, message)
        if (general_info /= marginality_spectrum_ok) then
            if (general_info == marginality_spectrum_invalid) then
                info = axisymmetric_spectrum_invalid_input
            else
                info = axisymmetric_spectrum_compute_error
            end if
            return
        end if
        call assign_result(general, toroidal_mode, poloidal_max, result)
        info = axisymmetric_spectrum_ok
        message = ""
    end subroutine compute_axisymmetric_spectrum

    subroutine build_axisymmetric_mode_table(toroidal_mode, poloidal_max, &
            mode_m, mode_n, normal_stored_power)
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
        do mode = 1, size(mode_m)
            if (mode_m(mode) > 0) normal_stored_power(mode) = &
                1.0_dp - 0.5_dp * real(mode_m(mode), dp)
        end do
    end subroutine build_axisymmetric_mode_table

    subroutine validate_input(equilibrium, toroidal_mode, poloidal_max, &
            degree, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: toroidal_mode, poloidal_max, degree
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        info = axisymmetric_spectrum_invalid_input
        if (toroidal_mode <= 0) then
            message = "toroidal mode must be positive"
        else if (poloidal_max < 1) then
            message = "poloidal maximum must be positive"
        else if (degree < 1 .or. degree > 4) then
            message = "FEEC degree must be between 1 and 4"
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

    subroutine assign_result(general, toroidal_mode, poloidal_max, result)
        type(marginality_spectrum_result_t), intent(in) :: general
        integer, intent(in) :: toroidal_mode, poloidal_max
        type(axisymmetric_spectrum_result_t), intent(out) :: result

        result%has_eigenpair = general%has_eigenpair
        result%field_periods = general%field_periods
        result%toroidal_mode = toroidal_mode
        result%poloidal_max = poloidal_max
        result%mode_count = general%mode_count
        result%radial_surfaces = general%radial_surfaces
        result%parity_class = general%parity_class
        result%degree = general%degree
        result%negative_count = general%negative_count
        result%lowest_eigenvalue = general%lowest_eigenvalue
        result%certificate = general%certificate
        result%eigenpair_residual = general%eigenpair_residual
        result%force_balance_residual = general%force_balance_residual
    end subroutine assign_result

end module axisymmetric_spectrum
