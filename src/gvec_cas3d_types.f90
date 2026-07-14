module gvec_cas3d_types
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: radial_grid_half = 1
    integer, parameter, public :: radial_grid_full = 2

    type, public :: harmonic_pair_t
        real(dp), allocatable :: cosine(:, :, :)
        real(dp), allocatable :: sine(:, :, :)
    end type harmonic_pair_t

    type, public :: gvec_cas3d_equilibrium_t
        integer :: schema_version = 0
        integer :: field_periods = 0
        integer :: winding = 0
        integer :: radial_grid = 0
        logical :: stellarator_symmetric = .false.
        logical :: has_chart_metric = .false.
        logical :: has_boozer_position_frame = .false.
        real(dp) :: beta_average = 0.0_dp
        integer, allocatable :: poloidal_modes(:)
        integer, allocatable :: toroidal_modes(:)
        real(dp), allocatable :: rho(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: pressure(:)
        real(dp), allocatable :: b_theta_average(:)
        real(dp), allocatable :: b_zeta_average(:)
        real(dp), allocatable :: toroidal_flux(:)
        real(dp), allocatable :: poloidal_flux(:)
        real(dp), allocatable :: rotational_transform(:)
        type(harmonic_pair_t) :: mod_b
        type(harmonic_pair_t) :: xhat
        type(harmonic_pair_t) :: yhat
        type(harmonic_pair_t) :: zhat
        type(harmonic_pair_t) :: jacobian
        type(harmonic_pair_t) :: g_tt
        type(harmonic_pair_t) :: g_tz
        type(harmonic_pair_t) :: g_zz
        type(harmonic_pair_t) :: g_st
        type(harmonic_pair_t) :: g_sz
        type(harmonic_pair_t) :: second_form_tt
        type(harmonic_pair_t) :: second_form_tz
        type(harmonic_pair_t) :: second_form_zz
        type(harmonic_pair_t) :: b_contravariant_theta
        type(harmonic_pair_t) :: b_contravariant_zeta
    end type gvec_cas3d_equilibrium_t

    public :: equilibrium_is_axisymmetric

contains

    pure function equilibrium_is_axisymmetric(equilibrium) result(axisymmetric)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        logical :: axisymmetric

        ! Cartesian positions rotate under toroidal symmetry and therefore
        ! contain n=+/-1 harmonics even for an axisymmetric torus.
        axisymmetric = harmonic_is_axisymmetric(equilibrium%jacobian, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%mod_b, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%g_tt, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%g_tz, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%g_zz, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%g_st, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%g_sz, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%second_form_tt, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%second_form_tz, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric(equilibrium%second_form_zz, &
            equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric( &
            equilibrium%b_contravariant_theta, equilibrium%toroidal_modes) &
            .and. harmonic_is_axisymmetric( &
            equilibrium%b_contravariant_zeta, &
            equilibrium%toroidal_modes)
    end function equilibrium_is_axisymmetric

    pure function harmonic_is_axisymmetric(harmonic, toroidal_modes) &
            result(axisymmetric)
        type(harmonic_pair_t), intent(in) :: harmonic
        integer, intent(in) :: toroidal_modes(:)
        logical :: axisymmetric
        real(dp) :: scale, tolerance
        integer :: i

        axisymmetric = .false.
        if (.not. allocated(harmonic%cosine)) return
        if (.not. allocated(harmonic%sine)) return
        if (size(harmonic%cosine, 3) /= size(toroidal_modes)) return
        if (any(shape(harmonic%cosine) /= shape(harmonic%sine))) return
        scale = max(maxval(abs(harmonic%cosine)), maxval(abs(harmonic%sine)))
        tolerance = 1.0e-8_dp * max(1.0_dp, scale)
        do i = 1, size(toroidal_modes)
            if (toroidal_modes(i) == 0) cycle
            if (maxval(abs(harmonic%cosine(:, :, i))) > tolerance) return
            if (maxval(abs(harmonic%sine(:, :, i))) > tolerance) return
        end do
        axisymmetric = .true.
    end function harmonic_is_axisymmetric

end module gvec_cas3d_types
