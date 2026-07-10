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
        integer :: field_periods = 0
        integer :: winding = 0
        integer :: radial_grid = 0
        logical :: stellarator_symmetric = .false.
        logical :: has_chart_metric = .false.
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

end module gvec_cas3d_types
