! The Pfirsch-Schlueter magnetic differential equation
!   P d(beta)/d(theta) + T d(beta)/d(zeta) = rhs - <rhs>,
!   rhs = sqrt(g) (mu0 p' + G' B^zeta + I' B^theta),
! must be solved with the forcing's own bandwidth.  A Jacobian cubed from a
! single cos(theta) shaping term has harmonics up to m=3 although the
! position table stops at m=1; truncating the solve to the position table
! left 36% of the Pfirsch-Schlueter current unresolved.
program test_magnetic_differential_equation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use export_surface_geometry, only: build_angular_grids, mercier_ok, mu0, &
        solve_beta_derivatives_modes, surface_data_t, two_pi
    implicit none
    integer, parameter :: n_theta = 32, n_zeta = 8
    real(dp), parameter :: poloidal_slope = 0.37_dp, toroidal_slope = 1.3_dp
    real(dp), parameter :: pressure_slope = -2.0e4_dp
    type(surface_data_t) :: surface
    real(dp), allocatable :: theta(:), zeta(:), beta(:, :), beta_theta(:, :)
    real(dp), allocatable :: beta_zeta(:, :), forcing(:, :), residual(:, :)
    integer :: j, k, info

    call build_angular_grids(n_theta, n_zeta, theta, zeta)
    allocate (surface%jacobian(n_theta, n_zeta), &
        surface%b_theta(n_theta, n_zeta), surface%b_zeta(n_theta, n_zeta))
    do k = 1, n_zeta
        do j = 1, n_theta
            surface%jacobian(j, k) = 2.0_dp * (1.0_dp + 0.3_dp &
                * cos(two_pi * theta(j)) + 0.1_dp &
                * cos(two_pi * (theta(j) - zeta(k))))**3
        end do
    end do
    surface%b_theta = poloidal_slope / surface%jacobian
    surface%b_zeta = toroidal_slope / surface%jacobian
    call solve_beta_derivatives_modes([0, 1], [0, 1, -1], surface, theta, &
        zeta, 0.0_dp, 0.0_dp, pressure_slope, poloidal_slope, &
        toroidal_slope, beta, beta_theta, beta_zeta, info=info)
    if (info /= mercier_ok) error stop 'magnetic differential equation failed'
    forcing = mu0 * pressure_slope * surface%jacobian
    forcing = forcing - sum(forcing) / real(size(forcing), dp)
    residual = poloidal_slope * beta_theta + toroidal_slope * beta_zeta &
        - forcing
    write (*, '(a,es10.3)') 'relative residual ', &
        maxval(abs(residual)) / maxval(abs(forcing))
    if (maxval(abs(residual)) > 1.0e-12_dp * maxval(abs(forcing))) &
        error stop 'magnetic differential equation is truncated'
end program test_magnetic_differential_equation
