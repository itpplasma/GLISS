! The Pfirsch-Schlueter magnetic differential equation
!   P d(beta)/d(theta) + T d(beta)/d(zeta) = rhs - <rhs>,
!   rhs = sqrt(g) (mu0 p' + G' B^zeta + I' B^theta),
! must be solved with the forcing's own bandwidth.  A Jacobian cubed from a
! single cos(theta) shaping term has harmonics up to m=3 although the
! position table stops at m=1; truncating the solve to the position table
! left 36% of the Pfirsch-Schlueter current unresolved.
!
! Near a rational surface the exact inverse amplifies any residual resonant
! forcing by 1/(m chi' - n Phi'). With the spread of chi'/Phi' over the
! radial data cell, the solve uses the cell-averaged inverse, which stays
! bounded as the resonance is approached and agrees with the exact inverse
! away from it.
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
    real(dp), parameter :: offsets(3) = [1.0e-3_dp, 1.0e-6_dp, 1.0e-9_dp]
    real(dp) :: exact_peak(3), regular_peak(3)
    integer :: j, k, info, case

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

    ! Away from resonance a small spread changes the solution by (w/D)^2.
    call solve_beta_derivatives_modes([0, 1], [0, 1, -1], surface, theta, &
        zeta, 0.0_dp, 0.0_dp, pressure_slope, poloidal_slope, &
        toroidal_slope, beta, beta_theta, beta_zeta, info=info, &
        iota_spread=1.0e-3_dp)
    if (info /= mercier_ok) error stop 'regularized solve failed'
    residual = poloidal_slope * beta_theta + toroidal_slope * beta_zeta &
        - forcing
    if (maxval(abs(residual)) > 1.0e-5_dp * maxval(abs(forcing))) &
        error stop 'regularization changes a nonresonant solution'

    ! Approach the m=1, n=1 resonance chi' = Phi' of the cos(theta - zeta)
    ! Jacobian harmonic, whose forcing does not vanish there.
    do case = 1, size(offsets)
        surface%b_theta = toroidal_slope * (1.0_dp + offsets(case)) &
            / surface%jacobian
        call solve_beta_derivatives_modes([0, 1], [0, 1, -1], surface, &
            theta, zeta, 0.0_dp, 0.0_dp, pressure_slope, &
            toroidal_slope * (1.0_dp + offsets(case)), toroidal_slope, &
            beta, beta_theta, beta_zeta, info=info, iota_spread=1.0e-2_dp)
        if (info /= mercier_ok) error stop 'near-resonant solve failed'
        regular_peak(case) = maxval(abs(beta))
        call solve_beta_derivatives_modes([0, 1], [0, 1, -1], surface, &
            theta, zeta, 0.0_dp, 0.0_dp, pressure_slope, &
            toroidal_slope * (1.0_dp + offsets(case)), toroidal_slope, &
            beta, beta_theta, beta_zeta, info=info)
        exact_peak(case) = huge(1.0_dp)
        if (info == mercier_ok) exact_peak(case) = maxval(abs(beta))
    end do
    write (*, '(a,3es10.3)') 'exact peak ', exact_peak
    write (*, '(a,3es10.3)') 'regularized peak ', regular_peak
    if (exact_peak(2) < 100.0_dp * exact_peak(1)) &
        error stop 'exact inverse is not resonant'
    ! The cell-averaged response D/(D^2 + w^2) is at most 1/(2w) and has a
    ! finite limit at D = 0, while the exact inverse grows like 1/D.
    if (maxval(regular_peak) > 1.0e-3_dp * exact_peak(3)) &
        error stop 'regularized inverse is not bounded at resonance'
    if (abs(regular_peak(3) - regular_peak(2)) > 1.0e-2_dp * regular_peak(2)) &
        error stop 'regularized inverse has no limit at resonance'

    ! At zero iota every n=0,m/=0 harmonic has D=0 but a positive width.
    ! A nonzero cos(theta) forcing has the finite regularized solution beta=0;
    ! only (m,n)=(0,0) is the gauge.  Arithmetic mode_scale=0 must not reject
    ! this physical resonance before the regularized inverse is considered.
    do k = 1, n_zeta
        do j = 1, n_theta
            surface%jacobian(j, k) = 2.0_dp + 0.3_dp * cos(two_pi * theta(j))
        end do
    end do
    surface%b_theta = 0.0_dp
    surface%b_zeta = toroidal_slope / surface%jacobian
    call solve_beta_derivatives_modes([0, 1], [0], surface, theta, zeta, &
        0.0_dp, 0.0_dp, pressure_slope, 0.0_dp, toroidal_slope, &
        beta, beta_theta, beta_zeta, info=info, iota_spread=1.0e-2_dp)
    if (info /= mercier_ok) error stop 'regularized zero-iota forcing rejected'
    if (maxval(abs(beta)) /= 0.0_dp) error stop 'zero-iota response is not zero'
end program test_magnetic_differential_equation
