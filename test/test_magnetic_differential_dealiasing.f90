program test_magnetic_differential_dealiasing
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use export_surface_geometry, only: build_angular_grids, mercier_ok, &
        mercier_angular_alias_error, magnetic_differential_grid_is_dealiased, &
        mercier_invalid_input, solve_beta_derivatives_modes, surface_data_t, &
        two_pi, mu0
    use gvec_cas3d_types, only: harmonic_pair_t
    implicit none

    call check(64, .false.)
    call check(256, .true.)
    call check_uniform_grid()
    ! Toroidal forcing has the same cubic bandwidth requirement.
    if (magnetic_differential_grid_is_dealiased([0], [0, 24], 8, 64)) &
        error stop 'toroidal aliasing grid was admitted'
    if (.not. magnetic_differential_grid_is_dealiased([0], [0, 24], 8, 256)) &
        error stop 'dealiased toroidal grid was rejected'
    print *, 'PASS'

contains

    subroutine check(n_theta, admitted)
        integer, intent(in) :: n_theta
        logical, intent(in) :: admitted
        real(dp), parameter :: amplitude = 0.3_dp, poloidal = 0.37_dp, toroidal = 1.3_dp
        type(surface_data_t) :: surface
        type(harmonic_pair_t) :: harmonics
        real(dp), allocatable :: theta(:), zeta(:), beta(:, :), bt(:, :), bz(:, :)
        real(dp) :: expected(3), measured(3)
        integer :: j, k, info

        call build_angular_grids(n_theta, 8, theta, zeta)
        allocate (surface%jacobian(n_theta, 8), surface%b_theta(n_theta, 8), &
            surface%b_zeta(n_theta, 8))
        do k = 1, 8
            do j = 1, n_theta
                surface%jacobian(j, k) = &
                    (1.0_dp + amplitude * cos(24.0_dp * two_pi * theta(j)))**3
            end do
        end do
        surface%b_theta = poloidal / surface%jacobian
        surface%b_zeta = toroidal / surface%jacobian
        call solve_beta_derivatives_modes([0, 24], [0], surface, theta, zeta, &
            0.0_dp, 0.0_dp, 1.0_dp / mu0, poloidal, toroidal, &
            beta, bt, bz, harmonics, info)
        if (.not. admitted) then
            if (info /= mercier_angular_alias_error) &
                error stop 'aliased cubic forcing was not rejected'
            return
        end if
        if (info /= mercier_ok) error stop 'dealiased cubic forcing rejected'
        ! Independent expansion of (1+a cos 24theta)^3: only 0,24,48,72
        ! occur. On 64 points its 48 and 72 harmonics fold to false 16 and 8.
        if (abs(harmonics%sine(1, 9, 1)) > 1.0e-13_dp .or. &
            abs(harmonics%sine(1, 17, 1)) > 1.0e-13_dp) &
            error stop 'cubic forcing contains false low harmonics'
        expected(1) = (3.0_dp * amplitude + 0.75_dp * amplitude**3) &
            / (two_pi * 24.0_dp * poloidal)
        expected(2) = 1.5_dp * amplitude**2 / (two_pi * 48.0_dp * poloidal)
        expected(3) = 0.25_dp * amplitude**3 / (two_pi * 72.0_dp * poloidal)
        measured(1) = harmonics%sine(1, 25, 1)
        measured(2) = harmonics%sine(1, 49, 1)
        measured(3) = harmonics%sine(1, 73, 1)
        if (maxval(abs(measured - expected)) > 1.0e-13_dp) &
            error stop 'dealiased response disagrees with exact cubic oracle'
    end subroutine check

    subroutine check_uniform_grid()
        type(surface_data_t) :: surface
        real(dp), allocatable :: theta(:), zeta(:), beta(:, :), bt(:, :), bz(:, :)
        integer :: j, k, info

        call build_angular_grids(16, 8, theta, zeta)
        theta = modulo(theta + 0.137_dp, 1.0_dp)
        allocate (surface%jacobian(16, 8), surface%b_theta(16, 8), &
            surface%b_zeta(16, 8))
        do k = 1, 8
            do j = 1, 16
                surface%jacobian(j, k) = 2.0_dp + 0.3_dp * cos(two_pi * theta(j))
            end do
        end do
        surface%b_theta = 0.37_dp / surface%jacobian
        surface%b_zeta = 1.3_dp / surface%jacobian
        call solve_beta_derivatives_modes([0, 1], [0], surface, theta, zeta, &
            0.0_dp, 0.0_dp, 1.0_dp / mu0, 0.37_dp, 1.3_dp, &
            beta, bt, bz, info=info)
        if (info /= mercier_ok) error stop 'shifted uniform angular grid rejected'
        do j = 1, 16
            if (maxval(abs(bt(j, :) - 0.3_dp / 0.37_dp &
                * cos(two_pi * theta(j)))) > 1.0e-13_dp) &
                error stop 'shifted-grid derivative disagrees with exact forcing'
        end do
        theta(2) = theta(1)
        call solve_beta_derivatives_modes([0, 1], [0], surface, theta, zeta, &
            0.0_dp, 0.0_dp, 1.0_dp / mu0, 0.37_dp, 1.3_dp, &
            beta, bt, bz, info=info)
        if (info /= mercier_invalid_input) error stop 'repeated angular nodes admitted'
        theta(2) = theta(1) + 1.0_dp / 16.0_dp + 0.01_dp
        call solve_beta_derivatives_modes([0, 1], [0], surface, theta, zeta, &
            0.0_dp, 0.0_dp, 1.0_dp / mu0, 0.37_dp, 1.3_dp, &
            beta, bt, bz, info=info)
        if (info /= mercier_invalid_input) error stop 'nonuniform angular nodes admitted'
    end subroutine check_uniform_grid

end program test_magnetic_differential_dealiasing
