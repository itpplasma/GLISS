! Toroidal regression for the axis-regular |m|=1 trial space.
!
! The fixtures are coarse (16 half-grid surfaces, M=8) exports of the GPEC
! analytic Solov'ev equilibrium (R0=1 m, a=0.33 m, elongation 1.6, F=1 T m)
! for q0=1.035 and q0=1.045, generated from public sources by
! benchmarks/solovev/public.  Independent GPEC/DCON Newcomb scans place the
! fixed-boundary n=1 marginal point at q0=1.0396 (archived bracket
! 1.039062-1.039843, converged 1.03956-1.03959).  Without the s^(-1/2)
! tangential axis factor the |m|=1 space was non-conforming and reported a
! spurious unstable mode whose eigenvalue scaled like the radial mesh width,
! so both fixtures counted one negative eigenvalue.
!
! The conforming space ties the leading |m|=1 coefficients of xi^s and eta
! at the axis (#35), so its eigenvalues bound the converged ones from above.
! On 16 surfaces degree two is too coarse for this near-marginal mode
! (lowest 1.7e2 at q0=1.035); degree three resolves it (-3.1e2, converged
! near -3.46e2 in benchmarks/solovev/public/convergence.sh).
program test_solovev_axis_regularity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use axisymmetric_spectrum, only: axisymmetric_spectrum_ok, &
        axisymmetric_spectrum_result_t, compute_axisymmetric_spectrum
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none
    character(len=1024) :: directory

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop 'supply the fixture directory'
    call check(trim(directory) // '/solovev_q1.035.nc', 1, 6, 3)
    call check(trim(directory) // '/solovev_q1.045.nc', 0, 6, 3)
    ! The wider table has near-null directions at the zero-shift inertia
    ! probe; a growth-free block factorization must still resolve it.
    call check(trim(directory) // '/solovev_q1.045.nc', 0, 10, 2)
    ! A configurable admitted grid preserves the independent DCON signs.
    call check(trim(directory) // '/solovev_q1.035.nc', 1, 6, 3, 256)
    call check(trim(directory) // '/solovev_q1.045.nc', 0, 6, 3, 256)
    write (*, '(a)') 'Solov''ev axis-regularity regression passed'

contains

    subroutine check(path, expected, poloidal_max, degree, angular_theta)
        character(len=*), intent(in) :: path
        integer, intent(in) :: expected, poloidal_max, degree
        integer, optional, intent(in) :: angular_theta
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        type(axisymmetric_spectrum_result_t) :: result
        character(len=256) :: message
        integer :: info, n_theta

        n_theta = 64
        if (present(angular_theta)) n_theta = angular_theta
        call read_gvec_cas3d_file(path, equilibrium, info)
        if (info /= reader_ok) error stop 'Solov''ev fixture read failed'
        call compute_axisymmetric_spectrum(equilibrium, 1, poloidal_max, &
            degree, &
            .true., &
            result, info, message, angular_theta=n_theta, angular_zeta=8)
        if (info /= axisymmetric_spectrum_ok) error stop 'spectrum failed'
        if (result%angular_theta /= n_theta .or. result%angular_zeta /= 8) &
            error stop 'spectrum reports the wrong angular grid'
        write (*, '(a,a,i0,a,es12.4,a,es10.2)') path, ' negative count ', &
            result%negative_count, ' lowest ', result%lowest_eigenvalue, &
            ' certificate ', result%certificate
        if (result%negative_count /= expected) &
            error stop 'n=1 stability disagrees with the DCON marginal point'
        if ((result%lowest_eigenvalue < 0.0_dp) .neqv. (expected > 0)) &
            error stop 'lowest eigenvalue sign disagrees with the inertia'
        if (result%certificate >= abs(result%lowest_eigenvalue)) &
            error stop 'certificate does not resolve the eigenvalue sign'
    end subroutine check

end program test_solovev_axis_regularity
