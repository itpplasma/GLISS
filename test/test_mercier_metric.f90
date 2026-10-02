! A truncated metric series of a regular chart can lose positive
! definiteness. Mercier used to clip det(g) to zero, so |grad psi| vanished
! and every flux-surface integral became 0/0 = NaN without an error. The
! Solov'ev fixture is regular; adding a uniform g_tz of the size of
! sqrt(g_tt g_zz) makes its tangential metric indefinite.
program test_mercier_metric
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: compute_mercier, mercier_metric_error, &
        mercier_ok, mercier_result_t
    implicit none
    character(len=1024) :: directory
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(mercier_result_t) :: result
    integer :: info, m0, n0

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop 'supply the fixture directory'
    call read_gvec_cas3d_file(trim(directory) // '/solovev_q1.035.nc', &
        equilibrium, info)
    if (info /= reader_ok) error stop 'cannot read the Solov''ev fixture'
    call compute_mercier(equilibrium, 64, 32, result, info)
    if (info /= mercier_ok) error stop 'regular fixture was rejected'
    if (.not. all(result%d_mercier == result%d_mercier)) &
        error stop 'regular fixture has NaN Mercier terms'

    m0 = findloc(equilibrium%poloidal_modes, 0, dim=1)
    n0 = findloc(equilibrium%toroidal_modes, 0, dim=1)
    equilibrium%g_tz%cosine(:, m0, n0) = equilibrium%g_tz%cosine(:, m0, n0) &
        + 2.0_dp * sqrt(abs(equilibrium%g_tt%cosine(:, m0, n0) &
        * equilibrium%g_zz%cosine(:, m0, n0)))
    call compute_mercier(equilibrium, 64, 32, result, info)
    if (info /= mercier_metric_error) &
        error stop 'indefinite metric was not reported'
    write (*, '(a)') 'Mercier metric positivity check passed'
end program test_mercier_metric
