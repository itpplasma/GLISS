! The two Fourier parity classes decouple only for stellarator-symmetric
! equilibria. An equilibrium that breaks the symmetry, here the Solov'ev
! fixture with a small up-down asymmetric zhat harmonic, must be refused
! with a dedicated status instead of being solved class by class.
program test_asymmetric_rejection
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_asymmetric, fixed_boundary_ok, fixed_boundary_problem_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none
    character(len=1024) :: directory
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(fixed_boundary_problem_t) :: problem
    integer :: info, m1, n0

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop 'supply the fixture directory'
    call read_gvec_cas3d_file(trim(directory) // '/solovev_q1.045.nc', &
        equilibrium, info)
    if (info /= reader_ok) error stop 'cannot read the Solov''ev fixture'
    call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, 1.0_dp, &
        1.0_dp, [1, 2], [1, 1], 1, problem, info, 32, 8)
    if (info /= fixed_boundary_ok) error stop 'symmetric fixture was rejected'

    m1 = findloc(equilibrium%poloidal_modes, 1, dim=1)
    n0 = findloc(equilibrium%toroidal_modes, 0, dim=1)
    ! zhat is odd (sine) under stellarator symmetry; a cosine part shifts
    ! the surfaces up-down asymmetrically.
    equilibrium%zhat%cosine(:, m1, n0) = 0.05_dp &
        * equilibrium%zhat%sine(:, m1, n0)
    equilibrium%stellarator_symmetric = .false.
    call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, 1.0_dp, &
        1.0_dp, [1, 2], [1, 1], 1, problem, info, 32, 8)
    if (info /= fixed_boundary_asymmetric) &
        error stop 'asymmetric equilibrium was not refused'
    write (*, '(a)') 'asymmetric equilibrium refused'
end program test_asymmetric_rejection
