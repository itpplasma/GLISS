! Coupled parity operator for equilibria without stellarator symmetry (#10).
!
! Oracle: shifting the poloidal angle origin of a symmetric export by
! theta0 leaves the equilibrium physically unchanged but stores
! every harmonic with both parities, so the file is no longer stellarator
! symmetric. The coupled operator (parity class 0: both parities of every
! mode) must reproduce the union of the two decoupled class spectra, on the
! symmetric file and on the shifted one; a grid-aligned shift keeps the
! angular quadrature invariant, so the agreement is at roundoff. The
! marginality counts add up the same way, Mercier is invariant, and the
! decoupled marginality classes refuse the shifted file by name. The
! axisymmetric Solov'ev export has degenerate classes; the three-dimensional
! QA export (Landreman-Paul 2021 low resolution, nfp = 2) does not.
program test_asymmetric_coupling
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_full_spectrum_t, fixed_boundary_is_coupled, &
        fixed_boundary_ok, fixed_boundary_problem_t, &
        solve_fixed_boundary_full_spectrum
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t
    use marginality_spectrum, only: compute_marginality_spectrum, &
        marginality_spectrum_asymmetric, marginality_spectrum_ok, &
        marginality_spectrum_result_t
    use mercier_diagnostic, only: compute_mercier, mercier_ok, &
        mercier_result_t
    implicit none

    character(len=1024) :: directory
    integer :: n_theta, n_zeta
    integer, allocatable :: mode_m(:), mode_n(:)

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop 'supply the fixture directory'
    n_theta = 32
    n_zeta = 8
    mode_m = [0, 1, 2]
    mode_n = [1, 1, 1]
    call run_case(trim(directory) // '/solovev_q1.035.nc', .true.)
    n_zeta = 16
    mode_m = [0, 1, 1, 2]
    mode_n = [1, 1, -1, 1]
    call run_case(trim(directory) // '/qa_lowres.nc', .false.)
    write (*, '(a)') 'coupled parity operator reproduces the class union'

contains

    subroutine run_case(path, with_mercier)
        character(len=*), intent(in) :: path
        logical, intent(in) :: with_mercier
        type(gvec_cas3d_equilibrium_t) :: symmetric, shifted
        real(dp), allocatable :: decoupled(:)
        integer :: status

        call read_gvec_cas3d_file(path, symmetric, status)
        call require(status == reader_ok, 'cannot read ' // path)
        shifted = symmetric
        call shift_poloidal_origin(shifted, 3.0_dp / real(n_theta, dp))
        shifted%stellarator_symmetric = .false.
        call class_union(symmetric, decoupled)
        call check_coupled(symmetric, .true., decoupled, 'symmetric file')
        call check_coupled(shifted, .false., decoupled, 'shifted file')
        ! A file flagged symmetric whose operator is not takes the coupled
        ! operator instead of dropping the coupling.
        shifted%stellarator_symmetric = .true.
        call check_coupled(shifted, .false., decoupled, 'mislabelled file')
        shifted%stellarator_symmetric = .false.
        call check_marginality(symmetric, shifted)
        if (with_mercier) call check_mercier(symmetric, shifted)
    end subroutine run_case

    subroutine shift_poloidal_origin(equilibrium, theta0)
        type(gvec_cas3d_equilibrium_t), intent(inout) :: equilibrium
        real(dp), intent(in) :: theta0

        ! f(theta + theta0) for f = c cos(phi) + s sin(phi) with the phase
        ! 2 pi m theta + ...: every stored field is a scalar under a shift
        ! of the angle origin, including the metric components.
        call shift_pair(equilibrium%mod_b, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%xhat, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%yhat, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%zhat, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%jacobian, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%g_tt, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%g_tz, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%g_zz, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%g_st, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%g_sz, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%second_form_tt, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%second_form_tz, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%second_form_zz, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%b_contravariant_theta, equilibrium%poloidal_modes, &
            theta0)
        call shift_pair(equilibrium%b_contravariant_zeta, equilibrium%poloidal_modes, &
            theta0)
    end subroutine shift_poloidal_origin

    subroutine shift_pair(pair, poloidal_modes, theta0)
        type(harmonic_pair_t), intent(inout) :: pair
        integer, intent(in) :: poloidal_modes(:)
        real(dp), intent(in) :: theta0
        real(dp) :: angle, cosine, sine
        integer :: m, n, surface

        if (.not. allocated(pair%cosine)) return
        do n = 1, size(pair%cosine, 3)
            do m = 1, size(pair%cosine, 2)
                angle = 2.0_dp * acos(-1.0_dp) &
                    * real(poloidal_modes(m), dp) * theta0
                do surface = 1, size(pair%cosine, 1)
                    cosine = pair%cosine(surface, m, n)
                    sine = pair%sine(surface, m, n)
                    pair%cosine(surface, m, n) = cosine * cos(angle) &
                        + sine * sin(angle)
                    pair%sine(surface, m, n) = sine * cos(angle) &
                        - cosine * sin(angle)
                end do
            end do
        end do
    end subroutine shift_pair

    subroutine class_union(equilibrium, eigenvalues)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), allocatable, intent(out) :: eigenvalues(:)
        type(fixed_boundary_problem_t) :: problem
        type(fixed_boundary_full_spectrum_t) :: first, second
        integer :: info

        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0_dp, mode_m, mode_n, 1, problem, info, n_theta, n_zeta)
        call require(info == fixed_boundary_ok, 'decoupled problem failed')
        call require(.not. fixed_boundary_is_coupled(problem), &
            'symmetric file was coupled')
        call solve_fixed_boundary_full_spectrum(problem, 1, first, info)
        call require(info == fixed_boundary_ok, 'class 1 spectrum failed')
        call solve_fixed_boundary_full_spectrum(problem, 2, second, info)
        call require(info == fixed_boundary_ok, 'class 2 spectrum failed')
        allocate (eigenvalues(size(first%eigenvalues) &
            + size(second%eigenvalues)))
        eigenvalues(:size(first%eigenvalues)) = first%eigenvalues
        eigenvalues(size(first%eigenvalues) + 1:) = second%eigenvalues
        call sort(eigenvalues)
    end subroutine class_union

    subroutine check_coupled(equilibrium, force, reference, label)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        logical, intent(in) :: force
        real(dp), intent(in) :: reference(:)
        character(len=*), intent(in) :: label
        type(fixed_boundary_problem_t) :: problem
        type(fixed_boundary_full_spectrum_t) :: coupled
        integer :: info

        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0_dp, mode_m, mode_n, 1, problem, info, n_theta, &
            n_zeta, coupled=force)
        call require(info == fixed_boundary_ok, label // ': build failed')
        call require(fixed_boundary_is_coupled(problem), &
            label // ': operator is not coupled')
        call solve_fixed_boundary_full_spectrum(problem, 1, coupled, info)
        call require(info /= fixed_boundary_ok, &
            label // ': coupled problem accepted parity class 1')
        call solve_fixed_boundary_full_spectrum(problem, 0, coupled, info)
        call require(info == fixed_boundary_ok, label // ': spectrum failed')
        call require(size(coupled%eigenvalues) == size(reference), &
            label // ': coupled dimension differs from the class union')
        call sort(coupled%eigenvalues)
        call require(maxval(abs(coupled%eigenvalues - reference)) &
            <= 1.0e-9_dp * maxval(abs(reference)), &
            label // ': coupled spectrum differs from the class union')
    end subroutine check_coupled

    subroutine check_marginality(symmetric, shifted)
        type(gvec_cas3d_equilibrium_t), intent(in) :: symmetric, shifted
        type(marginality_spectrum_result_t) :: first, second, coupled
        real(dp), allocatable :: powers(:)
        character(len=256) :: message
        integer :: info

        powers = merge(1.0_dp - 0.5_dp * real(mode_m, dp), 0.0_dp, mode_m > 0)
        call compute_marginality_spectrum(symmetric, mode_m, mode_n, powers, &
            1, 2, n_theta, n_zeta, .true., first, info, message)
        call require(info == marginality_spectrum_ok, 'class 1 marginality')
        call compute_marginality_spectrum(symmetric, mode_m, mode_n, powers, &
            2, 2, n_theta, n_zeta, .true., second, info, message)
        call require(info == marginality_spectrum_ok, 'class 2 marginality')
        call compute_marginality_spectrum(shifted, mode_m, mode_n, powers, &
            0, 2, n_theta, n_zeta, .true., coupled, info, message)
        call require(info == marginality_spectrum_ok, 'coupled marginality')
        call require(coupled%negative_count &
            == first%negative_count + second%negative_count, &
            'coupled marginality count differs from the class sum')
        call require(abs(coupled%lowest_eigenvalue &
            - min(first%lowest_eigenvalue, second%lowest_eigenvalue)) &
            <= 1.0e-8_dp * abs(first%lowest_eigenvalue), &
            'coupled marginality lowest eigenvalue differs')
        call compute_marginality_spectrum(shifted, mode_m, mode_n, powers, &
            1, 2, n_theta, n_zeta, .false., first, info, message)
        call require(info == marginality_spectrum_asymmetric .and. &
            index(message, 'stellarator symmetry') > 0, &
            'decoupled marginality accepted the shifted file')
    end subroutine check_marginality

    subroutine check_mercier(symmetric, shifted)
        type(gvec_cas3d_equilibrium_t), intent(in) :: symmetric, shifted
        type(mercier_result_t) :: reference, rotated
        integer :: info

        call compute_mercier(symmetric, n_theta, n_zeta, reference, info)
        call require(info == mercier_ok, 'symmetric Mercier failed')
        call compute_mercier(shifted, n_theta, n_zeta, rotated, info)
        call require(info == mercier_ok, 'shifted Mercier failed')
        call require(maxval(abs(rotated%d_mercier - reference%d_mercier)) &
            <= 1.0e-9_dp * maxval(abs(reference%d_mercier)), &
            'Mercier changes under a shift of the angle origin')
    end subroutine check_mercier

    pure subroutine sort(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: key
        integer :: i, j

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') 'FAIL: ' // message
            error stop 1
        end if
    end subroutine require

end program test_asymmetric_coupling
