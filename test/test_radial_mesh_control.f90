! Explicit radial meshes preserve the axis constraint and quadrature.
! DCON provides the independent signs of the two public Solov'ev fixtures;
! this fixed-export FEEC refinement does not qualify fresh export convergence.
program test_radial_mesh_control
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_problem_assembly_support, only: build_problem_breaks, &
        compatible_support_ok, compatible_support_invalid
    use axisymmetric_spectrum, only: compute_axisymmetric_spectrum, &
        axisymmetric_spectrum_result_t, axisymmetric_spectrum_ok, &
        axisymmetric_spectrum_invalid_input
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use marginality_spectrum, only: compute_marginality_spectrum, &
        marginality_spectrum_invalid, marginality_spectrum_result_t
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(axisymmetric_spectrum_result_t) :: result, uniform_result
    real(dp), allocatable :: breaks(:)
    real(dp) :: uniform(65), radius_mesh(65), bad(3)
    character(len=1024) :: directory, path
    character(len=256) :: message
    integer :: info, boundary, degree, fixture, expected

    call check_validation()
    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) directory = 'test/data'
    do boundary = 0, 64
        uniform(boundary + 1) = real(boundary, dp) / 64.0_dp
        radius_mesh(boundary + 1) = uniform(boundary + 1)**2
    end do
    do fixture = 1, 2
        if (fixture == 1) then
            path = trim(directory) // '/solovev_q1.035.nc'
            expected = 1
        else
            path = trim(directory) // '/solovev_q1.045.nc'
            expected = 0
        end if
        call read_gvec_cas3d_file(trim(path), equilibrium, info)
        call require(info == reader_ok, 'Solovev fixture read failed')
        if (fixture == 1) then
            call compute_axisymmetric_spectrum(equilibrium, 1, 6, 1, &
                .true., uniform_result, info, message, radial_cells=64)
            call require(info == axisymmetric_spectrum_ok, &
                'default uniform mesh failed')
            call compute_axisymmetric_spectrum(equilibrium, 1, 6, 1, &
                .true., result, info, message, radial_breaks=uniform)
            call require(info == axisymmetric_spectrum_ok, &
                'explicit uniform mesh failed')
            call require(result%negative_count == uniform_result%negative_count, &
                'explicit uniform mesh changed inertia')
            call require(abs(result%lowest_eigenvalue &
                - uniform_result%lowest_eigenvalue) < 1.0e-10_dp &
                * abs(uniform_result%lowest_eigenvalue), &
                'explicit uniform mesh changed the eigenvalue')
            bad = [0.0_dp, 0.4_dp, 0.9_dp]
            call compute_axisymmetric_spectrum(equilibrium, 1, 6, 1, &
                .false., result, info, message, radial_breaks=bad)
            call require(info == axisymmetric_spectrum_invalid_input, &
                'spectrum accepted a mesh with the wrong endpoint')
        end if
        do degree = 1, 4
            call compute_axisymmetric_spectrum(equilibrium, 1, 6, degree, &
                .true., result, info, message, radial_breaks=radius_mesh)
            call require(info == axisymmetric_spectrum_ok, &
                'radius-graded spectrum failed')
            call require(result%radial_surfaces == 64, &
                'explicit mesh cell count was not derived')
            call require(result%negative_count == expected, &
                'radius-graded inertia disagrees with DCON')
            call require((result%lowest_eigenvalue < 0.0_dp) &
                .eqv. (expected > 0), 'graded eigenvalue sign is wrong')
            call require(result%certificate < abs(result%lowest_eigenvalue), &
                'graded certificate does not resolve the sign')
            write(*,'(a,i0,a,i0,a,es12.4)') 'fixture ', fixture, &
                ' degree ', degree, ' lowest ', result%lowest_eigenvalue
        end do
    end do
    print *, 'PASS explicit radial mesh control'

contains

    subroutine check_validation()
        integer :: case_index
        real(dp) :: candidate(3), nan
        type(gvec_cas3d_equilibrium_t) :: empty
        type(marginality_spectrum_result_t) :: marginal

        call compute_marginality_spectrum(empty, [1], [1], [0.5_dp], &
            1, 2, 64, 8, .false., marginal, info, message)
        call require(info == marginality_spectrum_invalid, &
            'marginality accepted an unallocated equilibrium radial grid')

        call build_problem_breaks(4, breaks, info)
        call require(info == compatible_support_ok, 'uniform mesh failed')
        call require(size(breaks) == 5, 'uniform mesh has the wrong size')
        call require(all(breaks == [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]), &
            'default uniform breaks changed')
        candidate = [0.0_dp, 0.2_dp, 1.0_dp]
        call build_problem_breaks(16, breaks, info, radial_breaks=candidate)
        call require(info == compatible_support_ok, 'valid explicit mesh failed')
        call require(all(breaks == candidate), 'explicit breaks were changed')
        call build_problem_breaks(16, breaks, info, 2, candidate)
        call require(info == compatible_support_ok, 'matching cell count failed')
        call build_problem_breaks(16, breaks, info, 0, candidate)
        call require(info == compatible_support_ok, 'zero cell count failed')
        call build_problem_breaks(16, breaks, info, 3, candidate)
        call require(info == compatible_support_invalid, &
            'inconsistent cell count was accepted')
        call build_problem_breaks(16, breaks, info, -1, candidate)
        call require(info == compatible_support_invalid, &
            'negative cell count was accepted')
        call build_problem_breaks(16, breaks, info, radial_breaks=candidate(:2))
        call require(info == compatible_support_invalid, 'short mesh was accepted')
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        do case_index = 1, 6
            candidate = [0.0_dp, 0.2_dp, 1.0_dp]
            select case (case_index)
            case (1)
                candidate(1) = 0.01_dp
            case (2)
                candidate(3) = 0.99_dp
            case (3)
                candidate(2) = 0.0_dp
            case (4)
                candidate(2) = -0.1_dp
            case (5)
                candidate(2) = 1.1_dp
            case (6)
                candidate(2) = nan
            end select
            call build_problem_breaks(16, breaks, info, radial_breaks=candidate)
            call require(info == compatible_support_invalid, &
                'invalid explicit mesh was accepted')
        end do
    end subroutine check_validation

    subroutine require(condition, text)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: text
        if (.not. condition) then
            print *, text
            error stop 1
        end if
    end subroutine require
end program test_radial_mesh_control
