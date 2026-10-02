program test_plasma_vacuum_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        diagnose_fixed_boundary_energy, fixed_boundary_energy_terms_t, &
        fixed_boundary_invalid, fixed_boundary_is_free, fixed_boundary_ok, &
        fixed_boundary_problem_t, fixed_boundary_spectrum_result_t, &
        fixed_boundary_vacuum, fixed_boundary_vacuum_mesh, &
        fixed_boundary_wall, &
        solve_fixed_boundary_class
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use plasma_vacuum_boundary, only: build_plasma_boundary_surface, &
        plasma_vacuum_geometry_error, plasma_vacuum_model_t, &
        plasma_vacuum_ok, vacuum_wall_conformal, &
        vacuum_wall_none, vacuum_wall_surface
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    implicit none

    character(len=*), parameter :: fixture = "plasma_vacuum_boundary.nc"
    ! GPEC analytic Solov'ev family: edge toroidal flux and plasma extent.
    real(dp), parameter :: edge_flux = 0.7071679553076464_dp
    real(dp), parameter :: edge_r(2) = [0.5831_dp, 1.288_dp]
    real(dp), parameter :: edge_z = 0.5642_dp
    integer, parameter :: nu = 24, nv = 12
    integer, parameter :: mode_m(5) = [0, 1, 1, 2, 2]
    integer, parameter :: mode_n(5) = [1, -1, 1, -1, 1]
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    character(len=1024) :: directory
    integer :: info

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop "supply the fixture directory"
    call test_multiple_cover_rejected()
    call read_gvec_cas3d_file(trim(directory)//"/solovev_q1.035.nc", &
        equilibrium, info)
    call require(info == reader_ok, "fixture read failed")
    call test_edge_surface()
    call test_wall_ordering_and_energy()
    call test_rejections()
    write (*, "(a)") "PASS"

contains

    subroutine test_multiple_cover_rejected()
        type(gvec_cas3d_equilibrium_t) :: cylinder
        type(primitive_equilibrium_spline_t) :: spline
        type(fixed_boundary_problem_t) :: problem
        type(plasma_vacuum_model_t) :: model
        real(dp), allocatable :: surface(:, :, :), normal(:, :, :)
        real(dp) :: fp, ft

        ! The synthetic cylinder rotates its n=-1 frame harmonic with
        ! winding 1: the edge traces the torus twice and bounds no vacuum.
        call create_cylinder_fixture(fixture, surfaces=9)
        call read_gvec_cas3d_file(fixture, cylinder, info)
        call require(info == reader_ok, "cylinder fixture read failed")
        call fit_primitive_equilibrium(cylinder, spline, info)
        call require(info == primitive_equilibrium_ok, "spline fit failed")
        call build_plasma_boundary_surface(spline, 12, 24, surface, normal, &
            fp, ft, info)
        call require(info == plasma_vacuum_geometry_error, &
            "a doubly covered edge was accepted")
        model%nu = 12
        model%nv = 24
        call build_fixed_boundary_problem(cylinder, 5.0_dp / 3.0_dp, &
            2.0_dp, 1.0_dp, [1, 2], [1, 1], 1, problem, info, vacuum=model)
        call require(info == fixed_boundary_vacuum, &
            "a free boundary on a doubly covered edge was accepted")
    end subroutine test_multiple_cover_rejected

    subroutine test_edge_surface()
        type(primitive_equilibrium_spline_t) :: spline
        real(dp), allocatable :: surface(:, :, :), normal(:, :, :)
        real(dp) :: fp, ft, radius(nu), ring(3), outward
        integer :: i, k

        call fit_primitive_equilibrium(equilibrium, spline, info)
        call require(info == primitive_equilibrium_ok, "spline fit failed")
        call build_plasma_boundary_surface(spline, nu, nv, surface, normal, &
            fp, ft, info)
        call require(info == plasma_vacuum_ok, "edge surface failed")
        ! s is the normalized toroidal flux, so |ft| is the edge flux.
        call require(abs(abs(ft) / edge_flux - 1.0_dp) < 1.0e-6_dp, &
            "edge toroidal flux slope is wrong")
        call require(fp /= 0.0_dp, "edge poloidal flux slope vanishes")
        do k = 1, nv
            do i = 1, nu
                radius(i) = norm2(surface(1:2, i, k))
            end do
            call require(abs(minval(radius) - edge_r(1)) < 2.0e-3_dp &
                .and. abs(maxval(radius) - edge_r(2)) < 2.0e-3_dp &
                .and. abs(maxval(abs(surface(3, :, k))) - edge_z) &
                < 2.0e-2_dp, "edge is not the Solov'ev boundary")
            ! Normals point away from the geometric centre ring R = 1 m.
            do i = 1, nu
                ring(1:2) = surface(1:2, i, k) / radius(i)
                ring(3) = 0.0_dp
                outward = dot_product(normal(:, i, k), surface(:, i, k) - ring)
                call require(outward > 0.0_dp, "edge normal is not outward")
            end do
        end do
    end subroutine test_edge_surface

    subroutine test_wall_ordering_and_energy()
        real(dp), parameter :: distances(2) = [0.1_dp, 0.03_dp]
        type(fixed_boundary_spectrum_result_t) :: fixed, free, walled(2)
        type(fixed_boundary_energy_terms_t) :: energy
        type(plasma_vacuum_model_t) :: model
        integer :: wall

        call solve(model, fixed, .false.)
        model%nu = nu
        model%nv = nv
        model%wall_kind = vacuum_wall_none
        call solve(model, free, .true., energy)
        model%wall_kind = vacuum_wall_conformal
        do wall = 1, size(distances)
            model%wall_distance = distances(wall)
            call solve(model, walled(wall), .true.)
        end do
        ! The fixed-boundary space is the free space with a vanishing edge,
        ! and the vacuum energy is nonnegative and grows as the wall
        ! approaches: min-max orders the lowest eigenvalues.
        call require(free%lowest_eigenvalue <= walled(1)%lowest_eigenvalue, &
            "a wall lowered the free-boundary eigenvalue")
        do wall = 2, size(distances)
            call require(walled(wall - 1)%lowest_eigenvalue &
                <= walled(wall)%lowest_eigenvalue, &
                "a closer wall lowered the eigenvalue")
        end do
        call require(walled(size(distances))%lowest_eigenvalue &
            <= fixed%lowest_eigenvalue * (1.0_dp + 1.0e-10_dp), &
            "the free boundary exceeds the fixed boundary")
        call require(free%lowest_eigenvalue < fixed%lowest_eigenvalue, &
            "the edge displacement does not lower the eigenvalue")
        ! The external kink of this q0 = 1.035 edge needs a close wall.
        call require(free%negative_count > fixed%negative_count &
            .and. walled(2)%negative_count == fixed%negative_count, &
            "the wall does not stabilize the external kink")
        call require(free%unknowns > fixed%unknowns, &
            "the free boundary did not retain the edge coefficients")
        call require(energy%vacuum_energy > 0.0_dp, &
            "the free-boundary mode has no vacuum energy")
        call require(energy%closure_error <= energy%closure_tolerance, &
            "the energy decomposition does not close")
    end subroutine test_wall_ordering_and_energy

    subroutine test_rejections()
        type(fixed_boundary_problem_t) :: problem
        type(plasma_vacuum_model_t) :: model
        type(primitive_equilibrium_spline_t) :: spline
        real(dp), allocatable :: surface(:, :, :), normal(:, :, :)
        real(dp) :: fp, ft

        model%nu = 4
        model%nv = nv
        call build(model, problem, info)
        call require(info == fixed_boundary_vacuum_mesh, &
            "an edge mesh aliasing m=2 was accepted")
        model%nu = nu
        model%wall_kind = vacuum_wall_conformal
        model%wall_distance = -0.1_dp
        call build(model, problem, info)
        call require(info == fixed_boundary_invalid, &
            "a negative wall distance was accepted")
        model%wall_kind = 7
        call build(model, problem, info)
        call require(info == fixed_boundary_invalid, &
            "an unknown wall model was accepted")

        call fit_primitive_equilibrium(equilibrium, spline, info)
        call require(info == primitive_equilibrium_ok, "spline fit failed")
        call build_plasma_boundary_surface(spline, nu, nv, surface, normal, &
            fp, ft, info)
        call require(info == plasma_vacuum_ok, "edge surface failed")
        model%wall_kind = vacuum_wall_surface
        allocate (model%wall(3, nu, nv))
        model%wall = surface - 0.05_dp * normal
        call build(model, problem, info)
        call require(info == fixed_boundary_wall, &
            "a wall inside the plasma was accepted")
    end subroutine test_rejections

    subroutine solve(model, result, free, energy)
        type(plasma_vacuum_model_t), intent(in) :: model
        type(fixed_boundary_spectrum_result_t), intent(out) :: result
        logical, intent(in) :: free
        type(fixed_boundary_energy_terms_t), optional, intent(out) :: energy
        type(fixed_boundary_problem_t) :: problem

        if (free) then
            call build(model, problem, info)
        else
            call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
                1.0_dp, 1.0e-8_dp, mode_m, mode_n, 2, problem, info, 64, 32)
        end if
        call require(info == fixed_boundary_ok, "problem construction failed")
        call require(fixed_boundary_is_free(problem) .eqv. free, &
            "the boundary kind is not reported")
        call solve_fixed_boundary_class(problem, 1, result, info)
        call require(info == fixed_boundary_ok, "solve failed")
        if (present(energy)) then
            call diagnose_fixed_boundary_energy(problem, 1, &
                result%eigenvector, energy, info)
            call require(info == fixed_boundary_ok, &
                "energy decomposition failed")
        end if
    end subroutine solve

    subroutine build(model, problem, status)
        type(plasma_vacuum_model_t), intent(in) :: model
        type(fixed_boundary_problem_t), intent(out) :: problem
        integer, intent(out) :: status

        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0e-8_dp, mode_m, mode_n, 2, problem, status, 64, 32, &
            vacuum=model)
    end subroutine build

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: "//message
            error stop 1
        end if
    end subroutine require

end program test_plasma_vacuum_boundary
