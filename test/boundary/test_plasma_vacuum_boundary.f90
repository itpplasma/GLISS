program test_plasma_vacuum_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use boundary_mesh_geometry, only: periodic_surface_is_embedded, &
        triangle_pair_is_embedded
    use cylinder_fixture, only: create_cylinder_fixture
    use compatible_three_component_problem, only: &
        build_compatible_three_component_classes, &
        build_compatible_three_component_problem, build_compatible_vacuum_energy, &
        compatible_three_component_asymmetric, compatible_three_component_ok, &
        compatible_three_component_problem_t
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        diagnose_fixed_boundary_energy, fixed_boundary_energy_terms_t, &
        fixed_boundary_asymmetric, fixed_boundary_full_spectrum_t, &
        fixed_boundary_invalid, fixed_boundary_is_coupled, &
        fixed_boundary_is_free, fixed_boundary_ok, &
        fixed_boundary_problem_t, fixed_boundary_spectrum_result_t, &
        fixed_boundary_vacuum, fixed_boundary_vacuum_mesh, &
        fixed_boundary_wall, &
        solve_fixed_boundary_class, solve_fixed_boundary_full_spectrum
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
    call test_triangle_pair_embedding()
    call test_periodic_mesh_embedding()
    call test_multiple_cover_rejected()
    call read_gvec_cas3d_file(trim(directory)//"/solovev_q1.035.nc", &
        equilibrium, info)
    call require(info == reader_ok, "fixture read failed")
    call test_edge_surface()
    call test_wall_ordering_and_energy()
    call test_wall_parity_coupling()
    call test_rejections()
    write (*, "(a)") "PASS"

contains

    subroutine test_triangle_pair_embedding()
        real(dp) :: first(3, 3), second(3, 3)

        first(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        first(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        first(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        second = first
        second(:, 3) = [0.0_dp, -1.0_dp, 0.0_dp]
        call require(triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 2, 4]), "coplanar adjacent faces were rejected")
        second(:, 3) = [0.5_dp, 0.25_dp, 0.0_dp]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 2, 4]), "overlap beyond a shared edge was accepted")
        second(3, 3) = 1.0_dp
        call require(triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 2, 4]), "a noncoplanar crease was rejected")

        second(:, 2) = [-1.0_dp, 0.0_dp, 0.0_dp]
        second(:, 3) = [0.0_dp, -1.0_dp, 0.0_dp]
        call require(triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 4, 5]), "disjoint coplanar vertex cones were rejected")
        second(:, 2) = [1.0_dp, 1.0_dp, 0.0_dp]
        second(:, 3) = [-1.0_dp, 1.0_dp, 0.0_dp]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 4, 5]), "coplanar vertex-cone overlap was accepted")
        second(:, 2) = [1.0_dp, 1.0_dp, 1.0_dp]
        second(:, 3) = [1.0_dp, 1.0_dp, -1.0_dp]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 4, 5]), "noncoplanar vertex-cone overlap was accepted")
        second(1:2, 2:3) = -second(1:2, 2:3)
        call require(triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 4, 5]), "separated noncoplanar vertex cones were rejected")
        second(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp]
        second(:, 3) = [1.0_dp, 0.0_dp, -1.0_dp]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 4, 5]), "a vertex-cone boundary-ray contact was accepted")
        second(:, 2) = [1.0_dp, 1.0_dp, epsilon(1.0_dp)]
        second(:, 3) = [1.0_dp, -1.0_dp, -epsilon(1.0_dp)]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 4, 5]), "near-coplanar overlapping cones were accepted")

        second = first
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [4, 5, 6]), "coincident faces with distinct nodes were accepted")
        second(3, :) = 1.0_dp
        call require(triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [4, 5, 6]), "separated faces were rejected")
        second(:, 1) = [0.25_dp, 0.25_dp, -1.0_dp]
        second(:, 2) = [0.25_dp, 0.25_dp, 1.0_dp]
        second(:, 3) = [0.75_dp, 0.25_dp, 0.0_dp]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [4, 5, 6]), "transversely crossing faces were accepted")
        second = 0.1_dp + 0.1_dp * first
        second(3, :) = 0.0_dp
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [4, 5, 6]), "coplanar containment was accepted")
        second = -first
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [4, 5, 6]), "contact of distinct vertex IDs was accepted")

        ! Crossing does not depend on the ratio of triangle sizes. A short
        ! edge must not be discarded by a dimensional determinant cutoff.
        first(:, 1) = [-1.0_dp, -1.0_dp, 0.0_dp]
        first(:, 2) = [1.0_dp, -1.0_dp, 0.0_dp]
        first(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        second(:, 1) = [0.0_dp, -1.0e-14_dp, -1.0e-14_dp]
        second(:, 2) = [0.0_dp, 1.0e-14_dp, -1.0e-14_dp]
        second(:, 3) = [0.0_dp, 0.0_dp, 1.0e-14_dp]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [4, 5, 6]), "a tiny triangle crossing a large face was accepted")
        first(:, 1) = 0.0_dp
        first(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        first(:, 3) = [1.0_dp, 1.0e-8_dp, 0.0_dp]
        second(:, 1) = 0.0_dp
        second(:, 2) = [1.0_dp, 0.5e-8_dp, 1.0_dp]
        second(:, 3) = [1.0_dp, 0.5e-8_dp, -1.0_dp]
        call require(.not. triangle_pair_is_embedded(first, second, [1, 2, 3], &
            [1, 4, 5]), "a crossing at a skinny shared-vertex face was accepted")
    end subroutine test_triangle_pair_embedding

    subroutine test_periodic_mesh_embedding()
        integer, parameter :: mesh_u = 24, mesh_v = 48
        real(dp) :: surface(3, mesh_u, mesh_v), changed(3, mesh_u, mesh_v)
        real(dp) :: theta, phi, radius
        integer :: i, k, orientation, original_orientation, start, finish, rate

        do k = 1, mesh_v
            phi = 2.0_dp * acos(-1.0_dp) * real(k - 1, dp) / real(mesh_v, dp)
            do i = 1, mesh_u
                theta = 2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / real(mesh_u, dp)
                radius = 3.0_dp + 0.5_dp * cos(theta)
                surface(1, i, k) = radius * cos(phi)
                surface(2, i, k) = radius * sin(phi)
                surface(3, i, k) = 0.5_dp * sin(theta)
            end do
        end do
        call system_clock(start, rate)
        call require(periodic_surface_is_embedded(surface, original_orientation), &
            "an embedded torus mesh was rejected")
        call system_clock(finish)
        write (*, '(a,f10.6,a)') 'edge embedding 2304 triangles: ', &
            real(finish - start, dp) / real(rate, dp), ' seconds'
        call require(periodic_surface_is_embedded( &
            surface(:, 1:mesh_u:3, 1:mesh_v:12), orientation), &
            "a coarse 8-by-4 embedded torus mesh was rejected")
        changed = surface
        changed(1, :, :) = changed(1, :, :) + 10.0_dp
        call require(periodic_surface_is_embedded(changed, orientation), &
            "translation changed mesh admission")
        call require(orientation == original_orientation, &
            "translation changed mesh orientation")
        changed(1, :, :) = surface(3, :, :)
        changed(2, :, :) = -surface(2, :, :)
        changed(3, :, :) = surface(1, :, :)
        call require(periodic_surface_is_embedded(changed, orientation), &
            "rotation changed mesh admission")
        call require(orientation == original_orientation, &
            "rotation changed mesh orientation")
        changed = 1.0e-6_dp * surface
        call require(periodic_surface_is_embedded(changed, orientation), &
            "uniform scaling changed mesh admission")
        changed = surface(:, mesh_u:1:-1, :)
        call require(periodic_surface_is_embedded(changed, orientation), &
            "reversed mesh orientation was rejected")
        call require(orientation == -original_orientation, &
            "reversing a parameter did not reverse the orientation")

        ! A global polynomial diffeomorphism sends the torus to a smooth
        ! embedded surface with a figure-eight x-y projection. Its two lobes
        ! remain separated in z; projected winding cannot admit this surface.
        changed(1, :, :) = surface(1, :, :)
        changed(2, :, :) = surface(3, :, :) &
            + 2.0_dp * surface(1, :, :) * surface(2, :, :) / 3.0_dp
        changed(3, :, :) = surface(2, :, :)
        call require(periodic_surface_is_embedded(changed, orientation), &
            "an embedded mesh with a figure-eight projection was rejected")
        changed = surface
        changed(:, :, mesh_v / 2 + 1:) = surface(:, :, :mesh_v / 2)
        call require(.not. periodic_surface_is_embedded(changed, orientation), &
            "a duplicate half of the periodic mesh was accepted")
    end subroutine test_periodic_mesh_embedding

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
        call test_rigid_boundary_admission(cylinder)
    end subroutine test_multiple_cover_rejected

    subroutine test_rigid_boundary_admission(cylinder)
        type(gvec_cas3d_equilibrium_t), intent(in) :: cylinder
        type(gvec_cas3d_equilibrium_t) :: single, changed
        type(primitive_equilibrium_spline_t) :: spline
        real(dp), allocatable :: reference(:, :, :), reference_normal(:, :, :)
        real(dp), allocatable :: surface(:, :, :), normal(:, :, :)
        real(dp) :: fp, ft, reference_fp, reference_ft, ring(3), radius
        integer :: zero_m, zero_n, negative_n, i, k

        single = cylinder
        single%winding = 0
        call fit_primitive_equilibrium(single, spline, info)
        call require(info == primitive_equilibrium_ok, "single-cover fit failed")
        call build_plasma_boundary_surface(spline, 12, 24, reference, &
            reference_normal, reference_fp, reference_ft, info)
        call require(info == plasma_vacuum_ok, "single-cover torus was rejected")
        do k = 1, size(reference, 3)
            do i = 1, size(reference, 2)
                radius = norm2(reference(1:2, i, k))
                ring(1:2) = 3.0_dp * reference(1:2, i, k) / radius
                ring(3) = 0.0_dp
                call require(dot_product(reference_normal(:, i, k), &
                    reference(:, i, k) - ring) > 0.0_dp, &
                    "single-cover torus normals point inward")
            end do
        end do
        zero_m = findloc(single%poloidal_modes, 0, dim=1)
        zero_n = findloc(single%toroidal_modes, 0, dim=1)
        negative_n = findloc(single%toroidal_modes, -1, dim=1)
        changed = single
        changed%xhat%cosine(:, zero_m, zero_n) = &
            changed%xhat%cosine(:, zero_m, zero_n) + 10.0_dp
        call fit_primitive_equilibrium(changed, spline, info)
        call require(info == primitive_equilibrium_ok, "translated fit failed")
        call build_plasma_boundary_surface(spline, 12, 24, surface, normal, &
            fp, ft, info)
        call require(info == plasma_vacuum_ok, "translated torus was rejected")
        surface(1, :, :) = surface(1, :, :) - 10.0_dp
        call require(maxval(abs(surface - reference)) < 1.0e-11_dp &
            .and. maxval(abs(normal - reference_normal)) < 1.0e-11_dp, &
            "translation changed boundary geometry or outward normals")
        call require(abs(fp - reference_fp) + abs(ft - reference_ft) < 1.0e-12_dp, &
            "translation changed boundary flux slopes")

        ! A proper rigid rotation puts the torus axis along x rather than z.
        changed = single
        changed%xhat = single%zhat
        changed%yhat%cosine = -single%yhat%cosine
        changed%yhat%sine = -single%yhat%sine
        changed%zhat = single%xhat
        call fit_primitive_equilibrium(changed, spline, info)
        call require(info == primitive_equilibrium_ok, "rotated fit failed")
        call build_plasma_boundary_surface(spline, 12, 24, surface, normal, &
            fp, ft, info)
        call require(info == plasma_vacuum_ok, "rotated torus was rejected")
        call require(maxval(abs(surface(1, :, :) - reference(3, :, :))) &
            + maxval(abs(surface(2, :, :) + reference(2, :, :))) &
            + maxval(abs(surface(3, :, :) - reference(1, :, :))) < 1.0e-11_dp, &
            "rotation changed boundary positions")
        call require(maxval(abs(normal(1, :, :) - reference_normal(3, :, :))) &
            + maxval(abs(normal(2, :, :) + reference_normal(2, :, :))) &
            + maxval(abs(normal(3, :, :) - reference_normal(1, :, :))) &
            < 1.0e-11_dp, "rotation changed outward normals")

        ! A spindle torus is physically self-intersecting, although its mean
        ! projected ring winds once. Retain rejection beyond double covers.
        changed = single
        changed%xhat%cosine(:, zero_m, negative_n) = 0.2_dp
        changed%yhat%sine(:, zero_m, negative_n) = 0.2_dp
        call fit_primitive_equilibrium(changed, spline, info)
        call require(info == primitive_equilibrium_ok, "spindle fit failed")
        call build_plasma_boundary_surface(spline, 12, 24, surface, normal, &
            fp, ft, info)
        call require(info == plasma_vacuum_geometry_error, &
            "a self-intersecting spindle torus was accepted")
    end subroutine test_rigid_boundary_admission

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
        call build_plasma_boundary_surface(spline, 8, 4, surface, normal, &
            fp, ft, info)
        call require(info == plasma_vacuum_ok, &
            "a coarse 8-by-4 Solov'ev edge was rejected")
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

    subroutine test_wall_parity_coupling()
        type(plasma_vacuum_model_t) :: model
        type(fixed_boundary_problem_t) :: automatic, forced, reflected
        type(compatible_three_component_problem_t) :: full, first, second
        type(fixed_boundary_full_spectrum_t) :: automatic_spectrum, forced_spectrum
        type(fixed_boundary_full_spectrum_t) :: reflected_spectrum
        type(fixed_boundary_energy_terms_t) :: energy
        real(dp), allocatable :: centered(:, :), shifted(:, :), mirror(:, :)
        real(dp), allocatable :: vector(:)
        real(dp) :: theta, phi, radius, scale, powers(size(mode_m))
        real(dp) :: boundary(2 * size(mode_m)), expected, separated, quadratic
        integer :: i, k, modes, position(2), full_index, target

        model%nu = nu
        model%nv = nv
        model%wall_kind = vacuum_wall_surface
        allocate (model%wall(3, nu, nv))
        do k = 1, nv
            phi = 2.0_dp * acos(-1.0_dp) * real(k - 1, dp) / real(nv, dp)
            do i = 1, nu
                theta = 2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / real(nu, dp)
                radius = 0.935_dp + 0.75_dp * cos(theta)
                model%wall(1, i, k) = radius * cos(phi)
                model%wall(2, i, k) = radius * sin(phi)
                model%wall(3, i, k) = 0.75_dp * sin(theta)
            end do
        end do
        call build_compatible_vacuum_energy(equilibrium, mode_m, mode_n, &
            model, centered, info)
        call require(info == compatible_three_component_ok, &
            "centered shell vacuum assembly failed")
        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0e-8_dp, mode_m, mode_n, 1, automatic, info, 64, 32, &
            vacuum=model, radial_cells=2)
        call require(info == fixed_boundary_ok, "centered shell build failed")
        call require(.not. fixed_boundary_is_coupled(automatic), &
            "symmetric centered shell unnecessarily coupled plasma parity")
        model%wall(3, :, :) = model%wall(3, :, :) + 0.04_dp
        call build_compatible_vacuum_energy(equilibrium, mode_m, mode_n, &
            model, shifted, info)
        call require(info == compatible_three_component_ok, &
            "translated enclosing shell vacuum assembly failed")
        modes = size(mode_m)
        scale = maxval(abs(shifted))
        write (*, '(a,2es14.5)') 'vacuum cross-parity relative centered/shifted: ', &
            maxval(abs(centered(:modes, modes + 1:))) / maxval(abs(centered)), &
            maxval(abs(shifted(:modes, modes + 1:))) / scale
        call require(maxval(abs(shifted(:modes, modes + 1:))) > 1.0e-4_dp * scale, &
            "translated enclosing shell did not mix parity")
        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0e-8_dp, mode_m, mode_n, 1, automatic, info, 64, 32, &
            vacuum=model, radial_cells=2)
        call require(info == fixed_boundary_ok, "translated shell build failed")
        call require(fixed_boundary_is_coupled(automatic), &
            "translated enclosing wall lost its cross-parity energy")
        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0e-8_dp, mode_m, mode_n, 1, forced, info, 64, 32, &
            coupled=.false., vacuum=model, radial_cells=2)
        call require(info == fixed_boundary_asymmetric, &
            "explicitly separated classes accepted a parity-breaking wall")

        do i = 1, modes
            powers(i) = 0.0_dp
            if (mode_m(i) > 0) powers(i) = 1.0_dp - 0.5_dp * real(mode_m(i), dp)
        end do
        call build_compatible_three_component_classes(equilibrium, &
            5.0_dp / 3.0_dp, 1.0_dp, mode_m, mode_n, powers, 1, 64, 32, &
            first, second, info, vacuum=model, vacuum_energy=shifted, radial_cells=2)
        call require(info == compatible_three_component_asymmetric, &
            "paired native assembly discarded a supplied vacuum cross block")
        call build_compatible_three_component_classes(equilibrium, &
            5.0_dp / 3.0_dp, 1.0_dp, mode_m, mode_n, powers, 1, 64, 32, &
            first, second, info, vacuum=model, radial_cells=2)
        call require(info == compatible_three_component_asymmetric, &
            "paired native assembly discarded its computed vacuum cross block")

        ! Select the largest independently assembled cross entry. A boundary
        ! displacement using both parities has the full quadratic form;
        ! the separated principal blocks would lose exactly twice this entry.
        position = 1
        do k = 1, modes
            do i = 1, modes
                if (abs(shifted(i, modes + k)) &
                    > abs(shifted(position(1), modes + position(2)))) then
                    position(1) = i
                    position(2) = k
                end if
            end do
        end do
        boundary = 0.0_dp
        boundary(position(1)) = 1.0_dp
        boundary(modes + position(2)) = 1.0_dp
        expected = quadratic_form(shifted, boundary)
        separated = shifted(position(1), position(1)) &
            + shifted(modes + position(2), modes + position(2))
        call require(abs(expected - separated) > 1.0e-4_dp * scale, &
            "cross-parity oracle does not distinguish the separated model")
        call build_compatible_three_component_problem(equilibrium, &
            5.0_dp / 3.0_dp, 1.0_dp, mode_m, mode_n, powers, 0, 1, 64, 32, &
            full, info, vacuum=model, vacuum_energy=shifted, radial_cells=2)
        call require(info == compatible_three_component_ok, &
            "explicit coupled native assembly failed")
        allocate (vector(size(full%stiffness, 1)), source=0.0_dp)
        do i = 1, 2 * modes
            full_index = (full%h1_dofs - 1) * 2 * modes + i
            target = full%axis_tie%target(full_index)
            call require(target > 0, "edge trial was removed by an axis tie")
            vector(target) = boundary(i) / full%axis_tie%factor(full_index)
        end do
        quadratic = quadratic_form(full%vacuum, vector)
        call require(abs(quadratic - expected) < 1.0e-11_dp * scale, &
            "coupled vacuum scatter lost a cross term")
        call diagnose_fixed_boundary_energy(automatic, 0, vector, energy, info)
        call require(info == fixed_boundary_ok, "automatic coupled energy failed")
        call require(abs(energy%vacuum_energy - expected) < 1.0e-11_dp * scale, &
            "automatic coupled energy differs from the full vacuum quadratic")
        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0e-8_dp, mode_m, mode_n, 1, forced, info, 64, 32, &
            coupled=.true., vacuum=model, radial_cells=2)
        call require(info == fixed_boundary_ok, "forced coupled shell build failed")
        call solve_fixed_boundary_full_spectrum(automatic, 0, automatic_spectrum, info)
        call require(info == fixed_boundary_ok, "automatic coupled spectrum failed")
        call solve_fixed_boundary_full_spectrum(forced, 0, forced_spectrum, info)
        call require(info == fixed_boundary_ok, "forced coupled spectrum failed")
        call require(maxval(abs(automatic_spectrum%eigenvalues &
            - forced_spectrum%eigenvalues)) < 1.0e-10_dp &
            * maxval(abs(forced_spectrum%eigenvalues)), &
            "automatic coupling changed the explicit coupled spectrum")

        ! Reflecting the entire enclosing wall reverses the cross block and
        ! preserves both diagonal blocks. Its full spectrum is congruent.
        model%wall(3, :, :) = model%wall(3, :, :) - 0.08_dp
        call build_compatible_vacuum_energy(equilibrium, mode_m, mode_n, &
            model, mirror, info)
        call require(info == compatible_three_component_ok, "reflected wall failed")
        call require(maxval(abs(mirror(:modes, :modes) &
            - shifted(:modes, :modes))) < 1.0e-11_dp * scale, &
            "wall reflection changed the cosine diagonal vacuum block")
        call require(maxval(abs(mirror(modes + 1:, modes + 1:) &
            - shifted(modes + 1:, modes + 1:))) < 1.0e-11_dp * scale, &
            "wall reflection changed the sine diagonal vacuum block")
        call require(maxval(abs(mirror(:modes, modes + 1:) &
            + shifted(:modes, modes + 1:))) < 1.0e-11_dp * scale, &
            "wall reflection did not reverse the vacuum cross block")
        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            1.0_dp, 1.0e-8_dp, mode_m, mode_n, 1, reflected, info, 64, 32, &
            vacuum=model, radial_cells=2)
        call require(info == fixed_boundary_ok, "reflected shell build failed")
        call solve_fixed_boundary_full_spectrum(reflected, 0, reflected_spectrum, info)
        call require(info == fixed_boundary_ok, "reflected spectrum failed")
        call require(maxval(abs(reflected_spectrum%eigenvalues &
            - automatic_spectrum%eigenvalues)) < 1.0e-10_dp &
            * maxval(abs(automatic_spectrum%eigenvalues)), &
            "wall reflection changed the complete coupled spectrum")
    end subroutine test_wall_parity_coupling

    function quadratic_form(matrix, vector) result(value)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp) :: value
        integer :: row, column

        value = 0.0_dp
        do column = 1, size(vector)
            do row = 1, size(vector)
                value = value + vector(row) * matrix(row, column) * vector(column)
            end do
        end do
    end function quadratic_form

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
