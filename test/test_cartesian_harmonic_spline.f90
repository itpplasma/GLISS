program test_cartesian_harmonic_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cartesian_harmonic_spline, only: cartesian_harmonic_invalid, &
        cartesian_harmonic_ok, cartesian_harmonic_spline_t, &
        cartesian_jet_grid_t, evaluate_cartesian_harmonic_spline, &
        fit_cartesian_harmonic_spline
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_problem_invalid, &
        compatible_problem_ok, compatible_quadrature_cas3d_midpoint, &
        compatible_cell_trace_t, &
        compatible_two_component_problem_t
    use compatible_three_component_problem, only: &
        build_compatible_three_component_problem, &
        compatible_three_component_invalid, compatible_three_component_ok, &
        compatible_three_component_problem_t
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_half
    use primitive_equilibrium_spline, only: evaluate_primitive_equilibrium, &
        fit_primitive_equilibrium, primitive_equilibrium_invalid, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_geometry_grid, only: build_primitive_geometry_grid, &
        primitive_geometry_grid_invalid, primitive_geometry_grid_ok, &
        primitive_geometry_grid_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        radial_cubic_spline_grid_t, radial_cubic_spline_ok
    use symmetric_eigensolver, only: solve_symmetric_generalized, &
        symmetric_eigensolver_ok
    use variable_block_tridiagonal, only: variable_block_ok, &
        variable_block_to_dense
    implicit none

    real(dp), parameter :: nodes(6) = [0.04_dp, 0.17_dp, 0.36_dp, &
        0.58_dp, 0.79_dp, 0.96_dp]
    integer, parameter :: modes_m(2) = [0, 1]
    integer, parameter :: modes_n(3) = [-1, 0, 1]
    real(dp), parameter :: major_radius = 3.0_dp, minor_radius = 0.7_dp
    real(dp), parameter :: theta(2) = [0.17_dp, 0.31_dp]
    real(dp), parameter :: zeta(2) = [0.23_dp, 0.44_dp]
    real(dp), parameter :: query_s = 0.36_dp
    type(radial_cubic_spline_grid_t) :: grid
    type(harmonic_pair_t) :: x, y, z
    type(cartesian_harmonic_spline_t) :: spline
    type(cartesian_jet_grid_t) :: jet
    integer :: info

    call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
    call require(info == radial_cubic_spline_ok, "grid construction failed")
    call build_torus_coefficients(x, y, z)
    call fit_cartesian_harmonic_spline(grid, modes_m, modes_n, x, y, z, &
        spline, info)
    call require(info == cartesian_harmonic_ok, "Cartesian fit failed")
    call evaluate_cartesian_harmonic_spline(grid, spline, query_s, theta, &
        zeta, jet, info)
    call require(info == cartesian_harmonic_ok, &
        "Cartesian evaluation failed")
    call check_manufactured_torus(jet)
    call check_geometry_grid(jet)
    call check_geometry_failures(jet)
    call check_equilibrium_spline(jet, x, y, z)
    call check_invalid_inputs(grid, x, y, z, spline)
    write (*, "(a)") "PASS"

contains

    subroutine build_torus_coefficients(x_pair, y_pair, z_pair)
        type(harmonic_pair_t), intent(out) :: x_pair, y_pair, z_pair
        real(dp) :: radius
        integer :: surface

        allocate (x_pair%cosine(size(nodes), 2, 3), &
            x_pair%sine(size(nodes), 2, 3), &
            y_pair%cosine(size(nodes), 2, 3), &
            y_pair%sine(size(nodes), 2, 3), &
            z_pair%cosine(size(nodes), 2, 3), &
            z_pair%sine(size(nodes), 2, 3))
        x_pair%cosine = 0.0_dp
        x_pair%sine = 0.0_dp
        y_pair%cosine = 0.0_dp
        y_pair%sine = 0.0_dp
        z_pair%cosine = 0.0_dp
        z_pair%sine = 0.0_dp
        do surface = 1, size(nodes)
            radius = minor_radius * sqrt(nodes(surface))
            x_pair%cosine(surface, 1, 3) = major_radius
            x_pair%cosine(surface, 2, 1) = 0.5_dp * radius
            x_pair%cosine(surface, 2, 3) = 0.5_dp * radius
            y_pair%sine(surface, 1, 3) = -major_radius
            y_pair%sine(surface, 2, 1) = 0.5_dp * radius
            y_pair%sine(surface, 2, 3) = -0.5_dp * radius
            z_pair%sine(surface, 2, 2) = radius
        end do
    end subroutine build_torus_coefficients

    subroutine check_manufactured_torus(actual)
        type(cartesian_jet_grid_t), intent(in) :: actual
        real(dp) :: expected(3, 10)
        integer :: j, k

        do k = 1, size(zeta)
            do j = 1, size(theta)
                call torus_jet(query_s, theta(j), zeta(k), expected)
                call require(close_vector(actual%value(j, k, :), &
                    expected(:, 1)), "position differs")
                call require(close_vector(actual%radial(j, k, :), &
                    expected(:, 2)), "radial derivative differs")
                call require(close_vector(actual%poloidal(j, k, :), &
                    expected(:, 3)), "poloidal derivative differs")
                call require(close_vector(actual%toroidal(j, k, :), &
                    expected(:, 4)), "toroidal derivative differs")
                call require(close_vector(actual%radial_radial(j, k, :), &
                    expected(:, 5)), "second radial derivative differs")
                call require(close_vector(actual%radial_poloidal(j, k, :), &
                    expected(:, 6)), "radial-poloidal derivative differs")
                call require(close_vector(actual%radial_toroidal(j, k, :), &
                    expected(:, 7)), "radial-toroidal derivative differs")
                call require(close_vector(actual%poloidal_poloidal(j, k, :), &
                    expected(:, 8)), "second poloidal derivative differs")
                call require(close_vector(actual%poloidal_toroidal(j, k, :), &
                    expected(:, 9)), "mixed angular derivative differs")
                call require(close_vector(actual%toroidal_toroidal(j, k, :), &
                    expected(:, 10)), "second toroidal derivative differs")
            end do
        end do
    end subroutine check_manufactured_torus

    subroutine check_geometry_grid(torus_jet)
        type(cartesian_jet_grid_t), intent(in) :: torus_jet
        type(primitive_geometry_grid_t) :: geometry
        real(dp), parameter :: phi_slope = 7.0_dp, chi_slope = -3.0_dp
        integer, parameter :: periods = 5
        real(dp) :: b_theta, b_zeta, expected_j, expected_js, expected_jt
        real(dp) :: expected_metric(3, 3), expected_metric_s(3, 3)
        real(dp) :: expected_second(2, 2), expected_contravariant(2)
        real(dp) :: expected_covariant(2)
        real(dp) :: expected_b, expected_b_s(2), expected_covariant_s(2)
        real(dp) :: angle, radius, surface_radius, two_pi
        integer :: j, k, status

        call build_primitive_geometry_grid(torus_jet, periods, phi_slope, &
            chi_slope, geometry, status, 0.0_dp, 0.0_dp)
        call require(status == primitive_geometry_grid_ok, &
            "primitive geometry grid failed")
        two_pi = 2.0_dp * acos(-1.0_dp)
        radius = minor_radius * sqrt(query_s)
        do k = 1, size(zeta)
            do j = 1, size(theta)
                angle = two_pi * theta(j)
                surface_radius = major_radius + radius * cos(angle)
                expected_metric = 0.0_dp
                expected_metric(1, 1) = minor_radius**2 / (4.0_dp * query_s)
                expected_metric(2, 2) = (two_pi * radius)**2
                expected_metric(3, 3) = (two_pi * surface_radius)**2
                expected_metric_s = 0.0_dp
                expected_metric_s(1, 1) = &
                    -minor_radius**2 / (4.0_dp * query_s**2)
                expected_metric_s(2, 2) = (two_pi * minor_radius)**2
                expected_metric_s(3, 3) = 4.0_dp * acos(-1.0_dp)**2 &
                    * surface_radius * minor_radius * cos(angle) &
                    / sqrt(query_s)
                expected_j = -2.0_dp * acos(-1.0_dp)**2 &
                    * minor_radius**2 * surface_radius
                expected_js = -acos(-1.0_dp)**2 * minor_radius**3 &
                    * cos(angle) / sqrt(query_s)
                expected_jt = 4.0_dp * acos(-1.0_dp)**3 &
                    * minor_radius**3 * sqrt(query_s) * sin(angle)
                b_theta = -chi_slope / (real(periods, dp) * expected_j)
                b_zeta = -phi_slope / expected_j
                expected_contravariant(1) = b_theta
                expected_contravariant(2) = b_zeta
                expected_covariant(1) = expected_metric(2, 2) * b_theta
                expected_covariant(2) = expected_metric(3, 3) * b_zeta
                expected_b_s(1) = -b_theta * expected_js / expected_j
                expected_b_s(2) = -b_zeta * expected_js / expected_j
                expected_covariant_s(1) = expected_metric_s(2, 2) * b_theta &
                    + expected_metric(2, 2) * expected_b_s(1)
                expected_covariant_s(2) = expected_metric_s(3, 3) * b_zeta &
                    + expected_metric(3, 3) * expected_b_s(2)
                expected_b = sqrt(b_theta * expected_covariant(1) &
                    + b_zeta * expected_covariant(2))
                expected_second = 0.0_dp
                expected_second(1, 1) = -two_pi**2 * radius
                expected_second(2, 2) = -two_pi**2 * surface_radius &
                    * cos(angle)
                call require(close_matrix(geometry%metric(j, k, :, :), &
                    expected_metric), "torus metric differs")
                call require(close_matrix(geometry%metric_radial(j, k, :, :), &
                    expected_metric_s), "torus radial metric differs")
                call require(close_scalar(geometry%signed_jacobian(j, k), &
                    expected_j), "torus signed Jacobian differs")
                call require(close_scalar(geometry%jacobian_s(j, k), &
                    expected_js) .and. close_scalar( &
                    geometry%jacobian_theta(j, k), expected_jt) .and. &
                    close_scalar(geometry%jacobian_zeta(j, k), 0.0_dp), &
                    "torus Jacobian derivatives differ")
                call require(close_vector(geometry%b_contravariant(j, k, :), &
                    expected_contravariant), &
                    "torus contravariant field differs")
                call require(close_vector(geometry%b_covariant(j, k, :), &
                    expected_covariant), "torus covariant field differs")
                call require(close_vector( &
                    geometry%b_contravariant_radial(j, k, :), &
                    expected_b_s), "torus radial contravariant field differs")
                call require(close_vector( &
                    geometry%b_covariant_radial(j, k, :), &
                    expected_covariant_s), &
                    "torus radial covariant field differs")
                call require(close_scalar(geometry%mod_b(j, k), expected_b), &
                    "torus magnetic magnitude differs")
                call require(close_matrix(geometry%second_form(j, k, :, :), &
                    expected_second), "torus second form differs")
            end do
        end do
    end subroutine check_geometry_grid

    subroutine check_geometry_failures(valid_jet)
        type(cartesian_jet_grid_t), intent(in) :: valid_jet
        type(cartesian_jet_grid_t) :: invalid_jet
        type(primitive_geometry_grid_t) :: geometry
        integer :: status

        call build_primitive_geometry_grid(valid_jet, 0, 7.0_dp, -3.0_dp, &
            geometry, status)
        call require(status == primitive_geometry_grid_invalid, &
            "zero field periods were accepted by geometry grid")
        call build_primitive_geometry_grid(valid_jet, 5, 0.0_dp, 0.0_dp, &
            geometry, status)
        call require(status == primitive_geometry_grid_invalid, &
            "zero magnetic field was accepted by geometry grid")
        invalid_jet = valid_jet
        deallocate (invalid_jet%radial_toroidal)
        call build_primitive_geometry_grid(invalid_jet, 5, 7.0_dp, -3.0_dp, &
            geometry, status)
        call require(status == primitive_geometry_grid_invalid, &
            "incomplete Cartesian jet was accepted")
        call require(.not. allocated(geometry%metric), &
            "failed geometry grid retained output")
    end subroutine check_geometry_failures

    subroutine check_equilibrium_spline(reference_jet, x_pair, y_pair, z_pair)
        type(cartesian_jet_grid_t), intent(in) :: reference_jet
        type(harmonic_pair_t), intent(in) :: x_pair, y_pair, z_pair
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        type(primitive_equilibrium_spline_t) :: equilibrium_spline
        type(primitive_geometry_grid_t) :: actual, reference
        real(dp), allocatable :: kernel_fields(:, :, :), kernel_drive(:, :)
        real(dp) :: pressure, pressure_slope
        integer :: status

        equilibrium%radial_grid = radial_grid_half
        equilibrium%field_periods = 5
        equilibrium%s = nodes
        equilibrium%poloidal_modes = modes_m
        equilibrium%toroidal_modes = modes_n
        equilibrium%toroidal_flux = 2.0_dp + 7.0_dp * nodes
        equilibrium%poloidal_flux = -1.0_dp - 3.0_dp * nodes
        equilibrium%pressure = 1000.0_dp * (1.0_dp - nodes)**2
        equilibrium%xhat = x_pair
        equilibrium%yhat = y_pair
        equilibrium%zhat = z_pair
        call fit_primitive_equilibrium(equilibrium, equilibrium_spline, status)
        call require(status == primitive_equilibrium_ok, &
            "primitive equilibrium fit failed")
        call evaluate_primitive_equilibrium(equilibrium_spline, query_s, &
            theta, zeta, actual, pressure, pressure_slope, status)
        call require(status == primitive_equilibrium_ok, &
            "primitive equilibrium evaluation failed")
        call build_primitive_geometry_grid(reference_jet, 5, 7.0_dp, &
            -3.0_dp, reference, status, 0.0_dp, 0.0_dp)
        call require(close_rank4(actual%metric, reference%metric) &
            .and. close_rank2(actual%signed_jacobian, &
            reference%signed_jacobian) &
            .and. close_rank2(actual%jacobian_s, reference%jacobian_s) &
            .and. close_rank2(actual%jacobian_theta, &
            reference%jacobian_theta) &
            .and. close_rank2(actual%jacobian_zeta, &
            reference%jacobian_zeta) &
            .and. close_rank4(actual%metric_radial, &
            reference%metric_radial) &
            .and. close_rank3(actual%b_contravariant, &
            reference%b_contravariant) &
            .and. close_rank3(actual%b_contravariant_radial, &
            reference%b_contravariant_radial) &
            .and. close_rank3(actual%b_covariant, reference%b_covariant) &
            .and. close_rank3(actual%b_covariant_radial, &
            reference%b_covariant_radial) &
            .and. close_rank2(actual%mod_b, reference%mod_b) &
            .and. close_rank4(actual%second_form, reference%second_form), &
            "equilibrium composition changed primitive geometry")
        call require(close_scalar(pressure, 409.6_dp) &
            .and. close_scalar(pressure_slope, -1280.0_dp), &
            "equilibrium pressure jet differs")
        call evaluate_primitive_kernel_surface(equilibrium_spline, query_s, &
            theta, zeta, kernel_fields, kernel_drive, status)
        call require(status == primitive_kernel_ok, &
            "primitive kernel composition failed")
        call require(all(ieee_is_finite(kernel_fields)) &
            .and. all(ieee_is_finite(kernel_drive)), &
            "primitive kernel composition is nonfinite")
        call require(all(close_scalar(kernel_fields(:, :, 1), -7.0_dp)) &
            .and. all(close_scalar(kernel_fields(:, :, 2), 0.6_dp)) &
            .and. all(close_scalar(kernel_fields(:, :, 3), 0.0_dp)) &
            .and. all(close_scalar(kernel_fields(:, :, 4), 0.0_dp)), &
            "primitive kernel flux profiles differ")
        call require(close_rank2(kernel_fields(:, :, 7), &
            reference%signed_jacobian) &
            .and. close_rank2(kernel_fields(:, :, 8), reference%mod_b), &
            "primitive kernel geometry differs")
        call evaluate_primitive_equilibrium(equilibrium_spline, 1.0_dp, &
            theta, zeta, actual, pressure, pressure_slope, status)
        call require(status == primitive_equilibrium_ok, &
            "half-mesh edge reconstruction failed")
        call require(close_scalar(pressure, 0.0_dp) &
            .and. close_scalar(pressure_slope, 0.0_dp), &
            "half-mesh pressure edge differs")
        call evaluate_primitive_kernel_surface(equilibrium_spline, 1.0_dp, &
            theta, zeta, kernel_fields, kernel_drive, status)
        call require(status == primitive_kernel_ok, &
            "half-mesh edge kernel evaluation failed")
        call require(all(close_scalar(kernel_fields(:, :, 1), -7.0_dp)) &
            .and. all(close_scalar(kernel_fields(:, :, 2), 0.6_dp)), &
            "half-mesh flux edge slopes differ")
        call check_compatible_problem(equilibrium)
        call check_compatible_three_component_problem(equilibrium)
        equilibrium%winding = 1
        call fit_primitive_equilibrium(equilibrium, equilibrium_spline, status)
        call require(status == primitive_equilibrium_invalid, &
            "unidentified rotating position frame was accepted")
        equilibrium%winding = 0
        equilibrium%pressure(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        call fit_primitive_equilibrium(equilibrium, equilibrium_spline, status)
        call require(status == primitive_equilibrium_invalid, &
            "nonfinite equilibrium profile was accepted")
        call evaluate_primitive_equilibrium(equilibrium_spline, query_s, &
            theta, zeta, actual, pressure, pressure_slope, status)
        call require(status == primitive_equilibrium_invalid &
            .and. .not. allocated(actual%metric) &
            .and. pressure == 0.0_dp .and. pressure_slope == 0.0_dp, &
            "invalid equilibrium spline did not fail closed")
    end subroutine check_equilibrium_spline

    subroutine check_compatible_problem(equilibrium)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(compatible_two_component_problem_t) :: problem, sparse_problem
        type(compatible_cell_trace_t), allocatable :: traces(:)
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        real(dp), allocatable :: gauss_constraint(:, :)
        real(dp), allocatable :: reference_mass(:, :), reference_stiffness(:, :)
        real(dp), allocatable :: sparse_mass(:, :), sparse_stiffness(:, :)
        real(dp) :: scale
        integer :: degree, expected_h1, expected_l2, status, term

        do degree = 1, 4
            call build_compatible_two_component_problem(equilibrium, &
                [0, 1, 2], [0, 0, 0], [0.0_dp, 0.5_dp, 1.0_dp], 1, &
                degree, 16, 16, problem, status)
            call require(status == compatible_problem_ok, &
                "compatible problem assembly failed")
            call build_compatible_two_component_problem(equilibrium, &
                [0, 1, 2], [0, 0, 0], [0.0_dp, 0.5_dp, 1.0_dp], 1, &
                degree, 16, 16, sparse_problem, status, &
                sparse_storage=.true.)
            call require(status == compatible_problem_ok, &
                "sparse compatible problem assembly failed")
            call sparse_problem_to_canonical(sparse_problem, &
                sparse_stiffness, sparse_mass, status)
            call require(status == variable_block_ok, &
                "sparse compatible problem reconstruction failed")
            if (.not. allocated(sparse_stiffness)) then
                call require(.false., &
                    "sparse compatible stiffness was not allocated")
                return
            end if
            if (.not. allocated(sparse_mass)) then
                call require(.false., &
                    "sparse compatible mass was not allocated")
                return
            end if
            expected_h1 = size(equilibrium%s) * degree - 1
            expected_l2 = size(equilibrium%s) * degree
            call require(problem%h1_dofs == expected_h1 &
                .and. problem%l2_dofs == expected_l2, &
                "compatible problem dimensions differ")
            call require(problem%normal_unknowns == 3 * expected_h1 &
                .and. problem%eta_unknowns == 2 * expected_l2, &
                "compatible component dimensions differ")
            scale = max(1.0_dp, maxval(abs(problem%stiffness)))
            call require(maxval(abs(problem%stiffness &
                - transpose(problem%stiffness))) < 3.0e-14_dp * scale, &
                "compatible stiffness is not symmetric")
            call require(all(problem%mass == transpose(problem%mass)), &
                "compatible mass is not exactly symmetric")
            call require(maxval(abs(problem%stiffness &
                - sum(problem%stiffness_terms, dim=3))) &
                < 3.0e-14_dp * scale, "compatible energy terms do not sum")
            call require(maxval(abs(sparse_stiffness - problem%stiffness)) &
                < 5.0e-14_dp * scale, &
                "sparse and dense compatible stiffness differ")
            call require(maxval(abs(sparse_mass - problem%mass)) &
                < 5.0e-14_dp * max(1.0_dp, maxval(abs(problem%mass))), &
                "sparse and dense compatible mass differ")
            do term = 1, size(problem%stiffness_terms, 3)
                call require(all(problem%stiffness_terms(:, :, term) &
                    == transpose(problem%stiffness_terms(:, :, term))), &
                    "compatible energy term is not exactly symmetric")
            end do
            call solve_symmetric_generalized(problem%stiffness, problem%mass, &
                eigenvalues, eigenvectors, status)
            call require(status == symmetric_eigensolver_ok, &
                "compatible mass is not positive definite")
        end do
        reference_stiffness = problem%stiffness
        reference_mass = problem%mass
        call build_compatible_two_component_problem(equilibrium, &
            [0, 1, 2], [0, 0, 0], [0.0_dp, 0.5_dp, 1.0_dp], 1, 4, 16, 16, &
            problem, status, [1, 3, 6], traces)
        call require(status == compatible_problem_ok, &
            "traced compatible problem assembly failed")
        call require(all(problem%stiffness == reference_stiffness) &
            .and. all(problem%mass == reference_mass), &
            "operator tracing changed the assembled problem")
        call check_compatible_trace(traces)
        call build_compatible_two_component_problem(equilibrium, &
            [0, 1, 2], [0, 0, 0], [0.0_dp, 0.0_dp, 0.0_dp], 1, 1, &
            16, 16, problem, status)
        call require(status == compatible_problem_ok, &
            "degree-one Gauss reference problem failed")
        gauss_constraint = problem%stiffness_terms(:, :, 3)
        call build_compatible_two_component_problem(equilibrium, &
            [0, 1, 2], [0, 0, 0], [0.0_dp, 0.0_dp, 0.0_dp], 1, 1, &
            16, 16, problem, status, [3], traces, &
            radial_quadrature_policy=compatible_quadrature_cas3d_midpoint)
        call require(status == compatible_problem_ok, &
            "CAS3D midpoint compatible problem failed")
        call require(problem%quadrature_points == 1 &
            .and. problem%radial_quadrature_policy &
            == compatible_quadrature_cas3d_midpoint, &
            "CAS3D midpoint policy metadata differs")
        call require(size(traces) == 1, &
            "CAS3D midpoint trace count differs")
        call require(size(traces(1)%points) == 1, &
            "CAS3D midpoint trace did not contain one point")
        call require(abs(traces(1)%points(1)%coordinate &
            - 2.5_dp / real(size(equilibrium%s), dp)) &
            < 4.0_dp * epsilon(1.0_dp), &
            "CAS3D midpoint coordinate is off by one cell")
        call require(abs(traces(1)%points(1)%weight &
            - 1.0_dp / real(size(equilibrium%s), dp)) &
            < 4.0_dp * epsilon(1.0_dp), &
            "CAS3D midpoint radial weight differs")
        call require(all(traces(1)%points(1)%term_mask) &
            .and. traces(1)%points(1)%assembles_mass, &
            "CAS3D midpoint did not assemble every term and mass")
        scale = max(1.0_dp, maxval(abs(problem%stiffness)))
        call require(maxval(abs(problem%stiffness &
            - sum(problem%stiffness_terms, dim=3))) &
            < 3.0e-14_dp * scale, &
            "CAS3D midpoint energy terms do not close")
        call require(all(problem%stiffness_terms(:, :, 3) &
            == gauss_constraint), &
            "degree-one Gauss and midpoint K3 constraint terms differ")
        call build_compatible_two_component_problem(equilibrium, [1], [0], &
            [0.0_dp], 1, 2, 16, 16, problem, status, &
            radial_quadrature_policy=compatible_quadrature_cas3d_midpoint)
        call require(status == compatible_problem_invalid, &
            "degree-two CAS3D midpoint problem was accepted")
        call build_compatible_two_component_problem(equilibrium, [1], [0], &
            [0.0_dp], 1, 0, 16, 16, problem, status)
        call require(status == compatible_problem_invalid, &
            "degree-zero compatible problem was accepted")
    end subroutine check_compatible_problem

    subroutine sparse_problem_to_canonical(problem, stiffness, mass, info)
        type(compatible_two_component_problem_t), intent(in) :: problem
        real(dp), allocatable, intent(out) :: stiffness(:, :), mass(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: block_mass(:, :), block_stiffness(:, :)
        integer, allocatable :: offset(:), position(:)
        integer :: block, column, row, unknowns

        call variable_block_to_dense(problem%sparse_stiffness, &
            block_stiffness, info)
        if (info /= variable_block_ok) return
        call variable_block_to_dense(problem%sparse_mass, block_mass, info)
        if (info /= variable_block_ok) return
        unknowns = size(problem%sparse_block_index)
        allocate (offset(size(problem%sparse_stiffness%widths)), source=1)
        do block = 2, size(offset)
            offset(block) = offset(block - 1) &
                + problem%sparse_stiffness%widths(block - 1)
        end do
        allocate (position(unknowns))
        do row = 1, unknowns
            position(row) = offset(problem%sparse_block_index(row)) &
                + problem%sparse_local_index(row) - 1
        end do
        allocate (stiffness(unknowns, unknowns), mass(unknowns, unknowns))
        do column = 1, unknowns
            do row = 1, unknowns
                stiffness(row, column) = &
                    block_stiffness(position(row), position(column))
                mass(row, column) = &
                    block_mass(position(row), position(column))
            end do
        end do
        info = variable_block_ok
    end subroutine sparse_problem_to_canonical

    subroutine check_compatible_trace(traces)
        type(compatible_cell_trace_t), intent(in) :: traces(:)
        integer, parameter :: expected_cells(3) = [1, 3, 6]
        integer, parameter :: expected_map_size(3) = [24, 27, 24]
        integer :: cell, point

        call require(size(traces) == 3, "compatible trace count differs")
        do cell = 1, size(traces)
            call require(traces(cell)%cell == expected_cells(cell), &
                "compatible trace cell selection differs")
            call require(size(traces(cell)%points) == 9, &
                "compatible trace quadrature count differs")
            do point = 1, size(traces(cell)%points)
                call require(size(traces(cell)%points(point)%map) &
                    == expected_map_size(cell), &
                    "compatible trace map size differs")
                call require(all(shape(traces(cell)%points(point)%fields) &
                    == [16, 16, 13]), "compatible trace field shape differs")
                call require(all(ieee_is_finite( &
                    traces(cell)%points(point)%stiffness_terms)) &
                    .and. all(ieee_is_finite( &
                    traces(cell)%points(point)%mass)), &
                    "compatible trace matrix is nonfinite")
                call require(all(traces(cell)%points(point)%mass &
                    == transpose(traces(cell)%points(point)%mass)), &
                    "compatible trace mass is not symmetric")
            end do
        end do
    end subroutine check_compatible_trace

    subroutine check_compatible_three_component_problem(equilibrium)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(compatible_three_component_problem_t) :: problem
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        real(dp) :: scale
        integer :: degree, expected_h1, expected_l2, status, term

        do degree = 1, 4
            call build_compatible_three_component_problem(equilibrium, &
                5.0_dp / 3.0_dp, 2.0_dp, [0, 1, 2], [0, 0, 0], &
                [0.0_dp, 0.5_dp, 1.0_dp], 1, degree, 16, 16, problem, status)
            call require(status == compatible_three_component_ok, &
                "compatible three-component problem assembly failed")
            expected_h1 = size(equilibrium%s) * degree - 1
            expected_l2 = size(equilibrium%s) * degree
            call require(problem%h1_dofs == expected_h1 &
                .and. problem%l2_dofs == expected_l2, &
                "compatible three-component dimensions differ")
            call require(problem%normal_unknowns == 3 * expected_h1 &
                .and. problem%eta_unknowns == 2 * expected_l2 &
                .and. problem%mu_unknowns == 2 * expected_l2, &
                "compatible three-component activity differs")
            scale = max(1.0_dp, maxval(abs(problem%stiffness)))
            call require(maxval(abs(problem%stiffness &
                - transpose(problem%stiffness))) < 3.0e-14_dp * scale, &
                "compatible three-component stiffness is not symmetric")
            call require(all(problem%mass == transpose(problem%mass)), &
                "compatible three-component mass is not exactly symmetric")
            call require(maxval(abs(problem%stiffness &
                - sum(problem%stiffness_terms, dim=3))) &
                < 3.0e-14_dp * scale, "compatible energy terms do not sum")
            do term = 1, size(problem%stiffness_terms, 3)
                call require(all(problem%stiffness_terms(:, :, term) &
                    == transpose(problem%stiffness_terms(:, :, term))), &
                    "compatible energy term is not exactly symmetric")
            end do
            call solve_symmetric_generalized(problem%stiffness, problem%mass, &
                eigenvalues, eigenvectors, status)
            call require(status == symmetric_eigensolver_ok, &
                "compatible physical mass is not positive definite")
        end do
        call build_compatible_three_component_problem(equilibrium, 0.0_dp, &
            2.0_dp, [1], [0], [0.0_dp], 1, 2, 16, 16, problem, status)
        call require(status == compatible_three_component_invalid, &
            "zero adiabatic index was accepted")
    end subroutine check_compatible_three_component_problem

    subroutine torus_jet(s, theta_value, zeta_value, expected)
        real(dp), intent(in) :: s, theta_value, zeta_value
        real(dp), intent(out) :: expected(3, 10)
        real(dp) :: angle_t, angle_z, radius, radius_s, radius_ss
        real(dp) :: surface_radius, surface_theta
        real(dp) :: ct, st, cz, sz, two_pi

        two_pi = 2.0_dp * acos(-1.0_dp)
        angle_t = two_pi * theta_value
        angle_z = two_pi * zeta_value
        ct = cos(angle_t)
        st = sin(angle_t)
        cz = cos(angle_z)
        sz = sin(angle_z)
        radius = minor_radius * sqrt(s)
        radius_s = minor_radius / (2.0_dp * sqrt(s))
        radius_ss = -minor_radius / (4.0_dp * s**1.5_dp)
        surface_radius = major_radius + radius * ct
        surface_theta = -two_pi * radius * st
        expected(1, 1) = surface_radius * cz
        expected(2, 1) = surface_radius * sz
        expected(3, 1) = radius * st
        expected(1, 2) = radius_s * ct * cz
        expected(2, 2) = radius_s * ct * sz
        expected(3, 2) = radius_s * st
        expected(1, 3) = -two_pi * radius * st * cz
        expected(2, 3) = -two_pi * radius * st * sz
        expected(3, 3) = two_pi * radius * ct
        expected(1, 4) = -two_pi * surface_radius * sz
        expected(2, 4) = two_pi * surface_radius * cz
        expected(3, 4) = 0.0_dp
        expected(1, 5) = radius_ss * ct * cz
        expected(2, 5) = radius_ss * ct * sz
        expected(3, 5) = radius_ss * st
        expected(1, 6) = -two_pi * radius_s * st * cz
        expected(2, 6) = -two_pi * radius_s * st * sz
        expected(3, 6) = two_pi * radius_s * ct
        expected(1, 7) = -two_pi * radius_s * ct * sz
        expected(2, 7) = two_pi * radius_s * ct * cz
        expected(3, 7) = 0.0_dp
        expected(1, 8) = -two_pi**2 * radius * ct * cz
        expected(2, 8) = -two_pi**2 * radius * ct * sz
        expected(3, 8) = -two_pi**2 * radius * st
        expected(1, 9) = -two_pi * surface_theta * sz
        expected(2, 9) = two_pi * surface_theta * cz
        expected(3, 9) = 0.0_dp
        expected(1, 10) = -two_pi**2 * surface_radius * cz
        expected(2, 10) = -two_pi**2 * surface_radius * sz
        expected(3, 10) = 0.0_dp
    end subroutine torus_jet

    subroutine check_invalid_inputs(valid_grid, valid_x, valid_y, valid_z, &
            valid_spline)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        type(harmonic_pair_t), intent(in) :: valid_x, valid_y, valid_z
        type(cartesian_harmonic_spline_t), intent(in) :: valid_spline
        type(harmonic_pair_t) :: invalid_pair
        type(cartesian_harmonic_spline_t) :: invalid_spline
        type(cartesian_jet_grid_t) :: invalid_jet
        real(dp) :: invalid_theta(1)
        integer :: status

        allocate (invalid_pair%cosine(5, 2, 3), &
            invalid_pair%sine(5, 2, 3))
        invalid_pair%cosine = 0.0_dp
        invalid_pair%sine = 0.0_dp
        call fit_cartesian_harmonic_spline(valid_grid, modes_m, modes_n, &
            valid_x, invalid_pair, valid_z, invalid_spline, status)
        call require(status == cartesian_harmonic_invalid, &
            "mismatched position shape was accepted")
        invalid_pair = valid_y
        invalid_pair%cosine(2, 1, 1) = &
            ieee_value(0.0_dp, ieee_quiet_nan)
        call fit_cartesian_harmonic_spline(valid_grid, modes_m, modes_n, &
            valid_x, invalid_pair, valid_z, invalid_spline, status)
        call require(status == cartesian_harmonic_invalid, &
            "nonfinite position coefficient was accepted")
        call evaluate_cartesian_harmonic_spline(valid_grid, valid_spline, &
            0.0_dp, theta, zeta, invalid_jet, status)
        call require(status == cartesian_harmonic_invalid, &
            "singular axis jet was accepted")
        call require(.not. allocated(invalid_jet%value), &
            "failed axis evaluation retained output")
        invalid_theta = ieee_value(0.0_dp, ieee_quiet_nan)
        call evaluate_cartesian_harmonic_spline(valid_grid, valid_spline, &
            query_s, invalid_theta, zeta, invalid_jet, status)
        call require(status == cartesian_harmonic_invalid, &
            "nonfinite angle was accepted")
        call evaluate_cartesian_harmonic_spline(valid_grid, invalid_spline, &
            query_s, theta, zeta, invalid_jet, status)
        call require(status == cartesian_harmonic_invalid, &
            "invalid spline was accepted")
    end subroutine check_invalid_inputs

    function close_vector(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:), expected(:)
        logical :: matches

        matches = size(actual) == size(expected) &
            .and. all(abs(actual - expected) <= 3.0e-11_dp &
            * max(1.0_dp, abs(expected)))
    end function close_vector

    elemental function close_scalar(actual, expected) result(matches)
        real(dp), intent(in) :: actual, expected
        logical :: matches

        matches = abs(actual - expected) <= 3.0e-11_dp &
            * max(1.0_dp, abs(expected))
    end function close_scalar

    function close_matrix(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:, :), expected(:, :)
        logical :: matches

        matches = all(shape(actual) == shape(expected)) &
            .and. all(close_scalar(actual, expected))
    end function close_matrix

    function close_rank2(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:, :), expected(:, :)
        logical :: matches

        matches = all(shape(actual) == shape(expected)) &
            .and. all(close_scalar(actual, expected))
    end function close_rank2

    function close_rank3(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:, :, :), expected(:, :, :)
        logical :: matches

        matches = all(shape(actual) == shape(expected)) &
            .and. all(close_scalar(actual, expected))
    end function close_rank3

    function close_rank4(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:, :, :, :), expected(:, :, :, :)
        logical :: matches

        matches = all(shape(actual) == shape(expected)) &
            .and. all(close_scalar(actual, expected))
    end function close_rank4

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_cartesian_harmonic_spline
