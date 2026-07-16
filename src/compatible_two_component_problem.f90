module compatible_two_component_problem
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_family_point_assembly, only: &
        assemble_compatible_transformed_surface, &
        compatible_two_component_term_count
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_perpendicular_mass_surface
    use compatible_operator_trace_types, only: build_trace_radial_mass, &
        compatible_cell_trace_t, compatible_radial_point_trace_t
    use compatible_problem_assembly_support, only: apply_stored_power, &
        build_active_indices, build_uniform_breaks, &
        compatible_support_allocation, compatible_support_ok, &
        mode_table_is_unique, replicate_indexed_values, scale_matrix, &
        scale_tensor, scatter_matrix, sum_tensor, symmetrize_matrix, &
        symmetrize_tensor
    use compatible_radial_quadrature, only: accurate_nodes, &
        accurate_weights, build_constraint_quadrature, &
        compatible_quadrature_ok
    use export_surface_geometry, only: build_angular_grids
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    use radial_feec_complex, only: build_radial_feec_complex, &
        evaluate_radial_feec_complex, radial_feec_complex_t, radial_feec_ok
    use trial_space_topology, only: build_trial_space_topology, &
        trial_component_eta, trial_component_normal, trial_space_topology_t, &
        trial_topology_ok
    implicit none
    private

    integer, parameter, public :: compatible_problem_ok = 0
    integer, parameter, public :: compatible_problem_invalid = -1
    integer, parameter, public :: compatible_problem_assembly_error = -2
    integer, parameter, public :: compatible_problem_allocation_error = -3
    integer, parameter, public :: compatible_quadrature_gauss = 1
    integer, parameter, public :: compatible_quadrature_cas3d_midpoint = 2

    type, public :: compatible_two_component_problem_t
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        real(dp), allocatable :: stiffness_terms(:, :, :)
        integer :: degree = 0
        integer :: quadrature_points = 0
        integer :: h1_dofs = 0
        integer :: l2_dofs = 0
        integer :: normal_unknowns = 0
        integer :: eta_unknowns = 0
        logical :: has_physical_mass = .false.
        integer :: radial_quadrature_policy = 0
    end type compatible_two_component_problem_t

    public :: build_compatible_two_component_problem
    public :: compatible_cell_trace_t

    logical, parameter :: accurate_term(4) = [.true., .true., .false., .true.]
    logical, parameter :: constraint_term(4) = &
        [.false., .false., .true., .false.]
    logical, parameter :: all_terms(4) = [.true., .true., .true., .true.]

contains

    subroutine build_compatible_two_component_problem(equilibrium, mode_m, &
            mode_n, stored_power, parity_class, degree, n_theta, n_zeta, &
            problem, info, trace_cells, trace, density_kg_m3, &
            radial_quadrature_policy)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        type(compatible_two_component_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        integer, optional, intent(in) :: trace_cells(:)
        type(compatible_cell_trace_t), allocatable, optional, intent(out) :: trace(:)
        real(dp), optional, intent(in) :: density_kg_m3
        integer, optional, intent(in) :: radial_quadrature_policy
        type(primitive_equilibrium_spline_t) :: spline
        type(radial_feec_complex_t) :: complex
        type(trial_space_topology_t) :: topology
        real(dp), allocatable :: breaks(:), theta(:), zeta(:)
        integer, allocatable :: eta_rank(:), normal_rank(:), parity(:)
        integer :: allocation_status, intervals, local_info, quadrature_policy
        integer :: unknowns

        problem = compatible_two_component_problem_t()
        info = compatible_problem_invalid
        quadrature_policy = compatible_quadrature_gauss
        if (present(radial_quadrature_policy)) &
            quadrature_policy = radial_quadrature_policy
        if (.not. inputs_are_valid(equilibrium, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, n_zeta, quadrature_policy)) return
        if (present(density_kg_m3)) then
            if (.not. ieee_is_finite(density_kg_m3) &
                .or. density_kg_m3 <= 0.0_dp) return
        end if
        intervals = size(equilibrium%s)
        if (present(trace_cells) .neqv. present(trace)) return
        if (present(trace_cells)) then
            if (.not. trace_cells_are_valid(trace_cells, intervals)) return
            allocate (trace(size(trace_cells)), stat=allocation_status)
            if (allocation_status /= 0) then
                info = compatible_problem_allocation_error
                return
            end if
        end if
        allocate (breaks(intervals + 1), parity(size(mode_m)), &
            normal_rank(size(mode_m)), eta_rank(size(mode_m)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        call build_uniform_breaks(intervals, breaks, local_info)
        if (local_info /= compatible_support_ok) return
        parity = parity_class
        call build_radial_feec_complex(breaks, degree, .true., .true., &
            complex, local_info)
        if (local_info /= radial_feec_ok) return
        call build_trial_space_topology(mode_m, mode_n, parity, topology, &
            local_info)
        if (local_info /= trial_topology_ok) return
        call build_component_ranks(topology, normal_rank, eta_rank)
        problem%normal_unknowns = complex%h1_dofs &
            * count(topology%active(trial_component_normal, :))
        problem%eta_unknowns = complex%l2_dofs &
            * count(topology%active(trial_component_eta, :))
        unknowns = problem%normal_unknowns + problem%eta_unknowns
        if (unknowns < 1) return
        allocate (problem%stiffness(unknowns, unknowns), source=0.0_dp, &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        allocate (problem%mass(unknowns, unknowns), source=0.0_dp, &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        allocate (problem%stiffness_terms(unknowns, unknowns, &
            compatible_two_component_term_count), &
            source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        call fit_primitive_equilibrium(equilibrium, spline, local_info)
        if (local_info /= primitive_equilibrium_ok) then
            info = compatible_problem_assembly_error
            return
        end if
        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        if (present(trace)) then
            call assemble_problem(spline, complex, breaks, theta, zeta, &
                mode_m, mode_n, parity, stored_power, topology, normal_rank, &
                eta_rank, problem, info, trace_cells, trace, density_kg_m3, &
                quadrature_policy)
        else
            call assemble_problem(spline, complex, breaks, theta, zeta, &
                mode_m, mode_n, parity, stored_power, topology, normal_rank, &
                eta_rank, problem, info, quadrature_policy=quadrature_policy, &
                density_kg_m3=density_kg_m3)
        end if
        if (info /= compatible_problem_ok) return
        call symmetrize_matrix(problem%stiffness)
        call symmetrize_matrix(problem%mass)
        call symmetrize_tensor(problem%stiffness_terms)
        problem%degree = degree
        if (quadrature_policy == compatible_quadrature_gauss) then
            problem%quadrature_points = size(accurate_nodes)
        else
            problem%quadrature_points = 1
        end if
        problem%h1_dofs = complex%h1_dofs
        problem%l2_dofs = complex%l2_dofs
        problem%has_physical_mass = present(density_kg_m3)
        problem%radial_quadrature_policy = quadrature_policy
    end subroutine build_compatible_two_component_problem

    subroutine assemble_problem(spline, complex, breaks, theta, zeta, mode_m, &
            mode_n, parity, stored_power, topology, normal_rank, eta_rank, &
            problem, info, trace_cells, trace, density_kg_m3, quadrature_policy)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:), theta(:), zeta(:)
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: normal_rank(:), eta_rank(:)
        type(compatible_two_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        integer, optional, intent(in) :: trace_cells(:)
        type(compatible_cell_trace_t), optional, intent(inout) :: trace(:)
        real(dp), optional, intent(in) :: density_kg_m3
        integer, optional, intent(in) :: quadrature_policy
        real(dp), allocatable :: constraint_nodes(:), constraint_weights(:)
        real(dp) :: coordinate, half_width, midpoint, radial_weight
        integer :: cell, point, trace_index, trace_point
        integer :: policy

        info = compatible_problem_assembly_error
        policy = compatible_quadrature_gauss
        if (present(quadrature_policy)) policy = quadrature_policy
        call build_constraint_quadrature(complex%h1_degree, &
            constraint_nodes, constraint_weights, info)
        if (info /= compatible_quadrature_ok) return
        do cell = 1, size(breaks) - 1
            trace_index = 0
            if (present(trace_cells)) then
                trace_index = findloc(trace_cells, cell, dim=1)
                if (trace_index > 0) then
                    trace(trace_index)%cell = cell
                    if (policy == compatible_quadrature_gauss) then
                        allocate (trace(trace_index)%points( &
                            size(accurate_nodes) + size(constraint_nodes)))
                    else
                        allocate (trace(trace_index)%points(1))
                    end if
                end if
            end if
            midpoint = 0.5_dp * (breaks(cell) + breaks(cell + 1))
            half_width = 0.5_dp * (breaks(cell + 1) - breaks(cell))
            if (policy == compatible_quadrature_cas3d_midpoint) then
                radial_weight = 2.0_dp * half_width
                if (trace_index > 0) then
                    call assemble_radial_point(spline, complex, midpoint, &
                        radial_weight, theta, zeta, mode_m, mode_n, parity, &
                        stored_power, topology, normal_rank, eta_rank, &
                        problem, all_terms, .true., &
                        info, trace(trace_index)%points(1), density_kg_m3)
                else
                    call assemble_radial_point(spline, complex, midpoint, &
                        radial_weight, theta, zeta, mode_m, mode_n, parity, &
                        stored_power, topology, normal_rank, eta_rank, &
                        problem, all_terms, .true., &
                        info, density_kg_m3=density_kg_m3)
                end if
                if (info /= compatible_problem_ok) return
                cycle
            end if
            do point = 1, size(accurate_nodes)
                coordinate = midpoint + half_width * accurate_nodes(point)
                radial_weight = half_width * accurate_weights(point)
                if (trace_index > 0) then
                    call assemble_radial_point(spline, complex, coordinate, &
                        radial_weight, theta, zeta, mode_m, mode_n, parity, &
                        stored_power, topology, normal_rank, eta_rank, &
                        problem, accurate_term, .true., info, &
                        trace(trace_index)%points(point), density_kg_m3)
                else
                    call assemble_radial_point(spline, complex, coordinate, &
                        radial_weight, theta, zeta, mode_m, mode_n, parity, &
                        stored_power, topology, normal_rank, eta_rank, &
                        problem, accurate_term, .true., info, &
                        density_kg_m3=density_kg_m3)
                end if
                if (info /= compatible_problem_ok) return
            end do
            do point = 1, size(constraint_nodes)
                coordinate = midpoint + half_width * constraint_nodes(point)
                radial_weight = half_width * constraint_weights(point)
                trace_point = size(accurate_nodes) + point
                if (trace_index > 0) then
                    call assemble_radial_point(spline, complex, coordinate, &
                        radial_weight, theta, zeta, mode_m, mode_n, parity, &
                        stored_power, topology, normal_rank, eta_rank, &
                        problem, constraint_term, .false., info, &
                        trace(trace_index)%points(trace_point), density_kg_m3)
                else
                    call assemble_radial_point(spline, complex, coordinate, &
                        radial_weight, theta, zeta, mode_m, mode_n, parity, &
                        stored_power, topology, normal_rank, eta_rank, &
                        problem, constraint_term, .false., info, &
                        density_kg_m3=density_kg_m3)
                end if
                if (info /= compatible_problem_ok) return
            end do
        end do
        call sum_tensor(problem%stiffness_terms, problem%stiffness)
        info = compatible_problem_ok
    end subroutine assemble_problem

    subroutine assemble_radial_point(spline, complex, coordinate, weight, &
            theta, zeta, mode_m, mode_n, parity, stored_power, topology, &
            normal_rank, eta_rank, problem, term_mask, assemble_mass, info, &
            point_trace, density_kg_m3)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: coordinate, weight, theta(:), zeta(:)
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: normal_rank(:), eta_rank(:)
        type(compatible_two_component_problem_t), intent(inout) :: problem
        logical, intent(in) :: term_mask(:), assemble_mass
        integer, intent(out) :: info
        type(compatible_radial_point_trace_t), optional, intent(out) :: point_trace
        real(dp), optional, intent(in) :: density_kg_m3
        real(dp), allocatable :: fields(:, :, :), drive(:, :), h1(:), dh1(:)
        real(dp), allocatable :: l2(:), local(:, :), local_dh1(:, :)
        real(dp), allocatable :: local_mass(:, :), local_terms(:, :, :)
        real(dp), allocatable :: local_h1(:, :), local_l2(:, :)
        integer, allocatable :: h1_index(:), l2_index(:), map(:)
        real(dp) :: stiffness_scale
        integer :: local_info, term, trials

        info = compatible_problem_assembly_error
        call evaluate_radial_feec_complex(complex, coordinate, h1, dh1, l2, &
            local_info)
        if (local_info /= radial_feec_ok) return
        call build_active_indices(h1, h1_index, local_info, dh1)
        if (local_info == compatible_support_allocation) then
            info = compatible_problem_allocation_error
            return
        else if (local_info /= compatible_support_ok) then
            return
        end if
        call build_active_indices(l2, l2_index, local_info)
        if (local_info == compatible_support_allocation) then
            info = compatible_problem_allocation_error
            return
        else if (local_info /= compatible_support_ok) then
            return
        end if
        if (size(h1_index) < 1 .or. size(l2_index) < 1) return
        trials = size(mode_m)
        allocate (local_h1(size(h1_index), trials), &
            local_dh1(size(h1_index), trials), &
            local_l2(size(l2_index), trials), stat=local_info)
        if (local_info /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        call apply_stored_power(coordinate, stored_power, h1, dh1, h1_index, &
            local_h1, local_dh1, local_info)
        if (local_info /= compatible_support_ok) return
        call replicate_indexed_values(l2, l2_index, local_l2, local_info)
        if (local_info /= compatible_support_ok) return
        call evaluate_primitive_kernel_surface(spline, coordinate, theta, &
            zeta, fields, drive, local_info)
        if (local_info /= primitive_kernel_ok) return
        allocate (local(trials * (size(h1_index) + size(l2_index)), &
            trials * (size(h1_index) + size(l2_index))), source=0.0_dp, &
            stat=local_info)
        if (local_info /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        allocate (local_terms(size(local, 1), size(local, 2), &
            compatible_two_component_term_count), source=0.0_dp, &
            stat=local_info)
        if (local_info /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        call assemble_compatible_transformed_surface(fields, drive, mode_m, &
            mode_n, parity, spline%field_periods, local_h1, local_dh1, &
            local_l2, local, local_info, local_terms)
        if (local_info /= 0) return
        if (present(density_kg_m3) .and. assemble_mass) then
            allocate (local_mass(size(local, 1), size(local, 2)), &
                source=0.0_dp, stat=local_info)
            if (local_info /= 0) then
                info = compatible_problem_allocation_error
                return
            end if
            call assemble_compatible_perpendicular_mass_surface(fields, &
                density_kg_m3, mode_m, mode_n, parity, spline%field_periods, &
                local_h1, local_l2, weight, phase_assembly_transformed, &
                local_mass, local_info)
            if (local_info /= 0) return
        end if
        call build_local_map(complex, topology, normal_rank, eta_rank, &
            h1_index, l2_index, map)
        if (present(point_trace)) then
            point_trace%coordinate = coordinate
            point_trace%weight = weight
            point_trace%term_mask = term_mask
            point_trace%assembles_mass = assemble_mass
            point_trace%map = map
            point_trace%fields = fields
            point_trace%drive = drive
            point_trace%h1 = local_h1
            point_trace%dh1 = local_dh1
            point_trace%l2 = local_l2
            if (present(density_kg_m3)) then
                point_trace%stiffness_terms = local_terms
                call scale_tensor(point_trace%stiffness_terms, &
                    1.0_dp / vacuum_permeability)
                allocate (point_trace%mass(size(local, 1), size(local, 2)), &
                    source=0.0_dp)
                if (assemble_mass) then
                    point_trace%mass = local_mass
                    call scale_matrix(point_trace%mass, 1.0_dp / weight)
                end if
            else
                point_trace%stiffness_terms = local_terms
                call build_trace_radial_mass(local_h1, local_l2, &
                    point_trace%mass)
            end if
        end if
        stiffness_scale = weight
        if (present(density_kg_m3)) stiffness_scale = weight / vacuum_permeability
        do term = 1, compatible_two_component_term_count
            if (.not. term_mask(term)) cycle
            call scatter_matrix(map, local_terms(:, :, term), &
                stiffness_scale, problem%stiffness_terms(:, :, term))
        end do
        if (assemble_mass) then
            if (present(density_kg_m3)) then
                call scatter_matrix(map, local_mass, 1.0_dp, problem%mass)
            else
                call add_radial_mass(map, local_h1, local_l2, weight, &
                    problem%mass)
            end if
        end if
        info = compatible_problem_ok
    end subroutine assemble_radial_point

    subroutine build_local_map(complex, topology, normal_rank, eta_rank, &
            h1_index, l2_index, map)
        type(radial_feec_complex_t), intent(in) :: complex
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: normal_rank(:), eta_rank(:)
        integer, intent(in) :: h1_index(:), l2_index(:)
        integer, allocatable, intent(out) :: map(:)
        integer :: basis, column, eta_count, normal_count, trial, trials

        trials = size(normal_rank)
        normal_count = count(normal_rank > 0)
        eta_count = count(eta_rank > 0)
        allocate (map(trials * (size(h1_index) + size(l2_index))), source=0)
        do basis = 1, size(h1_index)
            do trial = 1, trials
                column = (basis - 1) * trials + trial
                if (topology%active(trial_component_normal, trial)) &
                    map(column) = (h1_index(basis) - 1) * normal_count &
                    + normal_rank(trial)
            end do
        end do
        do basis = 1, size(l2_index)
            do trial = 1, trials
                column = size(h1_index) * trials + (basis - 1) * trials + trial
                if (topology%active(trial_component_eta, trial)) &
                    map(column) = complex%h1_dofs * normal_count &
                    + (l2_index(basis) - 1) * eta_count + eta_rank(trial)
            end do
        end do
    end subroutine build_local_map

    subroutine add_radial_mass(map, h1, l2, weight, mass)
        integer, intent(in) :: map(:)
        real(dp), intent(in) :: h1(:, :), l2(:, :), weight
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: basis(size(map))
        integer :: a, b, basis_index, h1_columns, trial, trials

        trials = size(h1, 2)
        h1_columns = size(h1, 1) * trials
        do basis_index = 1, size(h1, 1)
            do trial = 1, trials
                basis((basis_index - 1) * trials + trial) = &
                    h1(basis_index, trial)
            end do
        end do
        do basis_index = 1, size(l2, 1)
            do trial = 1, trials
                basis(size(h1, 1) * trials + (basis_index - 1) * trials &
                    + trial) = l2(basis_index, trial)
            end do
        end do
        do b = 1, size(map)
            if (map(b) == 0) cycle
            do a = 1, size(map)
                if (map(a) == 0) cycle
                if (modulo(a - 1, trials) /= modulo(b - 1, trials)) cycle
                if ((a <= h1_columns) .neqv. (b <= h1_columns)) cycle
                mass(map(a), map(b)) = mass(map(a), map(b)) &
                    + weight * basis(a) * basis(b)
            end do
        end do
    end subroutine add_radial_mass

    pure subroutine build_component_ranks(topology, normal_rank, eta_rank)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(out) :: normal_rank(:), eta_rank(:)
        integer :: trial

        normal_rank = 0
        eta_rank = 0
        do trial = 1, size(normal_rank)
            if (topology%active(trial_component_normal, trial)) &
                normal_rank(trial) = count(normal_rank > 0) + 1
            if (topology%active(trial_component_eta, trial)) &
                eta_rank(trial) = count(eta_rank > 0) + 1
        end do
    end subroutine build_component_ranks

    function inputs_are_valid(equilibrium, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, n_zeta, quadrature_policy) &
            result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        integer, intent(in) :: quadrature_policy
        logical :: valid

        valid = size(equilibrium%s) >= 4 .and. size(mode_m) >= 1
        if (.not. valid) return
        valid = size(mode_n) == size(mode_m) &
            .and. size(stored_power) == size(mode_m)
        if (.not. valid) return
        valid = all(mode_m >= 0) .and. all(ieee_is_finite(stored_power))
        if (.not. valid) return
        if (any(mode_m == 0 .and. mode_n < 0)) then
            valid = .false.
            return
        end if
        valid = mode_table_is_unique(mode_m, mode_n)
        if (.not. valid) return
        valid = parity_class >= 1 .and. parity_class <= 2
        if (.not. valid) return
        valid = degree >= 1 .and. degree <= 4
        if (.not. valid) return
        valid = quadrature_policy == compatible_quadrature_gauss &
            .or. quadrature_policy == compatible_quadrature_cas3d_midpoint
        if (.not. valid) return
        if (quadrature_policy == compatible_quadrature_cas3d_midpoint) then
            valid = degree == 1
            if (.not. valid) return
        end if
        valid = n_theta >= 8 .and. n_zeta >= 8
    end function inputs_are_valid

    pure function trace_cells_are_valid(trace_cells, intervals) result(valid)
        integer, intent(in) :: trace_cells(:), intervals
        logical :: valid
        integer :: first, second

        valid = size(trace_cells) >= 1 .and. all(trace_cells >= 1) &
            .and. all(trace_cells <= intervals)
        if (.not. valid) return
        do first = 1, size(trace_cells)
            do second = 1, first - 1
                if (trace_cells(first) == trace_cells(second)) then
                    valid = .false.
                    return
                end if
            end do
        end do
    end function trace_cells_are_valid

end module compatible_two_component_problem
