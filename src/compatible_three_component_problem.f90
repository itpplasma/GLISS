module compatible_three_component_problem
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface, &
        compatible_stiffness_term_count
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use compatible_problem_assembly_support, only: apply_stored_power, &
        build_active_indices, build_uniform_breaks, &
        compatible_support_allocation, compatible_support_ok, &
        mode_table_is_unique, replicate_indexed_values, scatter_matrix, &
        sum_tensor, symmetrize_matrix, symmetrize_tensor
    use compatible_radial_quadrature, only: accurate_nodes, &
        accurate_weights, build_constraint_quadrature, &
        compatible_quadrature_ok
    use export_surface_geometry, only: build_angular_grids
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use phase_assembly_policy, only: phase_assembly_transformed
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
    use radial_feec_complex, only: build_radial_feec_complex, &
        evaluate_radial_feec_complex, radial_feec_complex_t, radial_feec_ok
    use trial_space_topology, only: build_trial_space_topology, &
        trial_component_eta, trial_component_mu, trial_component_normal, &
        trial_space_topology_t, trial_topology_ok
    implicit none
    private

    integer, parameter, public :: compatible_three_component_ok = 0
    integer, parameter, public :: compatible_three_component_invalid = -1
    integer, parameter, public :: compatible_three_component_assembly_error = -2
    integer, parameter, public :: compatible_three_component_allocation_error = -3

    type, public :: compatible_three_component_problem_t
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        real(dp), allocatable :: stiffness_terms(:, :, :)
        integer :: degree = 0
        integer :: quadrature_points = 0
        integer :: h1_dofs = 0
        integer :: l2_dofs = 0
        integer :: normal_unknowns = 0
        integer :: eta_unknowns = 0
        integer :: mu_unknowns = 0
    end type compatible_three_component_problem_t

    public :: build_compatible_three_component_problem

    logical, parameter :: accurate_term(5) = &
        [.true., .true., .false., .true., .false.]
    logical, parameter :: constraint_term(5) = &
        [.false., .false., .true., .false., .true.]

contains

    subroutine build_compatible_three_component_problem(equilibrium, &
            adiabatic_index, density_kg_m3, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, n_zeta, problem, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        type(compatible_three_component_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        type(primitive_equilibrium_spline_t) :: spline
        type(radial_feec_complex_t) :: complex
        type(trial_space_topology_t) :: topology
        real(dp), allocatable :: breaks(:), theta(:), zeta(:)
        integer, allocatable :: parity(:), ranks(:, :)
        integer :: allocation_status, intervals, local_info, unknowns

        problem = compatible_three_component_problem_t()
        info = compatible_three_component_invalid
        if (.not. inputs_are_valid(equilibrium, adiabatic_index, &
            density_kg_m3, mode_m, mode_n, stored_power, parity_class, &
            degree, n_theta, n_zeta)) return
        intervals = size(equilibrium%s)
        allocate (breaks(intervals + 1), parity(size(mode_m)), &
            ranks(3, size(mode_m)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
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
        call build_component_ranks(topology, ranks)
        problem%normal_unknowns = complex%h1_dofs &
            * count(topology%active(trial_component_normal, :))
        problem%eta_unknowns = complex%l2_dofs &
            * count(topology%active(trial_component_eta, :))
        problem%mu_unknowns = complex%l2_dofs &
            * count(topology%active(trial_component_mu, :))
        unknowns = problem%normal_unknowns + problem%eta_unknowns &
            + problem%mu_unknowns
        if (unknowns < 1) return
        allocate (problem%stiffness(unknowns, unknowns), &
            problem%mass(unknowns, unknowns), &
            problem%stiffness_terms(unknowns, unknowns, &
            compatible_stiffness_term_count), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        problem%stiffness = 0.0_dp
        problem%mass = 0.0_dp
        problem%stiffness_terms = 0.0_dp
        call fit_primitive_equilibrium(equilibrium, spline, local_info)
        if (local_info /= primitive_equilibrium_ok) then
            info = compatible_three_component_assembly_error
            return
        end if
        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        call assemble_problem(spline, complex, breaks, theta, zeta, &
            adiabatic_index, density_kg_m3, mode_m, mode_n, parity, &
            stored_power, topology, ranks, problem, info)
        if (info /= compatible_three_component_ok) return
        call symmetrize_matrix(problem%stiffness)
        call symmetrize_matrix(problem%mass)
        call symmetrize_tensor(problem%stiffness_terms)
        problem%degree = degree
        problem%quadrature_points = size(accurate_nodes)
        problem%h1_dofs = complex%h1_dofs
        problem%l2_dofs = complex%l2_dofs
    end subroutine build_compatible_three_component_problem

    subroutine assemble_problem(spline, complex, breaks, theta, zeta, &
            adiabatic_index, density, mode_m, mode_n, parity, stored_power, &
            topology, ranks, problem, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:), theta(:), zeta(:)
        real(dp), intent(in) :: adiabatic_index, density
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: ranks(:, :)
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        real(dp), allocatable :: constraint_nodes(:), constraint_weights(:)
        real(dp) :: coordinate, half_width, midpoint, radial_weight
        integer :: cell, point

        info = compatible_three_component_assembly_error
        call build_constraint_quadrature(complex%h1_degree, &
            constraint_nodes, constraint_weights, info)
        if (info /= compatible_quadrature_ok) return
        do cell = 1, size(breaks) - 1
            midpoint = 0.5_dp * (breaks(cell) + breaks(cell + 1))
            half_width = 0.5_dp * (breaks(cell + 1) - breaks(cell))
            do point = 1, size(accurate_nodes)
                coordinate = midpoint + half_width * accurate_nodes(point)
                radial_weight = half_width * accurate_weights(point)
                call assemble_radial_point(spline, complex, coordinate, &
                    radial_weight, theta, zeta, adiabatic_index, density, &
                    mode_m, mode_n, parity, stored_power, topology, ranks, &
                    problem, accurate_term, .true., info)
                if (info /= compatible_three_component_ok) return
            end do
            do point = 1, size(constraint_nodes)
                coordinate = midpoint + half_width * constraint_nodes(point)
                radial_weight = half_width * constraint_weights(point)
                call assemble_radial_point(spline, complex, coordinate, &
                    radial_weight, theta, zeta, adiabatic_index, density, &
                    mode_m, mode_n, parity, stored_power, topology, ranks, &
                    problem, constraint_term, .false., info)
                if (info /= compatible_three_component_ok) return
            end do
        end do
        call sum_tensor(problem%stiffness_terms, problem%stiffness)
        info = compatible_three_component_ok
    end subroutine assemble_problem

    subroutine assemble_radial_point(spline, complex, coordinate, weight, &
            theta, zeta, adiabatic_index, density, mode_m, mode_n, parity, &
            stored_power, topology, ranks, problem, term_mask, assemble_mass, &
            info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: coordinate, weight, theta(:), zeta(:)
        real(dp), intent(in) :: adiabatic_index, density
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: ranks(:, :)
        type(compatible_three_component_problem_t), intent(inout) :: problem
        logical, intent(in) :: term_mask(:), assemble_mass
        integer, intent(out) :: info
        real(dp), allocatable :: fields(:, :, :), drive(:, :), jacobian_s(:, :)
        real(dp), allocatable :: jacobian_t(:, :), jacobian_z(:, :), gamma_p(:, :)
        real(dp), allocatable :: h1(:), dh1(:), l2(:), local_h1(:, :)
        real(dp), allocatable :: local_dh1(:, :), local_l2(:, :)
        real(dp), allocatable :: local_k(:, :), local_m(:, :), local_terms(:, :, :)
        integer, allocatable :: h1_index(:), l2_index(:), map(:)
        real(dp) :: pressure
        integer :: allocation_status, local_info, trials

        info = compatible_three_component_assembly_error
        call evaluate_radial_feec_complex(complex, coordinate, h1, dh1, l2, &
            local_info)
        if (local_info /= radial_feec_ok) return
        call build_active_indices(h1, h1_index, local_info, dh1)
        if (local_info == compatible_support_allocation) then
            info = compatible_three_component_allocation_error
            return
        else if (local_info /= compatible_support_ok) then
            return
        end if
        call build_active_indices(l2, l2_index, local_info)
        if (local_info == compatible_support_allocation) then
            info = compatible_three_component_allocation_error
            return
        else if (local_info /= compatible_support_ok) then
            return
        end if
        if (size(h1_index) < 1 .or. size(l2_index) < 1) return
        trials = size(mode_m)
        allocate (local_h1(size(h1_index), trials), &
            local_dh1(size(h1_index), trials), &
            local_l2(size(l2_index), trials), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        call apply_stored_power(coordinate, stored_power, h1, dh1, h1_index, &
            local_h1, local_dh1, local_info)
        if (local_info /= compatible_support_ok) return
        call replicate_indexed_values(l2, l2_index, local_l2, local_info)
        if (local_info /= compatible_support_ok) return
        call evaluate_primitive_kernel_surface(spline, coordinate, theta, &
            zeta, fields, drive, local_info, jacobian_s, jacobian_t, &
            jacobian_z, pressure)
        if (local_info /= primitive_kernel_ok) return
        allocate (gamma_p(size(theta), size(zeta)), &
            source=adiabatic_index * pressure, stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        call allocate_local_matrices(trials, size(h1_index), size(l2_index), &
            local_k, local_m, local_terms, allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        call assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, &
            mode_n, parity, spline%field_periods, local_h1, local_dh1, &
            local_l2, weight, phase_assembly_transformed, local_k, local_info, &
            local_terms)
        if (local_info /= 0) return
        if (assemble_mass) then
            call assemble_compatible_physical_mass_surface(fields, density, &
                mode_m, mode_n, parity, spline%field_periods, local_h1, &
                local_l2, weight, phase_assembly_transformed, local_m, &
                local_info)
            if (local_info /= 0) return
        end if
        call build_local_map(complex, topology, ranks, h1_index, l2_index, map)
        if (assemble_mass) call scatter_matrix(map, local_m, 1.0_dp, &
            problem%mass)
        call scatter_terms(map, local_terms, term_mask, &
            problem%stiffness_terms)
        info = compatible_three_component_ok
    end subroutine assemble_radial_point

    subroutine allocate_local_matrices(trials, h1_count, l2_count, &
            stiffness, mass, terms, status)
        integer, intent(in) :: trials, h1_count, l2_count
        real(dp), allocatable, intent(out) :: stiffness(:, :), mass(:, :)
        real(dp), allocatable, intent(out) :: terms(:, :, :)
        integer, intent(out) :: status
        integer :: dimension

        dimension = trials * (h1_count + 2 * l2_count)
        allocate (stiffness(dimension, dimension), &
            mass(dimension, dimension), &
            terms(dimension, dimension, compatible_stiffness_term_count), &
            stat=status)
        if (status == 0) then
            stiffness = 0.0_dp
            mass = 0.0_dp
            terms = 0.0_dp
        end if
    end subroutine allocate_local_matrices

    subroutine build_local_map(complex, topology, ranks, h1_index, l2_index, map)
        type(radial_feec_complex_t), intent(in) :: complex
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: ranks(:, :), h1_index(:), l2_index(:)
        integer, allocatable, intent(out) :: map(:)
        integer :: basis, column, component, offset, trial, trials

        trials = size(ranks, 2)
        allocate (map(trials * (size(h1_index) + 2 * size(l2_index))), source=0)
        do basis = 1, size(h1_index)
            do trial = 1, trials
                column = (basis - 1) * trials + trial
                if (topology%active(trial_component_normal, trial)) &
                    map(column) = (h1_index(basis) - 1) &
                    * count(ranks(trial_component_normal, :) > 0) &
                    + ranks(trial_component_normal, trial)
            end do
        end do
        offset = complex%h1_dofs * count(ranks(trial_component_normal, :) > 0)
        do component = trial_component_eta, trial_component_mu
            do basis = 1, size(l2_index)
                do trial = 1, trials
                    column = size(h1_index) * trials &
                        + (component - trial_component_eta) &
                        * size(l2_index) * trials &
                        + (basis - 1) * trials + trial
                    if (topology%active(component, trial)) &
                        map(column) = offset + (l2_index(basis) - 1) &
                        * count(ranks(component, :) > 0) + ranks(component, trial)
                end do
            end do
            offset = offset + complex%l2_dofs &
                * count(ranks(component, :) > 0)
        end do
    end subroutine build_local_map

    subroutine scatter_terms(map, local, term_mask, global)
        integer, intent(in) :: map(:)
        real(dp), intent(in) :: local(:, :, :)
        logical, intent(in) :: term_mask(:)
        real(dp), intent(inout) :: global(:, :, :)
        integer :: term

        do term = 1, size(global, 3)
            if (.not. term_mask(term)) cycle
            call scatter_matrix(map, local(:, :, term), 1.0_dp, &
                global(:, :, term))
        end do
    end subroutine scatter_terms

    pure subroutine build_component_ranks(topology, ranks)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(out) :: ranks(:, :)
        integer :: component, trial

        ranks = 0
        do component = 1, 3
            do trial = 1, size(ranks, 2)
                if (topology%active(component, trial)) &
                    ranks(component, trial) = count(ranks(component, :) > 0) + 1
            end do
        end do
    end subroutine build_component_ranks

    function inputs_are_valid(equilibrium, adiabatic_index, density, mode_m, &
            mode_n, stored_power, parity_class, degree, n_theta, n_zeta) &
            result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density, stored_power(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        logical :: valid

        valid = size(equilibrium%s) >= 4 .and. equilibrium%field_periods >= 1
        valid = valid .and. size(mode_m) >= 1 .and. size(mode_n) == size(mode_m)
        valid = valid .and. size(stored_power) == size(mode_m)
        valid = valid .and. all(mode_m >= 0) &
            .and. all(ieee_is_finite(stored_power))
        valid = valid .and. ieee_is_finite(adiabatic_index) &
            .and. adiabatic_index > 0.0_dp
        valid = valid .and. ieee_is_finite(density) .and. density > 0.0_dp
        valid = valid .and. parity_class >= 1 .and. parity_class <= 2
        valid = valid .and. degree >= 1 .and. degree <= 4
        valid = valid .and. n_theta >= 8 .and. n_zeta >= 8
        if (.not. valid) return
        if (any(mode_m == 0 .and. mode_n < 0)) then
            valid = .false.
            return
        end if
        valid = mode_table_is_unique(mode_m, mode_n)
    end function inputs_are_valid

end module compatible_three_component_problem
