module compatible_three_component_problem
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_axis_regularity, only: axis_regularity_ok, axis_tie_t, &
        build_trial_axis_tie, tie_local_map
    use compatible_block_storage, only: allocate_compatible_blocks, &
        build_compatible_block_indices, compatible_block_allocation, &
        compatible_block_ok, scatter_symmetric_compatible_block, &
        symmetrize_compatible_blocks
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface, &
        compatible_stiffness_term_count
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use compatible_problem_assembly_support, only: apply_stored_power, &
        surface_preserves_parity, &
        apply_tangential_axis_weight, build_active_indices, &
        build_uniform_breaks, &
        compatible_support_allocation, compatible_support_ok, &
        mode_table_is_unique, replicate_indexed_values, scatter_matrix, &
        sum_tensor, symmetrize_matrix, symmetrize_tensor
    use compatible_radial_quadrature, only: accurate_nodes, &
        accurate_weights, axis_quadrature_points, build_axis_quadrature, &
        build_constraint_quadrature, compatible_quadrature_ok
    use export_surface_geometry, only: build_angular_grids
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use phase_assembly_policy, only: phase_assembly_transformed
    use phase_factor_topology, only: phase_cosine, phase_sine
    use plasma_vacuum_boundary, only: build_vacuum_edge_block, &
        plasma_vacuum_model_t, plasma_vacuum_ok, &
        plasma_vacuum_underresolved, plasma_vacuum_wall_not_nested
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
    use radial_feec_complex, only: build_radial_feec_complex, &
        evaluate_radial_feec_complex, radial_feec_complex_t, radial_feec_ok
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    use trial_space_topology, only: build_trial_space_topology, &
        trial_component_eta, trial_component_mu, trial_component_normal, &
        trial_space_topology_t, trial_topology_ok
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    private

    integer, parameter, public :: compatible_three_component_ok = 0
    integer, parameter, public :: compatible_three_component_invalid = -1
    integer, parameter, public :: compatible_three_component_assembly_error = -2
    integer, parameter, public :: compatible_three_component_allocation_error = -3
    ! The equilibrium breaks the (theta,zeta)->(-theta,-zeta) symmetry that
    ! decouples the two parity classes; their operator does not apply.
    integer, parameter, public :: compatible_three_component_asymmetric = -4
    ! Free-boundary vacuum failures: the edge mesh cannot resolve the mode
    ! table, the wall does not enclose the plasma, or the vacuum
    ! inductance is singular.
    integer, parameter, public :: compatible_three_component_vacuum_mesh = -5
    integer, parameter, public :: compatible_three_component_wall = -6
    integer, parameter, public :: compatible_three_component_vacuum = -7

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
        integer :: axis_quadrature_points = 0
        ! Parity class 0 couples both Fourier parities of every mode, the
        ! operator of an equilibrium without stellarator symmetry; its
        ! trials are the mode table with parity 1 followed by parity 2.
        logical :: coupled = .false.
        ! A free-boundary problem retains the edge normal coefficient and
        ! adds the vacuum stiffness, also stored alone for diagnostics.
        logical :: free_boundary = .false.
        real(dp), allocatable :: vacuum(:, :)
        ! Block-tridiagonal storage by radial basis group: normal, eta and
        ! mu unknowns of each group of degree consecutive functions. The
        ! dense arrays above stay unallocated when it is used.
        logical :: has_sparse_storage = .false.
        type(variable_block_tridiagonal_t) :: sparse_stiffness, sparse_mass
        type(variable_block_tridiagonal_t) :: sparse_vacuum
        type(variable_block_tridiagonal_t) :: &
            sparse_terms(compatible_stiffness_term_count)
        integer, allocatable :: sparse_block_index(:), sparse_local_index(:)
        type(axis_tie_t) :: axis_tie
    end type compatible_three_component_problem_t

    public :: build_compatible_three_component_problem
    public :: build_compatible_three_component_classes
    public :: build_compatible_vacuum_energy

    logical, parameter :: accurate_term(5) = &
        [.true., .true., .false., .true., .false.]
    logical, parameter :: constraint_term(5) = &
        [.false., .false., .true., .false., .true.]
    logical, parameter :: all_terms(5) = .true.
    ! Kinds of radial quadrature points.
    integer, parameter :: point_axis = 1, point_accurate = 2
    integer, parameter :: point_constraint = 3

    ! Local matrices of one radial point for one problem, tied to its axis
    ! constraint, and their map onto its global unknowns.
    type :: local_matrices_t
        integer, allocatable :: map(:)
        real(dp), allocatable :: mass(:, :), terms(:, :, :)
        logical :: has_mass = .false.
    end type local_matrices_t

    ! One radial point: the cosine-parity (or single) problem and, when
    ! assembled in pairs, the sine-parity problem.
    type :: radial_contribution_t
        type(local_matrices_t) :: cosine, sine
        integer :: orientation = 0
        integer :: info = compatible_three_component_assembly_error
    end type radial_contribution_t

    ! Per-problem data of an assembly pass.
    type :: class_setup_t
        type(trial_space_topology_t) :: topology
        integer, allocatable :: ranks(:, :)
        real(dp), allocatable :: block(:, :)
    end type class_setup_t

contains

    subroutine build_compatible_three_component_problem(equilibrium, &
            adiabatic_index, density_kg_m3, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, n_zeta, problem, info, vacuum, &
            sparse_storage, vacuum_energy)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        type(compatible_three_component_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        ! Present: free boundary with this vacuum and wall model.
        type(plasma_vacuum_model_t), optional, intent(in) :: vacuum
        ! True: assemble block-tridiagonal storage instead of dense arrays.
        logical, optional, intent(in) :: sparse_storage
        ! With vacuum: the trial-by-trial vacuum energy of this problem's
        ! trials (build_compatible_vacuum_energy), instead of solving it.
        real(dp), optional, intent(in) :: vacuum_energy(:, :)
        integer, allocatable :: parity(:), trial_m(:), trial_n(:)
        real(dp), allocatable :: trial_power(:)
        integer :: allocation_status, count

        info = compatible_three_component_invalid
        if (.not. inputs_are_valid(equilibrium, adiabatic_index, &
            density_kg_m3, mode_m, mode_n, stored_power, parity_class, &
            degree, n_theta, n_zeta)) return
        count = size(mode_m)
        if (parity_class == 0) count = 2 * count
        allocate (parity(count), trial_m(count), trial_n(count), &
            trial_power(count), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        if (parity_class == 0) then
            trial_m(:size(mode_m)) = mode_m
            trial_m(size(mode_m) + 1:) = mode_m
            trial_n(:size(mode_m)) = mode_n
            trial_n(size(mode_m) + 1:) = mode_n
            trial_power(:size(mode_m)) = stored_power
            trial_power(size(mode_m) + 1:) = stored_power
            parity(:size(mode_m)) = 1
            parity(size(mode_m) + 1:) = 2
        else
            trial_m = mode_m
            trial_n = mode_n
            trial_power = stored_power
            parity = parity_class
        end if
        problem%coupled = parity_class == 0
        problem%free_boundary = present(vacuum)
        if (present(sparse_storage)) problem%has_sparse_storage = sparse_storage
        call build_trials(equilibrium, adiabatic_index, density_kg_m3, &
            trial_m, trial_n, trial_power, parity, degree, n_theta, n_zeta, &
            problem, info, vacuum, vacuum_energy)
    end subroutine build_compatible_three_component_problem

    ! Both parity classes of a stellarator-symmetric problem in one pass:
    ! the sine-parity forms come from the angular products of the
    ! cosine-parity ones. vacuum_energy is the block of
    ! build_compatible_vacuum_energy (every mode with parity 1, then 2).
    ! An equilibrium that breaks the parity symmetry returns
    ! compatible_three_component_asymmetric.
    subroutine build_compatible_three_component_classes(equilibrium, &
            adiabatic_index, density_kg_m3, mode_m, mode_n, stored_power, &
            degree, n_theta, n_zeta, cosine_problem, sine_problem, info, &
            vacuum, sparse_storage, vacuum_energy)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: degree, n_theta, n_zeta
        type(compatible_three_component_problem_t), intent(out) :: &
            cosine_problem, sine_problem
        integer, intent(out) :: info
        type(plasma_vacuum_model_t), optional, intent(in) :: vacuum
        logical, optional, intent(in) :: sparse_storage
        real(dp), optional, intent(in) :: vacuum_energy(:, :)
        integer, allocatable :: parity(:)
        integer :: modes

        info = compatible_three_component_invalid
        if (.not. inputs_are_valid(equilibrium, adiabatic_index, &
            density_kg_m3, mode_m, mode_n, stored_power, 1, degree, n_theta, &
            n_zeta)) return
        modes = size(mode_m)
        if (present(vacuum) .and. present(vacuum_energy)) then
            if (size(vacuum_energy, 1) /= 2 * modes &
                .or. size(vacuum_energy, 2) /= 2 * modes) return
        end if
        allocate (parity(modes), source=phase_cosine)
        cosine_problem%free_boundary = present(vacuum)
        sine_problem%free_boundary = present(vacuum)
        if (present(sparse_storage)) then
            cosine_problem%has_sparse_storage = sparse_storage
            sine_problem%has_sparse_storage = sparse_storage
        end if
        if (present(vacuum) .and. present(vacuum_energy)) then
            call build_trials(equilibrium, adiabatic_index, density_kg_m3, &
                mode_m, mode_n, stored_power, parity, degree, n_theta, &
                n_zeta, cosine_problem, info, vacuum, &
                vacuum_energy(:modes, :modes), sine_problem, &
                vacuum_energy(modes + 1:, modes + 1:))
        else
            call build_trials(equilibrium, adiabatic_index, density_kg_m3, &
                mode_m, mode_n, stored_power, parity, degree, n_theta, &
                n_zeta, cosine_problem, info, vacuum, &
                sine_problem=sine_problem)
        end if
    end subroutine build_compatible_three_component_classes

    ! Vacuum energy of the normal trials of both parity classes in the
    ! coupled order: every mode with parity 1, then every mode with parity
    ! 2. A trial's edge flux depends only on its own mode and parity, so the
    ! problem of one parity class uses its principal sub-block, and the
    ! vacuum solve, which depends only on the edge and the wall, serves
    ! both classes.
    subroutine build_compatible_vacuum_energy(equilibrium, mode_m, mode_n, &
            vacuum, energy, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        type(plasma_vacuum_model_t), intent(in) :: vacuum
        real(dp), allocatable, intent(out) :: energy(:, :)
        integer, intent(out) :: info
        type(primitive_equilibrium_spline_t) :: spline
        type(trial_space_topology_t) :: topology
        integer, allocatable :: parity(:), trial_m(:), trial_n(:)
        integer :: local_info, modes

        info = compatible_three_component_invalid
        modes = size(mode_m)
        if (modes < 1 .or. size(mode_n) /= modes) return
        allocate (parity(2 * modes), trial_m(2 * modes), trial_n(2 * modes))
        trial_m(:modes) = mode_m
        trial_m(modes + 1:) = mode_m
        trial_n(:modes) = mode_n
        trial_n(modes + 1:) = mode_n
        parity(:modes) = 1
        parity(modes + 1:) = 2
        call build_trial_space_topology(trial_m, trial_n, parity, topology, &
            local_info)
        if (local_info /= trial_topology_ok) return
        call fit_primitive_equilibrium(equilibrium, spline, local_info)
        if (local_info /= primitive_equilibrium_ok) then
            info = compatible_three_component_assembly_error
            return
        end if
        call build_vacuum_block(spline, vacuum, topology, energy, info)
    end subroutine build_compatible_vacuum_energy

    ! sine_problem, with cosine-parity trials: the same modes with sine
    ! parity, assembled from the same angular products in the same pass.
    subroutine build_trials(equilibrium, adiabatic_index, density_kg_m3, &
            mode_m, mode_n, stored_power, parity, degree, n_theta, n_zeta, &
            problem, info, vacuum, vacuum_energy, sine_problem, &
            sine_vacuum_energy)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: degree, n_theta, n_zeta
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        type(plasma_vacuum_model_t), optional, intent(in) :: vacuum
        real(dp), optional, intent(in) :: vacuum_energy(:, :)
        type(compatible_three_component_problem_t), optional, &
            intent(inout) :: sine_problem
        real(dp), optional, intent(in) :: sine_vacuum_energy(:, :)
        type(primitive_equilibrium_spline_t) :: spline
        type(radial_feec_complex_t) :: complex
        type(class_setup_t) :: setup, sine_setup
        real(dp), allocatable :: breaks(:), theta(:), zeta(:)
        integer, allocatable :: sine_parity(:)
        integer :: allocation_status, intervals, local_info

        info = compatible_three_component_invalid
        intervals = size(equilibrium%s)
        allocate (breaks(intervals + 1), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        call build_uniform_breaks(intervals, breaks, local_info)
        if (local_info /= compatible_support_ok) return
        call build_radial_feec_complex(breaks, degree, .true., &
            .not. present(vacuum), complex, local_info)
        if (local_info /= radial_feec_ok) return
        if (present(sine_problem)) then
            if (any(parity /= phase_cosine)) return
        end if
        call fit_primitive_equilibrium(equilibrium, spline, local_info)
        if (local_info /= primitive_equilibrium_ok) then
            info = compatible_three_component_assembly_error
            return
        end if
        call prepare_class(spline, complex, mode_m, mode_n, parity, &
            stored_power, degree, problem, setup, info, vacuum, vacuum_energy)
        if (info /= compatible_three_component_ok) return
        if (present(sine_problem)) then
            allocate (sine_parity(size(parity)), source=phase_sine)
            call prepare_class(spline, complex, mode_m, mode_n, sine_parity, &
                stored_power, degree, sine_problem, sine_setup, info, vacuum, &
                sine_vacuum_energy)
            if (info /= compatible_three_component_ok) return
        end if
        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        if (present(sine_problem)) then
            call assemble_problem(spline, complex, breaks, theta, zeta, &
                adiabatic_index, density_kg_m3, mode_m, mode_n, parity, &
                stored_power, setup, problem, info, sine_setup, sine_problem)
        else
            call assemble_problem(spline, complex, breaks, theta, zeta, &
                adiabatic_index, density_kg_m3, mode_m, mode_n, parity, &
                stored_power, setup, problem, info)
        end if
        if (info /= compatible_three_component_ok) return
        call finish_class(complex, degree, setup, problem, info)
        if (info /= compatible_three_component_ok) return
        if (present(sine_problem)) call finish_class(complex, degree, &
            sine_setup, sine_problem, info)
    end subroutine build_trials

    ! Trial topology, unknown counts, vacuum block, axis tie and empty
    ! storage of one problem; the vacuum block is checked before the
    ! costly plasma assembly.
    subroutine prepare_class(spline, complex, mode_m, mode_n, parity, &
            stored_power, degree, problem, setup, info, vacuum, vacuum_energy)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:), degree
        real(dp), intent(in) :: stored_power(:)
        type(compatible_three_component_problem_t), intent(inout) :: problem
        type(class_setup_t), intent(out) :: setup
        integer, intent(out) :: info
        type(plasma_vacuum_model_t), optional, intent(in) :: vacuum
        real(dp), optional, intent(in) :: vacuum_energy(:, :)
        integer :: allocation_status, local_info, unknowns

        info = compatible_three_component_invalid
        call build_trial_space_topology(mode_m, mode_n, parity, &
            setup%topology, local_info)
        if (local_info /= trial_topology_ok) return
        allocate (setup%ranks(3, size(mode_m)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        call build_component_ranks(setup%topology, setup%ranks)
        problem%normal_unknowns = complex%h1_dofs &
            * count(setup%topology%active(trial_component_normal, :))
        problem%eta_unknowns = complex%l2_dofs &
            * count(setup%topology%active(trial_component_eta, :))
        problem%mu_unknowns = complex%l2_dofs &
            * count(setup%topology%active(trial_component_mu, :))
        unknowns = problem%normal_unknowns + problem%eta_unknowns &
            + problem%mu_unknowns
        if (unknowns < 1) return
        if (present(vacuum)) then
            if (present(vacuum_energy)) then
                if (size(vacuum_energy, 1) /= size(mode_m) &
                    .or. size(vacuum_energy, 2) /= size(mode_m)) return
                if (.not. all(ieee_is_finite(vacuum_energy))) return
                allocate (setup%block, source=vacuum_energy)
            else
                call build_vacuum_block(spline, vacuum, setup%topology, &
                    setup%block, info)
                if (info /= compatible_three_component_ok) return
                info = compatible_three_component_invalid
            end if
        end if
        call build_trial_axis_tie(spline, complex, mode_m, parity, &
            stored_power, setup%ranks(trial_component_normal, :), &
            setup%ranks(trial_component_eta, :), unknowns, .true., &
            problem%axis_tie, local_info)
        if (local_info /= axis_regularity_ok) then
            info = compatible_three_component_assembly_error
            return
        end if
        problem%eta_unknowns = problem%eta_unknowns &
            - problem%axis_tie%eliminated_count
        unknowns = problem%axis_tie%reduced_unknowns
        problem%axis_quadrature_points = axis_quadrature_points(degree, &
            maxval(mode_m))
        if (problem%has_sparse_storage) then
            call initialize_sparse_storage(complex, setup%topology, degree, &
                problem, info)
            if (info /= compatible_three_component_ok) return
        else
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
        end if
        info = compatible_three_component_ok
    end subroutine prepare_class

    subroutine finish_class(complex, degree, setup, problem, info)
        type(radial_feec_complex_t), intent(in) :: complex
        integer, intent(in) :: degree
        type(class_setup_t), intent(in) :: setup
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info

        ! Dense storage sums the terms before the vacuum is added; sparse
        ! storage adds the vacuum when it sums them.
        if (.not. problem%has_sparse_storage) &
            call sum_tensor(problem%stiffness_terms, problem%stiffness)
        if (allocated(setup%block)) then
            call add_vacuum_edge(setup%block, setup%ranks, complex, problem, &
                info)
            if (info /= compatible_three_component_ok) return
        end if
        if (problem%has_sparse_storage) then
            call finish_sparse_storage(problem, info)
            if (info /= compatible_three_component_ok) return
        else
            call symmetrize_matrix(problem%stiffness)
            call symmetrize_matrix(problem%mass)
            call symmetrize_tensor(problem%stiffness_terms)
        end if
        problem%degree = degree
        problem%quadrature_points = size(accurate_nodes)
        problem%h1_dofs = complex%h1_dofs
        problem%l2_dofs = complex%l2_dofs
        info = compatible_three_component_ok
    end subroutine finish_class

    subroutine build_vacuum_block(spline, vacuum, topology, block, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(plasma_vacuum_model_t), intent(in) :: vacuum
        type(trial_space_topology_t), intent(in) :: topology
        real(dp), allocatable, intent(out) :: block(:, :)
        integer, intent(out) :: info
        integer :: local_info

        call build_vacuum_edge_block(spline, vacuum, topology, block, &
            local_info)
        if (local_info == plasma_vacuum_ok) then
            info = compatible_three_component_ok
        else if (local_info == plasma_vacuum_underresolved) then
            info = compatible_three_component_vacuum_mesh
        else if (local_info == plasma_vacuum_wall_not_nested) then
            info = compatible_three_component_wall
        else
            info = compatible_three_component_vacuum
        end if
    end subroutine build_vacuum_block

    subroutine add_vacuum_edge(block, ranks, complex, problem, info)
        real(dp), intent(in) :: block(:, :)
        integer, intent(in) :: ranks(:, :)
        type(radial_feec_complex_t), intent(in) :: complex
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        real(dp), allocatable :: local(:, :), scale(:)
        integer, allocatable :: full_map(:), map(:)
        integer :: a, allocation_status, b, local_info, normals, trials

        trials = size(ranks, 2)
        normals = count(ranks(trial_component_normal, :) > 0)
        allocate (full_map(trials), map(trials), scale(trials), &
            local(trials, trials), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        ! The edge is the last H1 function, the only one nonzero at s = 1.
        do a = 1, trials
            full_map(a) = 0
            if (ranks(trial_component_normal, a) > 0) full_map(a) = &
                (complex%h1_dofs - 1) * normals &
                + ranks(trial_component_normal, a)
        end do
        call tie_local_map(problem%axis_tie, full_map, map, scale)
        do b = 1, trials
            do a = 1, trials
                local(a, b) = scale(a) * block(a, b) * scale(b)
            end do
        end do
        if (problem%has_sparse_storage) then
            call scatter_symmetric_compatible_block(map, local, 1.0_dp, &
                problem%sparse_block_index, problem%sparse_local_index, &
                problem%sparse_vacuum, local_info)
            info = compatible_three_component_assembly_error
            if (local_info /= compatible_block_ok) return
            info = compatible_three_component_ok
            return
        end if
        allocate (problem%vacuum(size(problem%stiffness, 1), &
            size(problem%stiffness, 2)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_three_component_allocation_error
            return
        end if
        problem%vacuum = 0.0_dp
        call scatter_matrix(map, local, 1.0_dp, problem%vacuum)
        problem%stiffness = problem%stiffness + problem%vacuum
        info = compatible_three_component_ok
    end subroutine add_vacuum_edge

    subroutine initialize_sparse_storage(complex, topology, degree, problem, &
            info)
        type(radial_feec_complex_t), intent(in) :: complex
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: degree
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        integer, allocatable :: widths(:)
        integer :: local_info, term

        info = compatible_three_component_assembly_error
        call build_compatible_block_indices(complex%h1_dofs, &
            complex%l2_dofs, count(topology%active(trial_component_normal, :)), &
            count(topology%active(trial_component_eta, :)), degree, widths, &
            problem%sparse_block_index, problem%sparse_local_index, &
            local_info, problem%axis_tie%eliminated, &
            count(topology%active(trial_component_mu, :)))
        if (local_info /= compatible_block_ok) return
        call allocate_sparse(widths, problem%sparse_stiffness, local_info)
        if (local_info /= compatible_block_ok) return
        call allocate_sparse(widths, problem%sparse_mass, local_info)
        if (local_info /= compatible_block_ok) return
        do term = 1, compatible_stiffness_term_count
            call allocate_sparse(widths, problem%sparse_terms(term), local_info)
            if (local_info /= compatible_block_ok) return
        end do
        call allocate_sparse(widths, problem%sparse_vacuum, local_info)
        if (local_info /= compatible_block_ok) return
        info = compatible_three_component_ok

    contains

        subroutine allocate_sparse(block_widths, blocks, status)
            integer, intent(in) :: block_widths(:)
            type(variable_block_tridiagonal_t), intent(out) :: blocks
            integer, intent(out) :: status

            call allocate_compatible_blocks(block_widths, blocks, status)
            if (status == compatible_block_allocation) &
                info = compatible_three_component_allocation_error
        end subroutine allocate_sparse
    end subroutine initialize_sparse_storage

    ! Stiffness = sum of the terms and the vacuum; every block pencil is
    ! symmetrized like the dense arrays.
    subroutine finish_sparse_storage(problem, info)
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        integer :: block, term

        do term = 1, compatible_stiffness_term_count
            call symmetrize_compatible_blocks(problem%sparse_terms(term))
        end do
        call symmetrize_compatible_blocks(problem%sparse_vacuum)
        call symmetrize_compatible_blocks(problem%sparse_mass)
        do block = 1, size(problem%sparse_stiffness%diagonal)
            problem%sparse_stiffness%diagonal(block)%values = &
                problem%sparse_vacuum%diagonal(block)%values
            do term = 1, compatible_stiffness_term_count
                problem%sparse_stiffness%diagonal(block)%values = &
                    problem%sparse_stiffness%diagonal(block)%values &
                    + problem%sparse_terms(term)%diagonal(block)%values
            end do
        end do
        do block = 1, size(problem%sparse_stiffness%lower)
            problem%sparse_stiffness%lower(block)%values = &
                problem%sparse_vacuum%lower(block)%values
            do term = 1, compatible_stiffness_term_count
                problem%sparse_stiffness%lower(block)%values = &
                    problem%sparse_stiffness%lower(block)%values &
                    + problem%sparse_terms(term)%lower(block)%values
            end do
        end do
        info = compatible_three_component_ok
    end subroutine finish_sparse_storage

    subroutine assemble_problem(spline, complex, breaks, theta, zeta, &
            adiabatic_index, density, mode_m, mode_n, parity, stored_power, &
            setup, problem, info, sine_setup, sine_problem)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:), theta(:), zeta(:)
        real(dp), intent(in) :: adiabatic_index, density
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(class_setup_t), intent(in) :: setup
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        type(class_setup_t), optional, intent(in) :: sine_setup
        type(compatible_three_component_problem_t), optional, &
            intent(inout) :: sine_problem
        type(radial_contribution_t), allocatable :: batch(:)
        real(dp), allocatable :: constraint_nodes(:), constraint_weights(:)
        real(dp), allocatable :: axis_nodes(:), axis_weights(:)
        real(dp), allocatable :: coordinates(:), weights(:)
        integer, allocatable :: kinds(:)
        real(dp) :: half_width, midpoint
        integer :: batch_size, cell, count, entry, first, last, orientation
        integer :: point, threads

        info = compatible_three_component_assembly_error
        call build_constraint_quadrature(complex%h1_degree, &
            constraint_nodes, constraint_weights, info)
        if (info /= compatible_quadrature_ok) return
        call build_axis_quadrature(breaks(2) - breaks(1), &
            problem%axis_quadrature_points, axis_nodes, axis_weights, info)
        if (info /= compatible_quadrature_ok) return
        ! Radial points in assembly order. The axis element takes every term
        ! at every point of the rule in sqrt(s), which integrates the regular
        ! trial products exactly; each later cell takes the accurate terms
        ! and the mass at the accurate nodes and the constraint terms at the
        ! constraint nodes.
        count = size(axis_nodes) + (size(breaks) - 2) &
            * (size(accurate_nodes) + size(constraint_nodes))
        allocate (coordinates(count), weights(count), kinds(count))
        coordinates(:size(axis_nodes)) = axis_nodes
        weights(:size(axis_nodes)) = axis_weights
        kinds(:size(axis_nodes)) = point_axis
        entry = size(axis_nodes)
        do cell = 2, size(breaks) - 1
            midpoint = 0.5_dp * (breaks(cell) + breaks(cell + 1))
            half_width = 0.5_dp * (breaks(cell + 1) - breaks(cell))
            do point = 1, size(accurate_nodes)
                entry = entry + 1
                coordinates(entry) = midpoint + half_width * accurate_nodes(point)
                weights(entry) = half_width * accurate_weights(point)
                kinds(entry) = point_accurate
            end do
            do point = 1, size(constraint_nodes)
                entry = entry + 1
                coordinates(entry) = midpoint &
                    + half_width * constraint_nodes(point)
                weights(entry) = half_width * constraint_weights(point)
                kinds(entry) = point_constraint
            end do
        end do
        ! The local matrices of a batch of points are computed in parallel
        ! and scattered in assembly order, so the sums are those of the
        ! serial loop. The first point fixes the chart orientation that the
        ! others check.
        threads = 1
        !$ threads = omp_get_max_threads()
        batch_size = 2 * threads
        allocate (batch(batch_size))
        orientation = 0
        first = 1
        do while (first <= count)
            last = min(count, first + batch_size - 1)
            if (first == 1) last = 1
            !$omp parallel do schedule(dynamic) firstprivate(orientation) &
            !$omp if (last > first)
            do entry = first, last
                if (present(sine_setup)) then
                    call compute_radial_point(spline, complex, &
                        coordinates(entry), weights(entry), theta, zeta, &
                        adiabatic_index, density, mode_m, mode_n, parity, &
                        stored_power, setup, problem, &
                        kinds(entry) /= point_constraint, orientation, &
                        batch(entry - first + 1), sine_setup, sine_problem)
                else
                    call compute_radial_point(spline, complex, &
                        coordinates(entry), weights(entry), theta, zeta, &
                        adiabatic_index, density, mode_m, mode_n, parity, &
                        stored_power, setup, problem, &
                        kinds(entry) /= point_constraint, orientation, &
                        batch(entry - first + 1))
                end if
            end do
            !$omp end parallel do
            if (first == 1) orientation = batch(1)%orientation
            do entry = first, last
                info = batch(entry - first + 1)%info
                if (info /= compatible_three_component_ok) return
                call scatter_point(batch(entry - first + 1)%cosine, &
                    kinds(entry), problem, info)
                if (info /= compatible_three_component_ok) return
                if (present(sine_problem)) then
                    call scatter_point(batch(entry - first + 1)%sine, &
                        kinds(entry), sine_problem, info)
                    if (info /= compatible_three_component_ok) return
                end if
            end do
            first = last + 1
        end do
        info = compatible_three_component_ok
    end subroutine assemble_problem

    ! Local matrices of one radial quadrature point, tied to the axis
    ! constraint, with the map onto the global unknowns; with sine_setup,
    ! also those of the sine-parity problem from the same products.
    subroutine compute_radial_point(spline, complex, coordinate, weight, &
            theta, zeta, adiabatic_index, density, mode_m, mode_n, parity, &
            stored_power, setup, problem, assemble_mass, orientation, &
            contribution, sine_setup, sine_problem)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: coordinate, weight, theta(:), zeta(:)
        real(dp), intent(in) :: adiabatic_index, density
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(class_setup_t), intent(in) :: setup
        type(compatible_three_component_problem_t), intent(in) :: problem
        logical, intent(in) :: assemble_mass
        integer, intent(inout) :: orientation
        type(radial_contribution_t), intent(out) :: contribution
        type(class_setup_t), optional, intent(in) :: sine_setup
        type(compatible_three_component_problem_t), optional, intent(in) :: &
            sine_problem
        real(dp), allocatable :: fields(:, :, :), drive(:, :)
        real(dp), allocatable :: jacobian_s(:, :), jacobian_t(:, :)
        real(dp), allocatable :: jacobian_z(:, :), gamma_p(:, :)
        real(dp), allocatable :: h1(:), dh1(:), l2(:), local_h1(:, :)
        real(dp), allocatable :: local_dh1(:, :), local_l2(:, :)
        real(dp), allocatable :: local_eta(:, :), local_k(:, :), sine_k(:, :)
        integer, allocatable :: h1_index(:), l2_index(:)
        real(dp) :: pressure
        integer :: allocation_status, local_info, trials

        contribution%info = compatible_three_component_assembly_error
        call evaluate_radial_feec_complex(complex, coordinate, h1, dh1, l2, &
            local_info)
        if (local_info /= radial_feec_ok) return
        call build_active_indices(h1, h1_index, local_info, dh1)
        if (local_info == compatible_support_allocation) then
            contribution%info = compatible_three_component_allocation_error
            return
        else if (local_info /= compatible_support_ok) then
            return
        end if
        call build_active_indices(l2, l2_index, local_info)
        if (local_info == compatible_support_allocation) then
            contribution%info = compatible_three_component_allocation_error
            return
        else if (local_info /= compatible_support_ok) then
            return
        end if
        if (size(h1_index) < 1 .or. size(l2_index) < 1) return
        trials = size(mode_m)
        allocate (local_h1(size(h1_index), trials), &
            local_dh1(size(h1_index), trials), &
            local_l2(size(l2_index), trials), &
            local_eta(size(l2_index), trials), stat=allocation_status)
        if (allocation_status /= 0) then
            contribution%info = compatible_three_component_allocation_error
            return
        end if
        call apply_stored_power(coordinate, stored_power, h1, dh1, h1_index, &
            local_h1, local_dh1, local_info)
        if (local_info /= compatible_support_ok) return
        call apply_tangential_axis_weight(coordinate, stored_power, l2, &
            l2_index, local_eta, local_info)
        if (local_info /= compatible_support_ok) return
        call replicate_indexed_values(l2, l2_index, local_l2, local_info)
        if (local_info /= compatible_support_ok) return
        call evaluate_primitive_kernel_surface(spline, coordinate, theta, &
            zeta, fields, drive, local_info, jacobian_s, jacobian_t, &
            jacobian_z, pressure, orientation=orientation)
        contribution%orientation = orientation
        if (local_info /= primitive_kernel_ok) return
        if (.not. problem%coupled) then
            if (.not. surface_preserves_parity(fields, drive, jacobian_s, &
                jacobian_t, jacobian_z)) then
                contribution%info = compatible_three_component_asymmetric
                return
            end if
        end if
        allocate (gamma_p(size(theta), size(zeta)), &
            source=adiabatic_index * pressure, stat=allocation_status)
        if (allocation_status /= 0) then
            contribution%info = compatible_three_component_allocation_error
            return
        end if
        call allocate_local_matrices(trials, size(h1_index), size(l2_index), &
            local_k, contribution%cosine%mass, contribution%cosine%terms, &
            allocation_status)
        if (allocation_status == 0 .and. present(sine_setup)) &
            call allocate_local_matrices(trials, size(h1_index), &
            size(l2_index), sine_k, contribution%sine%mass, &
            contribution%sine%terms, allocation_status)
        if (allocation_status /= 0) then
            contribution%info = compatible_three_component_allocation_error
            return
        end if
        if (present(sine_setup)) then
            call assemble_compatible_compressible_stiffness_surface(fields, &
                drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, &
                mode_n, parity, spline%field_periods, local_h1, local_dh1, &
                local_eta, local_l2, weight, phase_assembly_transformed, &
                local_k, local_info, contribution%cosine%terms, sine_k, &
                contribution%sine%terms)
        else
            call assemble_compatible_compressible_stiffness_surface(fields, &
                drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, &
                mode_n, parity, spline%field_periods, local_h1, local_dh1, &
                local_eta, local_l2, weight, phase_assembly_transformed, &
                local_k, local_info, contribution%cosine%terms)
        end if
        if (local_info /= 0) return
        if (assemble_mass) then
            if (present(sine_setup)) then
                call assemble_compatible_physical_mass_surface(fields, &
                    density, mode_m, mode_n, parity, spline%field_periods, &
                    local_h1, local_eta, local_l2, weight, &
                    phase_assembly_transformed, contribution%cosine%mass, &
                    local_info, contribution%sine%mass)
            else
                call assemble_compatible_physical_mass_surface(fields, &
                    density, mode_m, mode_n, parity, spline%field_periods, &
                    local_h1, local_eta, local_l2, weight, &
                    phase_assembly_transformed, contribution%cosine%mass, &
                    local_info)
            end if
            if (local_info /= 0) return
        end if
        call tie_local_matrices(complex, setup, problem%axis_tie, h1_index, &
            l2_index, assemble_mass, contribution%cosine)
        if (present(sine_setup)) call tie_local_matrices(complex, &
            sine_setup, sine_problem%axis_tie, h1_index, l2_index, &
            assemble_mass, contribution%sine)
        contribution%info = compatible_three_component_ok
    end subroutine compute_radial_point

    ! Map the local matrices onto the global unknowns of one problem and
    ! tie them to its axis constraint.
    subroutine tie_local_matrices(complex, setup, axis_tie, h1_index, &
            l2_index, assemble_mass, local)
        type(radial_feec_complex_t), intent(in) :: complex
        type(class_setup_t), intent(in) :: setup
        type(axis_tie_t), intent(in) :: axis_tie
        integer, intent(in) :: h1_index(:), l2_index(:)
        logical, intent(in) :: assemble_mass
        type(local_matrices_t), intent(inout) :: local
        real(dp), allocatable :: tie_scale(:)
        integer, allocatable :: full_map(:)

        call build_local_map(complex, setup%topology, setup%ranks, h1_index, &
            l2_index, full_map)
        allocate (local%map(size(full_map)), tie_scale(size(full_map)))
        call tie_local_map(axis_tie, full_map, local%map, tie_scale)
        call scale_local_tensor(tie_scale, local%terms)
        if (assemble_mass) call scale_local_matrix(tie_scale, local%mass)
        local%has_mass = assemble_mass
    end subroutine tie_local_matrices

    subroutine scatter_point(local, kind, problem, info)
        type(local_matrices_t), intent(in) :: local
        integer, intent(in) :: kind
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info

        select case (kind)
        case (point_axis)
            call scatter_local(local, all_terms, problem, info)
        case (point_accurate)
            call scatter_local(local, accurate_term, problem, info)
        case default
            call scatter_local(local, constraint_term, problem, info)
        end select
    end subroutine scatter_point

    subroutine scatter_local(local, term_mask, problem, info)
        type(local_matrices_t), intent(in) :: local
        logical, intent(in) :: term_mask(:)
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        integer :: local_info

        info = compatible_three_component_assembly_error
        if (problem%has_sparse_storage) then
            if (local%has_mass) then
                call scatter_symmetric_compatible_block(local%map, &
                    local%mass, 1.0_dp, problem%sparse_block_index, &
                    problem%sparse_local_index, problem%sparse_mass, &
                    local_info)
                if (local_info /= compatible_block_ok) return
            end if
            call scatter_sparse_terms(local%map, local%terms, term_mask, &
                problem, local_info)
            if (local_info /= compatible_block_ok) return
        else
            if (local%has_mass) call scatter_matrix(local%map, local%mass, &
                1.0_dp, problem%mass)
            call scatter_terms(local%map, local%terms, term_mask, &
                problem%stiffness_terms)
        end if
        info = compatible_three_component_ok
    end subroutine scatter_local

    subroutine scatter_sparse_terms(map, local, term_mask, problem, info)
        integer, intent(in) :: map(:)
        real(dp), intent(in) :: local(:, :, :)
        logical, intent(in) :: term_mask(:)
        type(compatible_three_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        integer :: term

        info = compatible_block_ok
        do term = 1, compatible_stiffness_term_count
            if (.not. term_mask(term)) cycle
            call scatter_symmetric_compatible_block(map, local(:, :, term), &
                1.0_dp, problem%sparse_block_index, &
                problem%sparse_local_index, problem%sparse_terms(term), info)
            if (info /= compatible_block_ok) return
        end do
    end subroutine scatter_sparse_terms

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

    pure subroutine scale_local_matrix(scale, matrix)
        real(dp), intent(in) :: scale(:)
        real(dp), intent(inout) :: matrix(:, :)
        integer :: a, b

        do b = 1, size(matrix, 2)
            do a = 1, size(matrix, 1)
                matrix(a, b) = scale(a) * scale(b) * matrix(a, b)
            end do
        end do
    end subroutine scale_local_matrix

    pure subroutine scale_local_tensor(scale, tensor)
        real(dp), intent(in) :: scale(:)
        real(dp), intent(inout) :: tensor(:, :, :)
        integer :: term

        do term = 1, size(tensor, 3)
            call scale_local_matrix(scale, tensor(:, :, term))
        end do
    end subroutine scale_local_tensor

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
        valid = valid .and. parity_class >= 0 .and. parity_class <= 2
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
