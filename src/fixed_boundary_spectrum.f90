module fixed_boundary_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_problem_assembly_support, only: angular_grid_aliases
    use compatible_three_component_problem, only: &
        build_compatible_three_component_classes, &
        build_compatible_three_component_problem, &
        build_compatible_vacuum_energy, &
        compatible_three_component_allocation_error, &
        compatible_three_component_asymmetric, &
        compatible_three_component_invalid, compatible_three_component_ok, &
        compatible_three_component_problem_t, &
        compatible_three_component_vacuum, &
        compatible_three_component_vacuum_mesh, &
        compatible_three_component_wall
    use dense_spectrum_support, only: certify_dense_spectrum_inertia, &
        certify_dense_spectrum_orthogonality, dense_spectrum_allocation, &
        dense_spectrum_is_certified, dense_spectrum_ok, &
        diagnose_dense_spectrum, refine_dense_spectrum, &
        unpermute_dense_vectors
    use fixed_boundary_energy, only: diagnose_fixed_boundary_energy_store, &
        fixed_boundary_energy_allocation, fixed_boundary_energy_invalid, &
        fixed_boundary_energy_ok, fixed_boundary_energy_store_t, &
        fixed_boundary_energy_term_count, fixed_boundary_energy_terms_t, &
        rayleigh_gradient_fixed_boundary_store
    use fixed_boundary_eigen_bracket, only: bracket_lowest_negative, &
        fixed_boundary_bracket_ok, bracket_lowest_positive, &
        certify_lowest_bracket
    use fixed_boundary_solver_controls, only: &
        fixed_boundary_solver_controls_t, valid_fixed_boundary_solver_controls
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use plasma_vacuum_boundary, only: plasma_vacuum_model_t, &
        valid_plasma_vacuum_model
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated, &
        symmetric_eigensolver_allocation, symmetric_eigensolver_ok
    use variable_block_tridiagonal, only: &
        variable_block_allocation, variable_block_ok, &
        variable_block_to_dense, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_eigenvalue_bound, &
        variable_generalized_ok
    use variable_spectrum_analysis, only: analyze_variable_spectrum, &
        variable_spectrum_ok, variable_spectrum_summary_t
    implicit none
    private

    integer, parameter, public :: fixed_boundary_ok = 0
    integer, parameter, public :: fixed_boundary_invalid = -1
    integer, parameter, public :: fixed_boundary_geometry_error = -2
    integer, parameter, public :: fixed_boundary_assembly_error = -3
    integer, parameter, public :: fixed_boundary_solver_error = -4
    integer, parameter, public :: fixed_boundary_allocation_error = -5
    integer, parameter, public :: fixed_boundary_asymmetric = -6
    ! Free boundary: the edge mesh does not resolve the mode table, the
    ! wall does not enclose the plasma, or the vacuum block is singular.
    integer, parameter, public :: fixed_boundary_vacuum_mesh = -7
    integer, parameter, public :: fixed_boundary_wall = -8
    integer, parameter, public :: fixed_boundary_vacuum = -9
    integer, parameter, public :: fixed_boundary_n_theta = 64
    integer, parameter, public :: fixed_boundary_n_zeta = 64

    type :: fixed_boundary_class_problem_t
        integer :: unknowns = 0
        integer :: normal_unknowns = 0
        integer :: eta_unknowns = 0
        integer :: mu_unknowns = 0
        integer, allocatable :: permutation(:)
        type(variable_block_tridiagonal_t) :: stiffness
        type(variable_block_tridiagonal_t) :: mass
        type(fixed_boundary_energy_store_t) :: energy
    end type fixed_boundary_class_problem_t

    type, public :: fixed_boundary_problem_t
        private
        logical :: ready = .false.
        logical :: has_chart_metric = .false.
        integer :: field_periods = 0
        integer :: degree = 0
        integer :: n_theta = fixed_boundary_n_theta
        integer :: n_zeta = fixed_boundary_n_zeta
        real(dp) :: adiabatic_index = 0.0_dp
        real(dp) :: density_kg_m3 = 0.0_dp
        real(dp) :: zero_floor = 0.0_dp
        type(fixed_boundary_solver_controls_t) :: solver_controls
        integer, allocatable :: mode_m(:)
        integer, allocatable :: mode_n(:)
        ! An equilibrium without stellarator symmetry couples both Fourier
        ! parities of every mode; its single class is parity class 0 and
        ! lives in classes(1).
        logical :: coupled = .false.
        ! The plasma edge couples to the vacuum (and wall) instead of
        ! being held fixed.
        logical :: free_boundary = .false.
        type(fixed_boundary_class_problem_t) :: classes(2)
    end type fixed_boundary_problem_t

    type, public :: fixed_boundary_spectrum_result_t
        logical :: has_chart_metric = .false.
        logical :: has_eigenvector = .false.
        integer :: field_periods = 0
        integer :: mode_count = 0
        integer :: parity_class = 0
        integer :: degree = 0
        integer :: angular_theta = 0
        integer :: angular_zeta = 0
        integer :: unknowns = 0
        integer :: normal_unknowns = 0
        integer :: eta_unknowns = 0
        integer :: mu_unknowns = 0
        integer :: negative_count = 0
        integer :: floor_count = 0
        real(dp) :: adiabatic_index = 0.0_dp
        real(dp) :: density_kg_m3 = 0.0_dp
        real(dp) :: zero_floor = 0.0_dp
        real(dp) :: lowest_eigenvalue = 0.0_dp
        real(dp) :: certificate = 0.0_dp
        real(dp) :: eigenpair_residual = 0.0_dp
        real(dp) :: eigenpair_resolution = 0.0_dp
        real(dp) :: inertia_interval = 0.0_dp
        real(dp), allocatable :: eigenvector(:)
    end type fixed_boundary_spectrum_result_t

    type, public :: fixed_boundary_full_spectrum_t
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: rayleigh_quotients(:)
        real(dp), allocatable :: residuals(:)
        real(dp), allocatable :: resolutions(:)
    end type fixed_boundary_full_spectrum_t

    public :: build_fixed_boundary_problem, diagnose_fixed_boundary_energy
    public :: fixed_boundary_is_coupled, fixed_boundary_is_free
    public :: fixed_boundary_energy_terms_t, fixed_boundary_unknown_count
    public :: fixed_boundary_rayleigh_gradient
    public :: set_fixed_boundary_solver_controls, solve_fixed_boundary_class
    public :: solve_fixed_boundary_full_spectrum

contains

    subroutine build_fixed_boundary_problem(equilibrium, adiabatic_index, &
            density_kg_m3, zero_floor, mode_m, mode_n, degree, problem, info, &
            angular_theta, angular_zeta, coupled, vacuum)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3, zero_floor
        integer, intent(in) :: mode_m(:), mode_n(:), degree
        type(fixed_boundary_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        integer, optional, intent(in) :: angular_theta, angular_zeta
        ! Force the coupled operator for a symmetric equilibrium; it then
        ! contains both parity classes as one problem.
        logical, optional, intent(in) :: coupled
        ! Present: the physical free-boundary problem with this vacuum model.
        type(plasma_vacuum_model_t), optional, intent(in) :: vacuum
        real(dp), allocatable :: stored_power(:)
        ! Default-initialized: holds no matrices.
        type(fixed_boundary_class_problem_t) :: empty_class
        real(dp), allocatable :: vacuum_energy(:, :)
        integer :: allocation_status, compatible_info, mode

        info = fixed_boundary_invalid
        if (present(angular_theta)) problem%n_theta = angular_theta
        if (present(angular_zeta)) problem%n_zeta = angular_zeta
        if (min(problem%n_theta, problem%n_zeta) < 1) return
        if (problem%n_theta > huge(1) / problem%n_zeta) return
        if (.not. valid_inputs(equilibrium, adiabatic_index, density_kg_m3, &
            zero_floor, mode_m, mode_n, degree, problem%n_theta, problem%n_zeta)) return
        if (present(vacuum)) then
            if (.not. valid_plasma_vacuum_model(vacuum)) return
        end if
        problem%free_boundary = present(vacuum)
        allocate (stored_power(size(mode_m)), problem%mode_m(size(mode_m)), &
            problem%mode_n(size(mode_n)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = fixed_boundary_allocation_error
            return
        end if
        do mode = 1, size(mode_m)
            stored_power(mode) = 0.0_dp
            if (mode_m(mode) > 0) stored_power(mode) = &
                1.0_dp - 0.5_dp * real(mode_m(mode), dp)
            problem%mode_m(mode) = mode_m(mode)
            problem%mode_n(mode) = mode_n(mode)
        end do
        ! An equilibrium whose reconstructed operator breaks the parity
        ! symmetry takes the coupled operator: the decoupled classes would
        ! silently drop the coupling. The storage flag of the file is not
        ! evidence; the admission test at every assembly point is.
        problem%coupled = .false.
        if (present(coupled)) problem%coupled = coupled
        ! One vacuum solve serves both parity classes and the coupled
        ! operator: each takes its principal sub-block.
        if (present(vacuum)) then
            call build_compatible_vacuum_energy(equilibrium, mode_m, mode_n, &
                vacuum, vacuum_energy, compatible_info)
            if (compatible_info /= compatible_three_component_ok) then
                info = class_status(compatible_info)
                return
            end if
        end if
        if (.not. problem%coupled) then
            ! Both parity classes come from one pass over the radial points.
            call assemble_classes(equilibrium, adiabatic_index, &
                density_kg_m3, mode_m, mode_n, stored_power, degree, &
                problem%classes, info, problem%n_theta, problem%n_zeta, &
                vacuum, vacuum_energy)
            if (info /= fixed_boundary_ok .and. &
                info /= fixed_boundary_asymmetric) return
            problem%coupled = info == fixed_boundary_asymmetric
        end if
        if (problem%coupled) then
            ! Discard the partial decoupled classes.
            problem%classes = empty_class
            if (present(vacuum)) then
                call assemble_class(equilibrium, adiabatic_index, &
                    density_kg_m3, mode_m, mode_n, stored_power, 0, degree, &
                    problem%classes(1), info, problem%n_theta, &
                    problem%n_zeta, vacuum, vacuum_energy)
            else
                call assemble_class(equilibrium, adiabatic_index, &
                    density_kg_m3, mode_m, mode_n, stored_power, 0, degree, &
                    problem%classes(1), info, problem%n_theta, problem%n_zeta)
            end if
            if (info /= fixed_boundary_ok) return
        end if
        problem%has_chart_metric = equilibrium%has_chart_metric
        problem%field_periods = equilibrium%field_periods
        problem%degree = degree
        problem%adiabatic_index = adiabatic_index
        problem%density_kg_m3 = density_kg_m3
        problem%zero_floor = zero_floor
        problem%ready = .true.
        info = fixed_boundary_ok
    end subroutine build_fixed_boundary_problem

    subroutine assemble_class(equilibrium, adiabatic_index, density_kg_m3, &
            mode_m, mode_n, stored_power, parity_class, degree, &
            class_problem, info, n_theta, n_zeta, vacuum, vacuum_energy)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:), parity_class, degree
        real(dp), intent(in) :: stored_power(:)
        type(fixed_boundary_class_problem_t), intent(out) :: class_problem
        integer, intent(out) :: info
        type(compatible_three_component_problem_t) :: compatible
        integer, intent(in) :: n_theta, n_zeta
        type(plasma_vacuum_model_t), optional, intent(in) :: vacuum
        real(dp), optional, intent(in) :: vacuum_energy(:, :)
        integer :: compatible_info

        call build_compatible_three_component_problem(equilibrium, &
            adiabatic_index, density_kg_m3, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, &
            n_zeta, compatible, compatible_info, vacuum, sparse_storage=.true., &
            vacuum_energy=vacuum_energy)
        if (compatible_info /= compatible_three_component_ok) then
            info = class_status(compatible_info)
            return
        end if
        call pack_class_problem(compatible, class_problem, info)
    end subroutine assemble_class

    subroutine assemble_classes(equilibrium, adiabatic_index, &
            density_kg_m3, mode_m, mode_n, stored_power, degree, classes, &
            info, n_theta, n_zeta, vacuum, vacuum_energy)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:), degree
        real(dp), intent(in) :: stored_power(:)
        type(fixed_boundary_class_problem_t), intent(inout) :: classes(2)
        integer, intent(out) :: info
        integer, intent(in) :: n_theta, n_zeta
        type(plasma_vacuum_model_t), optional, intent(in) :: vacuum
        real(dp), allocatable, intent(in) :: vacuum_energy(:, :)
        type(compatible_three_component_problem_t) :: cosine, sine
        integer :: compatible_info

        if (allocated(vacuum_energy)) then
            call build_compatible_three_component_classes(equilibrium, &
                adiabatic_index, density_kg_m3, mode_m, mode_n, &
                stored_power, degree, n_theta, n_zeta, cosine, sine, &
                compatible_info, vacuum, sparse_storage=.true., &
                vacuum_energy=vacuum_energy)
        else
            call build_compatible_three_component_classes(equilibrium, &
                adiabatic_index, density_kg_m3, mode_m, mode_n, &
                stored_power, degree, n_theta, n_zeta, cosine, sine, &
                compatible_info, vacuum, sparse_storage=.true.)
        end if
        info = class_status(compatible_info)
        if (info /= fixed_boundary_ok) return
        call pack_class_problem(cosine, classes(1), info)
        if (info /= fixed_boundary_ok) return
        call pack_class_problem(sine, classes(2), info)
    end subroutine assemble_classes

    pure integer function class_status(compatible_info) result(info)
        integer, intent(in) :: compatible_info

        select case (compatible_info)
        case (compatible_three_component_ok)
            info = fixed_boundary_ok
        case (compatible_three_component_allocation_error)
            info = fixed_boundary_allocation_error
        case (compatible_three_component_invalid)
            info = fixed_boundary_invalid
        case (compatible_three_component_asymmetric)
            info = fixed_boundary_asymmetric
        case (compatible_three_component_vacuum_mesh)
            info = fixed_boundary_vacuum_mesh
        case (compatible_three_component_wall)
            info = fixed_boundary_wall
        case (compatible_three_component_vacuum)
            info = fixed_boundary_vacuum
        case default
            info = fixed_boundary_assembly_error
        end select
    end function class_status

    ! The block pencil is adopted as assembled; permutation(p) is the
    ! assembly-order unknown stored at block-order position p.
    subroutine pack_class_problem(compatible, class_problem, info)
        type(compatible_three_component_problem_t), intent(in) :: compatible
        type(fixed_boundary_class_problem_t), intent(out) :: class_problem
        integer, intent(out) :: info
        integer, allocatable :: start(:)
        integer :: allocation_status, block, index, term, unknowns

        info = fixed_boundary_assembly_error
        if (.not. compatible%has_sparse_storage) return
        unknowns = size(compatible%sparse_block_index)
        allocate (class_problem%permutation(unknowns), &
            start(size(compatible%sparse_stiffness%widths)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = fixed_boundary_allocation_error
            return
        end if
        start(1) = 0
        do block = 2, size(start)
            start(block) = start(block - 1) &
                + compatible%sparse_stiffness%widths(block - 1)
        end do
        if (start(size(start)) + compatible%sparse_stiffness%widths( &
            size(start)) /= unknowns) return
        class_problem%permutation = 0
        do index = 1, unknowns
            class_problem%permutation(start(compatible%sparse_block_index( &
                index)) + compatible%sparse_local_index(index)) = index
        end do
        if (any(class_problem%permutation == 0)) return
        class_problem%stiffness = compatible%sparse_stiffness
        class_problem%mass = compatible%sparse_mass
        do term = 1, fixed_boundary_energy_term_count - 1
            class_problem%energy%terms(term) = compatible%sparse_terms(term)
        end do
        class_problem%energy%terms(fixed_boundary_energy_term_count) = &
            compatible%sparse_vacuum
        class_problem%unknowns = unknowns
        class_problem%normal_unknowns = compatible%normal_unknowns
        class_problem%eta_unknowns = compatible%eta_unknowns
        class_problem%mu_unknowns = compatible%mu_unknowns
        info = fixed_boundary_ok
    end subroutine pack_class_problem

    function valid_inputs(equilibrium, adiabatic_index, density_kg_m3, &
            zero_floor, mode_m, mode_n, degree, n_theta, n_zeta) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3, zero_floor
        integer, intent(in) :: mode_m(:), mode_n(:), degree
        integer, intent(in) :: n_theta, n_zeta
        logical :: valid
        integer :: first, second

        valid = .false.
        if (.not. ieee_is_finite(adiabatic_index) &
            .or. adiabatic_index <= 0.0_dp) return
        if (.not. ieee_is_finite(density_kg_m3) &
            .or. density_kg_m3 <= 0.0_dp) return
        if (.not. ieee_is_finite(zero_floor) .or. zero_floor <= 0.0_dp) return
        if (degree < 1 .or. degree > 4) return
        if (equilibrium%field_periods < 1) return
        if (size(mode_m) < 1 .or. size(mode_n) /= size(mode_m)) return
        do first = 1, size(mode_m)
            if (mode_m(first) < 0) return
            if (mode_m(first) == 0 .and. mode_n(first) < 0) return
            do second = 1, first - 1
                if (mode_m(first) == mode_m(second) &
                    .and. mode_n(first) == mode_n(second)) return
            end do
        end do
        if (angular_grid_aliases(equilibrium, mode_m, mode_n, &
            n_theta, n_zeta)) return
        valid = .true.
    end function valid_inputs

    subroutine diagnose_fixed_boundary_energy(problem, parity_class, vector, &
            result, info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        real(dp), intent(in) :: vector(:)
        type(fixed_boundary_energy_terms_t), intent(out) :: result
        integer, intent(out) :: info
        integer :: slot

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        slot = class_slot(problem, parity_class)
        if (slot == 0) return
        call diagnose_fixed_boundary_energy_store( &
            problem%classes(slot)%stiffness, &
            problem%classes(slot)%mass, &
            problem%classes(slot)%energy, &
            problem%classes(slot)%permutation, vector, result, info)
        info = map_fixed_boundary_energy_info(info)
    end subroutine diagnose_fixed_boundary_energy

    subroutine fixed_boundary_rayleigh_gradient(problem, parity_class, &
            vector, gradient, info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        real(dp), intent(in) :: vector(:)
        real(dp), allocatable, intent(out) :: gradient(:)
        integer, intent(out) :: info
        integer :: slot

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        slot = class_slot(problem, parity_class)
        if (slot == 0) return
        call rayleigh_gradient_fixed_boundary_store( &
            problem%classes(slot)%stiffness, &
            problem%classes(slot)%mass, &
            problem%classes(slot)%permutation, vector, gradient, info)
        info = map_fixed_boundary_energy_info(info)
    end subroutine fixed_boundary_rayleigh_gradient

    pure function map_fixed_boundary_energy_info(info) result(status)
        integer, intent(in) :: info
        integer :: status

        status = fixed_boundary_solver_error
        if (info == fixed_boundary_energy_ok) status = fixed_boundary_ok
        if (info == fixed_boundary_energy_allocation) &
            status = fixed_boundary_allocation_error
        if (info == fixed_boundary_energy_invalid) status = fixed_boundary_invalid
    end function map_fixed_boundary_energy_info

    subroutine solve_fixed_boundary_class(problem, parity_class, result, info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        type(fixed_boundary_spectrum_result_t), intent(out) :: result
        integer, intent(out) :: info
        type(variable_spectrum_summary_t) :: summary
        integer :: slot

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        slot = class_slot(problem, parity_class)
        if (slot == 0) return
        call initialize_result(problem, parity_class, result)
        call analyze_variable_spectrum( &
            problem%classes(slot)%stiffness, &
            problem%classes(slot)%mass, problem%zero_floor, summary, &
            info)
        if (info /= variable_spectrum_ok) then
            info = fixed_boundary_solver_error
            return
        end if
        result%negative_count = summary%negative_count
        result%floor_count = summary%zero_count
        call resolve_lowest(problem%classes(slot), summary, &
            problem%solver_controls, result, info)
        if (info /= fixed_boundary_ok) return
        result%certificate = result%inertia_interval &
            + result%eigenpair_residual + result%eigenpair_resolution
    end subroutine solve_fixed_boundary_class

    subroutine set_fixed_boundary_solver_controls(problem, controls, info)
        type(fixed_boundary_problem_t), intent(inout) :: problem
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        integer, intent(out) :: info

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        if (.not. valid_fixed_boundary_solver_controls(controls)) return
        problem%solver_controls = controls
        info = fixed_boundary_ok
    end subroutine set_fixed_boundary_solver_controls

    pure function class_slot(problem, parity_class) result(slot)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        integer :: slot

        ! Parity classes 1 and 2 of a decoupled problem, or class 0 of a
        ! coupled one; zero rejects the request.
        slot = 0
        if (problem%coupled) then
            if (parity_class == 0) slot = 1
        else if (parity_class == 1 .or. parity_class == 2) then
            slot = parity_class
        end if
    end function class_slot

    pure logical function fixed_boundary_is_free(problem) result(free)
        type(fixed_boundary_problem_t), intent(in) :: problem

        free = problem%free_boundary
    end function fixed_boundary_is_free

    pure logical function fixed_boundary_is_coupled(problem) result(coupled)
        type(fixed_boundary_problem_t), intent(in) :: problem

        coupled = problem%ready .and. problem%coupled
    end function fixed_boundary_is_coupled

    subroutine initialize_result(problem, parity_class, result)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        type(fixed_boundary_spectrum_result_t), intent(out) :: result
        integer :: slot

        slot = class_slot(problem, parity_class)
        result%has_chart_metric = problem%has_chart_metric
        result%field_periods = problem%field_periods
        result%mode_count = size(problem%mode_m)
        result%parity_class = parity_class
        result%degree = problem%degree
        result%angular_theta = problem%n_theta
        result%angular_zeta = problem%n_zeta
        result%adiabatic_index = problem%adiabatic_index
        result%density_kg_m3 = problem%density_kg_m3
        result%zero_floor = problem%zero_floor
        result%unknowns = problem%classes(slot)%unknowns
        result%normal_unknowns = &
            problem%classes(slot)%normal_unknowns
        result%eta_unknowns = problem%classes(slot)%eta_unknowns
        result%mu_unknowns = problem%classes(slot)%mu_unknowns
    end subroutine initialize_result

    subroutine fixed_boundary_unknown_count(problem, parity_class, unknowns, &
            info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        integer, intent(out) :: unknowns, info
        integer :: slot

        unknowns = 0
        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        slot = class_slot(problem, parity_class)
        if (slot == 0) return
        unknowns = problem%classes(slot)%unknowns
        info = fixed_boundary_ok
    end subroutine fixed_boundary_unknown_count

    subroutine solve_fixed_boundary_full_spectrum(problem, parity_class, &
            result, info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        type(fixed_boundary_full_spectrum_t), intent(out) :: result
        integer, intent(out) :: info
        type(fixed_boundary_spectrum_result_t) :: certified
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        integer :: slot

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        slot = class_slot(problem, parity_class)
        if (slot == 0) return
        call variable_block_to_dense(problem%classes(slot)%stiffness, &
            stiffness, info)
        if (info /= variable_block_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_assembly_error, info == variable_block_allocation)
            return
        end if
        call variable_block_to_dense(problem%classes(slot)%mass, &
            mass, info)
        if (info /= variable_block_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_assembly_error, info == variable_block_allocation)
            return
        end if
        call solve_symmetric_generalized_allocated(stiffness, mass, &
            result%eigenvalues, result%eigenvectors, info, equilibrate=.true.)
        if (info /= symmetric_eigensolver_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_solver_error, &
                info == symmetric_eigensolver_allocation)
            return
        end if
        call refine_dense_spectrum(problem%classes(slot)%stiffness, &
            problem%classes(slot)%mass, problem%solver_controls, &
            result%eigenvalues, result%eigenvectors, info)
        if (info == dense_spectrum_allocation) then
            info = fixed_boundary_allocation_error
            return
        else if (info /= dense_spectrum_ok) then
            info = fixed_boundary_solver_error
            return
        end if
        call certify_dense_spectrum_orthogonality( &
            problem%classes(slot)%mass, result%eigenvectors, info)
        if (info /= dense_spectrum_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_solver_error, info == dense_spectrum_allocation)
            return
        end if
        call diagnose_dense_spectrum(problem%classes(slot)%stiffness, &
            problem%classes(slot)%mass, result%eigenvalues, &
            result%eigenvectors, result%rayleigh_quotients, result%residuals, &
            result%resolutions, info)
        if (info /= dense_spectrum_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_solver_error, info == dense_spectrum_allocation)
            return
        end if
        call certify_dense_spectrum_inertia( &
            problem%classes(slot)%stiffness, &
            problem%classes(slot)%mass, result%eigenvalues, &
            result%residuals, result%resolutions, info)
        if (info /= dense_spectrum_ok) then
            info = fixed_boundary_solver_error
            return
        end if
        call unpermute_dense_vectors(result%eigenvectors, &
            problem%classes(slot)%permutation, info)
        if (info /= dense_spectrum_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_solver_error, info == dense_spectrum_allocation)
            return
        end if
        call solve_fixed_boundary_class(problem, parity_class, certified, info)
        if (info /= fixed_boundary_ok) return
        if (.not. dense_spectrum_is_certified(result%eigenvalues, &
            result%rayleigh_quotients, certified%zero_floor, &
            certified%negative_count, certified%floor_count, &
            certified%has_eigenvector, certified%lowest_eigenvalue, &
            certified%certificate)) then
            info = fixed_boundary_solver_error
            return
        end if
        info = fixed_boundary_ok
    end subroutine solve_fixed_boundary_full_spectrum

    ! The lowest eigenvalue above the floor band, or the lowest negative
    ! one. Bisection stops once the inertia bracket isolates it to
    ! isolation_relative; inverse iteration from the lower end converges to
    ! it, and two inertia probes around the iterate certify a bracket no
    ! wider than full bisection to the solver controls would give. Should the
    ! probes fail, the bracket is bisected to the controls instead.
    subroutine resolve_lowest(class_problem, summary, controls, result, info)
        type(fixed_boundary_class_problem_t), intent(in) :: class_problem
        type(variable_spectrum_summary_t), intent(in) :: summary
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        type(fixed_boundary_spectrum_result_t), intent(inout) :: result
        integer, intent(out) :: info
        real(dp), parameter :: isolation_relative = 1.0e-3_dp
        type(fixed_boundary_solver_controls_t) :: coarse
        real(dp), allocatable :: solver_vector(:)
        real(dp) :: interval, lower, margin, shift, upper
        integer :: allocation_status, base_count, index

        if (summary%negative_count == 0 .and. .not. summary%has_positive) then
            allocate (result%eigenvector(0), stat=allocation_status)
            if (allocation_status /= 0) then
                info = fixed_boundary_allocation_error
                return
            end if
            result%inertia_interval = summary%zero_floor
            info = fixed_boundary_ok
            return
        end if
        base_count = 0
        if (summary%negative_count == 0) base_count = summary%zero_count
        coarse = controls
        coarse%negative_bracket_relative = max(isolation_relative, &
            controls%negative_bracket_relative)
        call bracket(coarse, shift, interval, info)
        if (info /= fixed_boundary_ok) return
        call iterate(shift, info)
        if (info /= fixed_boundary_ok) return
        ! Half the width at which full bisection would stop.
        margin = 0.5_dp * controls%negative_bracket_relative &
            * abs(result%lowest_eigenvalue)
        if (base_count == 0) margin = margin + 0.5_dp &
            * controls%negative_bracket_floor * summary%zero_floor
        margin = max(margin, result%eigenpair_residual)
        call certify_lowest_bracket(class_problem%stiffness, &
            class_problem%mass, base_count, result%lowest_eigenvalue, margin, &
            lower, upper, info)
        if (info == fixed_boundary_bracket_ok) then
            result%inertia_interval = min(interval, upper - lower)
        else
            call bracket(controls, shift, result%inertia_interval, info)
            if (info /= fixed_boundary_ok) return
            call iterate(shift, info)
            if (info /= fixed_boundary_ok) return
        end if
        allocate (result%eigenvector(size(solver_vector)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = fixed_boundary_allocation_error
            return
        end if
        ! Report the eigenvector in assembly order.
        do index = 1, size(solver_vector)
            result%eigenvector(class_problem%permutation(index)) = &
                solver_vector(index)
        end do
        result%has_eigenvector = .true.
        info = fixed_boundary_ok
    contains
        subroutine bracket(stopping, lower_end, width, status)
            type(fixed_boundary_solver_controls_t), intent(in) :: stopping
            real(dp), intent(out) :: lower_end, width
            integer, intent(out) :: status

            if (summary%negative_count == 0) then
                call bracket_lowest_positive(class_problem%stiffness, &
                    class_problem%mass, summary%first_positive_lower, &
                    summary%first_positive_upper, lower_end, width, status, &
                    stopping)
            else
                call bracket_lowest_negative(class_problem%stiffness, &
                    class_problem%mass, summary%zero_floor, lower_end, &
                    width, status, stopping)
            end if
            status = merge(fixed_boundary_ok, fixed_boundary_solver_error, &
                status == fixed_boundary_bracket_ok)
        end subroutine bracket

        ! Inverse iteration from the lower end of an isolating bracket, then
        ! the rigorous M^-1-norm residual bound rather than the iteration's
        ! convergence metric.
        subroutine iterate(lower_end, status)
            real(dp), intent(in) :: lower_end
            integer, intent(out) :: status

            call iterate_variable_generalized_eigenvalue( &
                class_problem%stiffness, class_problem%mass, lower_end, &
                result%lowest_eigenvalue, solver_vector, &
                result%eigenpair_residual, result%eigenpair_resolution, &
                status, controls, validated=.true.)
            if (status /= variable_generalized_ok) then
                status = fixed_boundary_solver_error
                return
            end if
            call variable_eigenvalue_bound(class_problem%stiffness, &
                class_problem%mass, solver_vector, result%lowest_eigenvalue, &
                result%eigenpair_residual, status)
            status = merge(fixed_boundary_ok, fixed_boundary_solver_error, &
                status == variable_generalized_ok)
        end subroutine iterate
    end subroutine resolve_lowest

end module fixed_boundary_spectrum
