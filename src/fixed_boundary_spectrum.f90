module fixed_boundary_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_three_component_problem, only: &
        build_compatible_three_component_problem, &
        compatible_three_component_allocation_error, &
        compatible_three_component_invalid, compatible_three_component_ok, &
        compatible_three_component_problem_t
    use dense_spectrum_support, only: certify_dense_spectrum_inertia, &
        certify_dense_spectrum_orthogonality, dense_spectrum_allocation, &
        dense_spectrum_is_certified, dense_spectrum_ok, &
        diagnose_dense_spectrum, refine_dense_spectrum
    use fixed_boundary_energy, only: diagnose_fixed_boundary_energy_store, &
        fixed_boundary_energy_allocation, fixed_boundary_energy_invalid, &
        fixed_boundary_energy_ok, fixed_boundary_energy_store_t, &
        fixed_boundary_energy_terms_t, pack_fixed_boundary_energy_store, &
        rayleigh_gradient_fixed_boundary_store
    use fixed_boundary_eigen_bracket, only: bracket_lowest_negative, &
        fixed_boundary_bracket_ok
    use fixed_boundary_solver_controls, only: &
        fixed_boundary_solver_controls_t, valid_fixed_boundary_solver_controls
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated, &
        symmetric_eigensolver_allocation, symmetric_eigensolver_ok
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_allocation, variable_block_ok, &
        variable_block_to_dense, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_generalized_ok
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
        real(dp) :: adiabatic_index = 0.0_dp
        real(dp) :: density_kg_m3 = 0.0_dp
        real(dp) :: zero_floor = 0.0_dp
        type(fixed_boundary_solver_controls_t) :: solver_controls
        integer, allocatable :: mode_m(:)
        integer, allocatable :: mode_n(:)
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
    public :: fixed_boundary_energy_terms_t, fixed_boundary_unknown_count
    public :: fixed_boundary_rayleigh_gradient
    public :: set_fixed_boundary_solver_controls, solve_fixed_boundary_class
    public :: solve_fixed_boundary_full_spectrum

contains

    subroutine build_fixed_boundary_problem(equilibrium, adiabatic_index, &
            density_kg_m3, zero_floor, mode_m, mode_n, degree, problem, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3, zero_floor
        integer, intent(in) :: mode_m(:), mode_n(:), degree
        type(fixed_boundary_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        real(dp), allocatable :: stored_power(:)
        integer :: allocation_status, mode, parity_class

        info = fixed_boundary_invalid
        if (.not. valid_inputs(equilibrium, adiabatic_index, density_kg_m3, &
            zero_floor, mode_m, mode_n, degree)) return
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
        do parity_class = 1, 2
            call assemble_class(equilibrium, adiabatic_index, density_kg_m3, &
                mode_m, mode_n, stored_power, parity_class, degree, &
                problem%classes(parity_class), info)
            if (info /= fixed_boundary_ok) return
        end do
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
            class_problem, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:), parity_class, degree
        real(dp), intent(in) :: stored_power(:)
        type(fixed_boundary_class_problem_t), intent(out) :: class_problem
        integer, intent(out) :: info
        type(compatible_three_component_problem_t) :: compatible
        integer :: compatible_info

        call build_compatible_three_component_problem(equilibrium, &
            adiabatic_index, density_kg_m3, mode_m, mode_n, stored_power, &
            parity_class, degree, fixed_boundary_n_theta, &
            fixed_boundary_n_zeta, compatible, compatible_info)
        if (compatible_info /= compatible_three_component_ok) then
            if (compatible_info == compatible_three_component_allocation_error) then
                info = fixed_boundary_allocation_error
            else if (compatible_info == compatible_three_component_invalid) then
                info = fixed_boundary_invalid
            else
                info = fixed_boundary_assembly_error
            end if
            return
        end if
        call pack_class_problem(compatible, class_problem, info)
    end subroutine assemble_class

    subroutine pack_class_problem(compatible, class_problem, info)
        type(compatible_three_component_problem_t), intent(in) :: compatible
        type(fixed_boundary_class_problem_t), intent(out) :: class_problem
        integer, intent(out) :: info
        integer :: allocation_status, index, local_info, width(1)

        info = fixed_boundary_assembly_error
        width(1) = size(compatible%stiffness, 1)
        call pack_variable_blocks(compatible%stiffness, width, &
            class_problem%stiffness, local_info)
        if (local_info /= variable_block_ok) then
            if (local_info == variable_block_allocation) &
                info = fixed_boundary_allocation_error
            return
        end if
        call pack_variable_blocks(compatible%mass, width, class_problem%mass, &
            local_info)
        if (local_info /= variable_block_ok) then
            if (local_info == variable_block_allocation) &
                info = fixed_boundary_allocation_error
            return
        end if
        allocate (class_problem%permutation(width(1)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = fixed_boundary_allocation_error
            return
        end if
        do index = 1, width(1)
            class_problem%permutation(index) = index
        end do
        call pack_fixed_boundary_energy_store(compatible%stiffness_terms, &
            class_problem%permutation, width, class_problem%energy, local_info)
        if (local_info /= fixed_boundary_energy_ok) then
            if (local_info == fixed_boundary_energy_allocation) &
                info = fixed_boundary_allocation_error
            return
        end if
        class_problem%unknowns = width(1)
        class_problem%normal_unknowns = compatible%normal_unknowns
        class_problem%eta_unknowns = compatible%eta_unknowns
        class_problem%mu_unknowns = compatible%mu_unknowns
        info = fixed_boundary_ok
    end subroutine pack_class_problem

    function valid_inputs(equilibrium, adiabatic_index, density_kg_m3, &
            zero_floor, mode_m, mode_n, degree) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3, zero_floor
        integer, intent(in) :: mode_m(:), mode_n(:), degree
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
        valid = .true.
    end function valid_inputs

    subroutine diagnose_fixed_boundary_energy(problem, parity_class, vector, &
            result, info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        real(dp), intent(in) :: vector(:)
        type(fixed_boundary_energy_terms_t), intent(out) :: result
        integer, intent(out) :: info

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        if (parity_class < 1 .or. parity_class > 2) return
        call diagnose_fixed_boundary_energy_store( &
            problem%classes(parity_class)%stiffness, &
            problem%classes(parity_class)%mass, &
            problem%classes(parity_class)%energy, &
            problem%classes(parity_class)%permutation, vector, result, info)
        info = map_fixed_boundary_energy_info(info)
    end subroutine diagnose_fixed_boundary_energy

    subroutine fixed_boundary_rayleigh_gradient(problem, parity_class, &
            vector, gradient, info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        real(dp), intent(in) :: vector(:)
        real(dp), allocatable, intent(out) :: gradient(:)
        integer, intent(out) :: info

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        if (parity_class < 1 .or. parity_class > 2) return
        call rayleigh_gradient_fixed_boundary_store( &
            problem%classes(parity_class)%stiffness, &
            problem%classes(parity_class)%mass, &
            problem%classes(parity_class)%permutation, vector, gradient, info)
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

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        if (parity_class < 1 .or. parity_class > 2) return
        call initialize_result(problem, parity_class, result)
        call analyze_variable_spectrum( &
            problem%classes(parity_class)%stiffness, &
            problem%classes(parity_class)%mass, problem%zero_floor, summary, &
            info)
        if (info /= variable_spectrum_ok) then
            info = fixed_boundary_solver_error
            return
        end if
        result%negative_count = summary%negative_count
        result%floor_count = summary%zero_count
        call resolve_lowest(problem%classes(parity_class), summary, &
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

    subroutine initialize_result(problem, parity_class, result)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        type(fixed_boundary_spectrum_result_t), intent(out) :: result

        result%has_chart_metric = problem%has_chart_metric
        result%field_periods = problem%field_periods
        result%mode_count = size(problem%mode_m)
        result%parity_class = parity_class
        result%degree = problem%degree
        result%angular_theta = fixed_boundary_n_theta
        result%angular_zeta = fixed_boundary_n_zeta
        result%adiabatic_index = problem%adiabatic_index
        result%density_kg_m3 = problem%density_kg_m3
        result%zero_floor = problem%zero_floor
        result%unknowns = problem%classes(parity_class)%unknowns
        result%normal_unknowns = &
            problem%classes(parity_class)%normal_unknowns
        result%eta_unknowns = problem%classes(parity_class)%eta_unknowns
        result%mu_unknowns = problem%classes(parity_class)%mu_unknowns
    end subroutine initialize_result

    subroutine fixed_boundary_unknown_count(problem, parity_class, unknowns, &
            info)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        integer, intent(out) :: unknowns, info

        unknowns = 0
        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        if (parity_class < 1 .or. parity_class > 2) return
        unknowns = problem%classes(parity_class)%unknowns
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

        info = fixed_boundary_invalid
        if (.not. problem%ready) return
        if (parity_class < 1 .or. parity_class > 2) return
        call variable_block_to_dense(problem%classes(parity_class)%stiffness, &
            stiffness, info)
        if (info /= variable_block_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_assembly_error, info == variable_block_allocation)
            return
        end if
        call variable_block_to_dense(problem%classes(parity_class)%mass, &
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
        call refine_dense_spectrum(problem%classes(parity_class)%stiffness, &
            problem%classes(parity_class)%mass, problem%solver_controls, &
            result%eigenvalues, result%eigenvectors, info)
        if (info == dense_spectrum_allocation) then
            info = fixed_boundary_allocation_error
            return
        else if (info /= dense_spectrum_ok) then
            info = fixed_boundary_solver_error
            return
        end if
        call certify_dense_spectrum_orthogonality( &
            problem%classes(parity_class)%mass, result%eigenvectors, info)
        if (info /= dense_spectrum_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_solver_error, info == dense_spectrum_allocation)
            return
        end if
        call diagnose_dense_spectrum(problem%classes(parity_class)%stiffness, &
            problem%classes(parity_class)%mass, result%eigenvalues, &
            result%eigenvectors, result%rayleigh_quotients, result%residuals, &
            result%resolutions, info)
        if (info /= dense_spectrum_ok) then
            info = merge(fixed_boundary_allocation_error, &
                fixed_boundary_solver_error, info == dense_spectrum_allocation)
            return
        end if
        call certify_dense_spectrum_inertia( &
            problem%classes(parity_class)%stiffness, &
            problem%classes(parity_class)%mass, result%eigenvalues, &
            result%residuals, result%resolutions, info)
        if (info /= dense_spectrum_ok) then
            info = fixed_boundary_solver_error
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

    subroutine resolve_lowest(class_problem, summary, controls, result, info)
        type(fixed_boundary_class_problem_t), intent(in) :: class_problem
        type(variable_spectrum_summary_t), intent(in) :: summary
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        type(fixed_boundary_spectrum_result_t), intent(inout) :: result
        integer, intent(out) :: info
        real(dp), allocatable :: solver_vector(:)
        real(dp) :: shift
        integer :: allocation_status

        if (summary%negative_count == 0) then
            if (.not. summary%has_positive) then
                allocate (result%eigenvector(0), stat=allocation_status)
                if (allocation_status /= 0) then
                    info = fixed_boundary_allocation_error
                    return
                end if
                result%inertia_interval = summary%zero_floor
                info = fixed_boundary_ok
                return
            end if
            shift = 0.5_dp * (summary%first_positive_lower &
                + summary%first_positive_upper)
            result%inertia_interval = summary%first_positive_upper &
                - summary%first_positive_lower
        else
            call bracket_lowest_negative(class_problem%stiffness, &
                class_problem%mass, summary%zero_floor, shift, &
                result%inertia_interval, info, controls)
            if (info /= fixed_boundary_bracket_ok) then
                info = fixed_boundary_solver_error
                return
            end if
        end if
        call iterate_variable_generalized_eigenvalue( &
            class_problem%stiffness, class_problem%mass, shift, &
            result%lowest_eigenvalue, solver_vector, &
            result%eigenpair_residual, result%eigenpair_resolution, info, &
            controls)
        if (info /= variable_generalized_ok) then
            info = fixed_boundary_solver_error
            return
        end if
        allocate (result%eigenvector(size(solver_vector)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = fixed_boundary_allocation_error
            return
        end if
        result%eigenvector = solver_vector
        result%has_eigenvector = .true.
        info = fixed_boundary_ok
    end subroutine resolve_lowest

end module fixed_boundary_spectrum
