module fixed_boundary_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compressible_geometry, only: build_compressible_geometry, &
        compressible_geometry_ok
    use compressible_stiffness_family_assembly, only: &
        assemble_compressible_family_stiffness
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        dynamic_family_layout_t, dynamic_layout_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mass_density_policy, only: mass_density_ok, &
        mass_density_profile_t, validate_mass_density_profile
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_mass_family_assembly, only: assemble_physical_family_mass
    use radial_space_policy, only: radial_space_config_t
    use variable_block_tridiagonal, only: pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_generalized_ok, &
        variable_generalized_inertia
    use variable_spectrum_analysis, only: analyze_variable_spectrum, &
        variable_spectrum_ok, variable_spectrum_summary_t
    implicit none
    private

    integer, parameter, public :: fixed_boundary_ok = 0
    integer, parameter, public :: fixed_boundary_invalid = -1
    integer, parameter, public :: fixed_boundary_geometry_error = -2
    integer, parameter, public :: fixed_boundary_assembly_error = -3
    integer, parameter, public :: fixed_boundary_solver_error = -4
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
    end type fixed_boundary_class_problem_t

    type, public :: fixed_boundary_problem_t
        private
        logical :: ready = .false.
        logical :: has_chart_metric = .false.
        integer :: field_periods = 0
        integer :: radial_quadrature = 0
        real(dp) :: adiabatic_index = 0.0_dp
        real(dp) :: density_kg_m3 = 0.0_dp
        real(dp) :: zero_floor = 0.0_dp
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
        integer :: radial_quadrature = 0
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

    public :: build_fixed_boundary_problem
    public :: fixed_boundary_unknown_count
    public :: solve_fixed_boundary_class

contains

    subroutine build_fixed_boundary_problem(equilibrium, adiabatic_index, &
            density_kg_m3, zero_floor, mode_m, mode_n, radial_quadrature, &
            problem, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: adiabatic_index, density_kg_m3, zero_floor
        integer, intent(in) :: mode_m(:), mode_n(:), radial_quadrature
        type(fixed_boundary_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
        real(dp), allocatable :: jacobian_s(:, :, :), jacobian_t(:, :, :)
        real(dp), allocatable :: jacobian_z(:, :, :), gamma_p(:, :, :)
        integer :: parity_class

        info = fixed_boundary_invalid
        if (.not. valid_problem_inputs(adiabatic_index, density_kg_m3, &
            zero_floor, mode_m, mode_n, radial_quadrature)) return
        if (equilibrium%field_periods < 1) return
        call build_kernel_geometry(equilibrium, fixed_boundary_n_theta, &
            fixed_boundary_n_zeta, fields, drive, info)
        if (info /= mercier_ok) then
            info = fixed_boundary_geometry_error
            return
        end if
        call build_compressible_geometry(equilibrium, fixed_boundary_n_theta, &
            fixed_boundary_n_zeta, adiabatic_index, jacobian_s, jacobian_t, &
            jacobian_z, gamma_p, info)
        if (info /= compressible_geometry_ok) then
            info = fixed_boundary_geometry_error
            return
        end if
        do parity_class = 1, 2
            call assemble_class(equilibrium, fields, drive, jacobian_s, &
                jacobian_t, jacobian_z, gamma_p, density_kg_m3, mode_m, &
                mode_n, radial_quadrature, parity_class, &
                problem%classes(parity_class), info)
            if (info /= fixed_boundary_ok) return
        end do
        problem%has_chart_metric = equilibrium%has_chart_metric
        problem%field_periods = equilibrium%field_periods
        problem%radial_quadrature = radial_quadrature
        problem%adiabatic_index = adiabatic_index
        problem%density_kg_m3 = density_kg_m3
        problem%zero_floor = zero_floor
        allocate (problem%mode_m, source=mode_m)
        allocate (problem%mode_n, source=mode_n)
        problem%ready = .true.
        info = fixed_boundary_ok
    end subroutine build_fixed_boundary_problem

    function valid_problem_inputs(adiabatic_index, density_kg_m3, zero_floor, &
            mode_m, mode_n, radial_quadrature) result(valid)
        real(dp), intent(in) :: adiabatic_index, density_kg_m3, zero_floor
        integer, intent(in) :: mode_m(:), mode_n(:), radial_quadrature
        logical :: valid
        integer :: first, second

        valid = .false.
        if (.not. ieee_is_finite(adiabatic_index)) return
        if (.not. ieee_is_finite(density_kg_m3)) return
        if (.not. ieee_is_finite(zero_floor)) return
        if (adiabatic_index < 0.0_dp) return
        if (density_kg_m3 <= 0.0_dp .or. zero_floor <= 0.0_dp) return
        if (zero_floor > 0.125_dp * huge(zero_floor)) return
        if (size(mode_m) < 1 .or. size(mode_m) /= size(mode_n)) return
        if (radial_quadrature < 1 .or. radial_quadrature > 2) return
        do first = 1, size(mode_m)
            if (mode_m(first) < 0) return
            if (mode_m(first) == 0 .and. mode_n(first) < 0) return
            do second = 1, first - 1
                if (mode_m(first) == mode_m(second) &
                    .and. mode_n(first) == mode_n(second)) return
            end do
        end do
        valid = .true.
    end function valid_problem_inputs

    subroutine assemble_class(equilibrium, fields, drive, jacobian_s, &
            jacobian_t, jacobian_z, gamma_p, density_kg_m3, mode_m, mode_n, &
            radial_quadrature, parity_class, class_problem, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: fields(:, :, :, :), drive(:, :, :)
        real(dp), intent(in) :: jacobian_s(:, :, :), jacobian_t(:, :, :)
        real(dp), intent(in) :: jacobian_z(:, :, :), gamma_p(:, :, :)
        real(dp), intent(in) :: density_kg_m3
        integer, intent(in) :: mode_m(:), mode_n(:), radial_quadrature
        integer, intent(in) :: parity_class
        type(fixed_boundary_class_problem_t), intent(out) :: class_problem
        integer, intent(out) :: info
        type(dynamic_family_layout_t) :: layout, mass_layout
        type(mass_density_profile_t) :: density
        type(radial_space_config_t) :: radial_space
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        real(dp), allocatable :: stored_power(:)
        integer, allocatable :: trial_parity(:), widths(:)
        real(dp) :: radial_step

        info = fixed_boundary_assembly_error
        radial_space%quadrature_points = radial_quadrature
        allocate (trial_parity(size(mode_m)), source=parity_class)
        allocate (stored_power(size(mode_m)), source=0.0_dp)
        density%s = [0.0_dp, 1.0_dp]
        density%kilograms_per_cubic_metre = [density_kg_m3, density_kg_m3]
        call validate_mass_density_profile(density, info)
        if (info /= mass_density_ok) then
            info = fixed_boundary_invalid
            return
        end if
        radial_step = 1.0_dp / real(size(equilibrium%s), dp)
        call assemble_compressible_family_stiffness(fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, mode_n, &
            trial_parity, stored_power, equilibrium%field_periods, &
            radial_space, radial_step, phase_assembly_transformed, stiffness, &
            layout, info)
        if (info /= 0) then
            info = fixed_boundary_assembly_error
            return
        end if
        call assemble_physical_family_mass(fields, density, mode_m, mode_n, &
            trial_parity, stored_power, equilibrium%field_periods, &
            radial_space, radial_step, phase_assembly_transformed, mass, &
            mass_layout, info)
        if (info /= 0) then
            info = fixed_boundary_assembly_error
            return
        end if
        if (mass_layout%total_unknowns /= layout%total_unknowns) then
            info = fixed_boundary_assembly_error
            return
        end if
        call pack_class_matrices(stiffness, mass, layout, class_problem, &
            widths, info)
    end subroutine assemble_class

    subroutine pack_class_matrices(stiffness, mass, layout, class_problem, &
            widths, info)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :)
        type(dynamic_family_layout_t), intent(in) :: layout
        type(fixed_boundary_class_problem_t), intent(out) :: class_problem
        integer, allocatable, intent(out) :: widths(:)
        integer, intent(out) :: info

        call build_dynamic_block_permutation(layout, widths, &
            class_problem%permutation, info)
        if (info /= dynamic_layout_ok) then
            info = fixed_boundary_assembly_error
            return
        end if
        call pack_permuted_variable_blocks(stiffness, &
            class_problem%permutation, widths, class_problem%stiffness, info)
        if (info /= variable_block_ok) then
            info = fixed_boundary_assembly_error
            return
        end if
        call pack_permuted_variable_blocks(mass, class_problem%permutation, &
            widths, class_problem%mass, info)
        if (info /= variable_block_ok) then
            info = fixed_boundary_assembly_error
            return
        end if
        class_problem%unknowns = layout%total_unknowns
        class_problem%normal_unknowns = layout%normal_unknowns
        class_problem%eta_unknowns = layout%eta_unknowns
        class_problem%mu_unknowns = layout%mu_unknowns
        info = fixed_boundary_ok
    end subroutine pack_class_matrices

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
        call resolve_lowest(problem%classes(parity_class), summary, result, &
            info)
        if (info /= fixed_boundary_ok) return
        result%certificate = result%inertia_interval &
            + result%eigenpair_residual + result%eigenpair_resolution
    end subroutine solve_fixed_boundary_class

    subroutine initialize_result(problem, parity_class, result)
        type(fixed_boundary_problem_t), intent(in) :: problem
        integer, intent(in) :: parity_class
        type(fixed_boundary_spectrum_result_t), intent(out) :: result

        result%has_chart_metric = problem%has_chart_metric
        result%field_periods = problem%field_periods
        result%mode_count = size(problem%mode_m)
        result%parity_class = parity_class
        result%radial_quadrature = problem%radial_quadrature
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

    subroutine resolve_lowest(class_problem, summary, result, info)
        type(fixed_boundary_class_problem_t), intent(in) :: class_problem
        type(variable_spectrum_summary_t), intent(in) :: summary
        type(fixed_boundary_spectrum_result_t), intent(inout) :: result
        integer, intent(out) :: info
        real(dp), allocatable :: solver_vector(:)
        real(dp) :: shift

        if (summary%negative_count == 0) then
            if (.not. summary%has_positive) then
                allocate (result%eigenvector(0))
                result%inertia_interval = summary%zero_floor
                info = fixed_boundary_ok
                return
            end if
            shift = 0.5_dp * (summary%first_positive_lower &
                + summary%first_positive_upper)
            result%inertia_interval = summary%first_positive_upper &
                - summary%first_positive_lower
        else
            call bracket_lowest_negative(class_problem, summary%zero_floor, &
                shift, result%inertia_interval, info)
            if (info /= fixed_boundary_ok) return
        end if
        call iterate_variable_generalized_eigenvalue( &
            class_problem%stiffness, class_problem%mass, shift, &
            result%lowest_eigenvalue, solver_vector, &
            result%eigenpair_residual, result%eigenpair_resolution, info)
        if (info /= variable_generalized_ok) then
            info = fixed_boundary_solver_error
            return
        end if
        call unpermute_vector(solver_vector, class_problem%permutation, &
            result%eigenvector)
        result%has_eigenvector = .true.
        info = fixed_boundary_ok
    end subroutine resolve_lowest

    subroutine bracket_lowest_negative(class_problem, zero_floor, shift, &
            interval, info)
        type(fixed_boundary_class_problem_t), intent(in) :: class_problem
        real(dp), intent(in) :: zero_floor
        real(dp), intent(out) :: shift, interval
        integer, intent(out) :: info
        real(dp) :: lower, upper, middle
        integer :: count, iteration

        lower = -2.0_dp * zero_floor
        do iteration = 1, 200
            call variable_generalized_inertia(class_problem%stiffness, &
                class_problem%mass, lower, count, info)
            if (info /= variable_generalized_ok) then
                lower = lower * (1.0_dp + 1.0e-8_dp)
                cycle
            end if
            if (count == 0) exit
            lower = 2.0_dp * lower
        end do
        if (iteration > 200) then
            info = fixed_boundary_solver_error
            return
        end if
        upper = -zero_floor
        do iteration = 1, 200
            middle = 0.5_dp * (lower + upper)
            if (upper - lower <= 1.0e-9_dp * abs(middle) &
                + 1.0e-3_dp * zero_floor) exit
            call variable_generalized_inertia(class_problem%stiffness, &
                class_problem%mass, middle, count, info)
            if (info /= variable_generalized_ok) then
                middle = middle * (1.0_dp + 1.0e-8_dp)
                cycle
            end if
            if (count == 0) then
                lower = middle
            else
                upper = middle
            end if
        end do
        if (iteration > 200) then
            info = fixed_boundary_solver_error
            return
        end if
        shift = lower
        interval = upper - lower
        info = fixed_boundary_ok
    end subroutine bracket_lowest_negative

    subroutine unpermute_vector(permuted, permutation, original)
        real(dp), intent(in) :: permuted(:)
        integer, intent(in) :: permutation(:)
        real(dp), allocatable, intent(out) :: original(:)
        integer :: index

        allocate (original(size(permuted)))
        do index = 1, size(permuted)
            original(permutation(index)) = permuted(index)
        end do
    end subroutine unpermute_vector

end module fixed_boundary_spectrum
