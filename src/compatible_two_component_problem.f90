module compatible_two_component_problem
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_family_point_assembly, only: &
        assemble_compatible_transformed_surface
    use export_surface_geometry, only: build_angular_grids
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
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

    type, public :: compatible_two_component_problem_t
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        integer :: degree = 0
        integer :: quadrature_points = 0
        integer :: h1_dofs = 0
        integer :: l2_dofs = 0
        integer :: normal_unknowns = 0
        integer :: eta_unknowns = 0
    end type compatible_two_component_problem_t

    public :: build_compatible_two_component_problem

    real(dp), parameter :: gauss_nodes(5) = [-0.9061798459386640_dp, &
        -0.5384693101056831_dp, 0.0_dp, 0.5384693101056831_dp, &
        0.9061798459386640_dp]
    real(dp), parameter :: gauss_weights(5) = [0.2369268850561891_dp, &
        0.4786286704993665_dp, 0.5688888888888889_dp, &
        0.4786286704993665_dp, 0.2369268850561891_dp]

contains

    subroutine build_compatible_two_component_problem(equilibrium, mode_m, &
            mode_n, stored_power, parity_class, degree, n_theta, n_zeta, &
            problem, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        type(compatible_two_component_problem_t), intent(out) :: problem
        integer, intent(out) :: info
        type(primitive_equilibrium_spline_t) :: spline
        type(radial_feec_complex_t) :: complex
        type(trial_space_topology_t) :: topology
        real(dp), allocatable :: breaks(:), theta(:), zeta(:)
        integer, allocatable :: eta_rank(:), normal_rank(:), parity(:)
        integer :: allocation_status, intervals, local_info, unknowns

        problem = compatible_two_component_problem_t()
        info = compatible_problem_invalid
        if (.not. inputs_are_valid(equilibrium, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, n_zeta)) return
        intervals = size(equilibrium%s)
        allocate (breaks(intervals + 1), parity(size(mode_m)), &
            normal_rank(size(mode_m)), eta_rank(size(mode_m)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_problem_allocation_error
            return
        end if
        breaks = [(real(local_info, dp) / real(intervals, dp), &
            local_info=0, intervals)]
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
        call fit_primitive_equilibrium(equilibrium, spline, local_info)
        if (local_info /= primitive_equilibrium_ok) then
            info = compatible_problem_assembly_error
            return
        end if
        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        call assemble_problem(spline, complex, breaks, theta, zeta, mode_m, &
            mode_n, parity, stored_power, topology, normal_rank, eta_rank, &
            problem, info)
        if (info /= compatible_problem_ok) return
        problem%stiffness = 0.5_dp * (problem%stiffness &
            + transpose(problem%stiffness))
        problem%mass = 0.5_dp * (problem%mass + transpose(problem%mass))
        problem%degree = degree
        problem%quadrature_points = size(gauss_nodes)
        problem%h1_dofs = complex%h1_dofs
        problem%l2_dofs = complex%l2_dofs
    end subroutine build_compatible_two_component_problem

    subroutine assemble_problem(spline, complex, breaks, theta, zeta, mode_m, &
            mode_n, parity, stored_power, topology, normal_rank, eta_rank, &
            problem, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:), theta(:), zeta(:)
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: normal_rank(:), eta_rank(:)
        type(compatible_two_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        real(dp) :: coordinate, half_width, midpoint, radial_weight
        integer :: cell, point

        info = compatible_problem_assembly_error
        do cell = 1, size(breaks) - 1
            midpoint = 0.5_dp * (breaks(cell) + breaks(cell + 1))
            half_width = 0.5_dp * (breaks(cell + 1) - breaks(cell))
            do point = 1, size(gauss_nodes)
                coordinate = midpoint + half_width * gauss_nodes(point)
                radial_weight = half_width * gauss_weights(point)
                call assemble_radial_point(spline, complex, coordinate, &
                    radial_weight, theta, zeta, mode_m, mode_n, parity, &
                    stored_power, topology, normal_rank, eta_rank, problem, &
                    info)
                if (info /= compatible_problem_ok) return
            end do
        end do
        info = compatible_problem_ok
    end subroutine assemble_problem

    subroutine assemble_radial_point(spline, complex, coordinate, weight, &
            theta, zeta, mode_m, mode_n, parity, stored_power, topology, &
            normal_rank, eta_rank, problem, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: coordinate, weight, theta(:), zeta(:)
        integer, intent(in) :: mode_m(:), mode_n(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        type(trial_space_topology_t), intent(in) :: topology
        integer, intent(in) :: normal_rank(:), eta_rank(:)
        type(compatible_two_component_problem_t), intent(inout) :: problem
        integer, intent(out) :: info
        real(dp), allocatable :: fields(:, :, :), drive(:, :), h1(:), dh1(:)
        real(dp), allocatable :: l2(:), local(:, :), local_dh1(:, :)
        real(dp), allocatable :: local_h1(:, :), local_l2(:, :)
        integer, allocatable :: h1_index(:), l2_index(:), map(:)
        integer :: local_info, trials

        info = compatible_problem_assembly_error
        call evaluate_radial_feec_complex(complex, coordinate, h1, dh1, l2, &
            local_info)
        if (local_info /= radial_feec_ok) return
        h1_index = pack([(local_info, local_info=1, size(h1))], &
            h1 /= 0.0_dp .or. dh1 /= 0.0_dp)
        l2_index = pack([(local_info, local_info=1, size(l2))], l2 /= 0.0_dp)
        if (size(h1_index) < 1 .or. size(l2_index) < 1) return
        trials = size(mode_m)
        allocate (local_h1(size(h1_index), trials), &
            local_dh1(size(h1_index), trials), &
            local_l2(size(l2_index), trials))
        call apply_stored_power(coordinate, stored_power, h1(h1_index), &
            dh1(h1_index), local_h1, local_dh1, local_info)
        if (local_info /= compatible_problem_ok) return
        local_l2 = spread(l2(l2_index), 2, trials)
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
        call assemble_compatible_transformed_surface(fields, drive, mode_m, &
            mode_n, parity, spline%field_periods, local_h1, local_dh1, &
            local_l2, local, local_info)
        if (local_info /= 0) return
        call build_local_map(complex, topology, normal_rank, eta_rank, &
            h1_index, l2_index, map)
        call scatter_matrix(map, weight * local, problem%stiffness)
        call add_radial_mass(map, local_h1, local_l2, weight, problem%mass)
        info = compatible_problem_ok
    end subroutine assemble_radial_point

    subroutine apply_stored_power(coordinate, stored_power, h1, dh1, values, &
            derivatives, info)
        real(dp), intent(in) :: coordinate, stored_power(:), h1(:), dh1(:)
        real(dp), intent(out) :: values(:, :), derivatives(:, :)
        integer, intent(out) :: info
        real(dp) :: scale
        integer :: trial

        info = compatible_problem_invalid
        if (coordinate <= 0.0_dp .or. size(stored_power) /= size(values, 2)) &
            return
        do trial = 1, size(stored_power)
            scale = coordinate**(-stored_power(trial))
            if (.not. ieee_is_finite(scale)) return
            values(:, trial) = scale * h1
            derivatives(:, trial) = scale &
                * (dh1 - stored_power(trial) * h1 / coordinate)
        end do
        if (.not. all(ieee_is_finite(values)) &
            .or. .not. all(ieee_is_finite(derivatives))) return
        info = compatible_problem_ok
    end subroutine apply_stored_power

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

    subroutine scatter_matrix(map, local, global)
        integer, intent(in) :: map(:)
        real(dp), intent(in) :: local(:, :)
        real(dp), intent(inout) :: global(:, :)
        integer :: a, b

        do b = 1, size(map)
            if (map(b) == 0) cycle
            do a = 1, size(map)
                if (map(a) == 0) cycle
                global(map(a), map(b)) = global(map(a), map(b)) + local(a, b)
            end do
        end do
    end subroutine scatter_matrix

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
            parity_class, degree, n_theta, n_zeta) result(valid)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
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
        valid = n_theta >= 8 .and. n_zeta >= 8
    end function inputs_are_valid

    pure function mode_table_is_unique(mode_m, mode_n) result(unique)
        integer, intent(in) :: mode_m(:), mode_n(:)
        logical :: unique
        integer :: first, second

        unique = .false.
        do first = 1, size(mode_m)
            do second = 1, first - 1
                if (mode_m(first) == mode_m(second) &
                    .and. mode_n(first) == mode_n(second)) return
            end do
        end do
        unique = .true.
    end function mode_table_is_unique

end module compatible_two_component_problem
