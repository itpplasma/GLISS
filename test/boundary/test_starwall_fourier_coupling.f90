program test_starwall_fourier_coupling
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: build_resolved_dynamic_family_layout, &
        dynamic_family_layout_t, dynamic_layout_ok, normal_global_index
    use starwall_fourier_coupling, only: add_starwall_fourier_stiffness, &
        build_starwall_fourier_map, starwall_fourier_invalid_input, &
        starwall_fourier_invalid_layout, starwall_fourier_nonsymmetric, &
        starwall_fourier_ok, starwall_fourier_underresolved
    use trial_space_topology, only: build_trial_space_topology, &
        trial_phase_cosine, trial_phase_sine, trial_space_topology_t, &
        trial_topology_ok
    implicit none

    call test_phase_normalization_and_scatter()
    call test_invalid_inputs_are_atomic()
    write (*, "(a)") "PASS"

contains

    subroutine test_phase_normalization_and_scatter()
        integer, parameter :: nu = 7, nv = 9, intervals = 3, trials = 5
        type(dynamic_family_layout_t) :: layout
        type(trial_space_topology_t) :: topology
        real(dp), allocatable :: map(:, :), nodal(:, :), global(:, :)
        real(dp) :: expected, scale
        integer :: a, b, i, info, k, node

        call make_topology(topology)
        call build_resolved_dynamic_family_layout(topology, intervals, &
            layout, info, .true.)
        call require(info == dynamic_layout_ok, "free-boundary layout failed")
        call build_starwall_fourier_map(nu, nv, topology, map, info)
        call require(info == starwall_fourier_ok, "Fourier map failed")
        call require(all(shape(map) == [nu * nv, trials]), &
            "Fourier map shape is wrong")
        call require(maxval(abs(map(:, 1) - 1.0_dp)) == 0.0_dp, &
            "constant normal mode is not one")
        node = 2 + nu
        call require(abs(map(node, 3) - cos(2.0_dp * acos(-1.0_dp) &
            * (1.0_dp / real(nu, dp) - 1.0_dp / real(nv, dp)))) &
            < 8.0_dp * epsilon(1.0_dp), "cosine phase sign is wrong")
        call require(abs(map(1 + nu, 5) + sin(2.0_dp * acos(-1.0_dp) &
            / real(nv, dp))) < 8.0_dp * epsilon(1.0_dp), &
            "sine phase sign is wrong")
        call require(abs(dot_product(map(:, 1), map(:, 1)) &
            - real(nu * nv, dp)) < 1.0e-13_dp, &
            "constant Fourier normalization is wrong")
        do a = 2, trials
            call require(abs(dot_product(map(:, a), map(:, a)) &
                - 0.5_dp * real(nu * nv, dp)) < 1.0e-12_dp, &
                "nonconstant Fourier normalization is wrong")
        end do

        allocate (nodal(nu * nv, nu * nv), source=0.0_dp)
        do k = 1, nv
            do i = 1, nu
                node = i + nu * (k - 1)
                nodal(node, node) = 2.0_dp + 0.2_dp * sin(2.0_dp * acos(-1.0_dp) &
                    * (real(i - 1, dp) / real(nu, dp) &
                    + 2.0_dp * real(k - 1, dp) / real(nv, dp)))
            end do
        end do
        nodal = nodal + 0.03_dp * matmul(reshape(map(:, 3), [nu * nv, 1]), &
            reshape(map(:, 5), [1, nu * nv]))
        nodal = 0.5_dp * (nodal + transpose(nodal))
        allocate (global(layout%total_unknowns, layout%total_unknowns), &
            source=0.0_dp)
        do a = 1, layout%total_unknowns
            global(a, a) = 0.125_dp
        end do
        call add_starwall_fourier_stiffness(nu, nv, nodal, topology, &
            layout, global, info)
        call require(info == starwall_fourier_ok, "boundary scatter failed")
        do b = 1, trials
            do a = 1, trials
                expected = dot_product(map(:, a), matmul(nodal, map(:, b)))
                if (a == b) expected = expected + 0.125_dp
                call require(abs(global(normal_global_index(layout, intervals, a), &
                    normal_global_index(layout, intervals, b)) - expected) &
                    < 2.0e-12_dp * max(1.0_dp, abs(expected)), &
                    "Fourier stiffness was scattered with the wrong phase")
            end do
        end do
        scale = max(1.0_dp, maxval(abs(global)))
        call require(maxval(abs(global - transpose(global))) &
            < 64.0_dp * epsilon(1.0_dp) * scale, &
            "scattered stiffness is not symmetric")
        call require(global(normal_global_index(layout, intervals, 3), &
            normal_global_index(layout, intervals, 5)) > 0.5_dp, &
            "helical phase-sign gate did not couple the selected modes")
    end subroutine test_phase_normalization_and_scatter

    subroutine test_invalid_inputs_are_atomic()
        integer, parameter :: nu = 7, nv = 9, intervals = 3
        type(dynamic_family_layout_t) :: fixed_layout, free_layout
        type(trial_space_topology_t) :: duplicate, topology, toroidal_two
        real(dp), allocatable :: before(:, :), global(:, :), map(:, :), nodal(:, :)
        integer :: info

        call make_topology(topology)
        call build_resolved_dynamic_family_layout(topology, intervals, &
            free_layout, info, .true.)
        call require(info == dynamic_layout_ok, "free-boundary layout failed")
        call build_resolved_dynamic_family_layout(topology, intervals, &
            fixed_layout, info)
        call require(info == dynamic_layout_ok, "fixed-boundary layout failed")
        allocate (nodal(nu * nv, nu * nv), source=0.0_dp)
        nodal = identity(size(nodal, 1))
        allocate (global(free_layout%total_unknowns, free_layout%total_unknowns), &
            source=0.25_dp)
        allocate (before, source=global)

        call build_starwall_fourier_map(4, nv, topology, map, info)
        call require(info == starwall_fourier_underresolved, &
            "underresolved poloidal grid was accepted")
        call build_trial_space_topology([1], [2], [trial_phase_cosine], &
            toroidal_two, info)
        call require(info == trial_topology_ok, "toroidal fixture failed")
        call build_starwall_fourier_map(nu, 4, toroidal_two, map, info)
        call require(info == starwall_fourier_underresolved, &
            "underresolved toroidal grid was accepted")
        call build_trial_space_topology([1, 1], [1, 1], &
            [trial_phase_cosine, trial_phase_cosine], duplicate, info)
        call require(info == trial_topology_ok, "duplicate fixture failed")
        call build_starwall_fourier_map(nu, nv, duplicate, map, info)
        call require(info == starwall_fourier_invalid_input, &
            "duplicate Fourier trial was accepted")

        call add_starwall_fourier_stiffness(nu, nv, nodal, topology, &
            fixed_layout, global, info)
        call require(info == starwall_fourier_invalid_layout, &
            "fixed-boundary layout was accepted")
        call require(all(global == before), "failed layout call changed the matrix")
        nodal(1, 2) = 1.0_dp
        call add_starwall_fourier_stiffness(nu, nv, nodal, topology, &
            free_layout, global, info)
        call require(info == starwall_fourier_nonsymmetric, &
            "nonsymmetric nodal stiffness was accepted")
        call require(all(global == before), &
            "nonsymmetric input changed the matrix")
        nodal(1, 2) = 0.0_dp
        nodal(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call add_starwall_fourier_stiffness(nu, nv, nodal, topology, &
            free_layout, global, info)
        call require(info == starwall_fourier_invalid_input, &
            "nonfinite nodal stiffness was accepted")
        call require(all(global == before), "nonfinite input changed the matrix")
        nodal(1, 1) = 1.0_dp
        global = before
        global(1, 2) = 1.0_dp
        before = global
        call add_starwall_fourier_stiffness(nu, nv, nodal, topology, &
            free_layout, global, info)
        call require(info == starwall_fourier_nonsymmetric, &
            "nonsymmetric global stiffness was accepted")
        call require(all(global == before), &
            "nonsymmetric global stiffness was changed")
        global = 0.25_dp
        global(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        before = global
        call add_starwall_fourier_stiffness(nu, nv, nodal, topology, &
            free_layout, global, info)
        call require(info == starwall_fourier_invalid_input, &
            "nonfinite global stiffness was accepted")
        call require(same_values(global, before), &
            "nonfinite global stiffness was changed")
    end subroutine test_invalid_inputs_are_atomic

    subroutine make_topology(topology)
        type(trial_space_topology_t), intent(out) :: topology
        integer :: info

        call build_trial_space_topology([0, 1, 1, 2, 1], [0, 0, 1, -1, 1], &
            [trial_phase_cosine, trial_phase_cosine, trial_phase_cosine, &
            trial_phase_sine, trial_phase_sine], topology, info)
        call require(info == trial_topology_ok, "trial topology failed")
    end subroutine make_topology

    pure function identity(n) result(matrix)
        integer, intent(in) :: n
        real(dp) :: matrix(n, n)
        integer :: i

        matrix = 0.0_dp
        do i = 1, n
            matrix(i, i) = 1.0_dp
        end do
    end function identity

    pure function same_values(first, second) result(same)
        real(dp), intent(in) :: first(:, :), second(:, :)
        logical :: same

        same = all((first == second) .or. &
            (.not. ieee_is_finite(first) .and. .not. ieee_is_finite(second)))
    end function same_values

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_starwall_fourier_coupling
