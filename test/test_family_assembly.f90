program test_family_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use block_tridiagonal, only: block_tridiagonal_t
    use family_point_assembly, only: assemble_direct_surface, &
        assemble_transformed_surface
    use family_assembly, only: assemble_family_blocks, &
        assemble_family_stiffness, family_negative_count, &
        family_assembly_options_t, iterate_family_eigenvalue, &
        lowest_family_eigenvalue, phase_assembly_direct, &
        phase_assembly_transformed, surface_geometry_t
    use newcomb_limit, only: cylinder_profiles_t, &
        lowest_eigenvalue_single_mode
    use radial_space_policy, only: form_s_power_edge, radial_space_config_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: n_radial = 100, n_theta = 64, n_zeta = 32
    type(cylinder_profiles_t) :: profiles
    type(surface_geometry_t), allocatable :: geometry(:)
    real(dp) :: reference, family_value, pair_value, step
    real(dp) :: iterated
    integer :: i, info, negatives

    profiles%length = 6.0_dp * pi
    profiles%b_axial = 1.0_dp
    profiles%b_linear = 0.3_dp
    profiles%b_cubic = 0.4_dp
    step = 0.5_dp / real(n_radial, dp)

    allocate (geometry(n_radial))
    do i = 1, n_radial
        call cylinder_surface((real(i, dp) - 0.5_dp) * step, geometry(i))
    end do

    call lowest_eigenvalue_single_mode(profiles, 1, 1, 0.5_dp, &
        n_radial, reference)
    call lowest_family_eigenvalue(geometry, [1], [1], step, &
        family_value, info)
    call require(info == 0, "family assembly failed")
    call require(abs(family_value - reference) < 1.0e-6_dp * &
        abs(reference), &
        "single-mode family disagrees with the 1D assembly")

    call lowest_family_eigenvalue(geometry, [1, 2], [1, 1], step, &
        pair_value, info)
    call require(info == 0, "two-mode family assembly failed")
    call require(abs(pair_value - family_value) < 1.0e-10_dp * &
        abs(family_value), &
        "decoupled modes change the lowest eigenvalue")

    call check_block_assembly(geometry, step)
    call check_parity_classes(geometry, step)
    call check_symmetric_decoupling(step)
    call check_field_period_phase(geometry, step)
    call check_transformed_assembly()
    call check_radial_space_options(geometry, step)

    call iterate_family_eigenvalue(geometry, [1], [1], step, &
        1.05_dp * family_value, iterated, info)
    call require(info == 0, "inverse iteration failed")
    ! The dual-parity pair is degenerate up to quadrature round-off;
    ! the iterate may land anywhere in that split.
    call require(abs(iterated - family_value) < 1.0e-7_dp * &
        abs(family_value), &
        "inverse iteration disagrees with the dense solve")

    call family_negative_count(geometry, [1], [1], step, 0.0_dp, &
        negatives, info)
    call require(info == 0, "inertia count failed")
    call require(negatives >= 1, &
        "the unstable mode is invisible to the inertia count")
    call family_negative_count(geometry, [1], [1], step, &
        2.0_dp * family_value, negatives, info)
    call require(info == 0, "shifted inertia count failed")
    call require(negatives == 0, &
        "the inertia count below the lowest eigenvalue is not zero")
    write (*, "(a)") "PASS"

contains

    subroutine check_block_assembly(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        type(block_tridiagonal_t) :: blocks
        real(dp), allocatable :: stiffness(:, :), dense_block(:, :)
        real(dp) :: deviation
        integer :: trials, nodes, info, i

        call assemble_family_stiffness(geometry, [1], [1], step, &
            stiffness, info)
        call require(info == 0, "dense assembly failed")
        call assemble_family_blocks(geometry, [1], [1], step, blocks, &
            info)
        call require(info == 0, "block assembly failed")
        trials = size(blocks%diag, 1)
        nodes = size(blocks%diag, 3)
        deviation = 0.0_dp
        allocate (dense_block(trials, trials))
        do i = 1, nodes
            dense_block = stiffness(trials * (i - 1) + 1:trials * i, &
                trials * (i - 1) + 1:trials * i)
            deviation = max(deviation, maxval(abs(dense_block &
                - blocks%diag(:, :, i))))
            if (i < nodes) then
                dense_block = stiffness(trials * i + 1:trials * (i + 1), &
                    trials * (i - 1) + 1:trials * i)
                deviation = max(deviation, maxval(abs(dense_block &
                    - blocks%off(:, :, i))))
            end if
        end do
        call require(deviation == 0.0_dp, &
            "block assembly differs from the dense assembly")
    end subroutine check_block_assembly

    subroutine check_parity_classes(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        real(dp) :: both, first, second
        integer :: info
        type(family_assembly_options_t) :: options

        call lowest_family_eigenvalue(geometry, [1], [1], step, both, &
            info)
        call require(info == 0, "dual-class solve failed")
        options%parity_class = 1
        call lowest_family_eigenvalue(geometry, [1], [1], step, first, &
            info, options)
        call require(info == 0, "class-1 solve failed")
        options%parity_class = 2
        call lowest_family_eigenvalue(geometry, [1], [1], step, &
            second, info, options)
        call require(info == 0, "class-2 solve failed")
        call require(abs(first - second) < 1.0e-7_dp * abs(first), &
            "cylinder parity classes are not degenerate")
        call require(abs(min(first, second) - both) < 1.0e-10_dp &
            * abs(both), &
            "class split changes the lowest eigenvalue")
        call check_count_additivity(geometry, [1], [1], step, 0.0_dp)
        call check_count_additivity(geometry, [1], [1], step, &
            0.5_dp * both)
    end subroutine check_parity_classes

    subroutine check_symmetric_decoupling(step)
        real(dp), intent(in) :: step
        type(surface_geometry_t), allocatable :: shaped(:)
        real(dp), allocatable :: stiffness(:, :)
        real(dp) :: both, first, second, cross, scale
        integer :: info, i, row, column, trials
        type(family_assembly_options_t) :: options

        allocate (shaped(n_radial))
        do i = 1, n_radial
            call symmetric_surface((real(i, dp) - 0.5_dp) * step, &
                shaped(i))
        end do
        call assemble_family_stiffness(shaped, [1, 2], [1, 1], step, &
            stiffness, info)
        call require(info == 0, "symmetric assembly failed")
        trials = 4
        scale = maxval(abs(stiffness))
        cross = 0.0_dp
        do column = 1, size(stiffness, 2)
            do row = 1, size(stiffness, 1)
                if (mod(mod(row - 1, trials), 2) &
                    /= mod(mod(column - 1, trials), 2)) then
                    cross = max(cross, abs(stiffness(row, column)))
                end if
            end do
        end do
        call require(cross < 1.0e-10_dp * scale, &
            "stellarator-symmetric fields couple the parity classes")
        call lowest_family_eigenvalue(shaped, [1, 2], [1, 1], step, &
            both, info)
        call require(info == 0, "shaped dual solve failed")
        options%parity_class = 1
        call lowest_family_eigenvalue(shaped, [1, 2], [1, 1], step, &
            first, info, options)
        call require(info == 0, "shaped class-1 solve failed")
        options%parity_class = 2
        call lowest_family_eigenvalue(shaped, [1, 2], [1, 1], step, &
            second, info, options)
        call require(info == 0, "shaped class-2 solve failed")
        call require(abs(min(first, second) - both) < 1.0e-10_dp &
            * abs(both), &
            "shaped class split changes the lowest eigenvalue")
        call check_count_additivity(shaped, [1, 2], [1, 1], step, &
            0.0_dp)
        call check_count_additivity(shaped, [1, 2], [1, 1], step, &
            0.5_dp * both)
    end subroutine check_symmetric_decoupling

    subroutine check_count_additivity(geometry, mode_m, mode_n, step, &
            shift)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: step, shift
        integer :: info, count_both, count_first, count_second
        type(family_assembly_options_t) :: options

        call family_negative_count(geometry, mode_m, mode_n, step, &
            shift, count_both, info)
        call require(info == 0, "dual inertia count failed")
        options%parity_class = 1
        call family_negative_count(geometry, mode_m, mode_n, step, &
            shift, count_first, info, options)
        call require(info == 0, "class-1 inertia count failed")
        options%parity_class = 2
        call family_negative_count(geometry, mode_m, mode_n, step, &
            shift, count_second, info, options)
        call require(info == 0, "class-2 inertia count failed")
        call require(count_both == count_first + count_second, &
            "class inertia counts do not sum to the dual count")
    end subroutine check_count_additivity

    subroutine check_field_period_phase(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        type(family_assembly_options_t) :: one_period, three_periods
        real(dp), allocatable :: reference(:, :), repeated(:, :), wrong(:, :)
        real(dp) :: scale
        integer :: info

        one_period%parity_class = 1
        three_periods = one_period
        three_periods%field_periods = 3
        call assemble_family_stiffness(geometry, [1], [1], step, &
            reference, info, one_period)
        call require(info == 0, "one-period assembly failed")
        call assemble_family_stiffness(geometry, [1], [3], step, &
            repeated, info, three_periods)
        call require(info == 0, "three-period assembly failed")
        scale = max(1.0_dp, maxval(abs(reference)))
        call require(maxval(abs(repeated - reference)) < 1.0e-12_dp * scale, &
            "topological toroidal mode is not scaled by field periods")
        call assemble_family_stiffness(geometry, [1], [1], step, &
            wrong, info, three_periods)
        call require(info == 0, "fractional-wave assembly failed")
        call require(maxval(abs(wrong - reference)) > 1.0e-6_dp * scale, &
            "field-period scaling has no effect on the assembled operator")
        three_periods%field_periods = 0
        call assemble_family_stiffness(geometry, [1], [1], step, &
            wrong, info, three_periods)
        call require(info == -2, "invalid field-period count was accepted")
    end subroutine check_field_period_phase

    subroutine check_transformed_assembly()
        type(surface_geometry_t), allocatable :: shaped(:)
        integer, parameter :: phase_radial = 8
        real(dp) :: phase_step
        integer :: i

        phase_step = 0.5_dp / real(phase_radial, dp)
        allocate (shaped(phase_radial))
        do i = 1, size(shaped)
            call symmetric_surface((real(i, dp) - 0.5_dp) * phase_step, &
                shaped(i))
        end do
        call compare_phase_paths(shaped, phase_step, 1, [1, 2], [1, -1])
        call compare_phase_paths(shaped, phase_step, 3, [1, 2, 3, 4], &
            [1, -1, 4, -4])
        call compare_phase_paths(shaped, phase_step, 4, [1, 2, 3, 4], &
            [1, -1, 5, -5])
        call compare_phase_paths(shaped, phase_step, 4, [1, 2, 3, 4], &
            [2, -2, 6, -6])
        call check_uncondensed_phase_transform(shaped, phase_step)
        call check_phase_path_derivative(shaped, phase_step)
    end subroutine check_transformed_assembly

    subroutine compare_phase_paths(geometry, step, periods, mode_m, mode_n)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        integer, intent(in) :: periods, mode_m(:), mode_n(:)
        type(family_assembly_options_t) :: direct, transformed
        real(dp), allocatable :: direct_matrix(:, :), transformed_matrix(:, :)
        real(dp) :: scale
        integer :: info, direct_count, transformed_count

        direct%field_periods = periods
        direct%phase_assembly = phase_assembly_direct
        transformed = direct
        transformed%phase_assembly = phase_assembly_transformed
        call assemble_family_stiffness(geometry, mode_m, mode_n, step, &
            direct_matrix, info, direct)
        call require(info == 0, "direct phase assembly failed")
        call assemble_family_stiffness(geometry, mode_m, mode_n, step, &
            transformed_matrix, info, transformed)
        call require(info == 0, "transformed phase assembly failed")
        scale = max(1.0_dp, maxval(abs(direct_matrix)))
        call require(maxval(abs(direct_matrix - transformed_matrix)) &
            < 2.0e-11_dp * scale, &
            "one-period transform differs from the direct-period oracle")
        call require(maxval(abs(transformed_matrix &
            - transpose(transformed_matrix))) < 1.0e-12_dp * scale, &
            "one-period phase transform breaks symmetry")
        call family_negative_count(geometry, mode_m, mode_n, step, 0.0_dp, &
            direct_count, info, direct)
        call require(info == 0, "direct phase inertia failed")
        call family_negative_count(geometry, mode_m, mode_n, step, 0.0_dp, &
            transformed_count, info, transformed)
        call require(info == 0, "transformed phase inertia failed")
        call require(direct_count == transformed_count, &
            "phase transform changes matrix inertia")
    end subroutine compare_phase_paths

    subroutine check_uncondensed_phase_transform(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        integer, parameter :: periods = 3, trials = 8
        integer, parameter :: trial_m(trials) = [1, 1, 2, 2, 3, 3, 4, 4]
        integer, parameter :: trial_n(trials) = [1, 1, -1, -1, 4, 4, -4, -4]
        integer, parameter :: parity(trials) = [1, 2, 1, 2, 1, 2, 1, 2]
        type(radial_space_config_t) :: radial_space
        real(dp) :: fields(1, 1, 13), drive(1, 1)
        real(dp) :: short_fields(1, 1, 12)
        real(dp) :: direct(3 * trials, 3 * trials)
        real(dp) :: transformed(3 * trials, 3 * trials), scale
        integer :: bad_parity(trials)
        integer :: info, surface

        surface = size(geometry) / 2
        fields(1, 1, :) = geometry(surface)%fields(2, 3, :)
        drive(1, 1) = geometry(surface)%drive(2, 3)
        direct = 0.0_dp
        transformed = 0.0_dp
        call assemble_direct_surface(fields, drive, trial_m, trial_n, parity, &
            periods, radial_space, 0.25_dp, step, direct, info)
        call require(info == 0, "uncondensed direct assembly failed")
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            parity, periods, radial_space, 0.25_dp, step, transformed, info)
        call require(info == 0, "uncondensed transformed assembly failed")
        scale = max(1.0_dp, maxval(abs(direct)))
        call require(maxval(abs(direct - transformed)) < 2.0e-12_dp * scale, &
            "uncondensed phase transform differs from the direct oracle")
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            parity, 0, radial_space, 0.25_dp, step, transformed, info)
        call require(info == -1, "zero surface period count was accepted")
        bad_parity = parity
        bad_parity(1) = 0
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            bad_parity, periods, radial_space, 0.25_dp, step, transformed, &
            info)
        call require(info == -1, "invalid surface parity was accepted")
        short_fields = fields(:, :, 1:12)
        call assemble_transformed_surface(short_fields, drive, trial_m, &
            trial_n, parity, periods, radial_space, 0.25_dp, step, &
            transformed, info)
        call require(info == -1, "short surface field array was accepted")
    end subroutine check_uncondensed_phase_transform

    subroutine check_phase_path_derivative(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        type(surface_geometry_t), allocatable :: plus(:), minus(:)
        type(family_assembly_options_t) :: direct, transformed
        real(dp), allocatable :: direct_plus(:, :), direct_minus(:, :)
        real(dp), allocatable :: transformed_plus(:, :), transformed_minus(:, :)
        real(dp), allocatable :: direct_derivative(:, :)
        real(dp), allocatable :: transformed_derivative(:, :)
        real(dp), parameter :: perturbation = 1.0e-3_dp
        real(dp) :: angle, difference, direction, scale
        integer :: info, surface, j, l, field

        plus = geometry
        minus = geometry
        do surface = 1, size(geometry)
            do l = 1, size(geometry(surface)%fields, 2)
                do j = 1, size(geometry(surface)%fields, 1)
                    angle = 2.0_dp * pi * (real(j - 1, dp) &
                        / real(size(geometry(surface)%fields, 1), dp) &
                        - real(l - 1, dp) &
                        / real(size(geometry(surface)%fields, 2), dp))
                    do field = 1, 13
                        direction = 0.05_dp * real(field, dp) / 13.0_dp &
                            * max(1.0_dp, &
                            abs(geometry(surface)%fields(j, l, field))) &
                            * (1.0_dp + 0.1_dp * cos(angle))
                        plus(surface)%fields(j, l, field) = &
                            plus(surface)%fields(j, l, field) &
                            + perturbation * direction
                        minus(surface)%fields(j, l, field) = &
                            minus(surface)%fields(j, l, field) &
                            - perturbation * direction
                    end do
                    direction = max(1.0_dp, &
                        abs(geometry(surface)%drive(j, l))) &
                        * (0.03_dp + 0.01_dp * sin(angle))
                    plus(surface)%drive(j, l) = plus(surface)%drive(j, l) &
                        + perturbation * direction
                    minus(surface)%drive(j, l) = minus(surface)%drive(j, l) &
                        - perturbation * direction
                end do
            end do
        end do
        direct%field_periods = 3
        direct%phase_assembly = phase_assembly_direct
        transformed = direct
        transformed%phase_assembly = phase_assembly_transformed
        call assemble_family_stiffness(plus, [1, 2], [1, -1], step, &
            direct_plus, info, direct)
        call require(info == 0, "direct positive perturbation failed")
        call assemble_family_stiffness(minus, [1, 2], [1, -1], step, &
            direct_minus, info, direct)
        call require(info == 0, "direct negative perturbation failed")
        call assemble_family_stiffness(plus, [1, 2], [1, -1], step, &
            transformed_plus, info, transformed)
        call require(info == 0, "transformed positive perturbation failed")
        call assemble_family_stiffness(minus, [1, 2], [1, -1], step, &
            transformed_minus, info, transformed)
        call require(info == 0, "transformed negative perturbation failed")
        direct_derivative = (direct_plus - direct_minus) &
            / (2.0_dp * perturbation)
        transformed_derivative = (transformed_plus - transformed_minus) &
            / (2.0_dp * perturbation)
        scale = max(1.0_dp, maxval(abs(direct_derivative)))
        difference = maxval(abs(direct_derivative - transformed_derivative))
        if (difference >= 2.0e-8_dp * scale) then
            write (error_unit, "(a,2es24.16)") &
                "phase derivative difference and scale ", difference, scale
        end if
        call require(difference < 2.0e-8_dp * scale, &
            "phase transform changes the continuous-input derivative")
    end subroutine check_phase_path_derivative

    subroutine check_radial_space_options(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        type(family_assembly_options_t) :: options
        real(dp), allocatable :: reference(:, :), explicit_default(:, :)
        real(dp), allocatable :: weighted(:, :)
        real(dp) :: scale
        integer :: info

        call assemble_family_stiffness(geometry, [2], [1], step, &
            reference, info)
        call require(info == 0, "implicit radial-space assembly failed")
        call assemble_family_stiffness(geometry, [2], [1], step, &
            explicit_default, info, options)
        call require(info == 0, "explicit radial-space assembly failed")
        call require(all(explicit_default == reference), &
            "explicit default changes the assembled operator")
        options%radial_space%form_policy = form_s_power_edge
        call assemble_family_stiffness(geometry, [2], [1], step, &
            weighted, info, options)
        call require(info == 0, "weighted radial-space assembly failed")
        scale = max(1.0_dp, maxval(abs(weighted)))
        call require(maxval(abs(weighted - transpose(weighted))) &
            < 1.0e-12_dp * scale, "weighted radial space breaks symmetry")
        call require(maxval(abs(weighted - reference)) > 1.0e-6_dp * scale, &
            "form-function policy does not change the assembled operator")
        options%radial_space%normal_degree = 2
        call assemble_family_stiffness(geometry, [2], [1], step, &
            weighted, info, options)
        call require(info == -2, "unsupported radial space was accepted")
        options%radial_space%normal_degree = 1
        options%phase_assembly = 0
        call assemble_family_stiffness(geometry, [2], [1], step, &
            weighted, info, options)
        call require(info == -2, "unsupported phase assembly was accepted")
    end subroutine check_radial_space_options

    subroutine symmetric_surface(radius, surface)
        real(dp), intent(in) :: radius
        type(surface_geometry_t), intent(out) :: surface
        real(dp) :: angle, even
        integer :: j, l

        call cylinder_surface(radius, surface)
        do l = 1, n_zeta
            do j = 1, n_theta
                angle = 2.0_dp * pi &
                    * ((real(j, dp) - 1.0_dp) / real(n_theta, dp) &
                    - (real(l, dp) - 1.0_dp) / real(n_zeta, dp))
                even = 1.0_dp + 0.05_dp * cos(angle)
                surface%fields(j, l, 7) = surface%fields(j, l, 7) &
                    * even
                surface%fields(j, l, 8) = surface%fields(j, l, 8) &
                    * (1.0_dp + 0.03_dp * cos(angle))
                surface%fields(j, l, 9) = surface%fields(j, l, 9) &
                    * (1.0_dp + 0.04_dp * cos(angle))
                surface%fields(j, l, 10) = surface%fields(j, l, 10) &
                    * even
                surface%fields(j, l, 12) = 0.1_dp * sin(angle)
                surface%fields(j, l, 13) = 0.05_dp &
                    * surface%fields(j, l, 8) * sin(angle)
                surface%drive(j, l) = surface%drive(j, l) * even
            end do
        end do
    end subroutine symmetric_surface

    subroutine cylinder_surface(radius, surface)
        real(dp), intent(in) :: radius
        type(surface_geometry_t), intent(out) :: surface
        real(dp) :: b_theta, b_theta_slope, fields(13), drive
        real(dp) :: length

        length = profiles%length
        b_theta = profiles%b_linear * radius + profiles%b_cubic &
            * radius**3
        b_theta_slope = profiles%b_linear + 3.0_dp * profiles%b_cubic &
            * radius**2
        fields(1) = 2.0_dp * pi * radius * profiles%b_axial
        fields(2) = length * b_theta
        fields(3) = 2.0_dp * pi * profiles%b_axial
        fields(4) = length * b_theta_slope
        fields(5) = length * profiles%b_axial
        fields(6) = 2.0_dp * pi * radius * b_theta
        fields(7) = 2.0_dp * pi * length * radius
        fields(8) = sqrt(b_theta**2 + profiles%b_axial**2)
        fields(9) = 1.0_dp
        fields(10) = (b_theta_slope + b_theta / radius) &
            * profiles%b_axial
        fields(11) = -b_theta * (b_theta_slope + b_theta / radius)
        fields(12) = 0.0_dp
        fields(13) = 0.0_dp
        drive = 2.0_dp * b_theta * (b_theta_slope + b_theta / radius) &
            / radius
        allocate (surface%fields(n_theta, n_zeta, 13))
        allocate (surface%drive(n_theta, n_zeta))
        surface%fields = spread(spread(fields, 1, n_theta), 2, n_zeta)
        surface%drive = drive
    end subroutine cylinder_surface

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_family_assembly
