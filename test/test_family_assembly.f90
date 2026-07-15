program test_family_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use block_tridiagonal, only: block_tridiagonal_t
    use compatible_family_point_assembly, only: &
        assemble_compatible_direct_surface, &
        assemble_compatible_transformed_surface
    use family_point_assembly, only: assemble_direct_surface, &
        assemble_transformed_surface
    use family_assembly, only: assemble_family_blocks, &
        assemble_family_stiffness, condensed_surface_coefficients, &
        family_negative_count, &
        family_assembly_options_t, iterate_family_eigenvalue, &
        lowest_family_eigenvalue, phase_assembly_direct, &
        phase_assembly_transformed, surface_geometry_t
    use newcomb_limit, only: cylinder_profiles_t, &
        lowest_artificial_stiffness_level
    use mode_topology, only: build_mode_family, mode_family_t
    use radial_space_policy, only: form_s_power_edge, radial_space_config_t
    use terpsichore_topology, only: convert_terpsichore_mask, &
        terpsichore_mode_mask_t, terpsichore_mode_selection_t, &
        terpsichore_topology_config_t, terpsichore_topology_ok
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: n_radial = 100, n_theta = 64, n_zeta = 32
    type(cylinder_profiles_t) :: profiles
    type(surface_geometry_t), allocatable :: geometry(:)
    real(dp) :: reference, matched_family, family_value, pair_value, step
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

    call lowest_artificial_stiffness_level(profiles, 2, 1, 0.5_dp, &
        n_radial, reference, info)
    call require(info == 0, "one-dimensional reference solve failed")
    call lowest_family_eigenvalue(geometry, [2], [1], step, matched_family, info)
    call require(info == 0, "family assembly failed")
    call require(abs(reference - matched_family) < 1.0e-6_dp * &
        abs(reference), &
        "matched-axis family disagrees with the 1D assembly")

    call lowest_family_eigenvalue(geometry, [1], [1], step, &
        family_value, info)
    call require(info == 0, "restricted unit-harmonic assembly failed")

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
    call check_compatible_radial_columns(geometry(20))
    call check_terpsichore_selection_equivalence(geometry, step)
    call check_radial_space_options(geometry, step)
    call check_surface_coefficients(geometry(20), step)

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

    subroutine check_compatible_radial_columns(surface)
        type(surface_geometry_t), intent(in) :: surface
        integer, parameter :: modes(2) = [1, 2], toroidal(2) = [1, -1]
        integer, parameter :: parity(2) = [1, 2], periods = 3
        real(dp) :: fields(1, 1, 13), drive(1, 1)
        real(dp) :: h1_values(3, 2), h1_derivatives(3, 2)
        real(dp) :: l2_values(2, 2), bad_l2(2, 1)
        real(dp) :: direct(10, 10), transformed(10, 10), scale
        integer :: status

        fields(1, 1, :) = surface%fields(2, 3, :)
        drive(1, 1) = surface%drive(2, 3)
        h1_values = reshape([0.2_dp, 0.5_dp, 0.3_dp, &
            0.6_dp, 0.1_dp, 0.3_dp], shape(h1_values))
        h1_derivatives = reshape([-1.2_dp, 0.4_dp, 0.8_dp, &
            -0.7_dp, 1.1_dp, -0.4_dp], shape(h1_derivatives))
        l2_values = reshape([0.75_dp, 0.25_dp, 0.4_dp, 0.6_dp], &
            shape(l2_values))
        direct = 0.0_dp
        transformed = 0.0_dp
        call assemble_compatible_direct_surface(fields, drive, modes, &
            toroidal, parity, periods, h1_values, h1_derivatives, l2_values, &
            direct, status)
        call require(status == 0, "compatible direct assembly failed")
        call assemble_compatible_transformed_surface(fields, drive, modes, &
            toroidal, parity, periods, h1_values, h1_derivatives, l2_values, &
            transformed, status)
        call require(status == 0, "compatible transformed assembly failed")
        scale = max(1.0_dp, maxval(abs(direct)))
        call require(maxval(abs(direct - transformed)) &
            < 4.0e-13_dp * scale, &
            "compatible direct and transformed assemblies differ")
        call require(maxval(abs(transformed - transpose(transformed))) &
            < 2.0e-14_dp * scale, &
            "compatible transformed assembly is not symmetric")
        bad_l2 = l2_values(:, 1:1)
        call assemble_compatible_transformed_surface(fields, drive, modes, &
            toroidal, parity, periods, h1_values, h1_derivatives, bad_l2, &
            transformed, status)
        call require(status == -1, &
            "mis-sized compatible tangential basis was accepted")
    end subroutine check_compatible_radial_columns

    subroutine check_surface_coefficients(surface, step)
        type(surface_geometry_t), intent(in) :: surface
        real(dp), intent(in) :: step
        type(surface_geometry_t) :: repeated(2)
        type(family_assembly_options_t) :: options
        real(dp), allocatable :: f_matrix(:, :), g_matrix(:, :)
        real(dp), allocatable :: k_matrix(:, :), stiffness(:, :)
        real(dp) :: expected
        integer :: info

        call condensed_surface_coefficients(surface, [2], [1], [1], &
            0.2_dp, step, 1, f_matrix, k_matrix, g_matrix, info)
        call require(info == 0, "surface coefficient extraction failed")
        call require(maxval(abs(f_matrix - transpose(f_matrix))) < 1.0e-12_dp, &
            "surface F matrix is not symmetric")
        call require(maxval(abs(g_matrix - transpose(g_matrix))) < 1.0e-12_dp, &
            "surface G matrix is not symmetric")
        repeated(1) = surface
        repeated(2) = surface
        options%parity_class = 1
        call assemble_family_stiffness(repeated, [2], [1], step, &
            stiffness, info, options)
        call require(info == 0, "surface coefficient reconstruction failed")
        expected = 0.5_dp * step * g_matrix(1, 1) &
            + 2.0_dp * f_matrix(1, 1) / step
        call require(abs(stiffness(1, 1) - expected) < 1.0e-12_dp &
            * max(1.0_dp, abs(expected)), &
            "surface F/K/G matrices do not reconstruct the finite element")
    end subroutine check_surface_coefficients

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
        real(dp) :: both, eigen_tolerance, first, second
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
        call verify_parity_split(geometry, [1], [1], step, &
            "cylinder", eigen_tolerance)
        if (abs(min(first, second) - both) >= eigen_tolerance) then
            write (error_unit, "(3(a,es24.16))") &
                "dual=", both, " class_1=", first, " class_2=", second
        end if
        call require(abs(min(first, second) - both) < eigen_tolerance, &
            "class split changes the lowest eigenvalue")
        call check_count_additivity(geometry, [1], [1], step, 0.0_dp)
        call check_count_additivity(geometry, [1], [1], step, &
            0.5_dp * both)
    end subroutine check_parity_classes

    subroutine check_symmetric_decoupling(step)
        real(dp), intent(in) :: step
        type(surface_geometry_t), allocatable :: shaped(:)
        real(dp), allocatable :: stiffness(:, :)
        real(dp) :: both, eigen_tolerance, first, second, cross, scale
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
        call verify_parity_split(shaped, [1, 2], [1, 1], step, &
            "shaped", eigen_tolerance)
        if (abs(min(first, second) - both) >= eigen_tolerance) then
            write (error_unit, "(3(a,es24.16))") &
                "shaped dual=", both, " class_1=", first, &
                " class_2=", second
        end if
        call require(abs(min(first, second) - both) < eigen_tolerance, &
            "shaped class split changes the lowest eigenvalue")
        call check_count_additivity(shaped, [1, 2], [1, 1], step, &
            0.0_dp)
        call check_count_additivity(shaped, [1, 2], [1, 1], step, &
            0.5_dp * both)
    end subroutine check_symmetric_decoupling

    subroutine verify_parity_split(geometry, mode_m, mode_n, step, label, &
            eigen_tolerance)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: step
        character(len=*), intent(in) :: label
        real(dp), intent(out) :: eigen_tolerance
        type(family_assembly_options_t) :: options
        real(dp), allocatable :: dual(:, :), first(:, :), second(:, :)
        real(dp) :: block_difference, cross_block, matrix_scale
        integer :: column, dual_column, dual_row, info, row

        call assemble_family_stiffness(geometry, mode_m, mode_n, step, &
            dual, info)
        call require(info == 0, trim(label) // " dual matrix assembly failed")
        options%parity_class = 1
        call assemble_family_stiffness(geometry, mode_m, mode_n, step, &
            first, info, options)
        call require(info == 0, trim(label) // " class-1 matrix assembly failed")
        options%parity_class = 2
        call assemble_family_stiffness(geometry, mode_m, mode_n, step, &
            second, info, options)
        call require(info == 0, trim(label) // " class-2 matrix assembly failed")
        call require(size(dual, 1) == 2 * size(first, 1), &
            trim(label) // " dual matrix row count differs")
        call require(all(shape(first) == shape(second)), &
            trim(label) // " split matrix shapes differ")
        block_difference = 0.0_dp
        cross_block = 0.0_dp
        do column = 1, size(first, 2)
            do row = 1, size(first, 1)
                dual_row = 2 * row - 1
                dual_column = 2 * column - 1
                block_difference = max(block_difference, &
                    abs(dual(dual_row, dual_column) - first(row, column)), &
                    abs(dual(dual_row + 1, dual_column + 1) &
                    - second(row, column)))
                cross_block = max(cross_block, &
                    abs(dual(dual_row, dual_column + 1)), &
                    abs(dual(dual_row + 1, dual_column)))
            end do
        end do
        matrix_scale = max(1.0_dp, maxval(abs(dual)))
        if (block_difference > 128.0_dp * epsilon(1.0_dp) &
            * real(size(dual, 1), dp) * matrix_scale &
            .or. cross_block > 128.0_dp * epsilon(1.0_dp) &
            * real(size(dual, 1), dp) * matrix_scale) then
            write (error_unit, "(3(a,es24.16))") &
                trim(label) // " block_difference=", block_difference, &
                " cross_block=", cross_block, " matrix_scale=", matrix_scale
        end if
        call require(block_difference <= 128.0_dp * epsilon(1.0_dp) &
            * real(size(dual, 1), dp) * matrix_scale, &
            trim(label) // " split matrices differ from dual blocks")
        call require(cross_block <= 128.0_dp * epsilon(1.0_dp) &
            * real(size(dual, 1), dp) * matrix_scale, &
            trim(label) // " dual matrix couples parity blocks")
        eigen_tolerance = 256.0_dp * epsilon(1.0_dp) &
            * real(size(dual, 1), dp) * matrix_scale / step
    end subroutine verify_parity_split

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
        real(dp) :: zero_power_matrix(3 * trials, 3 * trials)
        real(dp) :: stored_power(trials)
        integer :: bad_parity(trials)
        integer :: info, surface

        surface = size(geometry) / 2
        fields(1, 1, :) = geometry(surface)%fields(2, 3, :)
        drive(1, 1) = geometry(surface)%drive(2, 3)
        stored_power = 0.0_dp
        stored_power(1:2) = 0.25_dp
        direct = 0.0_dp
        transformed = 0.0_dp
        call assemble_direct_surface(fields, drive, trial_m, trial_n, parity, &
            periods, radial_space, 0.25_dp, step, direct, info, stored_power)
        call require(info == 0, "uncondensed direct assembly failed")
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            parity, periods, radial_space, 0.25_dp, step, transformed, info, &
            stored_power)
        call require(info == 0, "uncondensed transformed assembly failed")
        scale = max(1.0_dp, maxval(abs(direct)))
        call require(maxval(abs(direct - transformed)) < 2.0e-12_dp * scale, &
            "uncondensed phase transform differs from the direct oracle")
        zero_power_matrix = 0.0_dp
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            parity, periods, radial_space, 0.25_dp, step, zero_power_matrix, &
            info)
        call require(info == 0, "zero-power surface assembly failed")
        call require(maxval(abs(transformed - zero_power_matrix)) &
            > 1.0e-8_dp * scale, &
            "stored normal power does not change the assembled operator")
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            parity, 0, radial_space, 0.25_dp, step, transformed, info, &
            stored_power)
        call require(info == -1, "zero surface period count was accepted")
        bad_parity = parity
        bad_parity(1) = 0
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            bad_parity, periods, radial_space, 0.25_dp, step, transformed, &
            info, stored_power)
        call require(info == -1, "invalid surface parity was accepted")
        call assemble_transformed_surface(fields, drive, trial_m, trial_n, &
            parity, periods, radial_space, 0.25_dp, step, transformed, info, &
            stored_power(1:trials - 1))
        call require(info == -1, "mis-sized stored-power table was accepted")
        short_fields = fields(:, :, 1:12)
        call assemble_transformed_surface(short_fields, drive, trial_m, &
            trial_n, parity, periods, radial_space, 0.25_dp, step, &
            transformed, info, stored_power)
        call require(info == -1, "short surface field array was accepted")
    end subroutine check_uncondensed_phase_transform

    subroutine check_terpsichore_selection_equivalence(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        integer, parameter :: modes = 5, surfaces = 8
        type(mode_family_t) :: family
        type(terpsichore_mode_selection_t) :: selection
        type(surface_geometry_t), allocatable :: shaped(:)
        real(dp) :: rectangular_power(modes)
        integer :: permutation(modes), i

        call build_comparison_selection(family, selection, rectangular_power, &
            permutation)
        allocate (shaped(surfaces))
        do i = 1, surfaces
            call coupled_surface((real(i, dp) - 0.5_dp) * step, shaped(i))
        end do
        call check_selection_stiffness(shaped, step, family, selection, &
            rectangular_power, permutation)
        call check_selection_blocks(shaped, step, family, selection, &
            rectangular_power, permutation)
        call check_selection_surface(shaped(4), step, family, selection, &
            rectangular_power, permutation)
        call check_family_input_validation(geometry, step)
    end subroutine check_terpsichore_selection_equivalence

    subroutine build_comparison_selection(family, selection, power, &
            permutation)
        type(mode_family_t), intent(out) :: family
        type(terpsichore_mode_selection_t), intent(out) :: selection
        real(dp), intent(out) :: power(:)
        integer, intent(out) :: permutation(:)
        type(terpsichore_topology_config_t) :: config
        type(terpsichore_mode_mask_t) :: mask
        integer :: info, i

        call build_mode_family(3, 1, 2, 1, family, info)
        call require(info == 0, "rectangular comparison family failed")
        config%equilibrium_periods = 3
        config%field_periods_per_stability_period = 3
        config%parfac = 0.0_dp
        config%qn = 0.25_dp
        mask%poloidal_min = 0
        mask%toroidal_min = -1
        allocate (mask%selected(3, 3), source=.false.)
        mask%selected(2:3, 1) = .true.
        mask%selected(1:3, 3) = .true.
        call convert_terpsichore_mask(config, mask, selection, info)
        call require(info == terpsichore_topology_ok, &
            "ragged comparison selection failed")
        call require(size(family%poloidal) == size(power), &
            "rectangular comparison family has the wrong size")
        call require(selection%field_periods == family%field_periods, &
            "ragged and rectangular period counts differ")
        power = 0.0_dp
        where (family%poloidal == 1) power = config%qn
        do i = 1, size(power)
            permutation(i) = find_mode(selection%poloidal, &
                selection%toroidal, family%poloidal(i), family%toroidal(i))
            call require(permutation(i) > 0, &
                "ragged selection is missing a rectangular mode")
        end do
        do i = 1, size(power)
            call require(count(permutation == i) == 1, &
                "ragged comparison mode is duplicated")
        end do
    end subroutine build_comparison_selection

    subroutine check_selection_stiffness(geometry, step, family, selection, &
            power, permutation)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step, power(:)
        type(mode_family_t), intent(in) :: family
        type(terpsichore_mode_selection_t), intent(in) :: selection
        integer, intent(in) :: permutation(:)
        type(family_assembly_options_t) :: options
        real(dp), allocatable :: rectangular(:, :), ragged(:, :), zero(:, :)
        real(dp), allocatable :: explicit_zero(:, :), negative(:, :)
        real(dp) :: negative_power(size(power)), scale
        integer :: info

        options%field_periods = family%field_periods
        options%parity_class = selection%parity_class
        call assemble_family_stiffness(geometry, family%poloidal, &
            family%toroidal, step, rectangular, info, options, power)
        call require(info == 0, "rectangular comparison assembly failed")
        call assemble_family_stiffness(geometry, selection%poloidal, &
            selection%toroidal, step, ragged, info, options, &
            selection%stored_variable_power)
        call require(info == 0, "ragged comparison assembly failed")
        scale = max(1.0_dp, maxval(abs(rectangular)))
        call require(permuted_matrix_difference(rectangular, ragged, &
            permutation) < 2.0e-12_dp * scale, &
            "ragged and rectangular condensed matrices differ")
        if (cross_mode_max(rectangular, size(power)) &
            <= 1.0e-10_dp * scale) then
            write (error_unit, "(a,2es24.16)") &
                "3D cross-mode maximum and scale ", &
                cross_mode_max(rectangular, size(power)), scale
        end if
        call require(cross_mode_max(rectangular, size(power)) &
            > 1.0e-10_dp * scale, "3D comparison has no mode coupling")
        call assemble_family_stiffness(geometry, family%poloidal, &
            family%toroidal, step, zero, info, options)
        call require(info == 0, "default stored-power assembly failed")
        call assemble_family_stiffness(geometry, family%poloidal, &
            family%toroidal, step, explicit_zero, info, options, &
            0.0_dp * power)
        call require(info == 0, "zero stored-power assembly failed")
        call require(all(zero == explicit_zero), &
            "zero stored power changes the assembled matrix")
        call require(maxval(abs(rectangular - zero)) > 1.0e-8_dp * scale, &
            "nonzero stored power does not change the condensed matrix")
        negative_power = -power
        call assemble_family_stiffness(geometry, family%poloidal, &
            family%toroidal, step, negative, info, options, negative_power)
        call require(info == 0, "negative stored-power assembly failed")
        scale = max(1.0_dp, maxval(abs(negative)))
        call require(maxval(abs(negative - transpose(negative))) &
            < 2.0e-12_dp * scale, &
            "negative stored power breaks matrix symmetry")
    end subroutine check_selection_stiffness

    subroutine check_selection_blocks(geometry, step, family, selection, &
            power, permutation)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step, power(:)
        type(mode_family_t), intent(in) :: family
        type(terpsichore_mode_selection_t), intent(in) :: selection
        integer, intent(in) :: permutation(:)
        type(family_assembly_options_t) :: options
        type(block_tridiagonal_t) :: rectangular, ragged, zero
        real(dp) :: difference, scale
        integer :: info, i

        options%field_periods = family%field_periods
        options%parity_class = selection%parity_class
        call assemble_family_blocks(geometry, family%poloidal, &
            family%toroidal, step, rectangular, info, options, power)
        call require(info == 0, "rectangular block assembly failed")
        call assemble_family_blocks(geometry, selection%poloidal, &
            selection%toroidal, step, ragged, info, options, &
            selection%stored_variable_power)
        call require(info == 0, "ragged block assembly failed")
        call assemble_family_blocks(geometry, family%poloidal, &
            family%toroidal, step, zero, info, options)
        call require(info == 0, "zero-power block assembly failed")
        difference = 0.0_dp
        do i = 1, size(rectangular%diag, 3)
            difference = max(difference, permuted_matrix_difference(&
                rectangular%diag(:, :, i), ragged%diag(:, :, i), permutation))
        end do
        do i = 1, size(rectangular%off, 3)
            difference = max(difference, permuted_matrix_difference(&
                rectangular%off(:, :, i), ragged%off(:, :, i), permutation))
        end do
        scale = max(1.0_dp, maxval(abs(rectangular%diag)))
        call require(difference < 2.0e-12_dp * scale, &
            "ragged and rectangular block matrices differ")
        call require(maxval(abs(rectangular%diag - zero%diag)) &
            > 1.0e-8_dp * scale, &
            "stored power does not change the production block path")
    end subroutine check_selection_blocks

    subroutine check_selection_surface(surface, step, family, selection, &
            power, permutation)
        type(surface_geometry_t), intent(in) :: surface
        real(dp), intent(in) :: step, power(:)
        type(mode_family_t), intent(in) :: family
        type(terpsichore_mode_selection_t), intent(in) :: selection
        integer, intent(in) :: permutation(:)
        type(radial_space_config_t) :: radial_space
        real(dp) :: rectangular(3 * size(power), 3 * size(power))
        real(dp) :: ragged(3 * size(power), 3 * size(power)), scale
        integer :: parity(size(power)), info

        parity = selection%parity_class
        rectangular = 0.0_dp
        ragged = 0.0_dp
        call assemble_transformed_surface(surface%fields, surface%drive, &
            family%poloidal, family%toroidal, parity, family%field_periods, &
            radial_space, 0.25_dp, step, rectangular, info, power)
        call require(info == 0, "rectangular surface comparison failed")
        call assemble_transformed_surface(surface%fields, surface%drive, &
            selection%poloidal, selection%toroidal, parity, &
            selection%field_periods, radial_space, 0.25_dp, step, ragged, &
            info, selection%stored_variable_power)
        call require(info == 0, "ragged surface comparison failed")
        scale = max(1.0_dp, maxval(abs(rectangular)))
        call require(permuted_matrix_difference(rectangular, ragged, &
            permutation) < 2.0e-12_dp * scale, &
            "ragged and rectangular uncondensed matrices differ")
    end subroutine check_selection_surface

    subroutine check_family_input_validation(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        real(dp), allocatable :: stiffness(:, :)
        integer :: info

        call assemble_family_stiffness(geometry, [1, 2], [1], step, &
            stiffness, info)
        call require(info == -2, "short toroidal mode table was accepted")
        call assemble_family_stiffness(geometry, [1], [1, 2], step, &
            stiffness, info)
        call require(info == -2, "long toroidal mode table was accepted")
        call assemble_family_stiffness(geometry, [1], [1], step, &
            stiffness, info, normal_stored_power=[0.0_dp, 0.0_dp])
        call require(info == -2, "mis-sized stored-power table was accepted")
    end subroutine check_family_input_validation

    pure function permuted_matrix_difference(first, second, permutation) &
            result(difference)
        real(dp), intent(in) :: first(:, :), second(:, :)
        integer, intent(in) :: permutation(:)
        real(dp) :: difference
        integer :: row, column

        difference = 0.0_dp
        do column = 1, size(first, 2)
            do row = 1, size(first, 1)
                difference = max(difference, abs(first(row, column) &
                    - second(permuted_trial_index(row, size(permutation), &
                    permutation), permuted_trial_index(column, &
                    size(permutation), permutation))))
            end do
        end do
    end function permuted_matrix_difference

    pure function cross_mode_max(matrix, modes) result(cross)
        real(dp), intent(in) :: matrix(:, :)
        integer, intent(in) :: modes
        real(dp) :: cross
        integer :: row, column

        cross = 0.0_dp
        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                if (modulo(row - 1, modes) == &
                    modulo(column - 1, modes)) cycle
                cross = max(cross, abs(matrix(row, column)))
            end do
        end do
    end function cross_mode_max

    pure function permuted_trial_index(index, trials, permutation) &
            result(permuted)
        integer, intent(in) :: index, trials, permutation(:)
        integer :: permuted, local

        local = modulo(index - 1, trials) + 1
        permuted = index - local + permutation(local)
    end function permuted_trial_index

    pure function find_mode(mode_m, mode_n, target_m, target_n) result(index)
        integer, intent(in) :: mode_m(:), mode_n(:), target_m, target_n
        integer :: index, i

        index = 0
        do i = 1, size(mode_m)
            if (mode_m(i) /= target_m) cycle
            if (mode_n(i) /= target_n) cycle
            index = i
            return
        end do
    end function find_mode

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

    subroutine coupled_surface(radius, surface)
        real(dp), intent(in) :: radius
        type(surface_geometry_t), intent(out) :: surface
        real(dp) :: factor, theta
        integer :: j, l

        call symmetric_surface(radius, surface)
        do l = 1, n_zeta
            do j = 1, n_theta
                theta = 2.0_dp * pi * (real(j, dp) - 1.0_dp) &
                    / real(n_theta, dp)
                factor = 1.0_dp + 0.02_dp * cos(theta)
                surface%fields(j, l, 7) = surface%fields(j, l, 7) * factor
                surface%fields(j, l, 8) = surface%fields(j, l, 8) * factor
                surface%drive(j, l) = surface%drive(j, l) * factor
            end do
        end do
    end subroutine coupled_surface

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
