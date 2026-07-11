program test_compressible_stiffness_family_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compressible_stiffness_assembly, only: &
        assemble_compressible_stiffness_surface_resolved
    use compressible_stiffness_family_assembly, only: &
        assemble_compressible_family_stiffness
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        dynamic_family_layout_t, eta_global_index, mu_global_index, &
        normal_global_index
    use mass_density_policy, only: mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use radial_space_policy, only: radial_space_config_t
    use physical_mass_family_assembly, only: assemble_physical_family_mass
    use symmetric_eigensolver, only: solve_symmetric_generalized
    use variable_block_tridiagonal, only: pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_to_dense, &
        variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, &
        variable_generalized_diagnostics, variable_generalized_ok
    implicit none

    integer, parameter :: intervals = 4, n_theta = 8, n_zeta = 8
    integer, parameter :: trial_m(2) = [1, 2]
    integer, parameter :: trial_n(2) = [1, 4]
    integer, parameter :: trial_parity(2) = [1, 2]
    real(dp), parameter :: stored_power(2) = 0.0_dp
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    real(dp), allocatable :: jacobian_radial(:, :, :)
    real(dp), allocatable :: jacobian_theta(:, :, :), jacobian_zeta(:, :, :)
    real(dp), allocatable :: gamma_pressure(:, :, :), direct(:, :)
    real(dp), allocatable :: transformed(:, :), zero(:, :), doubled(:, :)
    real(dp) :: probe(22), matrix_energy, element_energy
    type(dynamic_family_layout_t) :: layout, direct_layout
    type(radial_space_config_t) :: radial_space
    integer :: info

    call build_fixture(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, gamma_pressure)
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, gamma_pressure, phase_assembly_direct, direct, &
        direct_layout, info)
    call require(info == 0, "direct global stiffness assembly failed")
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, gamma_pressure, phase_assembly_transformed, &
        transformed, layout, info)
    call require(info == 0, "transformed global stiffness assembly failed")
    call require(layout%total_unknowns == 22, "global stiffness layout is wrong")
    call require(direct_layout%total_unknowns == layout%total_unknowns, &
        "phase backends produced different stiffness layouts")
    call require_close(direct, transformed, 1.0e-12_dp, &
        "global direct and transformed stiffness matrices differ")
    call require_close(transformed, transpose(transformed), 1.0e-13_dp, &
        "global compressible stiffness is not symmetric")
    call check_variable_block_structure(transformed, layout)
    call check_physical_generalized_problem(fields, transformed, layout)

    call build_probe(probe)
    matrix_energy = dot_product(probe, matmul(transformed, probe))
    call sum_element_energies(fields, drive, jacobian_radial, &
        jacobian_theta, jacobian_zeta, gamma_pressure, layout, probe, &
        element_energy, info)
    call require(info == 0, "global element-energy oracle failed")
    call require(abs(matrix_energy - element_energy) < 1.0e-12_dp &
        * max(1.0_dp, abs(matrix_energy), abs(element_energy)), &
        "global gather changes the sum of element energies")

    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, 0.0_dp * gamma_pressure, phase_assembly_transformed, &
        zero, layout, info)
    call require(info == 0, "zero-gamma global stiffness assembly failed")
    call require(maxval(abs(zero(layout%normal_unknowns &
        + layout%eta_unknowns + 1:, :))) < 1.0e-13_dp, &
        "zero-gamma global stiffness retained mu coupling")
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, 2.0_dp * gamma_pressure, phase_assembly_transformed, &
        doubled, layout, info)
    call require(info == 0, "scaled-gamma global stiffness assembly failed")
    call require_close(doubled(layout%normal_unknowns &
        + layout%eta_unknowns + 1:, :), &
        2.0_dp * transformed(layout%normal_unknowns &
        + layout%eta_unknowns + 1:, :), 1.0e-13_dp, &
        "global mu rows are not linear in gamma pressure")

    call assemble_compressible_family_stiffness(fields, drive, &
        jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure, &
        trial_m, trial_n, trial_parity, stored_power, 3, radial_space, &
        0.2_dp, phase_assembly_transformed, doubled, layout, info)
    call require(info /= 0, "inconsistent stiffness radial partition accepted")

    write (*, "(a)") "PASS"

contains

    subroutine assemble(local_fields, local_drive, local_jacobian_radial, &
            local_jacobian_theta, local_jacobian_zeta, local_gamma_pressure, &
            phase_assembly, stiffness, local_layout, local_info)
        real(dp), intent(in) :: local_fields(:, :, :, :), local_drive(:, :, :)
        real(dp), intent(in) :: local_jacobian_radial(:, :, :)
        real(dp), intent(in) :: local_jacobian_theta(:, :, :)
        real(dp), intent(in) :: local_jacobian_zeta(:, :, :)
        real(dp), intent(in) :: local_gamma_pressure(:, :, :)
        integer, intent(in) :: phase_assembly
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        type(dynamic_family_layout_t), intent(out) :: local_layout
        integer, intent(out) :: local_info

        call assemble_compressible_family_stiffness(local_fields, local_drive, &
            local_jacobian_radial, local_jacobian_theta, &
            local_jacobian_zeta, local_gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, 3, radial_space, 0.25_dp, &
            phase_assembly, stiffness, local_layout, local_info)
    end subroutine assemble

    subroutine sum_element_energies(local_fields, local_drive, &
            local_jacobian_radial, local_jacobian_theta, local_jacobian_zeta, &
            local_gamma_pressure, local_layout, global_vector, energy, info)
        real(dp), intent(in) :: local_fields(:, :, :, :), local_drive(:, :, :)
        real(dp), intent(in) :: local_jacobian_radial(:, :, :)
        real(dp), intent(in) :: local_jacobian_theta(:, :, :)
        real(dp), intent(in) :: local_jacobian_zeta(:, :, :)
        real(dp), intent(in) :: local_gamma_pressure(:, :, :), global_vector(:)
        type(dynamic_family_layout_t), intent(in) :: local_layout
        real(dp), intent(out) :: energy
        integer, intent(out) :: info
        real(dp) :: element(8, 8), local_vector(8), radial_coordinate
        integer :: interval

        energy = 0.0_dp
        do interval = 1, intervals
            radial_coordinate = (real(interval, dp) - 0.5_dp) / intervals
            call assemble_compressible_stiffness_surface_resolved( &
                local_fields(:, :, :, interval), local_drive(:, :, interval), &
                local_jacobian_radial(:, :, interval), &
                local_jacobian_theta(:, :, interval), &
                local_jacobian_zeta(:, :, interval), &
                local_gamma_pressure(:, :, interval), trial_m, trial_n, &
                trial_parity, stored_power, 3, radial_space, &
                radial_coordinate, 0.25_dp, phase_assembly_transformed, &
                element, info)
            if (info /= 0) return
            call extract_element_vector(local_layout, interval, global_vector, &
                local_vector)
            energy = energy + dot_product(local_vector, &
                matmul(element, local_vector))
        end do
        info = 0
    end subroutine sum_element_energies

    pure subroutine extract_element_vector(local_layout, interval, global, &
            local)
        type(dynamic_family_layout_t), intent(in) :: local_layout
        integer, intent(in) :: interval
        real(dp), intent(in) :: global(:)
        real(dp), intent(out) :: local(:)
        integer :: index, trial

        local = 0.0_dp
        do trial = 1, local_layout%trials
            index = normal_global_index(local_layout, interval - 1, trial)
            if (index > 0) local(trial) = global(index)
            index = normal_global_index(local_layout, interval, trial)
            if (index > 0) local(local_layout%trials + trial) = global(index)
            index = eta_global_index(local_layout, interval, trial)
            local(2 * local_layout%trials + trial) = global(index)
            index = mu_global_index(local_layout, interval, trial)
            local(3 * local_layout%trials + trial) = global(index)
        end do
    end subroutine extract_element_vector

    subroutine build_fixture(local_fields, local_drive, local_jacobian_radial, &
            local_jacobian_theta, local_jacobian_zeta, local_gamma_pressure)
        real(dp), allocatable, intent(out) :: local_fields(:, :, :, :)
        real(dp), allocatable, intent(out) :: local_drive(:, :, :)
        real(dp), allocatable, intent(out) :: local_jacobian_radial(:, :, :)
        real(dp), allocatable, intent(out) :: local_jacobian_theta(:, :, :)
        real(dp), allocatable, intent(out) :: local_jacobian_zeta(:, :, :)
        real(dp), allocatable, intent(out) :: local_gamma_pressure(:, :, :)
        real(dp) :: theta, zeta, s, two_pi
        integer :: i, j, k

        two_pi = 2.0_dp * acos(-1.0_dp)
        allocate (local_fields(n_theta, n_zeta, 13, intervals), source=0.0_dp)
        allocate (local_drive(n_theta, n_zeta, intervals))
        allocate (local_jacobian_radial(n_theta, n_zeta, intervals))
        allocate (local_jacobian_theta(n_theta, n_zeta, intervals))
        allocate (local_jacobian_zeta(n_theta, n_zeta, intervals))
        allocate (local_gamma_pressure(n_theta, n_zeta, intervals))
        do i = 1, intervals
            s = (real(i, dp) - 0.5_dp) / intervals
            do k = 1, n_zeta
                zeta = two_pi * real(k - 1, dp) / n_zeta
                do j = 1, n_theta
                    theta = two_pi * real(j - 1, dp) / n_theta
                    call set_point(local_fields(j, k, :, i), &
                        local_drive(j, k, i), local_jacobian_radial(j, k, i), &
                        local_jacobian_theta(j, k, i), &
                        local_jacobian_zeta(j, k, i), &
                        local_gamma_pressure(j, k, i), theta, zeta, s, two_pi)
                end do
            end do
        end do
    end subroutine build_fixture

    pure subroutine set_point(fields, local_drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, theta, zeta, s, &
            two_pi)
        real(dp), intent(out) :: fields(:), local_drive, jacobian_radial
        real(dp), intent(out) :: jacobian_theta, jacobian_zeta, gamma_pressure
        real(dp), intent(in) :: theta, zeta, s, two_pi

        fields = 0.0_dp
        fields(1:6) = [1.2_dp + 0.02_dp * s, 0.7_dp, 0.04_dp, -0.03_dp, &
            0.8_dp, 0.6_dp]
        fields(7) = -1.1_dp - 0.02_dp * cos(theta) - 0.03_dp * sin(zeta)
        fields(8:13) = [1.4_dp, 1.3_dp, 0.2_dp, -0.15_dp, 0.1_dp, -0.12_dp]
        local_drive = 0.05_dp + 0.01_dp * cos(theta - zeta)
        jacobian_radial = -0.08_dp + 0.01_dp * s
        jacobian_theta = 0.02_dp * two_pi * sin(theta)
        jacobian_zeta = -0.03_dp * two_pi * cos(zeta)
        gamma_pressure = 0.9_dp + 0.1_dp * cos(theta + zeta) + 0.02_dp * s
    end subroutine set_point

    pure subroutine build_probe(probe)
        real(dp), intent(out) :: probe(:)
        integer :: i

        do i = 1, size(probe)
            probe(i) = (-1.0_dp)**i * real(i, dp) / real(size(probe), dp)
        end do
    end subroutine build_probe

    subroutine check_variable_block_structure(matrix, local_layout)
        real(dp), intent(in) :: matrix(:, :)
        type(dynamic_family_layout_t), intent(in) :: local_layout
        type(variable_block_tridiagonal_t) :: blocks
        real(dp), allocatable :: blocked(:, :), expected(:, :)
        integer, allocatable :: widths(:), permutation(:)
        integer :: i, j, local_info

        call build_dynamic_block_permutation(local_layout, widths, &
            permutation, local_info)
        call require(local_info == 0, "stiffness block permutation failed")
        call pack_permuted_variable_blocks(matrix, permutation, widths, &
            blocks, local_info)
        call require(local_info == variable_block_ok, &
            "global stiffness is not variable-block tridiagonal")
        call variable_block_to_dense(blocks, blocked, local_info)
        allocate (expected(size(matrix, 1), size(matrix, 2)))
        do j = 1, size(matrix, 2)
            do i = 1, size(matrix, 1)
                expected(i, j) = matrix(permutation(i), permutation(j))
            end do
        end do
        call require_close(blocked, expected, 1.0e-14_dp, &
            "global stiffness variable-block reconstruction is wrong")
    end subroutine check_variable_block_structure

    subroutine check_physical_generalized_problem(local_fields, stiffness, &
            stiffness_layout)
        real(dp), intent(in) :: local_fields(:, :, :, :), stiffness(:, :)
        type(dynamic_family_layout_t), intent(in) :: stiffness_layout
        type(dynamic_family_layout_t) :: mass_layout
        type(mass_density_profile_t) :: density
        type(variable_block_tridiagonal_t) :: stiffness_blocks, mass_blocks
        real(dp), allocatable :: mass(:, :), blocked_stiffness(:, :)
        real(dp), allocatable :: blocked_mass(:, :), eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        integer, allocatable :: widths(:), permutation(:)
        integer :: local_info

        allocate (density%s(3), density%kilograms_per_cubic_metre(3))
        density%s = [0.0_dp, 0.5_dp, 1.0_dp]
        density%kilograms_per_cubic_metre = [2.0_dp, 3.0_dp, 5.0_dp]
        call assemble_physical_family_mass(local_fields, density, trial_m, &
            trial_n, trial_parity, stored_power, 3, radial_space, 0.25_dp, &
            phase_assembly_transformed, mass, mass_layout, local_info)
        call require(local_info == 0, "physical generalized mass failed")
        call require(mass_layout%total_unknowns &
            == stiffness_layout%total_unknowns, "physical K/M layouts differ")
        call build_dynamic_block_permutation(stiffness_layout, widths, &
            permutation, local_info)
        call pack_permuted_variable_blocks(stiffness, permutation, widths, &
            stiffness_blocks, local_info)
        call require(local_info == variable_block_ok, &
            "physical stiffness block packing failed")
        call pack_permuted_variable_blocks(mass, permutation, widths, &
            mass_blocks, local_info)
        call require(local_info == variable_block_ok, &
            "physical mass block packing failed")
        call variable_block_to_dense(stiffness_blocks, blocked_stiffness, &
            local_info)
        call require(local_info == variable_block_ok, &
            "physical stiffness block reconstruction failed")
        call variable_block_to_dense(mass_blocks, blocked_mass, local_info)
        call require(local_info == variable_block_ok, &
            "physical mass block reconstruction failed")
        call solve_symmetric_generalized(blocked_stiffness, blocked_mass, &
            eigenvalues, eigenvectors, local_info)
        call require(local_info == 0, "dense physical generalized solve failed")
        call check_physical_solution(stiffness_blocks, mass_blocks, &
            blocked_mass, eigenvalues, eigenvectors)
    end subroutine check_physical_generalized_problem

    subroutine check_physical_solution(stiffness, mass, dense_m, &
            eigenvalues, eigenvectors)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: dense_m(:, :)
        real(dp), intent(in) :: eigenvalues(:), eigenvectors(:, :)
        real(dp), allocatable :: vector(:)
        real(dp) :: eigenvalue, oracle_quotient, oracle_residual
        real(dp) :: oracle_resolution, overlap, residual, resolution
        real(dp) :: shift, scale, tolerance
        integer :: local_info

        scale = max(1.0_dp, abs(eigenvalues(1)))
        shift = eigenvalues(1) - 0.25_dp &
            * max(eigenvalues(2) - eigenvalues(1), 1.0e-6_dp * scale)
        call iterate_variable_generalized_eigenvalue(stiffness, mass, shift, &
            eigenvalue, vector, residual, resolution, &
            local_info)
        call require(local_info == variable_generalized_ok, &
            "physical variable generalized solve failed")
        call variable_generalized_diagnostics(stiffness, mass, &
            eigenvectors(:, 1), eigenvalues(1), oracle_quotient, &
            oracle_residual, oracle_resolution, local_info)
        call require(local_info == variable_generalized_ok, &
            "dense physical mode diagnostics failed")
        tolerance = 1.0e-12_dp * scale + resolution + oracle_resolution &
            + residual + oracle_residual
        overlap = dense_cluster_overlap(vector, dense_m, eigenvalues, &
            eigenvectors, tolerance)
        call require(abs(eigenvalue - oracle_quotient) <= tolerance, &
            "physical variable quotient disagrees with dense mode")
        if (abs(eigenvalue - eigenvalues(1)) > tolerance) then
            write (error_unit, "(a,5es24.16)") &
                "FAIL: physical eigenvalues and resolution bounds ", &
                eigenvalue, eigenvalues(1), resolution, oracle_resolution, &
                oracle_residual
            error stop 1
        end if
        call require(overlap > 1.0_dp - 1.0e-10_dp, &
            "physical variable mode misses the dense lowest eigenspace")
        call require(residual <= max(1.0e-12_dp * scale, resolution), &
            "physical variable residual exceeds its convergence certificate")
    end subroutine check_physical_solution

    function dense_cluster_overlap(vector, mass, eigenvalues, eigenvectors, &
            tolerance) result(overlap)
        real(dp), intent(in) :: vector(:), mass(:, :), eigenvalues(:)
        real(dp), intent(in) :: eigenvectors(:, :), tolerance
        real(dp) :: overlap, projection, mass_image(size(vector))
        integer :: i

        mass_image = matmul(mass, vector)
        overlap = 0.0_dp
        do i = 1, size(eigenvalues)
            if (eigenvalues(i) - eigenvalues(1) > tolerance) exit
            projection = dot_product(eigenvectors(:, i), mass_image)
            overlap = overlap + projection * projection
        end do
        overlap = sqrt(overlap)
    end function dense_cluster_overlap

    subroutine require_close(first, second, tolerance, message)
        real(dp), intent(in) :: first(:, :), second(:, :), tolerance
        character(len=*), intent(in) :: message
        real(dp) :: difference, scale

        difference = maxval(abs(first - second))
        scale = max(1.0_dp, maxval(abs(first)), maxval(abs(second)))
        if (difference > tolerance * scale) then
            write (error_unit, "(a,2es24.16)") "FAIL: " // message // &
                "; difference and scale: ", difference, scale
            error stop 1
        end if
    end subroutine require_close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_compressible_stiffness_family_assembly
