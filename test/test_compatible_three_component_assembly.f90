program test_compatible_three_component_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use compressible_stiffness_assembly, only: &
        assemble_compressible_stiffness_surface_resolved
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use physical_mass_assembly, only: &
        assemble_physical_mass_surface_resolved
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    implicit none

    integer, parameter :: n_theta = 12, n_zeta = 10, trials = 3
    integer, parameter :: mode_m(trials) = [1, 2, 1]
    integer, parameter :: mode_n(trials) = [1, 4, -2]
    integer, parameter :: parity(trials) = [1, 2, 1]
    real(dp), parameter :: stored_power(trials) = [0.0_dp, 0.25_dp, 0.5_dp]
    real(dp), allocatable :: fields(:, :, :), drive(:, :), jacobian_s(:, :)
    real(dp), allocatable :: jacobian_t(:, :), jacobian_z(:, :), gamma_p(:, :)
    real(dp) :: h1(2, trials), dh1(2, trials), l2(1, trials)
    real(dp) :: old_stiffness(12, 12), new_stiffness(12, 12)
    real(dp) :: old_mass(12, 12), new_mass(12, 12)
    type(radial_space_config_t) :: radial_space
    integer :: info, trial

    call build_fixture(fields, drive, jacobian_s, jacobian_t, jacobian_z, &
        gamma_p)
    do trial = 1, trials
        call evaluate_normal_basis(radial_space, mode_m(trial), 0.375_dp, &
            0.25_dp, radial_space%evaluation_coordinate, h1(:, trial), &
            dh1(:, trial), info, stored_power(trial))
        call require(info == radial_space_ok, "radial basis evaluation failed")
    end do
    l2 = 1.0_dp
    call compare_legacy_stiffness(phase_assembly_direct)
    call compare_legacy_stiffness(phase_assembly_transformed)
    call compare_legacy_mass(phase_assembly_direct)
    call compare_legacy_mass(phase_assembly_transformed)
    call check_arbitrary_compatible_columns()
    call check_invalid_shapes()
    write (*, "(a)") "PASS"

contains

    subroutine compare_legacy_stiffness(phase_assembly)
        integer, intent(in) :: phase_assembly

        call assemble_compressible_stiffness_surface_resolved(fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, mode_n, &
            parity, stored_power, 3, radial_space, 0.375_dp, 0.25_dp, &
            phase_assembly, old_stiffness, info)
        call require(info == 0, "legacy stiffness assembly failed")
        new_stiffness = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, &
            mode_n, parity, 3, h1, dh1, l2, 0.25_dp, phase_assembly, &
            new_stiffness, info)
        call require(info == 0, "compatible stiffness assembly failed")
        call require_close(old_stiffness, new_stiffness, 2.0e-12_dp, &
            "P1/P0/P0 stiffness specialization changed")
    end subroutine compare_legacy_stiffness

    subroutine compare_legacy_mass(phase_assembly)
        integer, intent(in) :: phase_assembly

        call assemble_physical_mass_surface_resolved(fields, 2.3_dp, mode_m, &
            mode_n, parity, stored_power, 3, radial_space, 0.375_dp, &
            0.25_dp, phase_assembly, old_mass, info)
        call require(info == 0, "legacy mass assembly failed")
        new_mass = 0.0_dp
        call assemble_compatible_physical_mass_surface(fields, 2.3_dp, &
            mode_m, mode_n, parity, 3, h1, l2, 0.25_dp, phase_assembly, &
            new_mass, info)
        call require(info == 0, "compatible mass assembly failed")
        call require_close(old_mass, new_mass, 2.0e-12_dp, &
            "P1/P0/P0 physical-mass specialization changed")
    end subroutine compare_legacy_mass

    subroutine check_arbitrary_compatible_columns()
        real(dp) :: general_h1(3, trials), general_dh1(3, trials)
        real(dp) :: general_l2(2, trials)
        real(dp) :: direct(21, 21), transformed(21, 21)
        real(dp) :: direct_mass(21, 21), transformed_mass(21, 21)

        general_h1 = reshape([(0.1_dp * real(trial, dp), trial=1, 9)], &
            shape(general_h1))
        general_dh1 = reshape([(-0.15_dp * real(trial, dp), trial=1, 9)], &
            shape(general_dh1))
        general_l2 = reshape([(0.2_dp * real(trial, dp), trial=1, 6)], &
            shape(general_l2))
        direct = 0.0_dp
        transformed = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, &
            mode_n, parity, 3, general_h1, general_dh1, general_l2, 0.2_dp, &
            phase_assembly_direct, direct, info)
        call require(info == 0, "arbitrary direct stiffness failed")
        call assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, &
            mode_n, parity, 3, general_h1, general_dh1, general_l2, 0.2_dp, &
            phase_assembly_transformed, transformed, info)
        call require(info == 0, "arbitrary transformed stiffness failed")
        call require_close(direct, transformed, 2.0e-12_dp, &
            "arbitrary stiffness phase paths differ")
        call require_close(transformed, transpose(transformed), 2.0e-13_dp, &
            "arbitrary stiffness is not symmetric")
        direct_mass = 0.0_dp
        transformed_mass = 0.0_dp
        call assemble_compatible_physical_mass_surface(fields, 2.3_dp, &
            mode_m, mode_n, parity, 3, general_h1, general_l2, 0.2_dp, &
            phase_assembly_direct, direct_mass, info)
        call require(info == 0, "arbitrary direct mass failed")
        call assemble_compatible_physical_mass_surface(fields, 2.3_dp, &
            mode_m, mode_n, parity, 3, general_h1, general_l2, 0.2_dp, &
            phase_assembly_transformed, transformed_mass, info)
        call require(info == 0, "arbitrary transformed mass failed")
        call require_close(direct_mass, transformed_mass, 2.0e-12_dp, &
            "arbitrary mass phase paths differ")
        call require_close(transformed_mass, transpose(transformed_mass), &
            2.0e-13_dp, "arbitrary physical mass is not symmetric")
    end subroutine check_arbitrary_compatible_columns

    subroutine check_invalid_shapes()
        real(dp) :: invalid(11, 11)

        invalid = 0.0_dp
        call assemble_compatible_physical_mass_surface(fields, 2.3_dp, &
            mode_m, mode_n, parity, 3, h1, l2, 0.25_dp, &
            phase_assembly_transformed, invalid, info)
        call require(info /= 0, "invalid compatible mass shape was accepted")
        call assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, &
            mode_n, parity, 3, h1, dh1, l2, 0.25_dp, &
            phase_assembly_transformed, invalid, info)
        call require(info /= 0, &
            "invalid compatible stiffness shape was accepted")
    end subroutine check_invalid_shapes

    subroutine build_fixture(local_fields, local_drive, local_jacobian_s, &
            local_jacobian_t, local_jacobian_z, local_gamma_p)
        real(dp), allocatable, intent(out) :: local_fields(:, :, :)
        real(dp), allocatable, intent(out) :: local_drive(:, :)
        real(dp), allocatable, intent(out) :: local_jacobian_s(:, :)
        real(dp), allocatable, intent(out) :: local_jacobian_t(:, :)
        real(dp), allocatable, intent(out) :: local_jacobian_z(:, :)
        real(dp), allocatable, intent(out) :: local_gamma_p(:, :)
        real(dp) :: theta, zeta
        integer :: j, k

        allocate (local_fields(n_theta, n_zeta, 13), source=0.0_dp)
        allocate (local_drive(n_theta, n_zeta), &
            local_jacobian_s(n_theta, n_zeta), &
            local_jacobian_t(n_theta, n_zeta), &
            local_jacobian_z(n_theta, n_zeta), &
            local_gamma_p(n_theta, n_zeta))
        do k = 1, n_zeta
            zeta = two_pi() * real(k - 1, dp) / real(n_zeta, dp)
            do j = 1, n_theta
                theta = two_pi() * real(j - 1, dp) / real(n_theta, dp)
                local_fields(j, k, 1:2) = [1.2_dp + 0.05_dp * cos(theta), &
                    0.7_dp + 0.03_dp * sin(zeta)]
                local_fields(j, k, 3:6) = [0.04_dp, -0.03_dp, &
                    0.8_dp + 0.02_dp * cos(zeta), &
                    0.6_dp - 0.01_dp * sin(theta)]
                local_fields(j, k, 7:9) = [-1.1_dp &
                    - 0.02_dp * cos(theta) - 0.03_dp * sin(zeta), &
                    1.4_dp + 0.04_dp * sin(theta + zeta), &
                    1.3_dp + 0.03_dp * cos(zeta)]
                local_fields(j, k, 10:13) = [0.2_dp - 0.01_dp * cos(theta), &
                    -0.15_dp, 0.1_dp - 0.01_dp * sin(theta), &
                    -0.12_dp + 0.02_dp * cos(zeta)]
                local_drive(j, k) = 0.05_dp + 0.01_dp * cos(theta - zeta)
                local_jacobian_s(j, k) = -0.08_dp &
                    + 0.01_dp * sin(theta + zeta)
                local_jacobian_t(j, k) = 0.02_dp * sin(theta) * two_pi()
                local_jacobian_z(j, k) = -0.03_dp * cos(zeta) * two_pi()
                local_gamma_p(j, k) = 0.9_dp + 0.1_dp * cos(theta + zeta)
            end do
        end do
    end subroutine build_fixture

    pure function two_pi() result(value)
        real(dp) :: value
        value = 2.0_dp * acos(-1.0_dp)
    end function two_pi

    subroutine require_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual(:, :), expected(:, :), tolerance
        character(len=*), intent(in) :: message
        call require(maxval(abs(actual - expected)) <= tolerance &
            * max(1.0_dp, maxval(abs(actual)), maxval(abs(expected))), message)
    end subroutine require_close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_compatible_three_component_assembly
