program test_physical_mass_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use mass_density_policy, only: mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use physical_mass_assembly, only: assemble_physical_mass_surface
    use radial_space_policy, only: radial_space_config_t
    implicit none

    integer, parameter :: n_theta = 12, n_zeta = 10
    integer, parameter :: trial_m(3) = [1, 2, 1]
    integer, parameter :: trial_n(3) = [1, 4, -2]
    integer, parameter :: trial_parity(3) = [1, 2, 1]
    real(dp), parameter :: stored_power(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), allocatable :: fields(:, :, :), direct(:, :), transformed(:, :)
    real(dp), allocatable :: scaled(:, :), step_scaled(:, :)
    type(mass_density_profile_t) :: density_profile, scaled_profile
    type(radial_space_config_t) :: radial_space
    real(dp) :: null_vector(12), probe_vector(12), energy
    integer :: i, info

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    call build_fields(fields)
    density_profile%s = [0.0_dp, 0.5_dp, 1.0_dp]
    density_profile%kilograms_per_cubic_metre = [2.0_dp, 3.0_dp, 5.0_dp]
    call assemble(fields, density_profile, radial_space, &
        phase_assembly_direct, direct, info)
    call require(info == 0, "direct physical mass assembly failed")
    call assemble(fields, density_profile, radial_space, &
        phase_assembly_transformed, transformed, info)
    call require(info == 0, "transformed physical mass assembly failed")
    call require(maxval(abs(direct - transformed)) < 1.0e-12_dp, &
        "one-period and all-period physical mass matrices differ")
    call require(maxval(abs(transformed - transpose(transformed))) &
        < 1.0e-13_dp, "physical mass element is not symmetric")

    null_vector = 0.0_dp
    null_vector(1:3) = 1.0_dp
    null_vector(4:6) = -1.0_dp
    energy = dot_product(null_vector, matmul(transformed, null_vector))
    call require(abs(energy) < 1.0e-12_dp, &
        "midpoint P1 element does not preserve its endpoint null direction")
    do i = 1, 12
        probe_vector(i) = real(i, dp)
    end do
    call require(dot_product(probe_vector, &
        matmul(transformed, probe_vector)) > 0.0_dp, &
        "physical mass element is not positive semidefinite")
    call require_positive_semidefinite(transformed)

    scaled_profile = density_profile
    scaled_profile%kilograms_per_cubic_metre = &
        2.0_dp * scaled_profile%kilograms_per_cubic_metre
    call assemble(fields, scaled_profile, radial_space, &
        phase_assembly_transformed, scaled, info)
    call require(info == 0, "scaled-density mass assembly failed")
    call require(maxval(abs(scaled - 2.0_dp * transformed)) < 1.0e-12_dp, &
        "assembled mass is not linear in density")
    call assemble_with_step(fields, density_profile, radial_space, 0.5_dp, &
        step_scaled, info)
    call require(info == 0, "scaled-step mass assembly failed")
    call require(maxval(abs(step_scaled - 2.0_dp * transformed)) &
        < 1.0e-12_dp, "assembled mass is not linear in radial step")

    fields(:, :, 8) = 0.0_dp
    call assemble(fields, density_profile, radial_space, &
        phase_assembly_transformed, scaled, info)
    call require(info /= 0, "zero magnetic field was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine assemble(fields, profile, radial_space, phase_assembly, &
            mass, info)
        real(dp), intent(in) :: fields(:, :, :)
        type(mass_density_profile_t), intent(in) :: profile
        type(radial_space_config_t), intent(in) :: radial_space
        integer, intent(in) :: phase_assembly
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: info

        call assemble_physical_mass_surface(fields, profile, trial_m, &
            trial_n, trial_parity, stored_power, 3, radial_space, 0.375_dp, &
            0.25_dp, phase_assembly, mass, info)
    end subroutine assemble

    subroutine assemble_with_step(fields, profile, radial_space, radial_step, &
            mass, info)
        real(dp), intent(in) :: fields(:, :, :), radial_step
        type(mass_density_profile_t), intent(in) :: profile
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: info

        call assemble_physical_mass_surface(fields, profile, trial_m, &
            trial_n, trial_parity, stored_power, 3, radial_space, 0.375_dp, &
            radial_step, phase_assembly_transformed, mass, info)
    end subroutine assemble_with_step

    subroutine require_positive_semidefinite(matrix)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: copy(size(matrix, 1), size(matrix, 2))
        real(dp) :: eigenvalues(size(matrix, 1)), work(8 * size(matrix, 1))
        integer :: info

        copy = matrix
        call dsyev("N", "U", size(matrix, 1), copy, size(matrix, 1), &
            eigenvalues, work, size(work), info)
        call require(info == 0, "mass semidefinite oracle failed")
        call require(eigenvalues(1) > -1.0e-12_dp, &
            "physical mass element has a negative eigenvalue")
    end subroutine require_positive_semidefinite

    subroutine build_fields(fields)
        real(dp), allocatable, intent(out) :: fields(:, :, :)
        real(dp) :: theta, zeta
        integer :: j, k

        allocate (fields(n_theta, n_zeta, 13), source=0.0_dp)
        do k = 1, n_zeta
            zeta = 2.0_dp * acos(-1.0_dp) * real(k - 1, dp) / n_zeta
            do j = 1, n_theta
                theta = 2.0_dp * acos(-1.0_dp) * real(j - 1, dp) / n_theta
                fields(j, k, 1) = 1.2_dp + 0.05_dp * cos(theta)
                fields(j, k, 2) = 0.7_dp + 0.03_dp * sin(zeta)
                fields(j, k, 5) = 0.8_dp + 0.02_dp * cos(zeta)
                fields(j, k, 6) = 0.6_dp - 0.01_dp * sin(theta)
                fields(j, k, 7) = -1.1_dp * (1.0_dp &
                    + 0.02_dp * cos(theta - zeta))
                fields(j, k, 8) = 1.4_dp + 0.04_dp * sin(theta + zeta)
                fields(j, k, 9) = 1.3_dp + 0.03_dp * cos(zeta)
                fields(j, k, 12) = 0.2_dp - 0.01_dp * sin(theta)
                fields(j, k, 13) = -0.15_dp + 0.03_dp * cos(zeta)
            end do
        end do
    end subroutine build_fields

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_physical_mass_assembly
