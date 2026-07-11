program test_physical_mass_family_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: dynamic_family_layout_t
    use mass_density_policy, only: mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use physical_mass_family_assembly, only: assemble_physical_family_mass
    use radial_space_policy, only: radial_space_config_t
    implicit none

    integer, parameter :: intervals = 4, n_theta = 8, n_zeta = 8
    integer, parameter :: trial_m(2) = [1, 2]
    integer, parameter :: trial_n(2) = [1, 4]
    integer, parameter :: trial_parity(2) = [1, 2]
    real(dp), parameter :: stored_power(2) = [0.0_dp, 0.0_dp]
    real(dp), allocatable :: fields(:, :, :, :), direct(:, :)
    real(dp), allocatable :: transformed(:, :), scaled(:, :)
    type(mass_density_profile_t) :: density_profile, scaled_profile
    type(dynamic_family_layout_t) :: layout, direct_layout
    type(radial_space_config_t) :: radial_space
    integer :: info

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *), work(*)
            real(dp), intent(out) :: w(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    call build_fields(fields)
    density_profile%s = [0.0_dp, 0.5_dp, 1.0_dp]
    density_profile%kilograms_per_cubic_metre = [2.0_dp, 3.0_dp, 5.0_dp]
    call assemble(fields, density_profile, phase_assembly_direct, direct, &
        direct_layout, info)
    call require(info == 0, "direct global mass assembly failed")
    call assemble(fields, density_profile, phase_assembly_transformed, &
        transformed, layout, info)
    call require(info == 0, "transformed global mass assembly failed")
    call require(layout%total_unknowns == 22, "global mass layout is wrong")
    call require(direct_layout%total_unknowns == layout%total_unknowns, &
        "phase backends produced different layouts")
    call require(maxval(abs(direct - transformed)) < 1.0e-12_dp, &
        "global direct and transformed mass matrices differ")
    call require(maxval(abs(transformed - transpose(transformed))) &
        < 1.0e-12_dp, "global physical mass is not symmetric")
    call require_positive_definite(transformed)
    call require(maxval(abs(transformed(1:layout%normal_unknowns, &
        layout%normal_unknowns + 1:))) > 1.0e-8_dp, &
        "global physical mass lost normal-tangential coupling")

    scaled_profile = density_profile
    scaled_profile%kilograms_per_cubic_metre = &
        2.0_dp * scaled_profile%kilograms_per_cubic_metre
    call assemble(fields, scaled_profile, phase_assembly_transformed, scaled, &
        layout, info)
    call require(info == 0, "scaled global mass assembly failed")
    call require(maxval(abs(scaled - 2.0_dp * transformed)) < 1.0e-12_dp, &
        "global mass is not linear in density")

    call assemble_physical_family_mass(fields, density_profile, trial_m, &
        trial_n, trial_parity, stored_power, 3, radial_space, 0.2_dp, &
        phase_assembly_transformed, scaled, layout, info)
    call require(info /= 0, "inconsistent radial partition was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine assemble(fields, profile, phase_assembly, mass, layout, info)
        real(dp), intent(in) :: fields(:, :, :, :)
        type(mass_density_profile_t), intent(in) :: profile
        integer, intent(in) :: phase_assembly
        real(dp), allocatable, intent(out) :: mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info

        call assemble_physical_family_mass(fields, profile, trial_m, trial_n, &
            trial_parity, stored_power, 3, radial_space, 0.25_dp, &
            phase_assembly, mass, layout, info)
    end subroutine assemble

    subroutine require_positive_definite(matrix)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: copy(size(matrix, 1), size(matrix, 2))
        real(dp) :: eigenvalues(size(matrix, 1)), work(8 * size(matrix, 1))
        integer :: info

        copy = matrix
        call dsyev("N", "U", size(matrix, 1), copy, size(matrix, 1), &
            eigenvalues, work, size(work), info)
        call require(info == 0, "global mass dense oracle failed")
        call require(eigenvalues(1) > 1.0e-10_dp, &
            "endpoint-constrained global mass is not positive definite")
    end subroutine require_positive_definite

    subroutine build_fields(fields)
        real(dp), allocatable, intent(out) :: fields(:, :, :, :)
        real(dp) :: theta, zeta, s
        integer :: i, j, k

        allocate (fields(n_theta, n_zeta, 13, intervals), source=0.0_dp)
        do i = 1, intervals
            s = (real(i, dp) - 0.5_dp) / real(intervals, dp)
            do k = 1, n_zeta
                zeta = 2.0_dp * acos(-1.0_dp) * real(k - 1, dp) / n_zeta
                do j = 1, n_theta
                    theta = 2.0_dp * acos(-1.0_dp) * real(j - 1, dp) / n_theta
                    call set_point_fields(fields(j, k, :, i), theta, zeta, s)
                end do
            end do
        end do
    end subroutine build_fields

    pure subroutine set_point_fields(fields, theta, zeta, s)
        real(dp), intent(out) :: fields(:)
        real(dp), intent(in) :: theta, zeta, s

        fields = 0.0_dp
        fields(1) = 1.2_dp + 0.05_dp * cos(theta) + 0.02_dp * s
        fields(2) = 0.7_dp + 0.03_dp * sin(zeta)
        fields(5) = 0.8_dp + 0.02_dp * cos(zeta)
        fields(6) = 0.6_dp - 0.01_dp * sin(theta)
        fields(7) = -1.1_dp * (1.0_dp + 0.02_dp * cos(theta - zeta))
        fields(8) = 1.4_dp + 0.04_dp * sin(theta + zeta) + 0.01_dp * s
        fields(9) = 1.3_dp + 0.03_dp * cos(zeta)
        fields(12) = 0.2_dp - 0.01_dp * sin(theta)
        fields(13) = -0.15_dp + 0.03_dp * cos(zeta)
    end subroutine set_point_fields

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_physical_mass_family_assembly
