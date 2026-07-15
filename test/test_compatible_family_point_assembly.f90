program test_compatible_family_point_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_family_point_assembly, only: &
        assemble_compatible_direct_surface, &
        assemble_compatible_transformed_surface, &
        compatible_two_component_term_count
    implicit none

    integer, parameter :: field_periods = 5
    integer, parameter :: n_theta = 8, n_zeta = 9
    integer, parameter :: trials = 15, h1_count = 2, l2_count = 2
    integer, parameter :: matrix_size = trials * (h1_count + l2_count)
    integer, parameter :: mode_m(trials) = [3, 4, 2, 4, 2, 5, 1, 5, 1, &
        6, 0, 6, 0, 3, 3]
    integer, parameter :: mode_n(trials) = [2, 2, 2, 7, -3, 2, 2, 12, &
        -8, 12, 8, 2, 2, 7, -3]
    real(dp) :: fields(n_theta, n_zeta, 13), drive(n_theta, n_zeta)
    real(dp) :: h1(h1_count, trials), dh1(h1_count, trials)
    real(dp) :: l2(l2_count, trials)
    integer :: parity(trials), trial

    call build_inputs(fields, drive, h1, dh1, l2)
    parity = 1
    call verify_equivalence(fields, drive, h1, dh1, l2, parity, "even")
    parity = 2
    call verify_equivalence(fields, drive, h1, dh1, l2, parity, "odd")
    parity = [(modulo(trial, 2) + 1, trial=1, trials)]
    call verify_equivalence(fields, drive, h1, dh1, l2, parity, "mixed")
    call verify_rejections(fields, drive, h1, dh1, l2)
    write (*, "(a)") "PASS"

contains

    subroutine build_inputs(surface_fields, surface_drive, h1_values, &
            h1_derivatives, l2_values)
        real(dp), intent(out) :: surface_fields(:, :, :), surface_drive(:, :)
        real(dp), intent(out) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(out) :: l2_values(:, :)
        real(dp) :: theta, zeta
        integer :: j, k, trial

        do k = 1, size(surface_fields, 2)
            zeta = real(k - 1, dp) / real(size(surface_fields, 2), dp)
            do j = 1, size(surface_fields, 1)
                theta = real(j - 1, dp) / real(size(surface_fields, 1), dp)
                surface_fields(j, k, 1) = 1.2_dp + 0.07_dp * cospi(2.0_dp * theta)
                surface_fields(j, k, 2) = 0.37_dp + 0.03_dp * sinpi(2.0_dp * zeta)
                surface_fields(j, k, 3) = -0.11_dp + 0.02_dp * theta
                surface_fields(j, k, 4) = 0.09_dp - 0.01_dp * zeta
                surface_fields(j, k, 5) = 0.43_dp + 0.05_dp * theta
                surface_fields(j, k, 6) = -0.28_dp + 0.04_dp * zeta
                surface_fields(j, k, 7) = 1.7_dp + 0.08_dp * cospi(2.0_dp * (theta - zeta))
                surface_fields(j, k, 8) = 2.3_dp + 0.06_dp * sinpi(2.0_dp * theta)
                surface_fields(j, k, 9) = 0.81_dp + 0.03_dp * cospi(2.0_dp * zeta)
                surface_fields(j, k, 10) = -0.17_dp + 0.02_dp * theta * zeta
                surface_fields(j, k, 11) = -0.26_dp + 0.01_dp * theta
                surface_fields(j, k, 12) = 0.14_dp - 0.02_dp * zeta
                surface_fields(j, k, 13) = -0.08_dp + 0.01_dp * theta
                surface_drive(j, k) = 0.21_dp + 0.04_dp * sinpi(2.0_dp * (theta + zeta))
            end do
        end do
        do trial = 1, trials
            h1_values(1, trial) = 0.31_dp + 0.013_dp * real(trial, dp)
            h1_values(2, trial) = -0.17_dp + 0.009_dp * real(trial, dp)
            h1_derivatives(1, trial) = -0.41_dp + 0.007_dp * real(trial, dp)
            h1_derivatives(2, trial) = 0.23_dp - 0.011_dp * real(trial, dp)
            l2_values(1, trial) = 0.19_dp - 0.005_dp * real(trial, dp)
            l2_values(2, trial) = -0.27_dp + 0.008_dp * real(trial, dp)
        end do
    end subroutine build_inputs

    subroutine verify_equivalence(surface_fields, surface_drive, h1_values, &
            h1_derivatives, l2_values, trial_parity, label)
        real(dp), intent(in) :: surface_fields(:, :, :), surface_drive(:, :)
        real(dp), intent(in) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(in) :: l2_values(:, :)
        integer, intent(in) :: trial_parity(:)
        character(len=*), intent(in) :: label
        real(dp) :: direct(matrix_size, matrix_size)
        real(dp) :: transformed(matrix_size, matrix_size)
        real(dp) :: direct_terms(matrix_size, matrix_size, &
            compatible_two_component_term_count)
        real(dp) :: transformed_terms(matrix_size, matrix_size, &
            compatible_two_component_term_count)
        real(dp) :: scale, term_scale
        integer :: info, term

        direct = 0.0_dp
        transformed = 0.0_dp
        direct_terms = 0.0_dp
        transformed_terms = 0.0_dp
        call assemble_compatible_direct_surface(surface_fields, surface_drive, &
            mode_m, mode_n, trial_parity, field_periods, h1_values, &
            h1_derivatives, l2_values, direct, info, direct_terms)
        call require(info == 0, label // " direct assembly failed")
        call assemble_compatible_transformed_surface(surface_fields, &
            surface_drive, mode_m, mode_n, trial_parity, field_periods, &
            h1_values, h1_derivatives, l2_values, transformed, info, &
            transformed_terms)
        call require(info == 0, label // " transformed assembly failed")
        scale = max(1.0_dp, maxval(abs(direct)))
        call require(maxval(abs(transformed - direct)) / scale < 8.0e-14_dp, &
            label // " direct and transformed matrices differ")
        do term = 1, compatible_two_component_term_count
            term_scale = max(1.0_dp, maxval(abs(direct_terms(:, :, term))))
            call require(maxval(abs(transformed_terms(:, :, term) &
                - direct_terms(:, :, term))) / term_scale < 8.0e-14_dp, &
                label // " direct and transformed energy terms differ")
        end do
        call require(maxval(abs(direct - sum(direct_terms, dim=3))) &
            / scale < 2.0e-15_dp, label // " direct term sum differs")
        call require(maxval(abs(transformed &
            - sum(transformed_terms, dim=3))) / scale < 2.0e-15_dp, &
            label // " transformed term sum differs")
        call require(maxval(abs(direct - transpose(direct))) / scale &
            < 2.0e-15_dp, label // " direct matrix is not symmetric")
        call require(maxval(abs(transformed - transpose(transformed))) &
            / scale < 2.0e-15_dp, &
            label // " transformed matrix is not symmetric")
    end subroutine verify_equivalence

    subroutine verify_rejections(surface_fields, surface_drive, h1_values, &
            h1_derivatives, l2_values)
        real(dp), intent(in) :: surface_fields(:, :, :), surface_drive(:, :)
        real(dp), intent(in) :: h1_values(:, :), h1_derivatives(:, :)
        real(dp), intent(in) :: l2_values(:, :)
        real(dp) :: matrix(matrix_size, matrix_size)
        real(dp) :: bad_fields(n_theta, n_zeta, 13)
        integer :: info, trial_parity(trials)

        matrix = 0.0_dp
        trial_parity = 1
        trial_parity(trials) = 0
        call assemble_compatible_transformed_surface(surface_fields, &
            surface_drive, mode_m, mode_n, trial_parity, field_periods, &
            h1_values, h1_derivatives, l2_values, matrix, info)
        call require(info /= 0, "invalid parity was accepted")
        trial_parity = 1
        bad_fields = surface_fields
        bad_fields(1, 1, 7) = ieee_value(0.0_dp, ieee_quiet_nan)
        call assemble_compatible_direct_surface(bad_fields, surface_drive, &
            mode_m, mode_n, trial_parity, field_periods, h1_values, &
            h1_derivatives, l2_values, matrix, info)
        call require(info /= 0, "nonfinite geometry field was accepted")
    end subroutine verify_rejections

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_compatible_family_point_assembly
