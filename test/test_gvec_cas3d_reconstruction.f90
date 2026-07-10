program test_gvec_cas3d_reconstruction
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_adapter, only: adapter_ok, reconstruct_fourier_scalar
    use gvec_cas3d_reconstruction, only: &
        periodic_sixth_order_derivatives, reconstruct_harmonic_grid, &
        reconstruction_invalid_spacing, reconstruction_invalid_surface, &
        reconstruction_nonfinite_input, reconstruction_ok, &
        reconstruction_shape_mismatch
    use gvec_cas3d_types, only: harmonic_pair_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tolerance = 2.0e-13_dp
    type(harmonic_pair_t) :: pair
    real(dp), allocatable :: values(:, :), derivative_theta(:, :)
    real(dp), allocatable :: derivative_zeta(:, :)
    integer :: info

    allocate (pair%cosine(2, 3, 3), pair%sine(2, 3, 3))
    pair%cosine = 0.0_dp
    pair%sine = 0.0_dp
    pair%cosine(1, 1, 1) = 2.0_dp
    pair%cosine(1, 2, 2) = 3.0_dp
    pair%sine(1, 2, 2) = 5.0_dp
    pair%cosine(1, 3, 3) = 4.0_dp
    pair%sine(1, 3, 3) = 6.0_dp

    call reconstruct_harmonic_grid(pair, 1, [0, 1, 2], [0, 1, -1], &
        [1.0_dp / 8.0_dp], [1.0_dp / 8.0_dp, 9.0_dp / 8.0_dp], &
        values, derivative_theta, derivative_zeta, info)
    call require(info == reconstruction_ok, "grid reconstruction failed")
    call require(maxval(abs(values - (5.0_dp + sqrt(2.0_dp)))) < tolerance, &
        "one-period Fourier value is wrong")
    call require(maxval(abs(derivative_theta - &
        (10.0_dp * pi - 20.0_dp * pi * sqrt(2.0_dp)))) < tolerance, &
        "poloidal derivative is wrong")
    call require(maxval(abs(derivative_zeta + &
        10.0_dp * pi * (1.0_dp + sqrt(2.0_dp)))) < tolerance, &
        "one-period toroidal derivative is wrong")

    call compare_scalar_oracle(pair)
    call test_invalid_inputs(pair)
    call test_periodic_derivatives()
    write (*, "(a)") "PASS"

contains

    subroutine compare_scalar_oracle(source_pair)
        type(harmonic_pair_t), intent(in) :: source_pair
        type(harmonic_pair_t) :: oracle_pair
        real(dp), parameter :: theta_grid(2) = [-0.13_dp, 0.27_dp]
        real(dp), parameter :: zeta_grid(2) = [0.02_dp, 0.61_dp]
        integer, parameter :: source_m(3) = [0, 1, 2]
        integer, parameter :: source_n(3) = [0, 1, -1]
        integer, allocatable :: paired_m(:), paired_n(:)
        real(dp), allocatable :: paired_cosine(:), paired_sine(:)
        real(dp) :: direct_value, direct_theta, direct_zeta_full
        integer :: theta_index, zeta_index

        call prepare_oracle_fixture(source_pair, source_m, source_n, &
            oracle_pair, paired_m, paired_n, paired_cosine, paired_sine)
        call reconstruct_harmonic_grid(oracle_pair, 2, source_m, source_n, &
            theta_grid, zeta_grid, values, derivative_theta, &
            derivative_zeta, info)
        call require(info == reconstruction_ok, "oracle grid reconstruction failed")
        do zeta_index = 1, 2
            do theta_index = 1, 2
                call reconstruct_fourier_scalar(theta_grid(theta_index), &
                    zeta_grid(zeta_index) / 5.0_dp, 5, paired_m, paired_n, &
                    paired_cosine, paired_sine, direct_value, direct_theta, &
                    direct_zeta_full, info)
                call require(info == adapter_ok, "scalar oracle failed")
                call require(abs(values(theta_index, zeta_index) - &
                    direct_value) < 2.0e-12_dp, "separable value changed")
                call require(abs(derivative_theta(theta_index, zeta_index) - &
                    direct_theta) < 2.0e-11_dp, &
                    "separable poloidal derivative changed")
                call require(abs(derivative_zeta(theta_index, zeta_index) - &
                    direct_zeta_full / 5.0_dp) < 2.0e-11_dp, &
                    "separable toroidal derivative changed")
            end do
        end do
    end subroutine compare_scalar_oracle

    subroutine prepare_oracle_fixture(source_pair, source_m, source_n, &
            oracle_pair, paired_m, paired_n, paired_cosine, paired_sine)
        type(harmonic_pair_t), intent(in) :: source_pair
        integer, intent(in) :: source_m(:), source_n(:)
        type(harmonic_pair_t), intent(out) :: oracle_pair
        integer, allocatable, intent(out) :: paired_m(:), paired_n(:)
        real(dp), allocatable, intent(out) :: paired_cosine(:), paired_sine(:)
        integer :: index, poloidal, toroidal

        oracle_pair = source_pair
        allocate (paired_m(9), paired_n(9), paired_cosine(9), paired_sine(9))
        index = 0
        do toroidal = 1, 3
            do poloidal = 1, 3
                oracle_pair%cosine(2, poloidal, toroidal) = &
                    0.1_dp * real(10 * poloidal + toroidal, dp)
                oracle_pair%sine(2, poloidal, toroidal) = &
                    -0.07_dp * real(poloidal + 10 * toroidal, dp)
                index = index + 1
                paired_m(index) = source_m(poloidal)
                paired_n(index) = source_n(toroidal)
                paired_cosine(index) = oracle_pair%cosine(2, poloidal, toroidal)
                paired_sine(index) = oracle_pair%sine(2, poloidal, toroidal)
            end do
        end do
    end subroutine prepare_oracle_fixture

    subroutine test_invalid_inputs(valid_pair)
        type(harmonic_pair_t), intent(in) :: valid_pair
        type(harmonic_pair_t) :: invalid_pair
        real(dp) :: nan

        call reconstruct_harmonic_grid(valid_pair, 0, [0, 1, 2], &
            [0, 1, -1], [0.0_dp], [0.0_dp], values, derivative_theta, &
            derivative_zeta, info)
        call require(info == reconstruction_invalid_surface, &
            "invalid radial surface was accepted")
        allocate (invalid_pair%cosine(1, 2, 1), invalid_pair%sine(1, 1, 1))
        call reconstruct_harmonic_grid(invalid_pair, 1, [0, 1], [0], &
            [0.0_dp], [0.0_dp], values, derivative_theta, derivative_zeta, &
            info)
        call require(info == reconstruction_shape_mismatch, &
            "mismatched harmonic arrays were accepted")
        invalid_pair = valid_pair
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        invalid_pair%cosine(1, 1, 1) = nan
        call reconstruct_harmonic_grid(invalid_pair, 1, [0, 1, 2], &
            [0, 1, -1], [0.0_dp], [0.0_dp], values, derivative_theta, &
            derivative_zeta, info)
        call require(info == reconstruction_nonfinite_input, &
            "nonfinite harmonic coefficient was accepted")
    end subroutine test_invalid_inputs

    subroutine test_periodic_derivatives()
        integer, parameter :: point_count = 65
        real(dp) :: field(point_count, point_count)
        real(dp) :: nan, theta, zeta
        real(dp), allocatable :: finite_theta(:, :), finite_zeta(:, :)
        integer :: theta_index, zeta_index

        do zeta_index = 1, point_count
            zeta = real(zeta_index - 1, dp) / real(point_count, dp)
            do theta_index = 1, point_count
                theta = real(theta_index - 1, dp) / real(point_count, dp)
                field(theta_index, zeta_index) = sin(2.0_dp * pi * theta) &
                    + 0.5_dp * cos(4.0_dp * pi * zeta)
            end do
        end do
        call periodic_sixth_order_derivatives(field, &
            1.0_dp / real(point_count, dp), 1.0_dp / real(point_count, dp), &
            finite_theta, finite_zeta, info)
        call require(info == reconstruction_ok, "periodic derivative failed")
        do zeta_index = 1, point_count
            zeta = real(zeta_index - 1, dp) / real(point_count, dp)
            do theta_index = 1, point_count
                theta = real(theta_index - 1, dp) / real(point_count, dp)
                call require(abs(finite_theta(theta_index, zeta_index) - &
                    2.0_dp * pi * cos(2.0_dp * pi * theta)) < 1.0e-5_dp, &
                    "theta derivative is wrong")
                call require(abs(finite_zeta(theta_index, zeta_index) + &
                    2.0_dp * pi * sin(4.0_dp * pi * zeta)) < 1.0e-5_dp, &
                    "zeta derivative is wrong")
            end do
        end do
        call periodic_sixth_order_derivatives(field, 0.0_dp, 1.0_dp, &
            finite_theta, finite_zeta, info)
        call require(info == reconstruction_invalid_spacing, &
            "zero derivative spacing was accepted")
        call periodic_sixth_order_derivatives(field(1:6, 1:7), 1.0_dp, &
            1.0_dp, finite_theta, finite_zeta, info)
        call require(info == reconstruction_shape_mismatch, &
            "short periodic grid was accepted")
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        field(1, 1) = nan
        call periodic_sixth_order_derivatives(field, 1.0_dp, 1.0_dp, &
            finite_theta, finite_zeta, info)
        call require(info == reconstruction_nonfinite_input, &
            "nonfinite periodic field was accepted")
    end subroutine test_periodic_derivatives

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_gvec_cas3d_reconstruction
