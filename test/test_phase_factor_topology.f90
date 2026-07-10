program test_phase_factor_topology
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use phase_factor_topology, only: build_phase_envelope_table, &
        evaluate_phase_envelope, phase_cosine, phase_factor_invalid, &
        phase_factor_ok, phase_envelope_table_t, phase_product_average, &
        phase_sine
    implicit none

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
    integer, parameter :: mode_m(4) = [1, 2, 3, 4]
    integer, parameter :: mode_n(4) = [1, -1, 4, -4]
    integer, parameter :: self_conjugate_n(2) = [2, -2]
    type(phase_envelope_table_t) :: table
    real(dp) :: theta, zeta, cosine, sine, phase
    integer :: info, mode

    call build_phase_envelope_table(3, 1, 1, mode_m, mode_n, table, info)
    call require(info == phase_factor_ok, "valid odd-period table failed")
    call require(all(table%envelope_poloidal == [0, 3, 2, 5]), &
        "envelope poloidal indices are wrong")
    call require(all(table%envelope_toroidal == [0, 0, 1, -1]), &
        "envelope toroidal indices are wrong")
    call require(all(table%base_sign == [1, -1, 1, -1]), &
        "base-phase signs are wrong")
    theta = 0.37_dp
    zeta = 0.61_dp
    do mode = 1, size(mode_m)
        call evaluate_phase_envelope(table, mode, theta, zeta, cosine, sine, &
            info)
        call require(info == phase_factor_ok, &
            "valid phase-envelope evaluation failed")
        phase = two_pi * (real(mode_m(mode), dp) * theta &
            - real(mode_n(mode), dp) * zeta / 3.0_dp)
        call require(abs(cosine - cos(phase)) < 1.0e-14_dp, &
            "cosine phase reconstruction failed")
        call require(abs(sine - sin(phase)) < 1.0e-14_dp, &
            "sine phase reconstruction failed")
    end do

    call build_phase_envelope_table(5, 2, 0, [0, 1], [2, -2], table, info)
    call require(info == phase_factor_ok, "valid odd topology was rejected")
    call build_phase_envelope_table(4, 2, 0, [0, 0], [2, -2], table, info)
    call require(info == phase_factor_ok, &
        "self-conjugate even topology was rejected")
    do mode = 1, 2
        call evaluate_phase_envelope(table, mode, theta, zeta, cosine, sine, &
            info)
        call require(info == phase_factor_ok, &
            "self-conjugate phase-envelope evaluation failed")
        phase = -two_pi * real(self_conjugate_n(mode), dp) * zeta / 4.0_dp
        call require(abs(cosine - cos(phase)) < 1.0e-14_dp, &
            "self-conjugate cosine reconstruction failed")
        call require(abs(sine - sin(phase)) < 1.0e-14_dp, &
            "self-conjugate sine reconstruction failed")
    end do
    call build_phase_envelope_table(1, 0, 0, [0], [0], table, info)
    call require(info == phase_factor_invalid, &
        "single-period envelope topology was accepted")
    call build_phase_envelope_table(3, 1, 0, [0], [0], table, info)
    call require(info == phase_factor_invalid, &
        "mode outside the family was accepted")
    call evaluate_phase_envelope(table, 1, theta, zeta, cosine, sine, info)
    call require(info == phase_factor_invalid, &
        "invalid phase-envelope table was evaluated")
    call require(ieee_is_nan(cosine), &
        "invalid cosine phase-envelope evaluation did not return NaN")
    call require(ieee_is_nan(sine), &
        "invalid sine phase-envelope evaluation did not return NaN")
    call check_product_table()
    call require(ieee_is_nan(phase_product_average(phase_cosine, &
        phase_cosine, 0.0_dp, 0.0_dp, 1, 1, 0)), &
        "zero field-period count did not return NaN")
    call require(ieee_is_nan(phase_product_average(0, phase_cosine, &
        0.0_dp, 0.0_dp, 1, 1, 3)), &
        "invalid phase kind did not return NaN")
    write (*, "(a)") "PASS"

contains

    subroutine check_product_table()
        real(dp) :: first_phase, second_phase, expected, actual
        integer :: periods, first_n, second_n, first_kind, second_kind, p

        first_phase = 0.37_dp
        second_phase = -0.91_dp
        do periods = 1, 7
            do first_n = -2 * periods, 2 * periods
                do second_n = -2 * periods, 2 * periods
                    do first_kind = phase_cosine, phase_sine
                        do second_kind = phase_cosine, phase_sine
                            expected = 0.0_dp
                            do p = 0, periods - 1
                                expected = expected + phase_value(first_kind, &
                                    first_phase - two_pi * real(first_n * p, &
                                    dp) / real(periods, dp)) &
                                    * phase_value(second_kind, &
                                    second_phase - two_pi &
                                    * real(second_n * p, dp) &
                                    / real(periods, dp))
                            end do
                            expected = expected / real(periods, dp)
                            actual = phase_product_average(first_kind, &
                                second_kind, first_phase, second_phase, &
                                first_n, second_n, periods)
                            call require(abs(actual - expected) < 2.0e-14_dp, &
                                "phase-product table disagrees with direct sum")
                        end do
                    end do
                end do
            end do
        end do
    end subroutine check_product_table

    pure function phase_value(kind, phase) result(value)
        integer, intent(in) :: kind
        real(dp), intent(in) :: phase
        real(dp) :: value

        value = cos(phase)
        if (kind == phase_sine) value = sin(phase)
    end function phase_value

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_phase_factor_topology
