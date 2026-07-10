program test_phase_factor_topology
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use phase_factor_topology, only: build_phase_factor_table, &
        evaluate_phase_factor, phase_factor_invalid, phase_factor_ok, &
        phase_factor_table_t
    implicit none

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
    integer, parameter :: mode_m(4) = [1, 2, 3, 4]
    integer, parameter :: mode_n(4) = [1, -1, 4, -4]
    type(phase_factor_table_t) :: table
    real(dp) :: theta, zeta, cosine, sine, phase
    integer :: info, mode

    call build_phase_factor_table(3, 1, 1, mode_m, mode_n, table, info)
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
        call evaluate_phase_factor(table, mode, theta, zeta, cosine, sine)
        phase = two_pi * (real(mode_m(mode), dp) * theta &
            - real(mode_n(mode), dp) * zeta / 3.0_dp)
        call require(abs(cosine - cos(phase)) < 1.0e-14_dp, &
            "cosine phase reconstruction failed")
        call require(abs(sine - sin(phase)) < 1.0e-14_dp, &
            "sine phase reconstruction failed")
    end do

    call build_phase_factor_table(5, 2, 0, [0, 1], [2, -2], table, info)
    call require(info == phase_factor_ok, "valid odd topology was rejected")
    call build_phase_factor_table(4, 2, 0, [0], [2], table, info)
    call require(info == phase_factor_invalid, &
        "half-integer even topology was accepted")
    call build_phase_factor_table(3, 1, 0, [0], [0], table, info)
    call require(info == phase_factor_invalid, &
        "mode outside the family was accepted")
    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_phase_factor_topology
