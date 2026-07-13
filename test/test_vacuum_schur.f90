program test_vacuum_schur
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use vacuum_schur, only: eliminate_vacuum, vacuum_schur_invalid_input, &
        vacuum_schur_not_spd, vacuum_schur_ok
    implicit none

    real(dp), parameter :: plasma(2, 2) = reshape([4.0_dp, 1.0_dp, &
        1.0_dp, 3.0_dp], [2, 2])
    real(dp), parameter :: vacuum(2, 2) = reshape([3.0_dp, 0.5_dp, &
        0.5_dp, 2.0_dp], [2, 2])
    real(dp), parameter :: coupling(2, 2) = reshape([0.4_dp, 0.1_dp, &
        -0.2_dp, 0.3_dp], [2, 2])

    call test_exact_elimination()
    call test_rectangular_coupling()
    call test_fixed_boundary_limit()
    call test_invalid_inputs()
    write (*, "(a)") "PASS"

contains

    subroutine test_exact_elimination()
        real(dp), allocatable :: effective(:, :), response(:, :)
        real(dp) :: inverse(2, 2), expected_response(2, 2)
        real(dp) :: expected_effective(2, 2), state(2), vacuum_state(2)
        real(dp) :: full_energy, reduced_energy
        integer :: info

        inverse = reshape([2.0_dp, -0.5_dp, -0.5_dp, 3.0_dp], [2, 2]) &
            / 5.75_dp
        expected_response = -matmul(inverse, transpose(coupling))
        expected_effective = plasma + matmul(coupling, expected_response)
        call eliminate_vacuum(plasma, vacuum, coupling, effective, response, info)
        call require(info == vacuum_schur_ok, "valid elimination failed")
        call require(maxval(abs(response - expected_response)) < 1.0e-14_dp, &
            "vacuum response disagrees with direct inverse")
        call require(maxval(abs(effective - expected_effective)) < 1.0e-14_dp, &
            "effective operator disagrees with direct Schur complement")
        call require(maxval(abs(matmul(vacuum, response) + &
            transpose(coupling))) < 1.0e-14_dp, "response is not stationary")
        call require(maxval(abs(effective - transpose(effective))) < &
            1.0e-14_dp, "effective operator is not symmetric")
        state = [0.7_dp, -0.4_dp]
        vacuum_state = matmul(response, state)
        full_energy = 0.5_dp * dot_product(state, matmul(plasma, state)) &
            + dot_product(state, matmul(coupling, vacuum_state)) &
            + 0.5_dp * dot_product(vacuum_state, matmul(vacuum, vacuum_state))
        reduced_energy = 0.5_dp * dot_product(state, matmul(effective, state))
        call require(abs(full_energy - reduced_energy) < 1.0e-14_dp, &
            "reduced energy does not equal stationary full energy")
    end subroutine test_exact_elimination

    subroutine test_rectangular_coupling()
        real(dp), parameter :: plasma_three(3, 3) = reshape([ &
            4.0_dp, 1.0_dp, 0.2_dp, 1.0_dp, 3.0_dp, -0.1_dp, &
            0.2_dp, -0.1_dp, 2.0_dp], [3, 3])
        real(dp), parameter :: coupling_three(3, 2) = reshape([ &
            0.4_dp, 0.1_dp, -0.3_dp, -0.2_dp, 0.3_dp, 0.25_dp], [3, 2])
        real(dp), allocatable :: effective(:, :), response(:, :)
        real(dp) :: inverse(2, 2), expected_response(2, 3)
        real(dp) :: expected_effective(3, 3)
        integer :: info

        inverse = reshape([2.0_dp, -0.5_dp, -0.5_dp, 3.0_dp], [2, 2]) &
            / 5.75_dp
        expected_response = -matmul(inverse, transpose(coupling_three))
        expected_effective = plasma_three + &
            matmul(coupling_three, expected_response)
        call eliminate_vacuum(plasma_three, vacuum, coupling_three, effective, &
            response, info)
        call require(info == vacuum_schur_ok, &
            "rectangular-coupling elimination failed")
        call require(all(shape(response) == [2, 3]), &
            "rectangular vacuum response has the wrong shape")
        call require(maxval(abs(response - expected_response)) < 1.0e-14_dp, &
            "rectangular vacuum response is transposed or incorrect")
        call require(maxval(abs(effective - expected_effective)) < 1.0e-14_dp, &
            "rectangular Schur complement is incorrect")
    end subroutine test_rectangular_coupling

    subroutine test_fixed_boundary_limit()
        real(dp), allocatable :: effective(:, :), response(:, :)
        real(dp) :: errors(3), scales(3)
        integer :: i, info

        scales = [1.0_dp, 10.0_dp, 100.0_dp]
        do i = 1, 3
            call eliminate_vacuum(plasma, scales(i) * vacuum, coupling, &
                effective, response, info)
            call require(info == vacuum_schur_ok, "scaled elimination failed")
            errors(i) = maxval(abs(effective - plasma))
        end do
        call require(errors(1) / errors(2) > 9.0_dp, &
            "first fixed-boundary contraction is too slow")
        call require(errors(2) / errors(3) > 9.0_dp, &
            "second fixed-boundary contraction is too slow")
    end subroutine test_fixed_boundary_limit

    subroutine test_invalid_inputs()
        real(dp), allocatable :: effective(:, :), response(:, :)
        real(dp) :: bad_plasma(2, 2), bad_vacuum(2, 2)
        integer :: info

        call eliminate_vacuum(plasma, vacuum, coupling(:, :1), effective, &
            response, info)
        call require(info == vacuum_schur_invalid_input, &
            "coupling shape mismatch was accepted")
        bad_plasma = plasma
        bad_plasma(1, 2) = 2.0_dp
        call eliminate_vacuum(bad_plasma, vacuum, coupling, effective, response, &
            info)
        call require(info == vacuum_schur_invalid_input, &
            "nonsymmetric plasma operator was accepted")
        bad_vacuum = vacuum
        bad_vacuum(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call eliminate_vacuum(plasma, bad_vacuum, coupling, effective, response, &
            info)
        call require(info == vacuum_schur_invalid_input, &
            "nonfinite vacuum operator was accepted")
        bad_vacuum = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
        call eliminate_vacuum(plasma, bad_vacuum, coupling, effective, response, &
            info)
        call require(info == vacuum_schur_not_spd, &
            "indefinite vacuum operator was accepted")
    end subroutine test_invalid_inputs

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_vacuum_schur
