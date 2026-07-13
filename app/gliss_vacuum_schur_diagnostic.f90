program gliss_vacuum_schur_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use vacuum_schur, only: eliminate_vacuum, vacuum_schur_ok
    implicit none

    real(dp), parameter :: plasma(2, 2) = reshape([4.0_dp, 1.0_dp, &
        1.0_dp, 3.0_dp], [2, 2])
    real(dp), parameter :: vacuum(2, 2) = reshape([3.0_dp, 0.5_dp, &
        0.5_dp, 2.0_dp], [2, 2])
    real(dp), parameter :: coupling(2, 2) = reshape([0.4_dp, 0.1_dp, &
        -0.2_dp, 0.3_dp], [2, 2])
    real(dp), parameter :: plasma_state(2) = [0.7_dp, -0.4_dp]
    real(dp), parameter :: scales(3) = [1.0_dp, 10.0_dp, 100.0_dp]
    integer :: i

    if (command_argument_count() /= 0) call fail_usage()
    write (*, "(a)") "vacuum_scale,effective_11,effective_12," // &
        "effective_21,effective_22,response_11,response_12,response_21," // &
        "response_22,response_residual,energy_closure,symmetry_defect," // &
        "fixed_boundary_error"
    do i = 1, size(scales)
        call report_scale(scales(i))
    end do

contains

    subroutine fail_usage()
        write (error_unit, "(a)") &
            "gliss_vacuum_schur_diagnostic: no arguments are accepted"
        flush (error_unit)
        stop 2
    end subroutine fail_usage

    subroutine report_scale(scale)
        real(dp), intent(in) :: scale
        real(dp), allocatable :: effective(:, :), response(:, :)
        real(dp) :: scaled_vacuum(2, 2), vacuum_state(2)
        real(dp) :: response_residual, full_energy, reduced_energy
        real(dp) :: energy_closure, symmetry_defect, fixed_boundary_error
        integer :: info

        scaled_vacuum = scale * vacuum
        call eliminate_vacuum(plasma, scaled_vacuum, coupling, effective, &
            response, info)
        if (info /= vacuum_schur_ok) then
            write (error_unit, "(a)") "vacuum Schur elimination failed"
            error stop 1
        end if
        vacuum_state = matmul(response, plasma_state)
        response_residual = maxval(abs(matmul(scaled_vacuum, response) &
            + transpose(coupling)))
        full_energy = 0.5_dp * dot_product(plasma_state, &
            matmul(plasma, plasma_state)) + dot_product(plasma_state, &
            matmul(coupling, vacuum_state)) + 0.5_dp * &
            dot_product(vacuum_state, matmul(scaled_vacuum, vacuum_state))
        reduced_energy = 0.5_dp * dot_product(plasma_state, &
            matmul(effective, plasma_state))
        energy_closure = abs(full_energy - reduced_energy)
        symmetry_defect = maxval(abs(effective - transpose(effective)))
        fixed_boundary_error = maxval(abs(effective - plasma))
        write (*, "(es24.16,12(',',es24.16))") scale, &
            effective(1, 1), effective(1, 2), effective(2, 1), &
            effective(2, 2), response(1, 1), response(1, 2), &
            response(2, 1), response(2, 2), response_residual, &
            energy_closure, symmetry_defect, fixed_boundary_error
    end subroutine report_scale

end program gliss_vacuum_schur_diagnostic
