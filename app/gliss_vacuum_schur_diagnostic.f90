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
        real(dp) :: coupling_image(2), effective_image(2), plasma_image(2)
        real(dp) :: scaled_vacuum(2, 2), vacuum_image(2), vacuum_state(2)
        real(dp) :: response_residual, full_energy, reduced_energy
        real(dp) :: energy_closure, symmetry_defect, fixed_boundary_error
        real(dp) :: residual_value
        integer :: column, info, inner, row

        scaled_vacuum = scale * vacuum
        call eliminate_vacuum(plasma, scaled_vacuum, coupling, effective, &
            response, info)
        if (info /= vacuum_schur_ok) then
            write (error_unit, "(a)") "vacuum Schur elimination failed"
            error stop 1
        end if
        call matrix_vector_product(response, plasma_state, vacuum_state)
        response_residual = 0.0_dp
        do column = 1, size(response, 2)
            do row = 1, size(response, 1)
                residual_value = coupling(column, row)
                do inner = 1, size(response, 1)
                    residual_value = residual_value &
                        + scaled_vacuum(row, inner) * response(inner, column)
                end do
                response_residual = max(response_residual, abs(residual_value))
            end do
        end do
        call matrix_vector_product(plasma, plasma_state, plasma_image)
        call matrix_vector_product(coupling, vacuum_state, coupling_image)
        call matrix_vector_product(scaled_vacuum, vacuum_state, vacuum_image)
        call matrix_vector_product(effective, plasma_state, effective_image)
        full_energy = 0.5_dp * dot_product(plasma_state, plasma_image) &
            + dot_product(plasma_state, coupling_image) + 0.5_dp &
            * dot_product(vacuum_state, vacuum_image)
        reduced_energy = 0.5_dp * dot_product(plasma_state, effective_image)
        energy_closure = abs(full_energy - reduced_energy)
        symmetry_defect = abs(effective(1, 2) - effective(2, 1))
        fixed_boundary_error = 0.0_dp
        do column = 1, 2
            do row = 1, 2
                fixed_boundary_error = max(fixed_boundary_error, &
                    abs(effective(row, column) - plasma(row, column)))
            end do
        end do
        write (*, "(es24.16,12(',',es24.16))") scale, &
            effective(1, 1), effective(1, 2), effective(2, 1), &
            effective(2, 2), response(1, 1), response(1, 2), &
            response(2, 1), response(2, 2), response_residual, &
            energy_closure, symmetry_defect, fixed_boundary_error
    end subroutine report_scale

    pure subroutine matrix_vector_product(matrix, vector, image)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp), intent(out) :: image(:)
        integer :: column, row

        do row = 1, size(matrix, 1)
            image(row) = 0.0_dp
            do column = 1, size(matrix, 2)
                image(row) = image(row) + matrix(row, column) * vector(column)
            end do
        end do
    end subroutine matrix_vector_product

end program gliss_vacuum_schur_diagnostic
