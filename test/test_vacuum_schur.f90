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
        real(dp) :: plasma_image(2), coupling_image(2), vacuum_image(2)
        real(dp) :: effective_image(2)
        real(dp) :: full_energy, reduced_energy
        integer :: info

        inverse = reshape([2.0_dp, -0.5_dp, -0.5_dp, 3.0_dp], [2, 2]) &
            / 5.75_dp
        call build_reference_elimination(plasma, inverse, coupling, &
            expected_response, expected_effective)
        call eliminate_vacuum(plasma, vacuum, coupling, effective, response, info)
        call require(info == vacuum_schur_ok, "valid elimination failed")
        call require(maxval(abs(response - expected_response)) < 1.0e-14_dp, &
            "vacuum response disagrees with direct inverse")
        call require(maxval(abs(effective - expected_effective)) < 1.0e-14_dp, &
            "effective operator disagrees with direct Schur complement")
        call require(stationarity_error(vacuum, response, coupling) &
            < 1.0e-14_dp, "response is not stationary")
        call require(symmetry_error(effective) < &
            1.0e-14_dp, "effective operator is not symmetric")
        state = [0.7_dp, -0.4_dp]
        call matrix_vector_product(response, state, vacuum_state)
        call matrix_vector_product(plasma, state, plasma_image)
        call matrix_vector_product(coupling, vacuum_state, coupling_image)
        call matrix_vector_product(vacuum, vacuum_state, vacuum_image)
        call matrix_vector_product(effective, state, effective_image)
        full_energy = 0.5_dp * dot_product(state, plasma_image) &
            + dot_product(state, coupling_image) &
            + 0.5_dp * dot_product(vacuum_state, vacuum_image)
        reduced_energy = 0.5_dp * dot_product(state, effective_image)
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
        call build_reference_elimination(plasma_three, inverse, &
            coupling_three, expected_response, expected_effective)
        call eliminate_vacuum(plasma_three, vacuum, coupling_three, effective, &
            response, info)
        call require(info == vacuum_schur_ok, &
            "rectangular-coupling elimination failed")
        call require(size(response, 1) == 2 .and. size(response, 2) == 3, &
            "rectangular vacuum response has the wrong shape")
        call require(maxval(abs(response - expected_response)) < 1.0e-14_dp, &
            "rectangular vacuum response is transposed or incorrect")
        call require(maxval(abs(effective - expected_effective)) < 1.0e-14_dp, &
            "rectangular Schur complement is incorrect")
    end subroutine test_rectangular_coupling

    subroutine test_fixed_boundary_limit()
        real(dp), allocatable :: effective(:, :), response(:, :)
        real(dp) :: errors(3), scaled_vacuum(2, 2), scales(3)
        integer :: i, info

        scales = [1.0_dp, 10.0_dp, 100.0_dp]
        do i = 1, 3
            scaled_vacuum = scales(i) * vacuum
            call eliminate_vacuum(plasma, scaled_vacuum, coupling, &
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

    pure subroutine build_reference_elimination(plasma_matrix, vacuum_inverse, &
            coupling_matrix, expected_response, expected_effective)
        real(dp), intent(in) :: plasma_matrix(:, :), vacuum_inverse(:, :)
        real(dp), intent(in) :: coupling_matrix(:, :)
        real(dp), intent(out) :: expected_response(:, :)
        real(dp), intent(out) :: expected_effective(:, :)
        integer :: column, inner, row

        do column = 1, size(coupling_matrix, 1)
            do row = 1, size(vacuum_inverse, 1)
                expected_response(row, column) = 0.0_dp
                do inner = 1, size(vacuum_inverse, 2)
                    expected_response(row, column) = &
                        expected_response(row, column) &
                        - vacuum_inverse(row, inner) &
                        * coupling_matrix(column, inner)
                end do
            end do
        end do
        do column = 1, size(plasma_matrix, 2)
            do row = 1, size(plasma_matrix, 1)
                expected_effective(row, column) = plasma_matrix(row, column)
                do inner = 1, size(coupling_matrix, 2)
                    expected_effective(row, column) = &
                        expected_effective(row, column) &
                        + coupling_matrix(row, inner) &
                        * expected_response(inner, column)
                end do
            end do
        end do
    end subroutine build_reference_elimination

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

    pure function stationarity_error(vacuum_matrix, response_matrix, &
            coupling_matrix) result(error)
        real(dp), intent(in) :: vacuum_matrix(:, :), response_matrix(:, :)
        real(dp), intent(in) :: coupling_matrix(:, :)
        real(dp) :: error, residual
        integer :: column, inner, row

        error = 0.0_dp
        do column = 1, size(response_matrix, 2)
            do row = 1, size(response_matrix, 1)
                residual = coupling_matrix(column, row)
                do inner = 1, size(vacuum_matrix, 2)
                    residual = residual + vacuum_matrix(row, inner) &
                        * response_matrix(inner, column)
                end do
                error = max(error, abs(residual))
            end do
        end do
    end function stationarity_error

    pure function symmetry_error(matrix) result(error)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: error
        integer :: column, row

        error = 0.0_dp
        do column = 1, size(matrix, 2)
            do row = 1, column - 1
                error = max(error, abs(matrix(row, column) &
                    - matrix(column, row)))
            end do
        end do
    end function symmetry_error

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_vacuum_schur
