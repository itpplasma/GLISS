program test_compatible_problem_assembly_support
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_problem_assembly_support, only: apply_l2_stored_power, &
        apply_stored_power, build_active_indices, build_uniform_breaks, &
        mode_table_is_unique, &
        evaluate_generalized_eigenpair, quadratic_form, &
        replicate_indexed_values, scatter_matrix, sum_tensor, &
        symmetrize_matrix, symmetrize_tensor
    implicit none

    call verify_stored_power_transform()
    call verify_local_index_helpers()
    call verify_mapped_scatter()
    call verify_matrix_reductions()
    call verify_eigenpair_diagnostics()
    call verify_mode_uniqueness()
    write (*, "(a)") "PASS"

contains

    subroutine verify_stored_power_transform()
        real(dp), parameter :: h1(2) = [2.0_dp, -3.0_dp]
        real(dp), parameter :: dh1(2) = [5.0_dp, 7.0_dp]
        real(dp), parameter :: powers(3) = [0.0_dp, 1.0_dp, -1.0_dp]
        integer, parameter :: indices(2) = [1, 2]
        real(dp) :: values(2, 3), derivatives(2, 3)
        integer :: info

        call apply_stored_power(0.25_dp, powers, h1, dh1, indices, values, &
            derivatives, info)
        call require(info == 0, "stored-power transform rejected valid input")
        call require(maxval(abs(values - reshape([2.0_dp, -3.0_dp, &
            8.0_dp, -12.0_dp, 0.5_dp, -0.75_dp], [2, 3]))) &
            < 1.0e-14_dp, "stored-power values differ")
        call require(maxval(abs(derivatives - reshape([5.0_dp, 7.0_dp, &
            -12.0_dp, 76.0_dp, 3.25_dp, -1.25_dp], [2, 3]))) &
            < 1.0e-14_dp, "stored-power product rule differs")
        call apply_stored_power(0.0_dp, powers, h1, dh1, indices, values, &
            derivatives, info)
        call require(info == -1, "axis stored-power evaluation was accepted")
        call apply_stored_power(0.25_dp, powers(:2), h1, dh1, indices, values, &
            derivatives, info)
        call require(info == -1, "mismatched stored-power table was accepted")
        call apply_l2_stored_power(0.25_dp, powers, h1, indices, values, info)
        call require(info == 0, "L2 stored-power transform rejected valid input")
        call require(maxval(abs(values - reshape([2.0_dp, -3.0_dp, &
            8.0_dp, -12.0_dp, 0.5_dp, -0.75_dp], [2, 3]))) &
            < 1.0e-14_dp, "L2 stored-power values differ")
        call apply_l2_stored_power(0.0_dp, powers, h1, indices, values, info)
        call require(info == -1, "axis L2 stored-power evaluation was accepted")
    end subroutine verify_stored_power_transform

    subroutine verify_local_index_helpers()
        real(dp), parameter :: values(5) = [2.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, -3.0_dp]
        real(dp), parameter :: derivatives(5) = [0.0_dp, 4.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp]
        real(dp) :: breaks(5), repeated(3, 2)
        integer, allocatable :: indices(:)
        integer :: info

        call build_active_indices(values, indices, info, derivatives)
        call require(info == 0 .and. all(indices == [1, 2, 5]), &
            "active first, derivative-only, or last index differs")
        call replicate_indexed_values(values, indices, repeated, info)
        call require(info == 0, "indexed replication rejected valid input")
        call require(maxval(abs(repeated - reshape([2.0_dp, 0.0_dp, &
            -3.0_dp, 2.0_dp, 0.0_dp, -3.0_dp], [3, 2]))) &
            < 1.0e-14_dp, "indexed replication differs")
        call build_uniform_breaks(4, breaks, info)
        call require(info == 0 .and. maxval(abs(breaks &
            - [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp])) &
            < 1.0e-14_dp, "uniform break endpoints or cell count differ")
        call build_uniform_breaks(5, breaks, info)
        call require(info == -1, "off-by-one break array was accepted")
    end subroutine verify_local_index_helpers

    subroutine verify_mapped_scatter()
        integer, parameter :: map(3) = [2, 0, 1]
        real(dp) :: global(2, 2), local(3, 3)

        local = reshape([11.0_dp, 21.0_dp, 31.0_dp, 12.0_dp, 22.0_dp, &
            32.0_dp, 13.0_dp, 23.0_dp, 33.0_dp], [3, 3])
        global = 0.0_dp
        call scatter_matrix(map, local, 2.0_dp, global)
        call require(maxval(abs(global - reshape([66.0_dp, 26.0_dp, &
            62.0_dp, 22.0_dp], [2, 2]))) < 1.0e-14_dp, &
            "mapped scatter has an inactive or off-by-one index")
    end subroutine verify_mapped_scatter

    subroutine verify_matrix_reductions()
        real(dp) :: matrix(2, 2), tensor(2, 2, 2), total(2, 2)

        matrix = reshape([1.0_dp, 5.0_dp, 3.0_dp, 7.0_dp], [2, 2])
        call symmetrize_matrix(matrix)
        call require(maxval(abs(matrix - reshape([1.0_dp, 4.0_dp, &
            4.0_dp, 7.0_dp], [2, 2]))) < 1.0e-14_dp, &
            "matrix symmetrization differs")
        tensor(:, :, 1) = reshape([1.0_dp, 6.0_dp, 2.0_dp, 4.0_dp], [2, 2])
        tensor(:, :, 2) = reshape([3.0_dp, 8.0_dp, 4.0_dp, 5.0_dp], [2, 2])
        call symmetrize_tensor(tensor)
        call sum_tensor(tensor, total)
        call require(maxval(abs(total - transpose(total))) < 1.0e-14_dp, &
            "tensor symmetrization is not symmetric")
        call require(maxval(abs(total - reshape([4.0_dp, 10.0_dp, &
            10.0_dp, 9.0_dp], [2, 2]))) < 1.0e-14_dp, &
            "tensor reduction differs")
    end subroutine verify_matrix_reductions

    subroutine verify_eigenpair_diagnostics()
        real(dp), parameter :: stiffness(2, 2) = reshape( &
            [2.0_dp, 0.0_dp, 0.0_dp, 6.0_dp], [2, 2])
        real(dp), parameter :: mass(2, 2) = reshape( &
            [1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2])
        real(dp), parameter :: vector(2) = [1.0_dp, 0.0_dp]
        real(dp) :: kinetic, potential, residual, value
        integer :: info

        call evaluate_generalized_eigenpair(stiffness, mass, vector, 2.0_dp, &
            kinetic, potential, residual, info)
        call require(info == 0 .and. abs(kinetic - 1.0_dp) < 1.0e-14_dp, &
            "kinetic-energy diagnostic differs")
        call require(abs(potential - 2.0_dp) < 1.0e-14_dp, &
            "potential-energy diagnostic differs")
        call require(residual < 1.0e-14_dp, "exact eigenpair residual is nonzero")
        call quadratic_form(stiffness, vector, value, info)
        call require(info == 0 .and. abs(value - potential) < 1.0e-14_dp, &
            "quadratic form differs from eigenpair diagnostic")
        call evaluate_generalized_eigenpair(stiffness, mass, vector, 3.0_dp, &
            kinetic, potential, residual, info)
        call require(info == 0 .and. residual > 0.0_dp, &
            "incorrect eigenvalue produced a zero residual")
    end subroutine verify_eigenpair_diagnostics

    subroutine verify_mode_uniqueness()
        call require(mode_table_is_unique([0, 1, 2], [1, -1, -1]), &
            "unique mode table was rejected")
        call require(.not. mode_table_is_unique([0, 1, 1], [1, -1, -1]), &
            "duplicate final mode was accepted")
        call require(.not. mode_table_is_unique([0, 1], [1]), &
            "mismatched mode arrays were accepted")
    end subroutine verify_mode_uniqueness

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_compatible_problem_assembly_support
