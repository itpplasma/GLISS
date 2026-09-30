program test_terpsichore_lowest
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use terpsichore_fixed_boundary_spectrum, only: &
        solve_terpsichore_fixed_boundary_file, solve_terpsichore_lowest, &
        terpsichore_fixed_boundary_result_t, terpsichore_fixed_spectrum_ok
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    implicit none

    call check_pencil([2.0_dp, 3.0_dp, 5.0_dp], 0, 2.0_dp)
    call check_pencil([3.0_dp, -1.0e-3_dp, 5.0_dp], 1, -1.0e-3_dp)
    call check_pencil([-4.0_dp, -2.0_dp, 5.0_dp], 2, -4.0_dp)
    call check_vacuum_file_rejected()
    print "(a)", "PASS"

contains

    subroutine check_pencil(diagonal, expected_count, expected_value)
        ! Stable pencils report the lowest nonnegative eigenpair instead of
        ! failing; unstable ones the lowest negative eigenpair.
        real(dp), intent(in) :: diagonal(3), expected_value
        integer, intent(in) :: expected_count
        type(variable_block_tridiagonal_t) :: stiffness, mass
        real(dp) :: dense(3, 3), identity(3, 3), eigenvalue, residual
        real(dp) :: resolution, certificate
        real(dp), allocatable :: vector(:)
        character(len=256) :: message
        integer :: count, i, info

        dense = 0.0_dp
        identity = 0.0_dp
        do i = 1, 3
            dense(i, i) = diagonal(i)
            identity(i, i) = 1.0_dp
        end do
        call pack_variable_blocks(dense, [1, 2], stiffness, info)
        call require(info == variable_block_ok, "stiffness packing failed")
        call pack_variable_blocks(identity, [1, 2], mass, info)
        call require(info == variable_block_ok, "mass packing failed")
        call solve_terpsichore_lowest(stiffness, mass, eigenvalue, vector, &
            residual, resolution, certificate, count, info, message)
        call require(info == terpsichore_fixed_spectrum_ok, trim(message))
        call require(count == expected_count, "wrong negative count")
        call require(abs(eigenvalue - expected_value) <= 1.0e-10_dp &
            * abs(expected_value), "wrong lowest eigenvalue")
    end subroutine check_pencil

    subroutine check_vacuum_file_rejected()
        ! An IVAC>0 FORT.23 has a longer radial record than the fixed layout.
        character(len=*), parameter :: path = "terpsichore_vacuum_probe.bin"
        type(terpsichore_fixed_boundary_result_t) :: result
        character(len=256) :: message
        real(dp) :: record(7 * 2 + 4 + 3)
        integer :: info, unit

        record = 0.5_dp
        open (newunit=unit, file=path, status="replace", action="write", &
            access="sequential", form="unformatted")
        write (unit) 2, 1, 1, 1, 1, 1
        write (unit) record
        close (unit)
        call solve_terpsichore_fixed_boundary_file(path, result, info, message)
        open (newunit=unit, file=path, status="old")
        close (unit, status="delete")
        call require(info /= terpsichore_fixed_spectrum_ok, &
            "vacuum FORT.23 was accepted")
        call require(index(message, "IVAC>0") > 0, &
            "vacuum FORT.23 error does not name IVAC")
    end subroutine check_vacuum_file_rejected

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_lowest
