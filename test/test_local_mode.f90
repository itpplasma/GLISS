program test_local_mode
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use local_mode_model, only: assemble_local_mode, vacuum_permeability
    use symmetric_eigensolver, only: solve_three_component_modes
    implicit none

    real(dp) :: stiffness(3, 3), mass(3, 3)
    real(dp) :: eigenvalues(3), eigenvectors(3, 3)
    real(dp), parameter :: wave_vector(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: magnetic_field(3) = [0.0_dp, 0.0_dp, 1.0e-3_dp]
    real(dp), parameter :: normal(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: alfven_squared, expected(3), residual(3)
    integer :: i, info

    alfven_squared = magnetic_field(3)**2 / &
        (vacuum_permeability * 2.0_dp)
    expected = [alfven_squared, alfven_squared, 5.0_dp / 6.0_dp]

    call assemble_local_mode(wave_vector, magnetic_field, 1.0_dp, 2.0_dp, &
        5.0_dp / 3.0_dp, 0.0_dp, normal, stiffness, mass)
    call require(maxval(abs(stiffness - transpose(stiffness))) < 1.0e-12_dp, &
        "stiffness is not symmetric")
    call solve_three_component_modes(stiffness, mass, eigenvalues, &
        eigenvectors, info)
    call require(info == 0, "stable eigenproblem failed")
    call require(maxval(abs(eigenvalues - expected)) < 1.0e-11_dp, &
        "stable spectrum is wrong")

    do i = 1, 3
        residual = matmul(stiffness, eigenvectors(:, i)) - &
            eigenvalues(i) * matmul(mass, eigenvectors(:, i))
        call require(maxval(abs(residual)) < 1.0e-11_dp, &
            "eigenpair residual is too large")
    end do

    call assemble_local_mode(wave_vector, magnetic_field, 1.0_dp, 2.0_dp, &
        5.0_dp / 3.0_dp, 2.0_dp * magnetic_field(3)**2 / &
        vacuum_permeability, normal, stiffness, mass)
    call solve_three_component_modes(stiffness, mass, eigenvalues, &
        eigenvectors, info)
    call require(info == 0, "unstable eigenproblem failed")
    call require(eigenvalues(1) < 0.0_dp, &
        "pressure-curvature drive did not destabilize the normal mode")

    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_local_mode
