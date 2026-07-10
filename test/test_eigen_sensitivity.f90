program test_eigen_sensitivity
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use local_mode_model, only: assemble_local_mode, vacuum_permeability
    use symmetric_eigensolver, only: solve_three_component_modes
    implicit none

    real(dp), parameter :: step = 1.0e-5_dp
    real(dp), parameter :: wave_vector(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: magnetic_field(3) = [0.0_dp, 0.0_dp, 1.0e-3_dp]
    real(dp), parameter :: normal(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: drive, scale, numerical, analytic

    scale = magnetic_field(3)**2 / vacuum_permeability
    drive = 0.5_dp * scale
    numerical = (lowest_eigenvalue(drive + step * scale) - &
        lowest_eigenvalue(drive - step * scale)) / (2.0_dp * step * scale)
    analytic = -0.5_dp
    if (abs(numerical - analytic) > 1.0e-8_dp) then
        write (error_unit, "(a,2es24.16)") &
            "FAIL eigenvalue sensitivity ", numerical, analytic
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    function lowest_eigenvalue(current_drive) result(value)
        real(dp), intent(in) :: current_drive
        real(dp) :: value
        real(dp) :: stiffness(3, 3), mass(3, 3)
        real(dp) :: eigenvalues(3), eigenvectors(3, 3)
        integer :: info

        call assemble_local_mode(wave_vector, magnetic_field, 1.0_dp, &
            2.0_dp, 5.0_dp / 3.0_dp, current_drive, normal, &
            stiffness, mass)
        call solve_three_component_modes(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        if (info /= 0) error stop "generalized eigensolver failed"
        value = eigenvalues(1)
    end function lowest_eigenvalue

end program test_eigen_sensitivity
