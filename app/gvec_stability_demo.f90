program gvec_stability_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use local_mode_model, only: assemble_local_mode, vacuum_permeability
    use symmetric_eigensolver, only: solve_three_component_modes
    implicit none

    real(dp) :: stiffness(3, 3), mass(3, 3)
    real(dp) :: eigenvalues(3), eigenvectors(3, 3)
    real(dp), parameter :: wave_vector(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: magnetic_field(3) = [0.0_dp, 0.0_dp, 1.0e-3_dp]
    real(dp), parameter :: normal(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: drive
    integer :: info

    drive = 2.0_dp * magnetic_field(3)**2 / vacuum_permeability
    call assemble_local_mode(wave_vector, magnetic_field, 1.0_dp, 2.0_dp, &
        5.0_dp / 3.0_dp, drive, normal, stiffness, mass)
    call solve_three_component_modes(stiffness, mass, eigenvalues, &
        eigenvectors, info)
    if (info /= 0) error stop "generalized eigensolver failed"

    write (*, "(a,3es16.8)") "omega^2 [s^-2]: ", eigenvalues
    if (eigenvalues(1) < 0.0_dp) then
        write (*, "(a,es16.8)") "unstable growth rate [s^-1]: ", &
            sqrt(-eigenvalues(1))
    else
        write (*, "(a)") "stable"
    end if
end program gvec_stability_demo
