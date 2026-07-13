program gliss_local_mode_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use local_mode_model, only: assemble_local_mode, vacuum_permeability
    use symmetric_eigensolver, only: solve_three_component_modes
    implicit none

    real(dp), parameter :: magnetic_field(3) = [0.0_dp, 0.0_dp, 1.0e-3_dp]
    real(dp), parameter :: normal(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: pressure = 1.0_dp
    real(dp), parameter :: density = 2.0_dp
    real(dp), parameter :: gamma_value = 5.0_dp / 3.0_dp

    if (command_argument_count() /= 0) call fail_usage()
    write (*, "(a)") "case,mode,kx_m_minus_1,ky_m_minus_1,kz_m_minus_1," // &
        "gamma,drive_pa,omega2_s_minus_2,abs_x,abs_y,abs_z,residual"
    call report_case("parallel", [0.0_dp, 0.0_dp, 1.0_dp], gamma_value, &
        0.0_dp)
    call report_case("oblique", [0.6_dp, 0.0_dp, 0.8_dp], gamma_value, &
        0.0_dp)
    call report_case("scaled_oblique", [1.5_dp, 0.0_dp, 2.0_dp], &
        gamma_value, 0.0_dp)
    call report_case("cold_oblique", [0.6_dp, 0.0_dp, 0.8_dp], 0.0_dp, &
        0.0_dp)
    call report_case("perpendicular", [1.0_dp, 0.0_dp, 0.0_dp], &
        gamma_value, 0.0_dp)
    call report_case("driven_parallel", [0.0_dp, 0.0_dp, 1.0_dp], &
        gamma_value, 2.0_dp)

contains

    subroutine fail_usage()
        write (error_unit, "(a)") &
            "gliss_local_mode_diagnostic: no arguments are accepted"
        flush (error_unit)
        stop 2
    end subroutine fail_usage

    subroutine report_case(name, wave_vector, gamma, drive_factor)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: wave_vector(3), gamma, drive_factor
        real(dp) :: stiffness(3, 3), mass(3, 3)
        real(dp) :: eigenvalues(3), eigenvectors(3, 3)
        real(dp) :: drive, norm, residual(3)
        integer :: i, info

        drive = drive_factor * dot_product(magnetic_field, magnetic_field) &
            / vacuum_permeability
        call assemble_local_mode(wave_vector, magnetic_field, pressure, &
            density, gamma, drive, normal, stiffness, mass)
        call solve_three_component_modes(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        if (info /= 0) then
            write (error_unit, "(a)") "local-mode eigensolve failed"
            error stop 1
        end if
        do i = 1, 3
            norm = norm2(eigenvectors(:, i))
            if (norm <= 0.0_dp) then
                write (error_unit, "(a)") "local-mode eigenvector is zero"
                error stop 1
            end if
            residual = matmul(stiffness, eigenvectors(:, i)) &
                - eigenvalues(i) * matmul(mass, eigenvectors(:, i))
            write (*, "(a,',',i0,10(',',es24.16))") trim(name), i, &
                wave_vector, gamma, drive, eigenvalues(i), &
                abs(eigenvectors(:, i)) / norm, maxval(abs(residual))
        end do
    end subroutine report_case

end program gliss_local_mode_diagnostic
