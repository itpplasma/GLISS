program test_local_mode
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use local_mode_model, only: assemble_local_mode, vacuum_permeability
    use symmetric_eigensolver, only: solve_three_component_modes
    implicit none

    real(dp), parameter :: magnetic_field(3) = [0.0_dp, 0.0_dp, 1.0e-3_dp]
    real(dp), parameter :: normal(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: pressure = 1.0_dp
    real(dp), parameter :: density = 2.0_dp
    real(dp), parameter :: adiabatic_index = 5.0_dp / 3.0_dp

    call test_parallel_dispersion()
    call test_oblique_dispersion()
    call test_wave_number_scaling()
    call test_cold_plasma_limit()
    call test_perpendicular_dispersion()
    call test_unstable_drive()

    write (*, "(a)") "PASS"

contains

    subroutine test_parallel_dispersion()
        real(dp), parameter :: wave_vector(3) = [0.0_dp, 0.0_dp, 1.0_dp]
        real(dp) :: eigenvalues(3), eigenvectors(3, 3)
        real(dp) :: alfven_squared, expected(3)

        alfven_squared = magnetic_field(3)**2 / &
            (vacuum_permeability * density)
        expected = [alfven_squared, alfven_squared, &
            adiabatic_index * pressure / density]
        call solve_modes(wave_vector, adiabatic_index, eigenvalues, eigenvectors)
        call require(maxval(abs(eigenvalues - expected)) < 1.0e-11_dp, &
            "parallel Alfven and sound spectrum is wrong")
    end subroutine test_parallel_dispersion

    subroutine test_oblique_dispersion()
        real(dp), parameter :: wave_vector(3) = [0.6_dp, 0.0_dp, 0.8_dp]
        real(dp) :: eigenvalues(3), eigenvectors(3, 3), expected(3)
        real(dp) :: polarization

        expected = ideal_mhd_spectrum(wave_vector, adiabatic_index)
        call solve_modes(wave_vector, adiabatic_index, eigenvalues, eigenvectors)
        call require(maxval(abs(eigenvalues - expected)) < 1.0e-11_dp, &
            "oblique shear-Alfven and magnetosonic spectrum is wrong")
        polarization = abs(eigenvectors(2, 2)) / &
            sqrt(dot_product(eigenvectors(:, 2), eigenvectors(:, 2)))
        call require(abs(polarization - 1.0_dp) < 1.0e-12_dp, &
            "oblique shear-Alfven polarization is wrong")
        call require(max(abs(eigenvectors(2, 1)) / &
            norm2(eigenvectors(:, 1)), abs(eigenvectors(2, 3)) / &
            norm2(eigenvectors(:, 3))) < 1.0e-12_dp, &
            "magnetosonic polarization leaves the wave-field plane")
    end subroutine test_oblique_dispersion

    subroutine test_wave_number_scaling()
        real(dp), parameter :: wave_vector(3) = [0.6_dp, 0.0_dp, 0.8_dp]
        real(dp), parameter :: scale = 2.5_dp
        real(dp) :: base_values(3), scaled_values(3)
        real(dp) :: base_vectors(3, 3), scaled_vectors(3, 3)
        real(dp) :: overlap
        integer :: i

        call solve_modes(wave_vector, adiabatic_index, base_values, base_vectors)
        call solve_modes(scale * wave_vector, adiabatic_index, scaled_values, &
            scaled_vectors)
        call require(maxval(abs(scaled_values - scale**2 * base_values)) < &
            1.0e-11_dp, "frequency squared does not scale with wave number squared")
        do i = 1, 3
            overlap = abs(dot_product(base_vectors(:, i), scaled_vectors(:, i))) / &
                (norm2(base_vectors(:, i)) * norm2(scaled_vectors(:, i)))
            call require(abs(overlap - 1.0_dp) < 1.0e-12_dp, &
                "wave-number scaling changes a polarization")
        end do
    end subroutine test_wave_number_scaling

    subroutine test_cold_plasma_limit()
        real(dp), parameter :: wave_vector(3) = [0.6_dp, 0.0_dp, 0.8_dp]
        real(dp) :: eigenvalues(3), eigenvectors(3, 3), expected(3)

        expected = ideal_mhd_spectrum(wave_vector, 0.0_dp)
        call solve_modes(wave_vector, 0.0_dp, eigenvalues, eigenvectors)
        call require(maxval(abs(eigenvalues - expected)) < 1.0e-11_dp, &
            "cold-plasma zero, shear-Alfven and fast spectrum is wrong")
    end subroutine test_cold_plasma_limit

    subroutine test_perpendicular_dispersion()
        real(dp), parameter :: wave_vector(3) = [1.0_dp, 0.0_dp, 0.0_dp]
        real(dp) :: eigenvalues(3), eigenvectors(3, 3), expected(3)

        expected = ideal_mhd_spectrum(wave_vector, adiabatic_index)
        call solve_modes(wave_vector, adiabatic_index, eigenvalues, eigenvectors)
        call require(maxval(abs(eigenvalues - expected)) < 1.0e-11_dp, &
            "perpendicular double-zero and fast spectrum is wrong")
    end subroutine test_perpendicular_dispersion

    subroutine test_unstable_drive()
        real(dp), parameter :: wave_vector(3) = [0.0_dp, 0.0_dp, 1.0_dp]
        real(dp) :: stiffness(3, 3), mass(3, 3)
        real(dp) :: eigenvalues(3), eigenvectors(3, 3), drive
        real(dp) :: expected(3), alfven_squared, polarization
        integer :: info

        alfven_squared = magnetic_field(3)**2 / &
            (vacuum_permeability * density)
        drive = 2.0_dp * magnetic_field(3)**2 / vacuum_permeability
        expected = [-alfven_squared, alfven_squared, &
            adiabatic_index * pressure / density]
        call assemble_local_mode(wave_vector, magnetic_field, pressure, density, &
            adiabatic_index, drive, normal, stiffness, mass)
        call solve_three_component_modes(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "unstable eigenproblem failed")
        call require(maxval(abs(eigenvalues - expected)) < 1.0e-11_dp, &
            "pressure-curvature drive spectrum is wrong")
        polarization = abs(eigenvectors(1, 1)) / norm2(eigenvectors(:, 1))
        call require(abs(polarization - 1.0_dp) < 1.0e-12_dp, &
            "unstable pressure-curvature mode is not normal")
    end subroutine test_unstable_drive

    subroutine solve_modes(wave_vector, gamma_value, eigenvalues, eigenvectors)
        real(dp), intent(in) :: wave_vector(3), gamma_value
        real(dp), intent(out) :: eigenvalues(3), eigenvectors(3, 3)
        real(dp) :: stiffness(3, 3), mass(3, 3), residual(3)
        integer :: i, info

        call assemble_local_mode(wave_vector, magnetic_field, pressure, density, &
            gamma_value, 0.0_dp, normal, stiffness, mass)
        call require(maxval(abs(stiffness - transpose(stiffness))) < 1.0e-12_dp, &
            "stiffness is not symmetric")
        call solve_three_component_modes(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "stable eigenproblem failed")
        do i = 1, 3
            residual = matmul(stiffness, eigenvectors(:, i)) - &
                eigenvalues(i) * matmul(mass, eigenvectors(:, i))
            call require(maxval(abs(residual)) < 1.0e-11_dp, &
                "eigenpair residual is too large")
        end do
    end subroutine solve_modes

    pure function ideal_mhd_spectrum(wave_vector, gamma_value) result(spectrum)
        real(dp), intent(in) :: wave_vector(3), gamma_value
        real(dp) :: spectrum(3), alfven_squared, sound_squared
        real(dp) :: wave_squared, parallel_squared, sum, discriminant

        alfven_squared = dot_product(magnetic_field, magnetic_field) / &
            (vacuum_permeability * density)
        sound_squared = gamma_value * pressure / density
        wave_squared = dot_product(wave_vector, wave_vector)
        parallel_squared = dot_product(wave_vector, magnetic_field)**2 / &
            dot_product(magnetic_field, magnetic_field)
        sum = wave_squared * (alfven_squared + sound_squared)
        discriminant = sqrt(max(0.0_dp, sum**2 - 4.0_dp * wave_squared * &
            parallel_squared * alfven_squared * sound_squared))
        spectrum = [0.5_dp * (sum - discriminant), &
            parallel_squared * alfven_squared, 0.5_dp * (sum + discriminant)]
    end function ideal_mhd_spectrum

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_local_mode
