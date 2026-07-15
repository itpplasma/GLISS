program test_newcomb_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use newcomb_limit, only: assemble_single_mode_stiffness, &
        cylinder_profiles_t, lowest_artificial_stiffness_level, &
        newcomb_axis_regular_power, newcomb_invalid_input, newcomb_ok, &
        newcomb_quadratic_coefficients
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    type(cylinder_profiles_t) :: profiles

    profiles%length = 6.0_dp * pi
    profiles%b_axial = 1.0_dp
    profiles%b_linear = 0.3_dp
    profiles%b_cubic = 0.4_dp

    call test_mode_dependent_axis_space()
    call test_published_line_bending_coefficient()
    call test_second_order_radial_coefficient()
    call test_full_operator_radial_convergence()
    call test_production_axis_asymptotics()
    call test_invalid_magnetic_field()
    call test_stability_signs()

    write (*, "(a)") "PASS"

contains

    subroutine test_mode_dependent_axis_space()
        real(dp), allocatable :: unit_mode(:, :), higher_mode(:, :)
        real(dp) :: a, b, c, expected, step
        integer :: info

        call assemble_single_mode_stiffness(profiles, 1, 1, 0.5_dp, 8, &
            unit_mode, info)
        call require(info == newcomb_ok, "unit-harmonic assembly failed")
        call assemble_single_mode_stiffness(profiles, 2, 1, 0.5_dp, 8, &
            higher_mode, info)
        call require(info == newcomb_ok, "higher-harmonic assembly failed")
        call require(size(unit_mode, 1) == 8, &
            "unit harmonic does not retain its regular axis value")
        call require(size(higher_mode, 1) == 7, &
            "higher harmonic does not vanish at the axis")
        call require(newcomb_axis_regular_power(0) == 1, &
            "axisymmetric regular power is wrong")
        call require(newcomb_axis_regular_power(1) == 0, &
            "unit-harmonic regular power is wrong")
        call require(newcomb_axis_regular_power(4) == 3, &
            "higher-harmonic regular power is wrong")
        call require(maxval(abs(unit_mode - transpose(unit_mode))) < 1.0e-13_dp, &
            "unit-harmonic stiffness is not symmetric")

        step = 0.5_dp / 8.0_dp
        call newcomb_quadratic_coefficients(profiles, 1, 1, 0.5_dp * step, &
            a, b, c, info)
        call require(info == newcomb_ok, "axis element coefficients failed")
        expected = step * (0.25_dp * a - b / step + c / step**2)
        call require(abs(unit_mode(1, 1) - expected) < 1.0e-13_dp &
            * max(1.0_dp, abs(expected)), &
            "unit-harmonic natural axis element is wrong")

        call assemble_single_mode_stiffness(profiles, 1, 1, 0.0_dp, 8, &
            unit_mode, info)
        call require(info == newcomb_invalid_input, &
            "zero cylinder radius was accepted")
    end subroutine test_mode_dependent_axis_space

    subroutine test_published_line_bending_coefficient()
        call check_line_bending_coefficient(2, 1)
        call check_line_bending_coefficient(-2, 1)
        call check_line_bending_coefficient(2, -3)
    end subroutine test_published_line_bending_coefficient

    subroutine check_line_bending_coefficient(mode_m, mode_n)
        integer, intent(in) :: mode_m, mode_n
        real(dp), parameter :: radius = 0.2_dp
        real(dp) :: a, b, c, b_theta, field_line_bending
        real(dp) :: k_parallel, k_squared, expected_f
        integer :: info

        call newcomb_quadratic_coefficients(profiles, mode_m, mode_n, radius, &
            a, b, c, info)
        call require(info == newcomb_ok, "Newcomb coefficient reduction failed")
        b_theta = profiles%b_linear * radius + profiles%b_cubic * radius**3
        k_parallel = two_pi() * real(mode_n, dp) / profiles%length
        field_line_bending = -real(mode_m, dp) * b_theta / radius &
            + k_parallel * profiles%b_axial
        k_squared = k_parallel**2 + real(mode_m, dp)**2 / radius**2
        expected_f = radius * field_line_bending**2 / k_squared
        call require(abs(c / (pi * profiles%length) - expected_f) &
            < 1.0e-12_dp * max(1.0_dp, abs(expected_f)), &
            "radial Newcomb coefficient disagrees with the published form")
    end subroutine check_line_bending_coefficient

    subroutine test_second_order_radial_coefficient()
        integer, parameter :: meshes(3) = [16, 32, 64]
        type(cylinder_profiles_t) :: theta_pinch
        real(dp) :: errors(3), exact, approximate
        real(dp) :: a, b, c, left, right, midpoint, slope, step
        integer :: i, info, level

        theta_pinch%length = two_pi()
        theta_pinch%b_axial = 1.0_dp
        exact = 4.0_dp * pi**2 * (0.5_dp - 1.0_dp + log(2.0_dp))
        do level = 1, size(meshes)
            step = 1.0_dp / real(meshes(level), dp)
            approximate = 0.0_dp
            do i = 1, meshes(level)
                left = real(i - 1, dp) * step
                right = real(i, dp) * step
                midpoint = 0.5_dp * (left + right)
                slope = ((1.0_dp - right**2) - (1.0_dp - left**2)) / step
                call newcomb_quadratic_coefficients(theta_pinch, 1, 1, &
                    midpoint, a, b, c, info)
                call require(info == newcomb_ok, &
                    "theta-pinch coefficient reduction failed")
                approximate = approximate + step * c * slope**2
            end do
            errors(level) = abs(approximate - exact)
        end do
        call require(all(errors(2:) < errors(:2)), &
            "radial coefficient does not converge monotonically")
        call require(all(errors(:2) / errors(2:) > 3.8_dp), &
            "radial coefficient is not second-order convergent")
    end subroutine test_second_order_radial_coefficient

    subroutine test_full_operator_radial_convergence()
        integer, parameter :: meshes(3) = [16, 32, 64]
        type(cylinder_profiles_t) :: theta_pinch
        real(dp), allocatable :: stiffness(:, :), displacement(:), image(:)
        real(dp) :: errors(3), reference, radius, step
        integer :: column, i, info, level

        theta_pinch%length = two_pi()
        theta_pinch%b_axial = 1.0_dp
        reference = reference_full_energy(theta_pinch, 16384)
        do level = 1, size(meshes)
            call assemble_single_mode_stiffness(theta_pinch, 1, 1, 1.0_dp, &
                meshes(level), stiffness, info)
            call require(info == newcomb_ok, "full theta-pinch assembly failed")
            allocate (displacement(size(stiffness, 1)), &
                image(size(stiffness, 1)))
            step = 1.0_dp / real(meshes(level), dp)
            do i = 1, size(displacement)
                radius = real(i - 1, dp) * step
                displacement(i) = 1.0_dp - radius**2
            end do
            do i = 1, size(image)
                image(i) = 0.0_dp
                do column = 1, size(displacement)
                    image(i) = image(i) &
                        + stiffness(i, column) * displacement(column)
                end do
            end do
            errors(level) = abs(dot_product(displacement, image) - reference)
            deallocate (stiffness, displacement, image)
        end do
        call require(all(errors(2:) < errors(:2)), &
            "full Newcomb operator does not converge monotonically")
        call require(all(errors(:2) / errors(2:) > 3.5_dp), &
            "full Newcomb operator is not second-order convergent")
    end subroutine test_full_operator_radial_convergence

    function reference_full_energy(theta_pinch, intervals) result(energy)
        type(cylinder_profiles_t), intent(in) :: theta_pinch
        integer, intent(in) :: intervals
        real(dp) :: energy, a, b, c, radius, step, value, slope
        integer :: i, info

        energy = 0.0_dp
        step = 1.0_dp / real(intervals, dp)
        do i = 1, intervals
            radius = (real(i, dp) - 0.5_dp) * step
            value = 1.0_dp - radius**2
            slope = -2.0_dp * radius
            call newcomb_quadratic_coefficients(theta_pinch, 1, 1, radius, &
                a, b, c, info)
            call require(info == newcomb_ok, "reference coefficient failed")
            energy = energy + step * (a * value**2 + 2.0_dp * b * value &
                * slope + c * slope**2)
        end do
    end function reference_full_energy

    subroutine test_production_axis_asymptotics()
        call check_axis_ratio(0, 1, 1.0_dp)
        call check_axis_ratio(1, 1, 0.0_dp)
        call check_axis_ratio(3, 1, 8.0_dp)
    end subroutine test_production_axis_asymptotics

    subroutine check_axis_ratio(mode_m, mode_n, expected)
        integer, intent(in) :: mode_m, mode_n
        real(dp), intent(in) :: expected
        real(dp), parameter :: radius = 1.0e-3_dp, delta = 1.0e-6_dp
        real(dp) :: a, b, c, b_left, b_right, dummy_a, dummy_c, ratio
        real(dp) :: polynomial(2), roots(2)
        integer :: info, magnitude

        call newcomb_quadratic_coefficients(profiles, mode_m, mode_n, radius, &
            a, b, c, info)
        call require(info == newcomb_ok, "axis coefficient reduction failed")
        call newcomb_quadratic_coefficients(profiles, mode_m, mode_n, &
            radius - delta, dummy_a, b_left, dummy_c, info)
        call require(info == newcomb_ok, "left axis coefficient failed")
        call newcomb_quadratic_coefficients(profiles, mode_m, mode_n, &
            radius + delta, dummy_a, b_right, dummy_c, info)
        call require(info == newcomb_ok, "right axis coefficient failed")
        ratio = radius**2 * (a - (b_right - b_left) / (2.0_dp * delta)) / c
        call require(abs(ratio - expected) < 2.0e-3_dp, &
            "production Newcomb axis ratio is wrong")
        magnitude = abs(mode_m)
        if (magnitude == 0) then
            roots(1) = -1.0_dp
            roots(2) = 1.0_dp
            polynomial = roots**2 - ratio
        else
            roots(1) = -1.0_dp - real(magnitude, dp)
            roots(2) = -1.0_dp + real(magnitude, dp)
            polynomial = roots * (roots + 2.0_dp) - ratio
        end if
        call require(maxval(abs(polynomial)) < 2.0e-3_dp, &
            "production Newcomb indicial polynomial is wrong")
    end subroutine check_axis_ratio

    subroutine test_invalid_magnetic_field()
        type(cylinder_profiles_t) :: zero_field
        real(dp) :: a, b, c, level
        integer :: info

        zero_field%length = 1.0_dp
        call newcomb_quadratic_coefficients(zero_field, 1, 1, 0.2_dp, &
            a, b, c, info)
        call require(info == newcomb_invalid_input, &
            "zero magnetic field was accepted")
        call lowest_artificial_stiffness_level(zero_field, 1, 1, 0.5_dp, &
            8, level, info)
        call require(info == newcomb_invalid_input .and. level == 0.0_dp, &
            "failed solve returned a physical stiffness level")
    end subroutine test_invalid_magnetic_field

    subroutine test_stability_signs()
        real(dp) :: unstable_coarse, unstable_fine, stable_mode
        integer :: info

        call lowest_artificial_stiffness_level(profiles, 1, 1, 0.5_dp, 100, &
            unstable_coarse, info)
        call require(info == newcomb_ok, "coarse unstable solve failed")
        call lowest_artificial_stiffness_level(profiles, 1, 1, 0.5_dp, 200, &
            unstable_fine, info)
        call require(info == newcomb_ok, "fine unstable solve failed")
        call lowest_artificial_stiffness_level(profiles, 2, 1, 0.5_dp, 100, &
            stable_mode, info)
        call require(info == newcomb_ok, "stable solve failed")
        call require(unstable_coarse < 0.0_dp, &
            "resonant restricted space is not unstable")
        call require(unstable_fine <= unstable_coarse, &
            "instability does not deepen under refinement")
        call require(stable_mode > 0.0_dp, "non-resonant mode is not stable")
    end subroutine test_stability_signs

    pure function two_pi() result(value)
        real(dp) :: value

        value = 2.0_dp * pi
    end function two_pi

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_newcomb_limit
