program test_two_component_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use two_component_kernel, only: two_component_components, &
        two_component_density
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: len = 6.0_dp * pi
    real(dp), parameter :: b1 = 0.3_dp, b3 = 0.4_dp, bz = 1.0_dp

    call check_reference(0.3_dp, 0.7_dp, -0.2_dp, 0.5_dp, -0.3_dp, &
        0.4_dp, 0.1_dp, 1.0822536130248883e-2_dp, &
        -4.9458926773973732e-1_dp, -2.3075713033433682_dp)
    call check_reference(0.5_dp, -0.4_dp, 0.3_dp, -0.25_dp, 0.2_dp, &
        -0.5_dp, 0.3_dp, -5.3051647697298445e-3_dp, &
        3.2939224238475452e-1_dp, 6.0633033498439228e-1_dp)
    call check_gradient()
    write (*, "(a)") "PASS"

contains

    pure function b_theta(radius) result(value)
        real(dp), intent(in) :: radius
        real(dp) :: value

        value = b1 * radius + b3 * radius**3
    end function b_theta

    subroutine kernel_inputs(radius, inputs)
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: inputs(12)
        real(dp) :: r

        r = radius
        inputs(1) = 2.0_dp * pi * r * bz
        inputs(2) = len * b_theta(r)
        inputs(3) = 2.0_dp * pi * bz
        inputs(4) = len * (b1 + 3.0_dp * b3 * r**2)
        inputs(5) = len * bz
        inputs(6) = 2.0_dp * pi * r * b_theta(r)
        inputs(7) = 2.0_dp * pi * len * r
        inputs(8) = sqrt(b_theta(r)**2 + bz**2)
        inputs(9) = 1.0_dp
        inputs(10) = (2.0_dp * b1 * r + 4.0_dp * b3 * r**3) * bz / r
        inputs(11) = -b_theta(r) * (2.0_dp * b1 * r &
            + 4.0_dp * b3 * r**3) / r
        inputs(12) = 0.0_dp
    end subroutine kernel_inputs

    subroutine check_reference(radius, xi_s, xi_s_s, xi_s_theta, &
            xi_s_zeta, eta_theta, eta_zeta, c1_ref, c2_ref, c3_ref)
        real(dp), intent(in) :: radius, xi_s, xi_s_s, xi_s_theta
        real(dp), intent(in) :: xi_s_zeta, eta_theta, eta_zeta
        real(dp), intent(in) :: c1_ref, c2_ref, c3_ref
        real(dp) :: inputs(12), c1, c2, c3

        call kernel_inputs(radius, inputs)
        call two_component_components(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta, &
            c1, c2, c3)
        call require(abs(c1 - c1_ref) < 1.0e-14_dp * max(1.0_dp, &
            abs(c1_ref)), "bending component is wrong")
        call require(abs(c2 - c2_ref) < 1.0e-14_dp * max(1.0_dp, &
            abs(c2_ref)), "shear component is wrong")
        call require(abs(c3 - c3_ref) < 1.0e-13_dp * max(1.0_dp, &
            abs(c3_ref)), "compression component is wrong")
    end subroutine check_reference

    subroutine check_gradient()
        real(dp) :: inputs(12), density_plus, density_minus, density
        real(dp) :: analytic, step, c1, c2, c3, drive

        call kernel_inputs(0.4_dp, inputs)
        drive = 2.0_dp * b_theta(0.4_dp) * (2.0_dp * b1 * 0.4_dp &
            + 4.0_dp * b3 * 0.4_dp**3) / 0.4_dp**2
        step = 1.0e-6_dp
        call two_component_density(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            drive, 0.6_dp + step, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, &
            -0.15_dp, density_plus)
        call two_component_density(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            drive, 0.6_dp - step, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, &
            -0.15_dp, density_minus)
        call two_component_density(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            drive, 0.6_dp, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, -0.15_dp, &
            density)
        call two_component_components(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), &
            inputs(9), inputs(10), inputs(11), inputs(12), inputs(12), &
            0.6_dp, 0.2_dp, -0.1_dp, 0.3_dp, 0.25_dp, -0.15_dp, &
            c1, c2, c3)
        analytic = (2.0_dp * c2 * shear_slope(inputs) &
            + 2.0_dp * c3 * compression_slope(inputs) &
            - 2.0_dp * drive * 0.6_dp) * abs(inputs(7))
        call require(abs((density_plus - density_minus) &
            / (2.0_dp * step) - analytic) < 1.0e-5_dp * max(1.0_dp, &
            abs(analytic)), "kernel gradient does not match")
    end subroutine check_gradient

    pure function shear_slope(inputs) result(value)
        real(dp), intent(in) :: inputs(12)
        real(dp) :: value

        value = -(sqrt(inputs(9)) / (inputs(8) * inputs(7))) &
            * (-(inputs(1) * inputs(4) - inputs(3) * inputs(2)) &
            + inputs(10) * inputs(7) / inputs(9))
    end function shear_slope

    pure function compression_slope(inputs) result(value)
        real(dp), intent(in) :: inputs(12)
        real(dp) :: value

        value = (1.0_dp / (inputs(8) * inputs(7))) &
            * (-(inputs(6) * inputs(4) + inputs(5) * inputs(3)) &
            - inputs(11) * inputs(7))
    end function compression_slope

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_two_component_kernel
