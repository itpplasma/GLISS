module enzyme_generated_geometry_fixture
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funloc, c_funptr
    use gliss_geometry_generated, only: generated_geometry_polynomial, &
        generated_geometry_field, generated_geometry_second_form
    implicit none

    real(c_double), parameter :: initial(3, 9) = reshape([ &
        1.1_c_double, 0.2_c_double, -0.3_c_double, &
        -0.4_c_double, 1.3_c_double, 0.2_c_double, &
        0.3_c_double, -0.1_c_double, 0.9_c_double, &
        0.2_c_double, -0.3_c_double, 0.1_c_double, &
        -0.1_c_double, 0.4_c_double, 0.3_c_double, &
        0.3_c_double, 0.2_c_double, -0.2_c_double, &
        -0.2_c_double, 0.3_c_double, 0.6_c_double, &
        0.4_c_double, -0.1_c_double, 0.2_c_double, &
        0.1_c_double, 0.5_c_double, -0.4_c_double], [3, 9])

    interface
        function enzyme_autodiff(function_pointer, alpha, chart, period, weight, axis) &
                result(gradient) bind(c, name="__enzyme_autodiff")
            import :: c_double, c_int, c_funptr
            type(c_funptr), value :: function_pointer
            real(c_double), value :: alpha
            integer(c_int), value :: chart, period, weight, axis
            real(c_double) :: gradient
        end function enzyme_autodiff
    end interface

contains

    subroutine case_inputs(family, chart, period, weight, axis, x, d, c, sign_j)
        integer, intent(in) :: family
        integer(c_int), intent(in) :: chart, period, weight, axis
        real(c_double), intent(out) :: x(27), d(27), c(13), sign_j
        real(c_double) :: jets(3, 9), polynomial(27), geometry(13)
        integer :: i, inputs, outputs

        jets = initial
        if (chart == 2) then
            jets(:, 2) = -jets(:, 2)
            jets(:, 5) = -jets(:, 5)
            jets(:, 8) = -jets(:, 8)
        end if
        polynomial = reshape(jets, [27])
        sign_j = real(3 - 2 * chart, c_double)
        x = 0.0_c_double
        select case (family)
        case (1)
            inputs = 27
            outputs = 13
            x = polynomial
        case (2)
            inputs = 6
            outputs = 5
            call generated_geometry_polynomial(polynomial, geometry)
            x(1:4) = geometry(4:7)
            x(5) = 0.67_c_double
            x(6) = -0.31_c_double
            x(7) = 1.0_c_double / real(period, c_double)
        case (3)
            inputs = 12
            outputs = 3
            call generated_geometry_polynomial(polynomial, geometry)
            x(1:3) = geometry(11:13)
            x(4:12) = reshape(jets(:, 7:9), [9])
        end select
        d = 0.0_c_double
        if (axis == 0) then
            do i = 1, inputs
                d(i) = sin(real(i + 3 * weight, c_double)) &
                    * (1.0_c_double + abs(x(i)))
            end do
        else
            d(axis) = 1.0_c_double
        end if
        c = 0.0_c_double
        do i = 1, outputs
            c(i) = cos(real(2 * i + weight, c_double))
        end do
    end subroutine case_inputs

    function polynomial_score(alpha, chart, period, weight, axis) bind(c) result(value)
        real(c_double), value :: alpha
        integer(c_int), value :: chart, period, weight, axis
        real(c_double) :: value, x(27), d(27), c(13), sign_j, y(13)

        call case_inputs(1, chart, period, weight, axis, x, d, c, sign_j)
        x = x + alpha * d
        call generated_geometry_polynomial(x(1:27), y)
        value = dot_product(c(1:13), y)
    end function polynomial_score

    function polynomial_gradient(chart, period, weight, axis) result(gradient)
        integer(c_int), intent(in) :: chart, period, weight, axis
        real(c_double) :: gradient

        gradient = enzyme_autodiff(c_funloc(polynomial_score), 0.0_c_double, &
            chart, period, weight, axis)
    end function polynomial_gradient

    function field_score(alpha, chart, period, weight, axis) bind(c) result(value)
        real(c_double), value :: alpha
        integer(c_int), value :: chart, period, weight, axis
        real(c_double) :: value, x(27), d(27), c(13), sign_j, y(5)

        call case_inputs(2, chart, period, weight, axis, x, d, c, sign_j)
        x = x + alpha * d
        call generated_geometry_field(x(1:7), y)
        value = dot_product(c(1:5), y)
    end function field_score

    function field_gradient(chart, period, weight, axis) result(gradient)
        integer(c_int), intent(in) :: chart, period, weight, axis
        real(c_double) :: gradient

        gradient = enzyme_autodiff(c_funloc(field_score), 0.0_c_double, &
            chart, period, weight, axis)
    end function field_gradient

    function second_form_score(alpha, chart, period, weight, axis) bind(c) result(value)
        real(c_double), value :: alpha
        integer(c_int), value :: chart, period, weight, axis
        real(c_double) :: value, x(27), d(27), c(13), sign_j, y(3)

        call case_inputs(3, chart, period, weight, axis, x, d, c, sign_j)
        x = x + alpha * d
        call generated_geometry_second_form(x(1:12), sign_j, y)
        value = dot_product(c(1:3), y)
    end function second_form_score

    function second_form_gradient(chart, period, weight, axis) result(gradient)
        integer(c_int), intent(in) :: chart, period, weight, axis
        real(c_double) :: gradient

        gradient = enzyme_autodiff(c_funloc(second_form_score), 0.0_c_double, &
            chart, period, weight, axis)
    end function second_form_gradient

end module enzyme_generated_geometry_fixture
