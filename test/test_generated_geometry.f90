program test_generated_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gliss_geometry_generated, only: generated_geometry_polynomial, &
        generated_geometry_polynomial_jvp, generated_geometry_polynomial_vjp, &
        generated_geometry_field, generated_geometry_field_jvp, &
        generated_geometry_field_vjp, generated_geometry_second_form, &
        generated_geometry_second_form_jvp, generated_geometry_second_form_vjp
    implicit none

    real(dp), parameter :: initial(3, 9) = reshape([ &
        1.2_dp, -0.4_dp, 0.7_dp, 0.3_dp, 1.1_dp, -0.2_dp, &
        -0.6_dp, 0.5_dp, 1.4_dp, 0.2_dp, -0.1_dp, 0.3_dp, &
        -0.4_dp, 0.3_dp, 0.2_dp, 0.1_dp, 0.2_dp, -0.3_dp, &
        0.6_dp, -0.2_dp, 0.1_dp, -0.1_dp, 0.4_dp, 0.2_dp, &
        0.3_dp, 0.1_dp, -0.5_dp], [3, 9])
    integer, parameter :: periods(3) = [1, 2, 5]
    real(dp) :: jets(3, 9), polynomial(27), field(7), second(12)
    real(dp) :: geometry(13), field_value(5), curvature(3), orientation
    real(dp) :: reference(13), changed(13), direction(27), tangent(13)
    real(dp) :: curvature_direction(12), base_field(5)
    real(dp) :: scale, expected(13), q(3), negative_q(3), coordinate_error(3), step
    integer :: chart, n, axis, level

    do chart = 1, 2
        jets = initial
        if (chart == 2) then
            ! Reverse theta in the polynomial map, including its mixed jets.
            jets(:, 2) = -jets(:, 2)
            jets(:, 5) = -jets(:, 5)
            jets(:, 8) = -jets(:, 8)
        end if
        polynomial = reshape(jets, [27])
        call generated_geometry_polynomial(polynomial, geometry)
        call oracle(1, polynomial, 1.0_dp, reference)
        call near(geometry, reference, 3.0e-14_dp, 'polynomial geometric oracle')
        orientation = sign(1.0_dp, geometry(7))
        call require(orientation * real(3 - 2 * chart, dp) > 0.0_dp, 'signed charts')
        call check_actions(1, polynomial, orientation, 13)
        do axis = 1, 3
            do level = 1, 3
                step = 2.0e-4_dp / real(3**(level - 1), dp)
                q = 0.0_dp
                q(axis) = step
                negative_q = -q
                coordinate_error(level) = abs((map_jacobian(jets, q) &
                    - map_jacobian(jets, negative_q)) / (2.0_dp * step) &
                    - geometry(7 + axis)) / max(1.0_dp, abs(geometry(7 + axis)))
            end do
            call require(maxval(coordinate_error) < 2.0e-7_dp &
                .and. minval(coordinate_error) < 2.0e-9_dp, &
                'polynomial coordinate Jacobian derivative')
        end do
        scale = 2.3_dp
        direction = scale * polynomial
        call generated_geometry_polynomial(direction, changed)
        expected(1:6) = scale**2 * geometry(1:6)
        expected(7:10) = scale**3 * geometry(7:10)
        expected(11:13) = scale**2 * geometry(11:13)
        call near(changed, expected, 5.0e-14_dp, 'geometry length units')
        direction = polynomial
        call generated_geometry_polynomial_jvp(polynomial, direction, changed, tangent)
        expected(1:6) = 2.0_dp * geometry(1:6)
        expected(7:10) = 3.0_dp * geometry(7:10)
        expected(11:13) = 2.0_dp * geometry(11:13)
        call near(tangent, expected, 5.0e-14_dp, 'geometry dilation tangent')
        base_field = 0.0_dp
        do n = 1, size(periods)
            field(1:4) = geometry(4:7)
            field(5) = 0.73_dp
            field(6) = -0.28_dp
            field(7) = 1.0_dp / real(periods(n), dp)
            call generated_geometry_field(field, field_value)
            call oracle(2, field, orientation, reference(1:5))
            call near(field_value, reference(1:5), 5.0e-14_dp, 'Cartesian field')
            if (n == 1) base_field = field_value
            call require(abs(field_value(1) * real(periods(n), dp) &
                - base_field(1)) < 3.0e-14_dp, 'field-period flux factor')
            call check_actions(2, field, orientation, 5)
            call check_field_scaling(field, field_value)
            call check_single_flux(field, orientation)
        end do
        second(1:3) = geometry(11:13)
        second(4:12) = reshape(jets(:, 7:9), [9])
        call generated_geometry_second_form(second, orientation, curvature)
        call oracle(3, second, orientation, reference(1:3))
        call near(curvature, reference(1:3), 5.0e-14_dp, 'outward second form')
        call check_actions(3, second, orientation, 3)
        curvature_direction = 0.0_dp
        curvature_direction(1:3) = second(1:3)
        call generated_geometry_second_form_jvp(second, curvature_direction, &
            orientation, changed(1:3), tangent(1:3))
        call near(tangent(1:3), [0.0_dp, 0.0_dp, 0.0_dp], &
            3.0e-14_dp, 'normal magnitude tangent')
        curvature_direction = second
        curvature_direction(1:3) = 2.0_dp * second(1:3)
        call generated_geometry_second_form_jvp(second, curvature_direction, &
            orientation, changed(1:3), tangent(1:3))
        call near(tangent(1:3), curvature, 5.0e-14_dp, 'second form length units')
        second(4:6) = second(4:6) + 0.4_dp * jets(:, 2) - 0.7_dp * jets(:, 3)
        second(7:9) = second(7:9) - 0.3_dp * jets(:, 2) + 0.2_dp * jets(:, 3)
        second(10:12) = second(10:12) + 0.8_dp * jets(:, 2) + 0.6_dp * jets(:, 3)
        call generated_geometry_second_form(second, orientation, changed(1:3))
        call near(changed(1:3), curvature, 5.0e-14_dp, 'tangential second jets')
    end do
    call check_field_dynamic_range()
    call check_curvature_dynamic_range()
    write (*, '(a)') 'PASS'

contains

    subroutine check_actions(family, x, sign_j, outputs)
        integer, intent(in) :: family, outputs
        real(dp), intent(in), contiguous :: x(:)
        real(dp), intent(in) :: sign_j
        real(dp) :: d(size(x)), c(outputs), values(outputs), forward(outputs)
        real(dp) :: reverse(size(x)), plus(outputs), minus(outputs)
        real(dp) :: trial(size(x)), fd(outputs), errors(3), scalar_errors(3), h
        integer :: k, i, level

        do k = 1, 3
            do i = 1, size(x)
                d(i) = sin(real(i + 2 * k, dp)) * (1.0_dp + abs(x(i)))
            end do
            do i = 1, outputs
                c(i) = cos(real(2 * i + k, dp))
            end do
            call actions(family, x, d, c, sign_j, values, forward, reverse)
            call oracle(family, x, sign_j, plus)
            call near(values, plus, 8.0e-14_dp, 'derivative primal output')
            call require(abs(dot_product(c, forward) - dot_product(d, reverse)) &
                < 2.0e-12_dp * max(1.0_dp, abs(dot_product(c, forward))), &
                'tangent/reverse duality')
            do level = 1, 3
                h = 2.0e-4_dp / real(3**(level - 1), dp)
                trial = x + h * d
                call oracle(family, trial, sign_j, plus)
                trial = x - h * d
                call oracle(family, trial, sign_j, minus)
                fd = (plus - minus) / (2.0_dp * h)
                errors(level) = maxval(abs(fd - forward) / max(1.0_dp, abs(forward)))
            end do
            call require(maxval(errors) < 2.0e-6_dp .and. minval(errors) < 2.0e-8_dp, &
                'independent three-step JVP')
            do i = 1, size(x)
                do level = 1, 3
                    h = 2.0e-4_dp * max(1.0_dp, abs(x(i))) &
                        / real(3**(level - 1), dp)
                    trial = x
                    trial(i) = x(i) + h
                    call oracle(family, trial, sign_j, plus)
                    trial(i) = x(i) - h
                    call oracle(family, trial, sign_j, minus)
                    scalar_errors(level) = abs(dot_product(c, plus - minus) &
                        / (2.0_dp * h) - reverse(i)) / max(1.0_dp, abs(reverse(i)))
                end do
                call require(maxval(scalar_errors) < 2.0e-6_dp &
                    .and. minval(scalar_errors) < 2.0e-8_dp, &
                    'independent three-step scalar VJP')
            end do
        end do
    end subroutine check_actions

    subroutine actions(family, x, d, c, sign_j, values, forward, reverse)
        integer, intent(in) :: family
        real(dp), intent(in), contiguous :: x(:), d(:), c(:)
        real(dp), intent(in) :: sign_j
        real(dp), intent(out), contiguous :: values(:), forward(:), reverse(:)
        real(dp) :: reverse_values(size(values))

        select case (family)
        case (1)
            call generated_geometry_polynomial_jvp(x, d, values, forward)
            call generated_geometry_polynomial_vjp(x, c, reverse_values, reverse)
        case (2)
            call generated_geometry_field_jvp(x, d, values, forward)
            call generated_geometry_field_vjp(x, c, reverse_values, reverse)
        case (3)
            call generated_geometry_second_form_jvp(x, d, sign_j, values, forward)
            call generated_geometry_second_form_vjp(x, c, sign_j, &
                reverse_values, reverse)
        end select
        call near(values, reverse_values, 3.0e-14_dp, 'forward/reverse primal outputs')
    end subroutine actions

    subroutine oracle(family, x, sign_j, y)
        integer, intent(in) :: family
        real(dp), intent(in) :: x(:), sign_j
        real(dp), intent(out) :: y(:)
        real(dp) :: basis(3, 3), altered(3, 3), jets_local(3, 9)
        real(dp) :: field_vector(3), normal(3), cross(3)
        integer, parameter :: derivative_jet(3, 3) = &
            reshape([4, 5, 6, 5, 7, 8, 6, 8, 9], [3, 3])
        integer :: first, last, k, axis

        select case (family)
        case (1)
            jets_local = reshape(x, [3, 9])
            basis = jets_local(:, 1:3)
            k = 0
            do first = 1, 3
                do last = first, 3
                    k = k + 1
                    y(k) = dot_product(basis(:, first), basis(:, last))
                end do
            end do
            y(7) = determinant(basis)
            do axis = 1, 3
                y(7 + axis) = 0.0_dp
                do k = 1, 3
                    altered = basis
                    altered(:, k) = jets_local(:, derivative_jet(k, axis))
                    y(7 + axis) = y(7 + axis) + determinant(altered)
                end do
            end do
            cross(1) = basis(2, 2) * basis(3, 3) - basis(3, 2) * basis(2, 3)
            cross(2) = basis(3, 2) * basis(1, 3) - basis(1, 2) * basis(3, 3)
            cross(3) = basis(1, 2) * basis(2, 3) - basis(2, 2) * basis(1, 3)
            y(11:13) = cross
        case (2)
            ! Realize the angular metric as Cartesian vectors, then form B.
            basis = 0.0_dp
            basis(1, 1) = sqrt(x(1))
            basis(1, 2) = x(2) / sqrt(x(1))
            basis(2, 2) = sqrt(x(3) - x(2)**2 / x(1))
            y(1) = -x(6) * x(7) / x(4)
            y(2) = -x(5) / x(4)
            field_vector = y(1) * basis(:, 1) + y(2) * basis(:, 2)
            y(3) = dot_product(field_vector, basis(:, 1))
            y(4) = dot_product(field_vector, basis(:, 2))
            y(5) = sqrt(dot_product(field_vector, field_vector))
        case (3)
            normal = sign_j * x(1:3) / sqrt(sum(x(1:3)**2))
            do k = 1, 3
                y(k) = dot_product(normal, x(3 * k + 1:3 * k + 3))
            end do
        end select
    end subroutine oracle

    function map_jacobian(jets_local, q) result(value)
        real(dp), intent(in) :: jets_local(3, 9), q(3)
        real(dp) :: value, basis(3, 3)

        basis(:, 1) = jets_local(:, 1) + jets_local(:, 4) * q(1) &
            + jets_local(:, 5) * q(2) + jets_local(:, 6) * q(3)
        basis(:, 2) = jets_local(:, 2) + jets_local(:, 5) * q(1) &
            + jets_local(:, 7) * q(2) + jets_local(:, 8) * q(3)
        basis(:, 3) = jets_local(:, 3) + jets_local(:, 6) * q(1) &
            + jets_local(:, 8) * q(2) + jets_local(:, 9) * q(3)
        value = determinant(basis)
    end function map_jacobian

    function determinant(b) result(value)
        real(dp), intent(in) :: b(3, 3)
        real(dp) :: value

        value = b(1, 1) * b(2, 2) * b(3, 3) + b(1, 2) * b(2, 3) * b(3, 1) &
            + b(1, 3) * b(2, 1) * b(3, 2) - b(1, 3) * b(2, 2) * b(3, 1) &
            - b(1, 2) * b(2, 1) * b(3, 3) - b(1, 1) * b(2, 3) * b(3, 2)
    end function determinant

    subroutine check_field_scaling(x, values)
        real(dp), intent(in) :: x(7), values(5)
        real(dp) :: trial(7), d(7), actual(5), tangent_field(5), target(5)
        real(dp), parameter :: a = 2.3_dp

        trial = x
        trial(1:3) = a**2 * x(1:3)
        trial(4) = a**3 * x(4)
        call generated_geometry_field(trial, actual)
        target(1:2) = values(1:2) / a**3
        target(3:4) = values(3:4) / a
        target(5) = values(5) / a**2
        call near(actual, target, 5.0e-14_dp, 'field length units')
        d = 0.0_dp
        d(1:3) = 2.0_dp * x(1:3)
        d(4) = 3.0_dp * x(4)
        call generated_geometry_field_jvp(x, d, actual, tangent_field)
        target(1:2) = -3.0_dp * values(1:2)
        target(3:4) = -values(3:4)
        target(5) = -2.0_dp * values(5)
        call near(tangent_field, target, 8.0e-14_dp, 'field dilation tangent')
        trial = x
        trial(5:6) = -2.0_dp * x(5:6)
        call generated_geometry_field(trial, actual)
        target(1:4) = -2.0_dp * values(1:4)
        target(5) = 2.0_dp * values(5)
        call near(actual, target, 5.0e-14_dp, 'signed flux units')
    end subroutine check_field_scaling

    subroutine check_single_flux(x, sign_j)
        real(dp), intent(in) :: x(7), sign_j
        real(dp) :: trial(7), actual(5), target(5)
        integer :: flux

        do flux = 5, 6
            trial = x
            trial(flux) = 0.0_dp
            call generated_geometry_field(trial, actual)
            call oracle(2, trial, sign_j, target)
            call near(actual, target, 5.0e-14_dp, 'individual zero flux')
            call check_actions(2, trial, sign_j, 5)
        end do
    end subroutine check_single_flux

    subroutine check_field_dynamic_range()
        real(dp) :: x(7)
        integer :: exponent, chart

        do exponent = -200, 200, 400
            x(1:3) = [2.0_dp, 0.3_dp, 1.5_dp]
            x(5) = 10.0_dp**real(exponent / 2, dp)
            x(6) = -2.0_dp * x(5)
            x(7) = 0.5_dp
            do chart = 1, 2
                x(4) = real(3 - 2 * chart, dp) * 10.0_dp**real(exponent, dp)
                call check_field_range_point(x)
            end do
        end do
        ! Cartesian basis scaled by 1e100: g~1e200, J~1e300, flux~1e200.
        ! The physical B is finite although metric*flux overflows.
        x(1:3) = [1.0e200_dp, 0.0_dp, 1.0e200_dp]
        x(5) = 1.0e200_dp
        x(6) = -2.0e200_dp
        x(7) = 0.5_dp
        do chart = 1, 2
            x(4) = real(3 - 2 * chart, dp) * 1.0e300_dp
            call check_field_range_point(x)
        end do
    end subroutine check_field_dynamic_range

    subroutine check_field_range_point(x)
        real(dp), intent(in) :: x(7)
        real(dp), parameter :: cotangent(5) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        real(dp) :: actual(5), target(5), d(7), forward(5), reverse(7)

        call generated_geometry_field(x, actual)
        call oracle(2, x, 1.0_dp, target)
        call require(actual(5) > 0.0_dp, 'nonzero field dynamic range')
        call require(all(abs(actual / target - 1.0_dp) < 5.0e-14_dp), &
            'signed field dynamic range relative oracle')
        d = 0.0_dp
        d(4) = x(4)
        call generated_geometry_field_jvp(x, d, actual, forward)
        call require(all(abs(forward / target + 1.0_dp) < 8.0e-14_dp), &
            'field dynamic range Jacobian tangent')
        call generated_geometry_field_vjp(x, cotangent, actual, reverse)
        call require(abs(reverse(4) * x(4) / target(5) + 1.0_dp) &
            < 8.0e-14_dp, 'field dynamic range magnitude reverse')
    end subroutine check_field_range_point

    subroutine check_curvature_dynamic_range()
        real(dp), parameter :: cotangent(3) = [1.0_dp, 0.0_dp, 0.0_dp]
        real(dp) :: x(12), d(12), actual(3), forward(3), reverse(12), sign_j, normal(3)
        integer :: chart

        x = 0.0_dp
        x(1) = 1.0e150_dp
        x(4:5) = 1.0e200_dp
        d = 0.0_dp
        d(2) = 1.0e100_dp
        do chart = 1, 2
            sign_j = real(3 - 2 * chart, dp)
            call generated_geometry_second_form(x, sign_j, actual)
            call require(abs(actual(1) / (sign_j * 1.0e200_dp) - 1.0_dp) &
                < 5.0e-14_dp, 'second form dynamic range value')
            call generated_geometry_second_form_jvp(x, d, sign_j, actual, forward)
            call require(abs(forward(1) / (sign_j * 1.0e150_dp) - 1.0_dp) &
                < 8.0e-14_dp, 'second form dynamic range tangent')
            call generated_geometry_second_form_vjp(x, cotangent, sign_j, &
                actual, reverse)
            call require(abs(reverse(2) / (sign_j * 1.0e50_dp) - 1.0_dp) &
                < 8.0e-14_dp, 'second form dynamic range cross pullback')
            normal = 0.0_dp
            normal(1) = sign_j
            call near(reverse(4:6), normal, &
                8.0e-14_dp, 'second form dynamic range second-jet pullback')
        end do
    end subroutine check_curvature_dynamic_range

    subroutine near(actual, target, tolerance, label)
        real(dp), intent(in) :: actual(:), target(:), tolerance
        character(*), intent(in) :: label

        call require(all(abs(actual - target) <= tolerance &
            * max(1.0_dp, abs(target))), label)
    end subroutine near

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        write (error_unit, '(a)') 'FAIL: ' // label
        error stop 1
    end subroutine require
end program test_generated_geometry
