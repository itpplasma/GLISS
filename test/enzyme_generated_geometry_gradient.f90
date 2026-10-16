program enzyme_generated_geometry_gradient
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use, intrinsic :: iso_fortran_env, only: error_unit
    use enzyme_generated_geometry_fixture, only: initial, polynomial_gradient, &
        field_gradient, second_form_gradient
    use gliss_geometry_generated, only: generated_geometry_polynomial, &
        generated_geometry_polynomial_jvp, generated_geometry_polynomial_vjp, &
        generated_geometry_field, generated_geometry_field_jvp, &
        generated_geometry_field_vjp, generated_geometry_second_form, &
        generated_geometry_second_form_jvp, generated_geometry_second_form_vjp
    implicit none

    integer, parameter :: periods(3) = [1, 2, 5]
    real(c_double) :: jets(3, 9), polynomial(27), geometry(13)
    real(c_double) :: field(7), second(12)
    real(c_double) :: direction(27), cotangent(13), orientation
    integer :: chart, n, checked

    checked = 0
    do chart = 1, 2
        jets = initial
        if (chart == 2) then
            jets(:, 2) = -jets(:, 2)
            jets(:, 5) = -jets(:, 5)
            jets(:, 8) = -jets(:, 8)
        end if
        polynomial = reshape(jets, [27])
        call generated_geometry_polynomial(polynomial, geometry)
        orientation = sign(1.0_c_double, geometry(7))
        call require(orientation * real(3 - 2 * chart, c_double) > 0.0_c_double, &
            'both nonzero signed charts')
        call check_family(1, polynomial, 13, chart, 1, checked)
        do n = 1, size(periods)
            field(1:4) = geometry(4:7)
            field(5) = 0.67_c_double
            field(6) = -0.31_c_double
            field(7) = 1.0_c_double / real(periods(n), c_double)
            call check_family(2, field, 5, chart, periods(n), checked)
        end do
        second(1:3) = geometry(11:13)
        second(4:12) = reshape(jets(:, 7:9), [9])
        call check_family(3, second, 3, chart, 1, checked)
    end do
    write (*, '(a,i0,a)') 'PASS ', checked, &
        ' Enzyme geometry contractions; chart orientation and NFP frozen'

contains

    subroutine check_family(family, x, outputs, chart_index, period, count)
        integer, intent(in) :: family, outputs, chart_index, period
        real(c_double), intent(in) :: x(:)
        integer, intent(inout) :: count
        real(c_double) :: values(outputs), reverse_values(outputs), tangent(outputs)
        real(c_double) :: reverse(size(x)), reference(outputs), plus(outputs)
        real(c_double) :: minus(outputs), trial(size(x)), errors(3), step
        real(c_double) :: ad, forward, backward, centered
        integer :: weight, axis, i, level, inputs

        select case (family)
        case (1)
            call generated_geometry_polynomial(x, values)
        case (2)
            call generated_geometry_field(x, values)
        case (3)
            call generated_geometry_second_form(x, orientation, values)
        end select
        call geometric_oracle(family, x, reference)
        call near(values, reference, 8.0e-14_c_double, 'value leaf geometric oracle')
        inputs = size(x)
        ! The final field input represents a discrete, fixed field period count.
        if (family == 2) inputs = inputs - 1
        do weight = 1, 3
            cotangent = 0.0_c_double
            do i = 1, outputs
                cotangent(i) = cos(real(2 * i + weight, c_double))
            end do
            do axis = 0, inputs
                direction = 0.0_c_double
                if (axis == 0) then
                    do i = 1, inputs
                        direction(i) = sin(real(i + 3 * weight, c_double)) &
                            * (1.0_c_double + abs(x(i)))
                    end do
                else
                    direction(axis) = 1.0_c_double
                end if
                call generated_actions(family, x, values, tangent, &
                    reverse_values, reverse)
                call geometric_oracle(family, x, reference)
                call near(values, reference, 8.0e-14_c_double, 'geometric value')
                call near(reverse_values, values, 5.0e-14_c_double, &
                    'VJP primal value')
                forward = dot_product(cotangent(1:outputs), tangent)
                backward = dot_product(direction(1:size(x)), reverse)
                select case (family)
                case (1)
                    ad = polynomial_gradient(int(chart_index, c_int), &
                        int(period, c_int), int(weight, c_int), int(axis, c_int))
                case (2)
                    ad = field_gradient(int(chart_index, c_int), &
                        int(period, c_int), int(weight, c_int), int(axis, c_int))
                case (3)
                    ad = second_form_gradient(int(chart_index, c_int), &
                        int(period, c_int), int(weight, c_int), int(axis, c_int))
                end select
                call near([ad], [forward], 4.0e-12_c_double, 'Enzyme versus JVP')
                call near([ad], [backward], 4.0e-12_c_double, 'Enzyme versus VJP')
                do level = 1, 3
                    step = 2.0e-4_c_double / real(3**(level - 1), c_double)
                    trial = x + step * direction(1:size(x))
                    call geometric_oracle(family, trial, plus)
                    trial = x - step * direction(1:size(x))
                    call geometric_oracle(family, trial, minus)
                    centered = dot_product(cotangent(1:outputs), plus - minus) &
                        / (2.0_c_double * step)
                    errors(level) = abs(centered - ad) &
                        / max(1.0_c_double, abs(ad))
                end do
                call require(maxval(errors) < 2.0e-6_c_double &
                    .and. minval(errors) < 2.0e-8_c_double, &
                    'independent geometric finite difference')
                count = count + 1
            end do
        end do
    end subroutine check_family

    subroutine generated_actions(family, x, values, tangent, reverse_values, reverse)
        integer, intent(in) :: family
        real(c_double), intent(in) :: x(:)
        real(c_double), intent(out) :: values(:), tangent(:), reverse_values(:)
        real(c_double), intent(out) :: reverse(:)

        select case (family)
        case (1)
            call generated_geometry_polynomial_jvp(x, direction, values, tangent)
            call generated_geometry_polynomial_vjp(x, cotangent, &
                reverse_values, reverse)
        case (2)
            call generated_geometry_field_jvp(x, direction(1:7), values, tangent)
            call generated_geometry_field_vjp(x, cotangent(1:5), &
                reverse_values, reverse)
        case (3)
            call generated_geometry_second_form_jvp(x, direction(1:12), &
                orientation, values, tangent)
            call generated_geometry_second_form_vjp(x, cotangent(1:3), &
                orientation, reverse_values, reverse)
        end select
    end subroutine generated_actions

    subroutine geometric_oracle(family, x, y)
        integer, intent(in) :: family
        real(c_double), intent(in) :: x(:)
        real(c_double), intent(out) :: y(:)
        integer, parameter :: second_jet(3, 3) = &
            reshape([4, 5, 6, 5, 7, 8, 6, 8, 9], [3, 3])
        real(c_double) :: local_jets(3, 9), basis(3, 3), altered(3, 3)
        real(c_double) :: magnetic(3), normal(3)
        integer :: first, last, k, coordinate

        select case (family)
        case (1)
            local_jets = reshape(x, [3, 9])
            basis = local_jets(:, 1:3)
            k = 0
            do first = 1, 3
                do last = first, 3
                    k = k + 1
                    y(k) = dot_product(basis(:, first), basis(:, last))
                end do
            end do
            y(7) = determinant(basis)
            do coordinate = 1, 3
                y(7 + coordinate) = 0.0_c_double
                do k = 1, 3
                    altered = basis
                    altered(:, k) = local_jets(:, second_jet(k, coordinate))
                    y(7 + coordinate) = y(7 + coordinate) + determinant(altered)
                end do
            end do
            y(11:13) = cross_product(basis(:, 2), basis(:, 3))
        case (2)
            ! Cartesian Gram realization supplies an independent field norm.
            basis = 0.0_c_double
            basis(1, 1) = sqrt(x(1))
            basis(1, 2) = x(2) / basis(1, 1)
            basis(2, 2) = sqrt(x(3) - basis(1, 2)**2)
            y(1:2) = [-x(6) * x(7), -x(5)] / x(4)
            magnetic = y(1) * basis(:, 1) + y(2) * basis(:, 2)
            y(3) = dot_product(magnetic, basis(:, 1))
            y(4) = dot_product(magnetic, basis(:, 2))
            y(5) = sqrt(dot_product(magnetic, magnetic))
        case (3)
            normal = orientation * x(1:3) / sqrt(dot_product(x(1:3), x(1:3)))
            do k = 1, 3
                y(k) = dot_product(normal, x(3 * k + 1:3 * k + 3))
            end do
        end select
    end subroutine geometric_oracle

    function cross_product(a, b) result(c)
        real(c_double), intent(in) :: a(3), b(3)
        real(c_double) :: c(3)

        c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), &
            a(1) * b(2) - a(2) * b(1)]
    end function cross_product

    function determinant(basis) result(value)
        real(c_double), intent(in) :: basis(3, 3)
        real(c_double) :: value

        value = dot_product(basis(:, 1), cross_product(basis(:, 2), basis(:, 3)))
    end function determinant

    subroutine near(actual, expected, tolerance, label)
        real(c_double), intent(in) :: actual(:), expected(:), tolerance
        character(*), intent(in) :: label

        call require(all(ieee_is_finite(actual)), label // ' finite')
        call require(all(ieee_is_finite(expected)), label // ' oracle finite')
        if (maxval(abs(actual - expected) / max(1.0_c_double, abs(expected))) &
                > tolerance) then
            write (error_unit, '(a,1x,a,1x,es24.16)') 'FAIL', label, &
                maxval(abs(actual - expected))
            error stop 1
        end if
    end subroutine near

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (.not. condition) then
            write (error_unit, '(a,1x,a)') 'FAIL', label
            error stop 1
        end if
    end subroutine require

end program enzyme_generated_geometry_gradient
