program test_primitive_cartesian_field
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use physical_constants, only: vacuum_permeability
    use primitive_cartesian_field, only: primitive_cartesian_field_point_t, &
        build_primitive_cartesian_field_point, primitive_cartesian_field_ok, &
        primitive_cartesian_field_invalid
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), length = 1.3_dp
    integer, parameter :: periods(2) = [1, 5]
    real(dp), parameter :: points(3, 2) = reshape([ &
        0.37_dp, 0.29_dp, -0.41_dp, 0.63_dp, -0.21_dp, 0.17_dp], [3, 2])
    type(primitive_cartesian_field_point_t) :: actual, plus, minus
    real(dp) :: p(3), trial(3), jets(3, 9), b(3), gradient(3, 3), curl(3)
    real(dp) :: fd(3, 3), errors(3), fd_maximum(3), sign_j, h, phi, s, theta, zeta
    integer :: chart, family, n, point, axis, level, info

    fd_maximum = 0.0_dp
    do family = 1, 2
        do chart = 1, 2
            sign_j = real(3 - 2 * chart, dp)
            do n = 1, size(periods)
                do point = 1, size(points, 2)
                    p = points(:, point)
                    call manufactured(p, sign_j, periods(n), family, actual, info)
                    call require(info == primitive_cartesian_field_ok, 'curved chart')
                    call manufactured_oracle(p, family, b, gradient)
                    call check_point(actual, p, b, gradient, 3.0e-13_dp)
                    call require(actual%signed_jacobian * sign_j > 0.0_dp, &
                        'curved chart orientation')
                    do level = 1, 3
                        h = 2.0e-4_dp / real(5**(level - 1), dp)
                        do axis = 1, 3
                            trial = p
                            trial(axis) = p(axis) + h
                            call manufactured(trial, sign_j, periods(n), &
                                family, plus, info)
                            call require(info == primitive_cartesian_field_ok, &
                                'FD plus')
                            trial(axis) = p(axis) - h
                            call manufactured(trial, sign_j, periods(n), &
                                family, minus, info)
                            call require(info == primitive_cartesian_field_ok, &
                                'FD minus')
                            fd(:, axis) = (plus%b - minus%b) / (2.0_dp * h)
                        end do
                        errors(level) = maxval(abs(fd - actual%gradient))
                        fd_maximum(level) = max(fd_maximum(level), errors(level))
                        call curl_of(fd, curl)
                        call near(actual%curl, curl, 3.0e-8_dp, 'Cartesian FD curl')
                    end do
                    call require(maxval(errors) < 3.0e-8_dp, &
                        'three-step Cartesian FD gradient')
                end do
            end do
        end do
    end do
    write (*, '(a,3es12.4)') 'Cartesian FD gradient max errors: ', fd_maximum
    do chart = 1, 2
        sign_j = real(3 - 2 * chart, dp)
        do n = 1, size(periods)
            do point = 1, 3
                s = 0.11_dp + 0.23_dp * real(point, dp)
                theta = 0.037_dp + 0.127_dp * real(point, dp)
                zeta = 0.021_dp + 0.083_dp * real(point, dp)
                call vacuum_torus(s, theta, zeta, periods(n), sign_j, p, jets, phi)
                call submit(p, jets, periods(n), phi, 0.0_dp, 0.0_dp, &
                    0.0_dp, actual, info)
                call require(info == primitive_cartesian_field_ok, 'vacuum torus')
                call vacuum_oracle(p, b, gradient)
                call check_point(actual, p, b, gradient, 2.0e-11_dp)
                call near(actual%curl, [0.0_dp, 0.0_dp, 0.0_dp], &
                    2.0e-11_dp, 'vacuum torus has no current')
            end do
            call check_screw_pinch(periods(n), sign_j)
        end do
    end do
    call check_field_dynamic_range()
    call check_magnetic_null()
    call check_invalid()
    write (*, '(a)') 'PASS'

contains

    subroutine manufactured(position, orientation, nfp, kind, field, status)
        real(dp), intent(in) :: position(3), orientation
        integer, intent(in) :: nfp, kind
        type(primitive_cartesian_field_point_t), intent(out) :: field
        integer, intent(out) :: status
        real(dp), parameter :: k = 0.7_dp, l = -0.4_dp, a = 0.6_dp
        real(dp) :: j(3, 9), c, d, cs, ds, scale

        c = 0.8_dp + 0.3_dp * position(1)
        cs = 0.3_dp
        d = -0.2_dp + 0.4_dp * position(1) + 0.1_dp * position(1)**2
        ds = 0.4_dp + 0.2_dp * position(1)
        j = 0.0_dp
        j(1, 1) = 1.0_dp
        j(3, 3) = 1.0_dp
        if (kind == 1) then
            j(3, 1) = k * position(2)
            j(2, 2) = orientation
            j(3, 2) = orientation * (k * position(1) + l * position(2))
            j(3, 5) = orientation * k
            j(3, 7) = l
        else
            scale = exp(a * position(3))
            j(2, 2) = orientation * scale
            j(2, 3) = a * position(2)
            j(2, 8) = orientation * a * scale
            j(2, 9) = a**2 * position(2)
        end if
        call submit(position, j, nfp, -orientation * d, -real(nfp, dp) * c, &
            -orientation * ds, -real(nfp, dp) * cs, field, status)
    end subroutine manufactured

    subroutine manufactured_oracle(position, kind, field, g)
        real(dp), intent(in) :: position(3)
        integer, intent(in) :: kind
        real(dp), intent(out) :: field(3), g(3, 3)
        real(dp), parameter :: k = 0.7_dp, l = -0.4_dp, a = 0.6_dp
        real(dp) :: c, d, ds, factor

        c = 0.8_dp + 0.3_dp * position(1)
        d = -0.2_dp + 0.4_dp * position(1) + 0.1_dp * position(1)**2
        ds = 0.4_dp + 0.2_dp * position(1)
        field = 0.0_dp
        g = 0.0_dp
        if (kind == 1) then
            factor = k * position(1) + l * position(2)
            field(2) = c
            field(3) = c * factor + d
            g(2, 1) = 0.3_dp
            g(3, 1) = 0.3_dp * factor + k * c + ds
            g(3, 2) = l * c
        else
            factor = exp(-a * position(3))
            field(2) = c + a * d * position(2) * factor
            field(3) = d * factor
            g(2, 1) = 0.3_dp + a * ds * position(2) * factor
            g(2, 2) = a * d * factor
            g(2, 3) = -a**2 * d * position(2) * factor
            g(3, 1) = ds * factor
            g(3, 3) = -a * d * factor
        end if
    end subroutine manufactured_oracle

    subroutine vacuum_torus(s0, t, z, nfp, orientation, position, j, flux)
        real(dp), intent(in) :: s0, t, z, orientation
        integer, intent(in) :: nfp
        real(dp), intent(out) :: position(3), j(3, 9), flux
        real(dp), parameter :: r0 = 3.0_dp, c = r0 - sqrt(r0**2 - 1.0_dp)
        complex(dp) :: v, omega, f, fv, fvv, w, ws, wt, wss, wst, wtt
        real(dp) :: d, q, logarithmic, logarithmic_s, rate, cosine, sine

        d = r0 - c * s0
        q = sqrt(c * s0 / (2.0_dp * r0 - c * s0))
        logarithmic = r0 / (s0 * (2.0_dp * r0 - c * s0))
        logarithmic_s = -2.0_dp * r0 * d / (s0 * (2.0_dp * r0 - c * s0))**2
        omega = cmplx(0.0_dp, -orientation * 2.0_dp * pi, dp)
        v = q * exp(omega * t)
        f = (1.0_dp + v) / (1.0_dp - v)
        fv = 2.0_dp / (1.0_dp - v)**2
        fvv = 4.0_dp / (1.0_dp - v)**3
        w = d * f
        ws = -c * f + d * fv * logarithmic * v
        wt = d * fv * omega * v
        wss = -2.0_dp * c * fv * logarithmic * v &
            + d * (fvv * (logarithmic * v)**2 &
            + fv * (logarithmic**2 + logarithmic_s) * v)
        wst = -c * fv * omega * v &
            + d * (fvv * logarithmic * omega * v**2 + fv * logarithmic * omega * v)
        wtt = d * (fvv * (omega * v)**2 + fv * omega**2 * v)
        rate = 2.0_dp * pi / real(nfp, dp)
        cosine = cos(rate * z)
        sine = sin(rate * z)
        call embed(w, cosine, sine, position)
        call embed(ws, cosine, sine, j(:, 1))
        call embed(wt, cosine, sine, j(:, 2))
        call embed(wss, cosine, sine, j(:, 4))
        call embed(wst, cosine, sine, j(:, 5))
        call embed(wtt, cosine, sine, j(:, 7))
        j(1, 3) = -rate * real(w, dp) * sine
        j(2, 3) = rate * real(w, dp) * cosine
        j(3, 3) = 0.0_dp
        j(1, 6) = -rate * real(ws, dp) * sine
        j(2, 6) = rate * real(ws, dp) * cosine
        j(3, 6) = 0.0_dp
        j(1, 8) = -rate * real(wt, dp) * sine
        j(2, 8) = rate * real(wt, dp) * cosine
        j(3, 8) = 0.0_dp
        j(1, 9) = -rate**2 * real(w, dp) * cosine
        j(2, 9) = -rate**2 * real(w, dp) * sine
        j(3, 9) = 0.0_dp
        flux = -orientation * 2.0_dp * pi * c * 3.0_dp
    end subroutine vacuum_torus

    subroutine embed(w, cosine, sine, vector)
        complex(dp), intent(in) :: w
        real(dp), intent(in) :: cosine, sine
        real(dp), intent(out) :: vector(3)

        vector(1) = real(w, dp) * cosine
        vector(2) = real(w, dp) * sine
        vector(3) = aimag(w)
    end subroutine embed

    subroutine vacuum_oracle(position, field, g)
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: field(3), g(3, 3)
        real(dp) :: radius2

        radius2 = sum(position(1:2)**2)
        field(1) = -3.0_dp * position(2) / radius2
        field(2) = 3.0_dp * position(1) / radius2
        field(3) = 0.0_dp
        g = 0.0_dp
        g(1, 1) = 6.0_dp * product(position(1:2)) / radius2**2
        g(2, 2) = -g(1, 1)
        g(1, 2) = 3.0_dp * (position(2)**2 - position(1)**2) / radius2**2
        g(2, 1) = g(1, 2)
    end subroutine vacuum_oracle

    subroutine check_screw_pinch(nfp, orientation)
        integer, intent(in) :: nfp
        real(dp), intent(in) :: orientation
        real(dp), parameter :: s0 = 0.42_dp, a = 0.02_dp, b0 = 0.7_dp
        real(dp) :: position(3), j(3, 9), field(3), g(3, 3), e_r(3), e_phi(3)
        real(dp) :: radius, angle, flux, current(3), force(3), pressure_gradient(3)
        type(primitive_cartesian_field_point_t) :: result
        integer :: status

        radius = sqrt(2.0_dp * s0)
        angle = orientation * 2.0_dp * pi * 0.173_dp
        e_r = 0.0_dp
        e_phi = 0.0_dp
        e_r(1) = cos(angle)
        e_r(2) = sin(angle)
        e_phi(1) = -sin(angle)
        e_phi(2) = cos(angle)
        position = radius * e_r
        position(3) = -0.21_dp * length
        j = 0.0_dp
        j(:, 1) = e_r / radius
        j(:, 2) = orientation * 2.0_dp * pi * radius * e_phi
        j(3, 3) = length
        j(:, 4) = -e_r / radius**3
        j(:, 5) = orientation * 2.0_dp * pi * e_phi / radius
        j(:, 7) = -(2.0_dp * pi)**2 * radius * e_r
        flux = -orientation * 2.0_dp * pi * b0
        call submit(position, j, nfp, flux, -real(nfp, dp) * length * a, &
            0.0_dp, 0.0_dp, result, status)
        call require(status == primitive_cartesian_field_ok, 'screw pinch')
        field(1) = -a * position(2)
        field(2) = a * position(1)
        field(3) = b0
        g = 0.0_dp
        g(1, 2) = -a
        g(2, 1) = a
        call check_point(result, position, field, g, 3.0e-13_dp)
        current = result%curl / vacuum_permeability
        call require(abs(current(3) / (2.0_dp * a / vacuum_permeability) - 1.0_dp) &
            < 3.0e-13_dp, 'screw pinch current SI units')
        call cross(current, result%b, force)
        pressure_gradient = -2.0_dp * a**2 / vacuum_permeability * radius * e_r
        call near(force, pressure_gradient, 3.0e-11_dp, &
            'finite-pressure screw pinch balance')
        call require(2000.0_dp - 2.0_dp * a**2 * s0 / vacuum_permeability > 0.0_dp, &
            'screw pinch has positive finite pressure')
    end subroutine check_screw_pinch

    subroutine submit(position, j, nfp, phi_s, chi_s, phi_ss, chi_ss, field, status)
        real(dp), intent(in) :: position(3), j(3, 9), phi_s, chi_s, phi_ss, chi_ss
        integer, intent(in) :: nfp
        type(primitive_cartesian_field_point_t), intent(out) :: field
        integer, intent(out) :: status

        call build_primitive_cartesian_field_point(position, j(:, 1), j(:, 2), &
            j(:, 3), j(:, 4), j(:, 5), j(:, 6), j(:, 7), j(:, 8), j(:, 9), &
            nfp, phi_s, chi_s, phi_ss, chi_ss, field, status)
    end subroutine submit

    subroutine check_field_dynamic_range()
        real(dp), parameter :: expected_b(3) = [0.0_dp, 0.0_dp, -1.0_dp]
        real(dp), parameter :: expected_curl(3) = [0.0_dp, -1.0_dp, 0.0_dp]
        type(primitive_cartesian_field_point_t) :: result
        real(dp) :: position(3), j(3, 9), scaled(3), gradient_expected(3, 3)
        integer :: status

        position = 0.0_dp
        j = 0.0_dp
        j(1, 1) = 1.0e100_dp
        j(2, 2) = 1.0e100_dp
        j(3, 3) = 1.0e100_dp
        j(1, 4) = 1.0e108_dp
        call submit(position, j, 1, 1.0e308_dp, 0.0_dp, 0.0_dp, 0.0_dp, result, status)
        call require(status == primitive_cartesian_field_ok, 'finite field dynamic range')
        scaled = result%b / 1.0e108_dp
        call near(scaled, expected_b, 3.0e-13_dp, 'field dynamic range value')
        gradient_expected = 0.0_dp
        gradient_expected(3, 1) = 1.0_dp
        call require(all(abs(result%gradient / 1.0e16_dp - gradient_expected) &
            < 3.0e-13_dp), 'field dynamic range gradient')
        scaled = result%curl / 1.0e16_dp
        call near(scaled, expected_curl, 3.0e-13_dp, 'field dynamic range curl')
        call require(abs(result%divergence / 1.0e16_dp) < 3.0e-13_dp, &
            'field dynamic range divergence')
    end subroutine check_field_dynamic_range

    subroutine check_magnetic_null()
        type(primitive_cartesian_field_point_t) :: result
        real(dp) :: j(3, 9), zeros(3)
        integer :: status

        j = 0.0_dp
        j(1, 1) = 1.0_dp
        j(2, 2) = 1.0_dp
        j(3, 3) = 1.0_dp
        zeros = 0.0_dp
        call submit(points(:, 1), j, 3, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, result, status)
        call require(status == primitive_cartesian_field_ok, 'regular magnetic null')
        call near(result%b, zeros, 0.0_dp, 'null magnetic field')
        call near(result%curl, zeros, 0.0_dp, 'null curl')
        call require(all(abs(result%gradient) <= tiny(1.0_dp)) &
            .and. abs(result%divergence) <= tiny(1.0_dp), 'null field derivatives')
    end subroutine check_magnetic_null

    subroutine check_invalid()
        type(primitive_cartesian_field_point_t) :: result
        real(dp) :: position(3), j(3, 9), zeros(3)
        integer :: status, kind, nfp

        zeros = 0.0_dp
        do kind = 1, 4
            position = points(:, 1)
            j = 0.0_dp
            j(1, 1) = 1.0_dp
            j(2, 2) = 1.0_dp
            j(3, 3) = 1.0_dp
            nfp = 2
            select case (kind)
            case (1)
                nfp = 0
            case (2)
                j(:, 3) = j(:, 2)
            case (3)
                position(2) = ieee_value(0.0_dp, ieee_quiet_nan)
            case (4)
                j(1, 9) = ieee_value(0.0_dp, ieee_quiet_nan)
            end select
            call submit(position, j, nfp, 0.3_dp, -0.7_dp, &
                0.1_dp, 0.2_dp, result, status)
            call require(status == primitive_cartesian_field_invalid, &
                'invalid field input')
            call near(result%position, zeros, 0.0_dp, 'failed point has no position')
            call near(result%b, zeros, 0.0_dp, 'failed point has no field')
            call near(result%curl, zeros, 0.0_dp, 'failed point has no curl')
            call require(all(abs(result%gradient) <= tiny(1.0_dp)) &
                .and. abs(result%divergence) <= tiny(1.0_dp) &
                .and. abs(result%signed_jacobian) <= tiny(1.0_dp), &
                'no partial field output')
        end do
    end subroutine check_invalid

    subroutine check_point(field, position, b_expected, g_expected, tolerance)
        type(primitive_cartesian_field_point_t), intent(in) :: field
        real(dp), intent(in) :: position(3), b_expected(3), g_expected(3, 3), tolerance
        real(dp) :: curl_expected(3)

        call near(field%position, position, tolerance, 'Cartesian position')
        call near(field%b, b_expected, tolerance, 'Cartesian magnetic field')
        call require(all(abs(field%gradient - g_expected) <= tolerance &
            * max(1.0_dp, abs(g_expected))), 'Cartesian magnetic gradient')
        call curl_of(g_expected, curl_expected)
        call near(field%curl, curl_expected, tolerance, 'Cartesian curl')
        call require(abs(field%divergence) < tolerance, 'solenoidal Cartesian field')
    end subroutine check_point

    subroutine curl_of(g, vector)
        real(dp), intent(in) :: g(3, 3)
        real(dp), intent(out) :: vector(3)

        vector(1) = g(3, 2) - g(2, 3)
        vector(2) = g(1, 3) - g(3, 1)
        vector(3) = g(2, 1) - g(1, 2)
    end subroutine curl_of

    subroutine cross(left, right, vector)
        real(dp), intent(in) :: left(3), right(3)
        real(dp), intent(out) :: vector(3)

        vector(1) = left(2) * right(3) - left(3) * right(2)
        vector(2) = left(3) * right(1) - left(1) * right(3)
        vector(3) = left(1) * right(2) - left(2) * right(1)
    end subroutine cross

    subroutine near(actual_value, expected_value, tolerance, label)
        real(dp), intent(in) :: actual_value(:), expected_value(:), tolerance
        character(*), intent(in) :: label

        call require(all(abs(actual_value - expected_value) <= tolerance &
            * max(1.0_dp, abs(expected_value))), label)
    end subroutine near

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        write (error_unit, '(a)') 'FAIL: ' // label
        error stop 1
    end subroutine require
end program test_primitive_cartesian_field
