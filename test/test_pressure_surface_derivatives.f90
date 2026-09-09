program test_pressure_surface_derivatives
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use export_surface_geometry, only: build_surface_kernel_fields, mercier_ok, &
        surface_data_t, surface_profiles_t, two_pi, mu0
    use pressure_surface_derivatives, only: build_pressure_surface_response, &
        pressure_derivative_ok, pressure_surface_jvp, pressure_surface_response_t, &
        pressure_surface_vjp
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        evaluate_radial_cubic_spline, fit_radial_cubic_spline, &
        radial_cubic_spline_grid_t, radial_cubic_spline_ok, radial_cubic_spline_t
    implicit none

    integer, parameter :: m(3) = [0, 1, 2], n(3) = [0, 1, -1]
    real(dp), parameter :: coordinate = 0.43_dp, gamma = 5.0_dp / 3.0_dp
    real(dp), parameter :: steps(3) = [1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp]
    type(surface_data_t) :: surface
    type(surface_profiles_t) :: profiles
    type(pressure_surface_response_t) :: response, invalid
    real(dp), allocatable :: theta(:), zeta(:), jacobian_slope(:, :)
    real(dp), allocatable :: nodes(:), pressure(:), tangent(:), gradient(:)
    real(dp), allocatable :: perturbed_pressure(:)
    real(dp), allocatable :: fields(:, :, :), drive(:, :), fcot(:, :, :), dcot(:, :)
    real(dp), allocatable :: plus(:, :, :), minus(:, :, :), dpplus(:, :), dpminus(:, :)
    real(dp) :: gp, gpplus, gpminus, lhs, rhs, gcot, scale, expected, nan
    integer :: resolution, count, chart_index, step, status
    logical :: chart

    do resolution = 16, 32, 16
        call make_surface(resolution)
        do count = 5, 9, 4
            call make_samples(count)
            do chart_index = 0, 1
                chart = chart_index == 1
                call build_pressure_surface_response(m, n, chart, surface, profiles, &
                    jacobian_slope, theta, zeta, nodes, pressure, coordinate, &
                    gamma, response, status)
                call require(status == pressure_derivative_ok, "response failed")
                call pressure_surface_jvp(response, tangent, fields, drive, gp, status)
                call require(status == pressure_derivative_ok, "JVP failed")
                expected = gamma * 1.0e5_dp * (1.0_dp + 2.0_dp * coordinate &
                    + 3.0_dp * coordinate**2 + 4.0_dp * coordinate**3)
                call require(abs(gp - expected) < 1.0e-12_dp * abs(expected), &
                    "cubic pressure interpolation derivative is inexact")
                expected = mu0 * 1.0e5_dp * (2.0_dp + 6.0_dp * coordinate &
                    + 12.0_dp * coordinate**2)
                call require(maxval(abs(fields(:, :, 11) - expected)) &
                    < 1.0e-12_dp * abs(expected), "cubic pressure slope is inexact")
                call primal(pressure, chart, plus, dpplus, gpplus)
                call require(maxval(abs(response%fields - plus)) < 1.0e-12_dp, &
                    "response primal fields differ from production")
                call require(maxval(abs(response%drive - dpplus)) < 1.0e-12_dp, &
                    "response primal drive differs from production")
                call require(abs(response%gamma_pressure - gpplus) < 1.0e-8_dp, &
                    "response primal gamma*p differs from production")
                fcot = cos(response%fields)
                dcot = sin(response%drive)
                gcot = 0.31_dp
                call pressure_surface_vjp(response, fcot, dcot, gcot, gradient, status)
                call require(status == pressure_derivative_ok, "VJP failed")
                lhs = sum(fields * fcot) + sum(drive * dcot) + gp * gcot
                rhs = dot_product(tangent, gradient)
                call require(abs(lhs - rhs) < 2.0e-12_dp * max(1.0_dp, abs(lhs)), &
                    "pressure sample JVP/VJP duality failed")
                do step = 1, size(steps)
                    perturbed_pressure = pressure + steps(step) * tangent
                    call primal(perturbed_pressure, chart, plus, dpplus, gpplus)
                    perturbed_pressure = pressure - steps(step) * tangent
                    call primal(perturbed_pressure, chart, minus, dpminus, gpminus)
                    scale = max(1.0_dp, maxval(abs(fields)))
                    call require(maxval(abs((plus - minus) / (2.0_dp * steps(step)) &
                        - fields)) < 1.0e-8_dp * scale, "field FD plateau failed")
                    scale = max(1.0_dp, maxval(abs(drive)))
                    call require(maxval(abs((dpplus - dpminus) &
                        / (2.0_dp * steps(step)) - drive)) < 1.0e-8_dp * scale, &
                        "drive FD plateau failed")
                    call require(abs((gpplus - gpminus) / (2.0_dp * steps(step)) &
                        - gp) < 1.0e-8_dp * abs(gp), "gamma*p FD plateau failed")
                end do
            end do
        end do
    end do
    call check_cardinal_directions()
    call check_invalid()
    call check_resonance()
    write (*, '(a)') "PASS pressure samples -> spline -> production surface JVP/VJP"

contains

    subroutine make_samples(size_samples)
        integer, intent(in) :: size_samples
        integer :: i

        if (allocated(nodes)) deallocate (nodes, pressure, tangent)
        allocate (nodes(size_samples), pressure(size_samples), tangent(size_samples))
        do i = 1, size_samples
            nodes(i) = (real(i, dp) - 0.5_dp) / real(size_samples, dp)
        end do
        pressure = 1.0e6_dp * (4.0_dp + 2.0_dp * nodes &
            - 0.4_dp * nodes**2 + 0.3_dp * nodes**3)
        tangent = 1.0e5_dp * (1.0_dp + 2.0_dp * nodes &
            + 3.0_dp * nodes**2 + 4.0_dp * nodes**3)
    end subroutine make_samples

    subroutine make_surface(points)
        integer, intent(in) :: points
        integer :: i, j
        real(dp) :: phase

        surface = surface_data_t()
        if (allocated(theta)) deallocate (theta, zeta, jacobian_slope)
        allocate (theta(points), zeta(points), jacobian_slope(points, points))
        do i = 1, points
            theta(i) = real(i - 1, dp) / real(points, dp)
            zeta(i) = theta(i)
        end do
        allocate (surface%jacobian(points, points), surface%g_tt(points, points), &
            surface%g_tz(points, points), surface%g_zz(points, points), &
            surface%b_theta(points, points), surface%b_zeta(points, points), &
            surface%g_st(points, points), surface%g_sz(points, points), &
            surface%mod_b(points, points))
        do j = 1, points
            do i = 1, points
                phase = two_pi * (theta(i) - zeta(j))
                surface%jacobian(i, j) = 1.0_dp + 0.15_dp * cos(phase)
                surface%g_tt(i, j) = 2.0_dp + 0.1_dp * sin(phase)
                surface%g_tz(i, j) = 0.2_dp
                surface%g_zz(i, j) = 3.0_dp
                surface%b_theta(i, j) = 0.7_dp + 0.04_dp * cos(phase)
                surface%b_zeta(i, j) = 1.3_dp
                surface%g_st(i, j) = 0.07_dp * cos(phase)
                surface%g_sz(i, j) = 0.05_dp * sin(phase)
                surface%mod_b(i, j) = 2.0_dp
                jacobian_slope(i, j) = 0.13_dp + 0.02_dp * cos(phase)
            end do
        end do
        profiles = surface_profiles_t(1.3_dp, 0.71_dp, 0.12_dp, -0.17_dp, &
            0.7_dp, 1.3_dp, 0.15_dp, -0.3_dp, 0.0_dp)
    end subroutine make_surface

    subroutine primal(samples, use_chart, field_values, drive_values, gamma_value)
        real(dp), intent(in) :: samples(:)
        logical, intent(in) :: use_chart
        real(dp), allocatable, intent(out) :: field_values(:, :, :), drive_values(:, :)
        real(dp), intent(out) :: gamma_value
        type(radial_cubic_spline_grid_t) :: grid
        type(radial_cubic_spline_t) :: spline
        type(surface_profiles_t) :: active
        real(dp) :: value, slope
        integer :: info

        ! Independent production primal: the scalar spline API and surface
        ! builder contain no calls to the new differentiated implementation.
        call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
        call require(info == radial_cubic_spline_ok, "oracle grid failed")
        call fit_radial_cubic_spline(grid, samples, spline, info)
        call require(info == radial_cubic_spline_ok, "oracle fit failed")
        call evaluate_radial_cubic_spline(grid, spline, coordinate, value, slope, info)
        call require(info == radial_cubic_spline_ok, "oracle interpolation failed")
        active = profiles
        active%pressure_slope = slope
        allocate (field_values(size(theta), size(zeta), 13), &
            drive_values(size(theta), size(zeta)))
        call build_surface_kernel_fields(m, n, use_chart, surface, active, &
            jacobian_slope, theta, zeta, field_values, drive_values, info)
        call require(info == mercier_ok, "oracle surface failed")
        call require(all(ieee_is_finite(field_values)), "nonfinite field oracle")
        call require(all(ieee_is_finite(drive_values)), "nonfinite drive oracle")
        gamma_value = gamma * value
    end subroutine primal

    subroutine check_cardinal_directions()
        integer :: sample, local_chart, local_step

        ! Discretization coverage above is not a convergence claim. Here each
        ! coordinate direction tests the full nine-dimensional nodal map on a
        ! nonuniform grid, beyond the global cubic reproduction subspace.
        call make_samples(9)
        nodes = [0.03_dp, 0.09_dp, 0.18_dp, 0.27_dp, 0.42_dp, &
            0.55_dp, 0.70_dp, 0.83_dp, 0.96_dp]
        pressure = 4.0e6_dp + 1.0e6_dp * sin(5.0_dp * nodes) &
            + 0.2e6_dp * cos(13.0_dp * nodes)
        do local_chart = 0, 1
            chart = local_chart == 1
            call build_pressure_surface_response(m, n, chart, surface, profiles, &
                jacobian_slope, theta, zeta, nodes, pressure, coordinate, &
                gamma, response, status)
            call require(status == pressure_derivative_ok, "nodal response failed")
            fcot = cos(response%fields)
            dcot = sin(response%drive)
            gcot = 0.31_dp
            call pressure_surface_vjp(response, fcot, dcot, gcot, gradient, status)
            call require(status == pressure_derivative_ok, "nodal VJP failed")
            do sample = 1, size(nodes)
                tangent = 0.0_dp
                ! h is a relative perturbation of this nodal pressure, which
                ! keeps the scalar-spline oracle above subtraction roundoff.
                tangent(sample) = pressure(sample)
                call pressure_surface_jvp(response, tangent, fields, drive, gp, status)
                call require(status == pressure_derivative_ok, "nodal JVP failed")
                lhs = sum(fields * fcot) + sum(drive * dcot) + gp * gcot
                rhs = dot_product(tangent, gradient)
                call require(abs(lhs - rhs) < 2.0e-12_dp * max(1.0_dp, abs(lhs)), &
                    "coordinate-basis JVP/VJP duality failed")
                do local_step = 1, size(steps)
                    perturbed_pressure = pressure + steps(local_step) * tangent
                    call primal(perturbed_pressure, chart, plus, dpplus, gpplus)
                    perturbed_pressure = pressure - steps(local_step) * tangent
                    call primal(perturbed_pressure, chart, minus, dpminus, gpminus)
                    scale = max(1.0_dp, maxval(abs(fields)))
                    call require(maxval(abs((plus - minus) &
                        / (2.0_dp * steps(local_step)) - fields)) &
                        < 1.0e-8_dp * scale, "nodal field FD plateau failed")
                    scale = max(1.0_dp, maxval(abs(drive)))
                    call require(maxval(abs((dpplus - dpminus) &
                        / (2.0_dp * steps(local_step)) - drive)) &
                        < 1.0e-8_dp * scale, "nodal drive FD plateau failed")
                    call require(abs((gpplus - gpminus) &
                        / (2.0_dp * steps(local_step)) - gp) &
                        < 1.0e-8_dp * max(1.0_dp, abs(gp)), &
                        "nodal gamma*p FD plateau failed")
                end do
            end do
        end do
    end subroutine check_cardinal_directions

    subroutine check_invalid()
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call pressure_surface_jvp(response, tangent(:2), fields, drive, gp, status)
        call require(status /= pressure_derivative_ok, "wrong tangent size accepted")
        tangent(1) = nan
        call pressure_surface_jvp(response, tangent, fields, drive, gp, status)
        call require(status /= pressure_derivative_ok, "NaN tangent accepted")
        call pressure_surface_vjp(response, fcot, dcot, nan, gradient, status)
        call require(status /= pressure_derivative_ok, "NaN cotangent accepted")
        call pressure_surface_vjp(response, fcot(:, :, :2), dcot, gcot, gradient, status)
        call require(status /= pressure_derivative_ok, "wrong cotangent shape accepted")
        call build_pressure_surface_response(m, n, .true., surface, profiles, &
            jacobian_slope, theta, zeta, nodes, pressure, -0.1_dp, gamma, &
            invalid, status)
        call require(status /= pressure_derivative_ok, "outside coordinate accepted")
        call build_pressure_surface_response(m, n, .true., surface, profiles, &
            jacobian_slope, theta, zeta, nodes, pressure, coordinate, 0.0_dp, &
            invalid, status)
        call require(status /= pressure_derivative_ok, "zero gamma accepted")
        pressure(1) = -1.0_dp
        call build_pressure_surface_response(m, n, .true., surface, profiles, &
            jacobian_slope, theta, zeta, nodes, pressure, coordinate, gamma, &
            invalid, status)
        call require(status /= pressure_derivative_ok, "negative pressure accepted")
        pressure(1) = nan
        call build_pressure_surface_response(m, n, .true., surface, profiles, &
            jacobian_slope, theta, zeta, nodes, pressure, coordinate, gamma, &
            invalid, status)
        call require(status /= pressure_derivative_ok, "NaN pressure accepted")
        call require(.not. invalid%ready, "invalid response marked ready")
    end subroutine check_invalid

    subroutine check_resonance()
        type(surface_profiles_t) :: active

        call make_samples(5)
        pressure = 4.0e6_dp
        profiles%flux_slope = 1.0_dp
        profiles%poloidal_slope = 1.0_dp
        profiles%covariant_theta_slope = 0.0_dp
        profiles%covariant_zeta_slope = 0.0_dp
        active = profiles
        active%pressure_slope = 0.0_dp
        ! Primal constant pressure has zero beta forcing, but a varying pressure
        ! tangent excites the (1,1) Jac harmonic exactly at resonance.
        call primal(pressure, .true., plus, dpplus, gpplus)
        call build_pressure_surface_response(m, n, .true., surface, active, &
            jacobian_slope, theta, zeta, nodes, pressure, coordinate, gamma, &
            invalid, status)
        call require(status /= pressure_derivative_ok, &
            "incompatible resonant pressure tangent accepted")
        surface%jacobian = 1.0_dp
        call build_pressure_surface_response(m, n, .true., surface, active, &
            jacobian_slope, theta, zeta, nodes, pressure, coordinate, gamma, &
            response, status)
        call require(status == pressure_derivative_ok, &
            "compatible mean-only pressure tangent rejected")
    end subroutine check_resonance

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') message
            error stop 1
        end if
    end subroutine require
end program test_pressure_surface_derivatives
