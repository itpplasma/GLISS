! These controls qualify generated leaves and expose a space obstruction.
! They do not replace the production tangential space or close Solov'ev #35.
program test_weighted_feec_compatibility
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use radial_feec_complex, only: radial_feec_complex_t, &
        build_radial_feec_complex, evaluate_radial_feec_complex, radial_feec_ok
    use compatible_problem_assembly_support, only: apply_stored_power, &
        apply_tangential_axis_weight, compatible_support_ok
    use gliss_weighted_feec_generated, only: generated_weighted_feec_jet, &
        generated_weighted_feec_jvp, generated_weighted_feec_vjp
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use phase_assembly_policy, only: phase_assembly_direct, phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)

    call check_generated_jets()
    call check_commuting_maps()
    call check_cylinder_witness()
    write (*, '(a)') 'weighted FEEC compatibility controls: PASS; Solovev gate unchanged'

contains

    subroutine check_generated_jets()
        real(dp), parameter :: powers(7) = [0.0_dp, 0.5_dp, -0.5_dp, &
            -1.0_dp, 1.0_dp, 2.0_dp, 0.25_dp]
        real(dp), parameter :: points(3) = [0.125_dp, 0.375_dp, 0.875_dp]
        real(dp), parameter :: seed(5) = [0.3_dp, -0.2_dp, 0.4_dp, -0.1_dp, 0.7_dp]
        real(dp), parameter :: cotangent(3) = [0.6_dp, -0.3_dp, 0.2_dp]
        real(dp) :: coefficients(0:4), local(0:4), njet(0:4), exact(0:4)
        real(dp) :: jet(3), tangent(3), gradient(5), jacobian(3, 5), alpha_jet(3)
        real(dp) :: reference_gradient(5)
        real(dp) :: s, alpha, exponent, term, actual(3), reference(3)
        integer :: degree, power, point, j, order, k

        do degree = 0, 4
            coefficients = 0.0_dp
            ! Independent monomial expansion of (1-s)^degree.
            do j = 0, degree
                coefficients(j) = real((-1)**j * binomial(degree, j), dp)
            end do
            do power = 1, size(powers)
                alpha = powers(power)
                do point = 1, size(points)
                    s = points(point)
                    call monomial_jet(coefficients, 0.0_dp, s, njet)
                    call monomial_jet(coefficients, alpha, s, exact)
                    call generated_weighted_feec_jet(s, alpha, njet(0), &
                        njet(1), njet(2), actual(1), actual(2), actual(3))
                    reference = exact(0:2)
                    call near_vector(actual, reference, 3.0e-12_dp, 'generated jet')
                    alpha_jet = 0.0_dp
                    do j = 0, degree
                        exponent = real(j, dp) - alpha
                        term = coefficients(j) * s**exponent
                        alpha_jet(1) = alpha_jet(1) - log(s) * term
                        alpha_jet(2) = alpha_jet(2) &
                            - (1.0_dp + exponent * log(s)) * term / s
                        alpha_jet(3) = alpha_jet(3) - (2.0_dp * exponent - 1.0_dp &
                            + exponent * (exponent - 1.0_dp) * log(s)) * term / s**2
                    end do
                    ! Three centered polynomials independently vary N, N_s, N_ss.
                    do k = 0, 2
                        local = 0.0_dp
                        do j = 0, k
                            local(j) = real(binomial(k, j), dp) * (-s)**(k - j)
                        end do
                        if (k == 2) local = 0.5_dp * local
                        call monomial_jet(local, alpha, s, exact)
                        jacobian(:, k + 3) = exact(0:2)
                    end do
                    call monomial_jet(coefficients, alpha, s, exact)
                    jacobian(:, 2) = alpha_jet
                    do order = 0, 2
                        jacobian(order + 1, 1) = exact(order + 1)
                        do k = 0, 2
                            jacobian(order + 1, 1) = jacobian(order + 1, 1) &
                                - jacobian(order + 1, k + 3) * njet(k + 1)
                        end do
                    end do
                    call generated_weighted_feec_jvp(s, alpha, njet(0), njet(1), &
                        njet(2), seed(1), seed(2), seed(3), seed(4), seed(5), &
                        jet, tangent)
                    reference = matmul(jacobian, seed)
                    call near_vector(jet, actual, 3.0e-12_dp, 'JVP primal')
                    call near_vector(tangent, reference, 8.0e-12_dp, 'independent JVP')
                    call generated_weighted_feec_vjp(s, alpha, njet(0), njet(1), &
                        njet(2), cotangent(1), cotangent(2), cotangent(3), jet, gradient)
                    call near_vector(jet, actual, 3.0e-12_dp, 'VJP primal')
                    reference_gradient = matmul(transpose(jacobian), cotangent)
                    call near_vector(gradient, reference_gradient, &
                        8.0e-12_dp, 'independent VJP')
                end do
            end do
        end do
        ! Accepted integer custom powers carry genuine derivative kernels.
        do degree = 1, 2
            coefficients = 0.0_dp
            coefficients(degree) = 1.0_dp
            call monomial_jet(coefficients, 0.0_dp, 0.375_dp, njet)
            call generated_weighted_feec_jet(0.375_dp, real(degree, dp), &
                njet(0), njet(1), njet(2), actual(1), actual(2), actual(3))
            reference = [1.0_dp, 0.0_dp, 0.0_dp]
            call near_vector(actual, reference, 3.0e-13_dp, 'custom-power kernel')
        end do
    end subroutine check_generated_jets

    subroutine check_commuting_maps()
        real(dp), parameter :: points(4) = [0.125_dp, 0.375_dp, 0.625_dp, 0.875_dp]
        real(dp), parameter :: powers(3) = [0.0_dp, 0.5_dp, 2.0_dp]
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: h1(:), dh1(:), l2(:)
        real(dp) :: values(1, 3), derivatives(1, 3), eta(2, 3), ordinary, s
        integer :: info, point, power

        call build_radial_feec_complex([0.0_dp, 0.5_dp, 1.0_dp], 1, &
            .true., .true., complex, info)
        call require(info == radial_feec_ok, 'P1 complex build')
        do point = 1, size(points)
            s = points(point)
            call evaluate_radial_feec_complex(complex, s, h1, dh1, l2, info)
            call require(info == radial_feec_ok, 'P1 complex evaluation')
            ordinary = dot_product(complex%derivative(:, 1), l2)
            call near(ordinary, merge(2.0_dp, -2.0_dp, s < 0.5_dp), &
                2.0e-14_dp, 'independent P1 derivative')
            call near(dh1(1), ordinary, 2.0e-14_dp, 'ordinary commuting map')
            call apply_stored_power(s, powers, h1, dh1, [1], values, derivatives, info)
            call require(info == compatible_support_ok, 'production weighted normal')
            call apply_tangential_axis_weight(s, powers, l2, [1, 2], eta, info)
            call require(info == compatible_support_ok, 'production weighted eta')
            do power = 1, size(powers)
                call near(values(1, power), s**(-powers(power)) * h1(1), &
                    2.0e-14_dp, 'weighted value')
                call near(derivatives(1, power) &
                    - dot_product(complex%derivative(:, 1), eta(:, power)), &
                    -powers(power) * s**(-powers(power) - 1.0_dp) * h1(1), &
                    3.0e-14_dp, 'weighted commutator')
            end do
        end do
        call evaluate_radial_feec_complex(complex, 1.0_dp, h1, dh1, l2, info)
        call require(info == radial_feec_ok, 'wall trace evaluation')
        call near(h1(1), 0.0_dp, 0.0_dp, 'normal wall trace')
        ! A left-only primitive must keep its wall generator: its integral is one.
        call build_radial_feec_complex([0.0_dp, 0.5_dp, 1.0_dp], 1, &
            .true., .false., complex, info)
        call require(info == radial_feec_ok, 'left-trace primitive build')
        call evaluate_radial_feec_complex(complex, 1.0_dp, h1, dh1, l2, info)
        call require(info == radial_feec_ok, 'primitive wall evaluation')
        call near(h1(size(h1)), 1.0_dp, 0.0_dp, 'retained primitive wall trace')
    end subroutine check_commuting_maps

    subroutine check_cylinder_witness()
        real(dp), parameter :: points(4) = [0.125_dp, 0.375_dp, 0.625_dp, 0.875_dp]
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: h1(:), dh1(:), l2(:)
        real(dp) :: values(1, 1), derivatives(1, 1), s, eta, coefficient(4)
        real(dp) :: a, b, x, y, cartesian_divergence
        real(dp) :: direct_energy, transformed_energy
        integer :: info, point

        call build_radial_feec_complex([0.0_dp, 0.5_dp, 1.0_dp], 1, &
            .true., .true., complex, info)
        call require(info == radial_feec_ok, 'cylinder P1 build')
        do point = 1, size(points)
            s = points(point)
            call evaluate_radial_feec_complex(complex, s, h1, dh1, l2, info)
            call require(info == radial_feec_ok, 'cylinder P1 evaluation')
            call apply_stored_power(s, [0.5_dp], h1, dh1, [1], &
                values, derivatives, info)
            call require(info == compatible_support_ok, 'cylinder production jet')
            eta = -0.5_dp * derivatives(1, 1) ! T=pi, m=1.
            a = values(1, 1) / (2.0_dp * sqrt(s))
            b = 2.0_dp * sqrt(s) * eta
            if (s < 0.5_dp) then
                call near(a, 1.0_dp, 2.0e-14_dp, 'axis Cartesian radial amplitude')
                call near(b, -1.0_dp, 2.0e-14_dp, 'axis Cartesian tangent amplitude')
                cartesian_divergence = 0.0_dp ! xi=(1,0) on the first cell.
            else
                call near(a, 1.0_dp / s - 1.0_dp, 2.0e-14_dp, 'outer radial oracle')
                call near(b, 1.0_dp / s + 1.0_dp, 2.0e-14_dp, 'outer tangent oracle')
                x = sqrt(s) * 0.6_dp
                y = sqrt(s) * 0.8_dp
                ! xi=((X^2-Y^2)/s^2-1,2XY/s^2), s=X^2+Y^2.
                cartesian_divergence = 4.0_dp * x / s**2 &
                    - 4.0_dp * x * (x**2 + y**2) / s**3
            end if
            call near(cartesian_divergence, 0.0_dp, 2.0e-13_dp, 'Cartesian divergence')
            ! Uniform Bz=1: curl(xi cross B)=(0,0,-div xi), hence raw deltaW=0.
            call cylinder_surface(s, values(1, 1), derivatives(1, 1), eta, &
                phase_assembly_direct, direct_energy)
            call cylinder_surface(s, values(1, 1), derivatives(1, 1), eta, &
                phase_assembly_transformed, transformed_energy)
            call near(direct_energy, transformed_energy, 2.0e-8_dp, 'phase energy routes')
            coefficient(point) = eta * sqrt(s)
        end do
        ! On the second cell N=2-2s, eta*sqrt(s)=(1+1/s)/2.
        ! Two exact samples differ, so no single old weighted-P0 coefficient fits.
        call near(coefficient(3) - coefficient(4), 8.0_dp / 35.0_dp, &
            3.0e-14_dp, 'positive old-space obstruction')
        call require(abs(coefficient(3) - coefficient(4)) > 0.2_dp, &
            'interior weighted-P0 mismatch unexpectedly disappeared')
    end subroutine check_cylinder_witness

    subroutine cylinder_surface(s, normal, derivative, eta, route, energy)
        real(dp), intent(in) :: s, normal, derivative, eta
        integer, intent(in) :: route
        real(dp), intent(out) :: energy
        real(dp) :: fields(8, 2, 13), zeros(8, 2), gamma_pressure(8, 2)
        real(dp) :: h(1, 1), dh(1, 1), e(1, 1), parallel(1, 1)
        real(dp) :: stiffness(3, 3), mass(3, 3), expected(2, 2), response(2)
        real(dp) :: image(3), vector(3), kinetic, scale, amplitude_a, amplitude_b
        real(dp) :: actual_entries(4), expected_entries(4)
        integer :: info, row, column

        fields = 0.0_dp
        fields(:, :, 1) = pi
        fields(:, :, 5) = 1.0_dp
        fields(:, :, 7) = pi
        fields(:, :, 8) = 1.0_dp
        fields(:, :, 9) = 4.0_dp * s
        zeros = 0.0_dp
        gamma_pressure = 1.0_dp
        h = normal
        dh = derivative
        e = eta
        parallel = 1.0_dp
        stiffness = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, zeros, &
            zeros, zeros, zeros, gamma_pressure, [1], [0], [1], 1, h, dh, &
            e, parallel, 1.0_dp, route, stiffness, info)
        call require(info == 0, 'cylinder stiffness assembly')
        response(1) = derivative
        response(2) = 2.0_dp * eta
        scale = 0.5_dp * pi * (1.0_dp / vacuum_permeability + 1.0_dp)
        do column = 1, 2
            do row = 1, 2
                expected(row, column) = scale * response(row) * response(column)
            end do
        end do
        do column = 1, 2
            do row = 1, 2
                actual_entries(row + 2 * (column - 1)) = stiffness(row, column)
                expected_entries(row + 2 * (column - 1)) = expected(row, column)
            end do
        end do
        call near_vector(actual_entries, expected_entries, &
            2.0e-13_dp, 'independent raw cylinder stiffness')
        vector = [1.0_dp, 1.0_dp, 0.0_dp]
        image = matmul(stiffness, vector)
        energy = dot_product(vector, image)
        call require(abs(energy) <= 128.0_dp * epsilon(1.0_dp) &
            * max(1.0_dp, sum(abs(stiffness))), 'zero-curl energy roundoff control')
        mass = 0.0_dp
        call assemble_compatible_physical_mass_surface(fields, 1.0_dp, &
            [1], [0], [1], 1, h, e, parallel, 1.0_dp, route, mass, info)
        call require(info == 0, 'cylinder physical mass assembly')
        image = matmul(mass, vector)
        kinetic = dot_product(vector, image)
        amplitude_a = normal / (2.0_dp * sqrt(s))
        amplitude_b = 2.0_dp * sqrt(s) * eta
        call near(kinetic, 0.5_dp * pi * (amplitude_a**2 + amplitude_b**2), &
            3.0e-13_dp, 'finite Cartesian physical mass')
        call require(kinetic > 0.0_dp, 'zero-curl trial has nonpositive physical mass')
    end subroutine cylinder_surface

    subroutine monomial_jet(coefficients, alpha, s, jet)
        real(dp), intent(in) :: coefficients(0:4), alpha, s
        real(dp), intent(out) :: jet(0:4)
        real(dp) :: factor, exponent
        integer :: j, order

        jet = 0.0_dp
        do j = 0, 4
            exponent = real(j, dp) - alpha
            factor = coefficients(j)
            do order = 0, 4
                jet(order) = jet(order) + factor * s**(exponent - real(order, dp))
                factor = factor * (exponent - real(order, dp))
            end do
        end do
    end subroutine monomial_jet

    integer function binomial(n, k) result(value)
        integer, intent(in) :: n, k
        integer :: j

        value = 1
        do j = 1, k
            value = value * (n - j + 1) / j
        end do
    end function binomial

    subroutine near(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected))) return
        write (error_unit, '(a,2(1x,es24.16))') message, actual, expected
        error stop 1
    end subroutine near

    subroutine near_vector(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual(:), expected(:), tolerance
        character(len=*), intent(in) :: message
        integer :: index

        call require(size(actual) == size(expected), 'oracle vector shape')
        do index = 1, size(actual)
            call near(actual(index), expected(index), tolerance, message)
        end do
    end subroutine near_vector

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, '(a)') message
        error stop 1
    end subroutine require

end program test_weighted_feec_compatibility
