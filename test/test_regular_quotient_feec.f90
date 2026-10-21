! Independent cylinder controls for the regular-quotient displacement space.
program test_regular_quotient_feec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use radial_feec_complex, only: radial_feec_complex_t, &
        build_radial_feec_complex, evaluate_radial_feec_complex
    use compatible_problem_assembly_support, only: apply_stored_power, &
        apply_tangential_axis_weight, apply_regular_normal_lift
    use compatible_family_point_assembly, only: &
        assemble_compatible_direct_surface, assemble_compatible_transformed_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_perpendicular_mass_surface
    use compatible_two_component_problem, only: evaluate_compatible_two_component_vector
    use compatible_radial_quadrature, only: accurate_nodes, accurate_weights
    use symmetric_eigensolver, only: solve_symmetric_generalized
    use gliss_weighted_feec_generated, only: generated_regular_eta_lift, &
        generated_regular_eta_lift_jvp, generated_regular_eta_lift_vjp
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: errors(3), value, rate, root, lo, hi
    integer :: degree, level, iteration, m, parity

    call check_lift_derivatives()
    do degree = 1, 4
        do m = 1, 3
            do parity = 1, 2
                call check_zero_curl(degree, m, parity)
            end do
        end do
    end do
    ! The axisymmetric fast cylinder satisfies the Bessel J1 equation;
    ! xi_r(1)=0 selects its first positive root, bracketed independently.
    lo = 3.0_dp
    hi = 4.0_dp
    do iteration = 1, 55
        root = 0.5_dp * (lo + hi)
        if (bessel_j1(root) > 0.0_dp) then
            lo = root
        else
            hi = root
        end if
    end do
    root = 0.5_dp * (lo + hi)
    do degree = 1, 4
        do level = 1, 3
            call cylinder_eigenvalue(degree, 2**level, value)
            errors(level) = abs(value - root**2)
        end do
        rate = log(errors(1) / errors(2)) / log(2.0_dp)
        write (*, '(a,i0,a,3es13.5,a,f8.3)') 'cylinder p=', degree, &
            ' errors ', errors, ' rate ', rate
        ! Use the coarse pair, before high-degree roundoff dominates.
        if (rate < 2.0_dp * degree - 0.6_dp) &
            error stop 'cylinder eigenvalue lost its expected convergence order'
        if (errors(3) > errors(2)) error stop 'cylinder refinement regressed'
    end do
    write (*, '(a)') 'regular quotient FEEC: PASS'

contains

    subroutine check_lift_derivatives()
        real(dp) :: primal, tangent, reverse, gradient(5), seed(5), expected(5)
        real(dp), parameter :: s = 0.37_dp, flux = 1.2_dp, normal = 0.4_dp
        real(dp), parameter :: sign = -1.0_dp, cotangent = 0.7_dp

        seed = [0.3_dp, -0.2_dp, 0.1_dp, 0.5_dp, -0.1_dp]
        call generated_regular_eta_lift(s, flux, sign, normal, pi, primal)
        expected = [sign*flux*normal/(4*pi*s**2), &
            -sign*normal/(4*pi*s), -flux*normal/(4*pi*s), &
            -sign*flux/(4*pi*s), sign*flux*normal/(4*pi**2*s)]
        call generated_regular_eta_lift_jvp(s, flux, sign, normal, pi, &
            seed(1), seed(2), seed(3), seed(4), seed(5), reverse, tangent)
        call near(reverse, primal, 1.0e-14_dp, 'lift JVP primal')
        call near(tangent, dot_product(expected, seed), 1.0e-14_dp, 'lift JVP')
        call generated_regular_eta_lift_vjp(s, flux, sign, normal, pi, cotangent, &
            reverse, gradient(1), gradient(2), gradient(3), gradient(4), gradient(5))
        if (maxval(abs(gradient - cotangent*expected)) > 1.0e-14_dp) &
            error stop 'lift VJP differs from independent rational derivative'
    end subroutine check_lift_derivatives

    subroutine check_zero_curl(degree, m, parity)
        integer, intent(in) :: degree, m, parity
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: h(:), dh(:), l(:), hv(:, :), dhv(:, :), lv(:, :)
        real(dp), allocatable :: lift(:, :), vector(:), coef(:), companion(:)
        real(dp), allocatable :: direct(:, :), transformed(:, :), mass(:, :), image(:)
        real(dp), allocatable :: normal_samples(:, :), eta_samples(:, :)
        real(dp) :: coordinates(7), flux_slopes(7), powers(1), stored_power(1)
        integer :: mode(1), phases(1)
        integer, allocatable :: indices(:)
        real(dp) :: fields(16, 2, 13), drive(16, 2), s, sigma, energy, kinetic
        real(dp) :: x, eta, a, b, expected
        integer :: info, i, n, point

        call build_radial_feec_complex([0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp], &
            degree, .false., .true., complex, info)
        if (info /= 0) error stop 'regular quotient complex failed'
        n = complex%h1_dofs
        allocate (coef(n), indices(n), hv(n, 1), dhv(n, 1), lv(n, 1), lift(n, 1))
        allocate (image(2*n))
        allocate (vector(2*n), direct(2*n, 2*n), transformed(2*n, 2*n), mass(2*n, 2*n))
        do i = 1, n
            indices(i) = i
            coef(i) = sin(real(i, dp))
        end do
        sigma = 1.0_dp
        if (parity == 2) sigma = -1.0_dp
        ! Ordinary FEEC D gives E=-sigma*T*U'/(2*pi*m), T=pi.
        companion = -sigma * matmul(complex%derivative, coef) / (2.0_dp*m)
        mode(1) = m
        phases(1) = parity
        powers(1) = -0.5_dp*m
        stored_power(1) = 1.0_dp - 0.5_dp*m
        flux_slopes = pi
        do point = 1, size(coordinates)
            coordinates(point) = 0.017_dp + 0.13_dp*point
        end do
        vector(:n) = coef
        vector(n + 1:) = companion
        call evaluate_compatible_two_component_vector(3, mode, [0], stored_power, &
            stored_power, parity, degree, coordinates, vector, normal_samples, &
            eta_samples, info, radial_breaks=[0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp], &
            regular_quotient=.true., flux_slope=flux_slopes)
        if (info /= 0) error stop 'public quotient vector evaluation failed'
        do point = 1, 7
            s = coordinates(point)
            call evaluate_radial_feec_complex(complex, s, h, dh, l, info)
            if (info /= 0) error stop 'regular quotient evaluation failed'
            call apply_stored_power(s, powers, h, dh, indices, hv, dhv, info)
            if (info /= 0) error stop 'regular quotient normal failed'
            call apply_tangential_axis_weight(s, powers, l, indices, lv, info)
            if (info /= 0) error stop 'regular quotient companion failed'
            call apply_regular_normal_lift(s, pi, mode, phases, hv, lift, info)
            if (info /= 0) error stop 'regular quotient lift failed'
            call cylinder_fields(s, fields)
            drive = 0.0_dp
            direct = 0.0_dp
            transformed = 0.0_dp
            call assemble_compatible_direct_surface(fields, drive, mode, [0], &
                phases, 1, hv, dhv, lv, direct, info, normal_eta=lift)
            if (info /= 0) error stop 'regular quotient direct energy failed'
            call assemble_compatible_transformed_surface(fields, drive, mode, [0], &
                phases, 1, hv, dhv, lv, transformed, info, normal_eta=lift)
            if (info /= 0) error stop 'regular quotient transformed energy failed'
            if (maxval(abs(direct - transformed)) > &
                2.0e-12_dp * max(1.0_dp, maxval(abs(direct)))) &
                error stop 'lifted phase routes disagree'
            image = matmul(direct, vector)
            energy = dot_product(vector, image)
            if (abs(energy) > 256*epsilon(1.0_dp)*sum(abs(direct))*sum(vector**2)) &
                error stop 'regular quotient zero-curl trial has nonzero energy'
            mass = 0.0_dp
            call assemble_compatible_perpendicular_mass_surface(fields, 1.0_dp, &
                mode, [0], phases, 1, hv, lv, 1.0_dp, 1, mass, info, lift)
            if (info /= 0) error stop 'regular quotient physical mass failed'
            image = matmul(mass, vector)
            kinetic = dot_product(vector, image)
            x = s**(0.5_dp*m)*dot_product(h, coef)
            eta = -sigma*0.5_dp/m*s**(0.5_dp*m) &
                * (dot_product(dh, coef) + 0.5_dp*m*dot_product(h, coef)/s)
            call near(normal_samples(point, 1), x, 2.0e-13_dp, 'evaluated normal')
            call near(eta_samples(point, 1), eta, 2.0e-13_dp, 'evaluated eta')
            a = x/(2*sqrt(s))
            b = 2*sqrt(s)*eta
            expected = 0.5_dp*pi*(a*a + b*b)
            call near(kinetic, expected, 3.0e-12_dp, 'Cartesian quotient mass')
            mass = 0.0_dp
            call assemble_compatible_perpendicular_mass_surface(fields, 1.0_dp, &
                mode, [0], phases, 1, hv, lv, 1.0_dp, 2, mass, info, lift)
            if (info /= 0) error stop 'regular quotient transformed mass failed'
            image = matmul(mass, vector)
            call near(dot_product(vector, image), expected, &
                3.0e-12_dp, 'transformed Cartesian quotient mass')
        end do
    end subroutine check_zero_curl

    subroutine cylinder_eigenvalue(degree, cells, eigenvalue)
        integer, intent(in) :: degree, cells
        real(dp), intent(out) :: eigenvalue
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: breaks(:), h(:), dh(:), l(:), hv(:, :), dhv(:, :)
        real(dp), allocatable :: lv(:, :), stiffness(:, :), mass(:, :), k(:, :), b(:, :)
        real(dp), allocatable :: eigenvalues(:), vectors(:, :)
        integer, allocatable :: indices(:)
        real(dp) :: fields(16, 2, 13), drive(16, 2), s, weight
        integer :: i, j, n, info

        allocate (breaks(cells + 1))
        do i = 0, cells
            breaks(i + 1) = real(i, dp)/cells
        end do
        call build_radial_feec_complex(breaks, degree, .false., .true., complex, info)
        if (info /= 0) error stop 'cylinder complex failed'
        n = complex%h1_dofs
        allocate (indices(n), hv(n, 1), dhv(n, 1), lv(n, 1))
        do i = 1, n
            indices(i) = i
        end do
        allocate (stiffness(n, n), mass(n, n), source=0.0_dp)
        allocate (k(2*n, 2*n), b(2*n, 2*n))
        drive = 0.0_dp
        do i = 1, cells
            do j = 1, size(accurate_nodes)
                s = (real(i, dp) - 0.5_dp + 0.5_dp*accurate_nodes(j))/cells
                weight = 0.5_dp*accurate_weights(j)/cells
                call evaluate_radial_feec_complex(complex, s, h, dh, l, info)
                if (info /= 0) error stop 'cylinder basis failed'
                call apply_stored_power(s, [-1.0_dp], h, dh, indices, hv, dhv, info)
                if (info /= 0) error stop 'cylinder normal failed'
                lv(:, 1) = l
                call cylinder_fields(s, fields)
                k = 0.0_dp
                b = 0.0_dp
                call assemble_compatible_transformed_surface(fields, drive, [0], [0], &
                    [1], 1, hv, dhv, lv, k, info)
                if (info /= 0) error stop 'cylinder stiffness failed'
                call assemble_compatible_perpendicular_mass_surface(fields, 1.0_dp, &
                    [0], [0], [1], 1, hv, lv, 1.0_dp, 2, b, info)
                if (info /= 0) error stop 'cylinder mass failed'
                stiffness = stiffness + weight*k(:n, :n)
                mass = mass + weight*b(:n, :n)
            end do
        end do
        call solve_symmetric_generalized(stiffness, mass, eigenvalues, vectors, info)
        if (info /= 0) error stop 'cylinder eigensolve failed'
        eigenvalue = eigenvalues(1)
    end subroutine cylinder_eigenvalue

    subroutine cylinder_fields(s, fields)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: fields(:, :, :)

        fields = 0.0_dp
        fields(:, :, 1) = pi
        fields(:, :, 5) = 1.0_dp
        fields(:, :, 7) = pi
        fields(:, :, 8) = 1.0_dp
        fields(:, :, 9) = 4.0_dp*s
    end subroutine cylinder_fields

    subroutine near(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: label

        if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) then
            write (*, *) label, actual, expected
            error stop 'regular quotient oracle differs'
        end if
    end subroutine near
end program test_regular_quotient_feec
