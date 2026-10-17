! Local Newcomb coefficients at a simple rational surface of an exact pinch.
! This is not a global finite-m stability or physical-mass eigenmode test.
program test_sheared_pinch_newcomb
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use phase_assembly_policy, only: phase_assembly_direct, phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), radius = 1.1_dp
    real(dp), parameter :: pressure_edge = 17000.0_dp, mass_budget = 1.0e-11_dp
    real(dp), parameter :: steps(3) = [1.0_dp / 16, 1.0_dp / 32, 1.0_dp / 64]
    real(dp), parameter :: profiles(3) = &
        [-1.8_dp, 4.0_dp - 4.0_dp * sqrt(2.0_dp), -1.0_dp]
    integer, parameter :: phases(2) = [phase_assembly_direct, phase_assembly_transformed]
    real(dp), parameter :: gammas(2) = [1.0_dp, 5.0_dp / 3.0_dp]
    real(dp), parameter :: matrix_budget = 1.0e-10_dp, coefficient_budget = 1.0e-8_dp
    real(dp), parameter :: condition_cap = 1.0e8_dp, local_budget = 1.0e-4_dp
    real(dp), parameter :: discriminant_budget = 1.0e-3_dp
    real(dp) :: plus(3), minus(3), alpha(3), beta(3), exact_alpha, exact_beta
    real(dp) :: c, h, alpha_r, beta_r, delta_r, exact_delta, ratio, gamma
    real(dp) :: worst_matrix, worst_coefficient, worst_condition, worst_local
    real(dp) :: worst_delta, worst_mass
    integer :: profile, level, parity, phase, handedness, gamma_index, cases, phase_index

    worst_matrix = 0.0_dp
    worst_coefficient = 0.0_dp
    worst_condition = 0.0_dp
    worst_local = 0.0_dp
    worst_delta = 0.0_dp
    worst_mass = 0.0_dp
    cases = 0
    do profile = 1, size(profiles)
        c = profiles(profile)
        exact_alpha = 0.5_dp * c**2
        exact_beta = -(2.0_dp + c)
        exact_delta = 1.0_dp + 4.0_dp * exact_beta / exact_alpha
        do gamma_index = 1, size(gammas)
            gamma = gammas(gamma_index)
            do handedness = -1, 1, 2
                do phase_index = 1, size(phases)
                    phase = phases(phase_index)
                    do parity = 1, 2
                        do level = 1, size(steps)
                            h = steps(level)
                            call point_coefficients(1.0_dp + h, c, gamma, handedness, &
                                phase, parity, plus)
                            call point_coefficients(1.0_dp - h, c, gamma, handedness, &
                                phase, parity, minus)
                            alpha(level) = ((1.0_dp + h) * plus(1) &
                                + (1.0_dp - h) * minus(1)) / (2.0_dp * h**2)
                            beta(level) = ((1.0_dp + h) * plus(3) &
                                + (1.0_dp - h) * minus(3)) / 2.0_dp &
                                - ((1.0_dp + h) * plus(2) &
                                - (1.0_dp - h) * minus(2)) / (2.0_dp * h)
                        end do
                        do level = 1, 2
                            ratio = abs((alpha(level) - exact_alpha) &
                                / (alpha(level + 1) - exact_alpha))
                            call require(ratio > 3.5_dp .and. ratio < 4.5_dp, &
                                'alpha central stencil lost second-order convergence')
                            ratio = abs((beta(level) - exact_beta) &
                                / (beta(level + 1) - exact_beta))
                            call require(ratio > 3.5_dp .and. ratio < 4.5_dp, &
                                'beta central stencil lost second-order convergence')
                        end do
                        alpha_r = (4.0_dp * alpha(3) - alpha(2)) / 3.0_dp
                        beta_r = (4.0_dp * beta(3) - beta(2)) / 3.0_dp
                        delta_r = 1.0_dp + 4.0_dp * beta_r / alpha_r
                        worst_local = max(worst_local, &
                            abs(alpha_r / exact_alpha - 1.0_dp), &
                            abs(beta_r - exact_beta) / max(1.0_dp, abs(exact_beta)))
                        worst_delta = max(worst_delta, abs(delta_r - exact_delta))
                        call require(abs(alpha_r / exact_alpha - 1.0_dp) &
                            < local_budget, &
                            'local Newcomb alpha differs')
                        call require(abs(beta_r - exact_beta) &
                            < local_budget * max(1.0_dp, abs(exact_beta)), &
                            'local Newcomb beta differs')
                        call require(abs(delta_r - exact_delta) < discriminant_budget, &
                            'local Frobenius discriminant differs')
                        if (profile == 1) call require(delta_r > 0.0_dp, &
                            'positive local discriminant control failed')
                        if (profile == 3) call require(delta_r < 0.0_dp, &
                            'negative local discriminant control failed')
                        cases = cases + 1
                        write (*, '(a,5i3,3es25.16)') 'LOCAL ', profile, gamma_index, &
                            handedness, phase, parity, alpha_r, beta_r, delta_r
                    end do
                end do
            end do
        end do
    end do
    write (*, '(a,i0)') 'local stencil cases: ', cases
    write (*, '(a,6es14.5)') 'worst matrix/coefficient/condition/local/delta/mass: ', &
        worst_matrix, worst_coefficient, worst_condition, worst_local, worst_delta, &
        worst_mass
    call require(cases == 48, 'local stencil case selection is empty or incomplete')
    write (*, '(a)') 'sheared pinch local Newcomb threshold: PASS'

contains

    subroutine point_coefficients(r, c, gamma, handedness, phase, parity, result)
        real(dp), intent(in) :: r, c, gamma
        integer, intent(in) :: handedness, phase, parity
        real(dp), intent(out) :: result(3)
        real(dp) :: fields(8, 8, 13), drive(8, 8), jr(8, 8), zero(8, 8)
        real(dp) :: gamma_p(8, 8), normal(2, 1), radial(2, 1), tangent(1, 1)
        real(dp) :: stiffness(4, 4), physical(4, 4), transform(4, 4), reference(4, 4)
        real(dp) :: mass(4, 4), physical_mass(4, 4), identity(4, 4)
        real(dp) :: unit(4), other(4), combined(4), diagonal(4), schur(2, 2), y(2, 2)
        real(dp) :: z, zp, b, b2, length, jacobian, pressure, flux_norm, sign_phase
        real(dp) :: determinant, largest, condition, l11, l21, l22, fieldline, wave2
        real(dp) :: expected(3), cross_field, error, scale, raw, numerical
        integer :: modes(1), toroidal(1), kinds(1), status, i, j, k, l

        z = sqrt(1.0_dp + c * (r**2 - 1.0_dp))
        zp = c * r / z
        b = r
        b2 = b**2 + z**2
        length = 2.0_dp * pi
        jacobian = real(handedness, dp) * 2.0_dp * pi * length * r
        pressure = pressure_edge + (2.0_dp + c) * (radius**2 - r**2) &
            / (2.0_dp * vacuum_permeability)
        call require(pressure > 0.0_dp .and. z > 0.0_dp, 'equilibrium domain failed')
        fields(:, :, 1) = 2.0_dp * pi * r * z
        fields(:, :, 2) = real(handedness, dp) * length * b
        fields(:, :, 3) = 2.0_dp * pi * (z + r * zp)
        fields(:, :, 4) = real(handedness, dp) * length
        fields(:, :, 5) = real(handedness, dp) * length * z
        fields(:, :, 6) = 2.0_dp * pi * r * b
        fields(:, :, 7) = jacobian
        fields(:, :, 8) = sqrt(b2)
        fields(:, :, 9) = 1.0_dp
        fields(:, :, 10) = 2.0_dp * z - b * zp
        fields(:, :, 11) = -(2.0_dp + c) * r
        fields(:, :, 12:13) = 0.0_dp
        drive = 4.0_dp
        jr = real(handedness, dp) * 2.0_dp * pi * length
        zero = 0.0_dp
        gamma_p = gamma * pressure
        normal(:, 1) = [1.0_dp, 0.0_dp]
        radial(:, 1) = [0.0_dp, 1.0_dp]
        tangent = 1.0_dp
        modes = 2
        toroidal = 2 * handedness
        kinds = parity
        stiffness = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, drive, jr, &
            zero, zero, gamma_p, modes, toroidal, kinds, 1, normal, radial, tangent, &
            tangent, 1.0_dp, phase, stiffness, status)
        call require(status == 0, 'production surface assembly failed')
        sign_phase = -1.0_dp
        if (parity == 2) sign_phase = 1.0_dp
        flux_norm = fields(1, 1, 1)**2 + fields(1, 1, 2)**2
        transform = 0.0_dp
        transform(1, 1) = 1.0_dp
        transform(2, 2) = 1.0_dp
        transform(3, 3) = sign_phase * z
        transform(3, 4) = -sign_phase * b
        transform(4, 4) = sign_phase * flux_norm / z
        physical = 0.0_dp
        do j = 1, 4
            do i = 1, 4
                do l = 1, 4
                    do k = 1, 4
                        physical(i, j) = physical(i, j) + transform(k, i) &
                            * stiffness(k, l) * transform(l, j)
                    end do
                end do
            end do
        end do
        physical = physical * (2.0_dp * vacuum_permeability / abs(jacobian))
        mass = 0.0_dp
        call assemble_compatible_physical_mass_surface(fields, 1.0_dp, modes, &
            toroidal, kinds, 1, normal, tangent, tangent, 1.0_dp, phase, mass, status)
        call require(status == 0, 'production point mass assembly failed')
        physical_mass = 0.0_dp
        do j = 1, 4
            do i = 1, 4
                do l = 1, 4
                    do k = 1, 4
                        physical_mass(i, j) = physical_mass(i, j) &
                            + transform(k, i) * mass(k, l) * transform(l, j)
                    end do
                end do
            end do
        end do
        identity = 0.0_dp
        identity(1, 1) = 1.0_dp
        identity(3, 3) = 1.0_dp
        identity(4, 4) = 1.0_dp
        error = maxval(abs(2.0_dp * physical_mass / abs(jacobian) - identity))
        worst_mass = max(worst_mass, error)
        call require(error < mass_budget, 'physical point mass transformation differs')
        do i = 1, 4
            unit = 0.0_dp
            unit(i) = 1.0_dp
            diagonal(i) = bernstein(unit, r, z, zp, c, gamma * pressure)
        end do
        do j = 1, 4
            do i = 1, 4
                unit = 0.0_dp
                other = 0.0_dp
                unit(i) = 1.0_dp
                other(j) = 1.0_dp
                combined = unit + other
                raw = bernstein(combined, r, z, zp, c, gamma * pressure)
                reference(i, j) = 0.5_dp * (raw - diagonal(i) - diagonal(j))
            end do
        end do
        scale = max(1.0_dp, maxval(abs(reference)))
        error = maxval(abs(physical - reference)) / scale
        worst_matrix = max(worst_matrix, error)
        call require(error < matrix_budget, 'full raw Bernstein matrix differs')
        determinant = physical(3, 3) * physical(4, 4) - physical(3, 4)**2
        call require(physical(3, 3) > 0.0_dp .and. determinant > 0.0_dp, &
            'physical tangential block is not positive definite')
        largest = 0.5_dp * (physical(3, 3) + physical(4, 4) &
            + sqrt((physical(3, 3) - physical(4, 4))**2 + 4.0_dp * physical(3, 4)**2))
        condition = largest**2 / determinant
        worst_condition = max(worst_condition, condition)
        call require(condition < condition_cap, &
            'physical tangential block is ill conditioned')
        l11 = sqrt(physical(3, 3))
        l21 = physical(4, 3) / l11
        l22 = sqrt(physical(4, 4) - l21**2)
        do j = 1, 2
            y(1, j) = physical(3, j) / l11
            y(2, j) = (physical(4, j) - l21 * y(1, j)) / l22
        end do
        do j = 1, 2
            do i = 1, 2
                schur(i, j) = physical(i, j) - y(1, i) * y(1, j) - y(2, i) * y(2, j)
            end do
        end do
        result(1) = schur(2, 2)
        result(2) = 0.5_dp * (schur(1, 2) + schur(2, 1))
        result(3) = schur(1, 1)
        fieldline = 2.0_dp * (1.0_dp - z)
        wave2 = 4.0_dp * (1.0_dp + 1.0_dp / r**2)
        cross_field = 2.0_dp * z / r - 2.0_dp * b
        expected(1) = fieldline**2 / wave2
        expected(2) = (4.0_dp * z**2 - (2.0_dp * b / r)**2) / (r * wave2)
        expected(3) = fieldline**2 + (z**2 - b**2) / r**2 - 2.0_dp * b / r &
            - cross_field**2 / (r**2 * wave2)
        numerical = maxval(abs(result - expected))
        worst_coefficient = max(worst_coefficient, numerical)
        call require(numerical < coefficient_budget, &
            'static reduced coefficients differ')
    end subroutine point_coefficients

    function bernstein(x, r, z, zp, c, gamma_pressure) result(energy)
        real(dp), intent(in) :: x(4), r, z, zp, c, gamma_pressure
        real(dp) :: energy, f, df, g, h, tangent, divergence, fieldline
        real(dp) :: qr, qt, qz, current_work

        f = x(1)
        df = x(2)
        g = x(3)
        h = x(4)
        tangent = z * g - r * h
        fieldline = 2.0_dp - 2.0_dp * z
        divergence = df + f / r - 2.0_dp * g / r + 2.0_dp * h
        qr = fieldline * f
        qt = 2.0_dp * tangent - f - r * df
        qz = -z * df - (zp + z / r) * f + 2.0_dp * tangent / r
        current_work = -zp * f * (qz - fieldline * h) &
            + 2.0_dp * f * (fieldline * g - qt)
        energy = qr**2 + qt**2 + qz**2 - current_work &
            - (2.0_dp + c) * r * f * divergence &
            + vacuum_permeability * gamma_pressure * divergence**2
    end function bernstein

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, '(a)') trim(message)
        error stop 1
    end subroutine require

end program test_sheared_pinch_newcomb
