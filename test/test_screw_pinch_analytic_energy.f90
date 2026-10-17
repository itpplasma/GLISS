program test_screw_pinch_analytic_energy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_family_point_assembly, only: &
        assemble_compatible_transformed_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface, &
        assemble_compatible_perpendicular_mass_surface
    use compatible_radial_quadrature, only: build_gauss_legendre, &
        compatible_quadrature_ok
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    use physical_mass_kernel, only: physical_mass_energy
    use three_component_kernel, only: three_component_density_value
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tolerance = 4.0e-12_dp
    real(dp), allocatable :: nodes(:), weights(:)
    real(dp) :: radius, linear_field, axial_field, length, edge_pressure
    real(dp) :: gamma, density, worst_error
    integer, parameter :: phases(2) = [phase_assembly_direct, &
        phase_assembly_transformed]
    integer :: fixture, orientation, phase_index, status

    ! Exact equilibrium: B_theta=A*r, B_z=Z,
    ! mu0*p=mu0*p_edge+A^2*(a^2-r^2). The fixed-wall displacement is
    ! xi_r=f=r*(1-r^2/a^2); it is Cartesian regular at the axis.
    ! Bernstein et al. (1958), eq. 3.13 (no gravity), gives
    ! 2*dW=int[(Q^2-jhat.dot.(Q cross xi)+mu0*p'*f*divxi)/mu0
    !         +gamma*p*(divxi)^2]dV, Q=curl(xi cross B), jhat=curl(B).
    ! Here Q_theta=-A*(f+r*f'), Q_z=-Z*(f'+f/r), jhat_z=2*A.
    ! Direct polynomial integration gives the independent absolute oracles:
    ! 2*dW=2*pi*L*((A^2*a^4+2*Z^2*a^2)/(3*mu0)
    !              +gamma*(2*p_edge*a^2+A^2*a^4/mu0)/3),
    ! int(rho*xi^2)dV=pi*L*rho*a^4/12.
    ! Schwab (1991), eqs. 3.2.1--3.2.3, instead completes the squares;
    ! its drive A_drive=4*A^2 is independently known for this equilibrium.
    ! This checks that formulation through production surface assemblies.
    ! It does not exercise the toroidal primitive-import/axis FEEC path.
    call build_gauss_legendre(16, nodes, weights, status)
    call require(status == compatible_quadrature_ok, 'radial quadrature failed')
    nodes = 0.5_dp * (nodes + 1.0_dp)
    weights = 0.5_dp * weights
    worst_error = 0.0_dp
    do fixture = 1, 6
        radius = 0.7_dp
        linear_field = 0.8_dp
        axial_field = 1.2_dp
        length = 3.3_dp
        edge_pressure = 17000.0_dp
        gamma = 5.0_dp / 3.0_dp
        density = 1.7e-7_dp
        select case (fixture)
        case (2)
            linear_field = -1.1_dp
            axial_field = -0.6_dp
        case (3)
            radius = 1.4_dp
            gamma = 0.0_dp
        case (4)
            linear_field = 0.0_dp
            axial_field = -2.1_dp
        case (5)
            edge_pressure = 0.0_dp
            length = 1.1_dp
            density = 6.8e-7_dp
        case (6)
            linear_field = 8.0e-5_dp
            axial_field = 1.2e-4_dp
            edge_pressure = 1.7e-4_dp
        end select
        do orientation = -1, 1, 2
            do phase_index = 1, size(phases)
                call check_equilibrium(radius, linear_field, axial_field, length, &
                    edge_pressure, gamma, density, orientation, &
                    phases(phase_index), worst_error)
            end do
        end do
    end do
    call check_zero_toroidal_flux()
    write (*, '(a,es12.4)') 'screw-pinch absolute energy/mass worst relative error: ', &
        worst_error
    write (*, '(a)') 'PASS'

contains

    subroutine check_equilibrium(a, blinear, bz, period_length, p_edge, &
            adiabatic_index, rho, handedness, phase_assembly, maximum_error)
        real(dp), intent(in) :: a, blinear, bz, period_length, p_edge
        real(dp), intent(in) :: adiabatic_index, rho
        integer, intent(in) :: handedness, phase_assembly
        real(dp), intent(inout) :: maximum_error
        integer, parameter :: angular_points = 4
        integer, parameter :: mode_m(1) = [0], mode_n(1) = [0], parity(1) = [1]
        real(dp) :: fields(angular_points, angular_points, 13)
        real(dp) :: drive(angular_points, angular_points)
        real(dp) :: jacobian_s(angular_points, angular_points)
        real(dp) :: angular_zero(angular_points, angular_points)
        real(dp) :: gamma_pressure(angular_points, angular_points)
        real(dp) :: h1(1, 1), dh1(1, 1), l2(1, 1)
        real(dp) :: two(2, 2), two_terms(2, 2, 4), full(3, 3)
        real(dp) :: terms(3, 3, 5), mass(3, 3), term_energy(5)
        real(dp) :: r, f, df, bmag, jacobian, p, radial_weight
        real(dp) :: reduced_energy, compressible_energy, kinetic_norm
        real(dp) :: expected_magnetic, expected_fluid, expected_norm
        real(dp) :: raw_energy, raw_magnetic, raw_fluid, raw_norm
        real(dp) :: q_theta, q_z, divergence, current_work, pressure_work
        real(dp) :: pressure_slope, error
        integer :: point, info

        expected_magnetic = 2.0_dp * pi * period_length &
            * (blinear**2 * a**4 + 2.0_dp * bz**2 * a**2) &
            / (3.0_dp * vacuum_permeability)
        expected_fluid = 2.0_dp * pi * period_length * adiabatic_index &
            * (2.0_dp * p_edge * a**2 &
            + blinear**2 * a**4 / vacuum_permeability) / 3.0_dp
        expected_norm = pi * period_length * rho * a**4 / 12.0_dp
        reduced_energy = 0.0_dp
        compressible_energy = 0.0_dp
        kinetic_norm = 0.0_dp
        raw_energy = 0.0_dp
        term_energy = 0.0_dp
        angular_zero = 0.0_dp
        l2 = 0.0_dp
        do point = 1, size(nodes)
            r = a * nodes(point)
            radial_weight = a * weights(point)
            f = r * (1.0_dp - r**2 / a**2)
            df = 1.0_dp - 3.0_dp * r**2 / a**2
            bmag = sqrt((blinear * r)**2 + bz**2)
            jacobian = real(handedness, dp) * 2.0_dp * pi * period_length * r
            p = p_edge + blinear**2 * (a**2 - r**2) / vacuum_permeability
            pressure_slope = -2.0_dp * blinear**2 * r
            ! Chart (r,theta,zeta): theta and zeta are unit turns;
            ! e_zeta=handedness*L*e_z. All currents below are covariant B,
            ! except field 10=curl(B).dot.B; field 11=mu0*dp_Pa/dr.
            fields(:, :, 1) = 2.0_dp * pi * r * bz
            fields(:, :, 2) = real(handedness, dp) * period_length * blinear * r
            fields(:, :, 3) = 2.0_dp * pi * bz
            fields(:, :, 4) = real(handedness, dp) * period_length * blinear
            fields(:, :, 5) = real(handedness, dp) * period_length * bz
            fields(:, :, 6) = 2.0_dp * pi * blinear * r**2
            fields(:, :, 7) = jacobian
            fields(:, :, 8) = bmag
            fields(:, :, 9) = 1.0_dp
            fields(:, :, 10) = 2.0_dp * blinear * bz
            fields(:, :, 11) = pressure_slope
            fields(:, :, 12:13) = 0.0_dp
            drive = 4.0_dp * blinear**2
            jacobian_s = real(handedness, dp) * 2.0_dp * pi * period_length
            gamma_pressure = adiabatic_index * p
            h1 = f
            dh1 = df
            two = 0.0_dp
            two_terms = 0.0_dp
            full = 0.0_dp
            terms = 0.0_dp
            mass = 0.0_dp
            call assemble_compatible_transformed_surface(fields, drive, mode_m, &
                mode_n, parity, 1, h1, dh1, l2, two, info, two_terms)
            call require(info == 0, 'two-component surface assembly failed')
            call assemble_compatible_compressible_stiffness_surface(fields, drive, &
                jacobian_s, angular_zero, angular_zero, gamma_pressure, mode_m, &
                mode_n, parity, 1, h1, dh1, l2, l2, radial_weight, &
                phase_assembly, full, info, terms)
            call require(info == 0, 'compressible surface assembly failed')
            call assemble_compatible_physical_mass_surface(fields, rho, mode_m, &
                mode_n, parity, 1, h1, l2, l2, radial_weight, phase_assembly, &
                mass, info)
            call require(info == 0, 'physical mass surface assembly failed')
            reduced_energy = reduced_energy &
                + radial_weight * two(1, 1) / vacuum_permeability
            compressible_energy = compressible_energy + full(1, 1)
            kinetic_norm = kinetic_norm + mass(1, 1)
            term_energy = term_energy + terms(1, 1, :)
            ! Independent cylindrical vector calculus, using the raw energy
            ! rather than production completed-square coefficients or drive.
            divergence = df + f / r
            q_theta = -blinear * (f + r * df)
            q_z = -bz * divergence
            current_work = 2.0_dp * blinear * f * q_theta
            pressure_work = pressure_slope * f * divergence
            raw_magnetic = (q_theta**2 + q_z**2 + current_work &
                + pressure_work) / vacuum_permeability
            raw_fluid = adiabatic_index * p * divergence**2
            raw_norm = rho * f**2
            raw_energy = raw_energy + radial_weight * abs(jacobian) &
                * (raw_magnetic + raw_fluid)
            error = relative_error(full(1, 1), radial_weight * abs(jacobian) &
                * (raw_magnetic + raw_fluid), expected_magnetic + expected_fluid)
            call require(error < tolerance, &
                'point energy differs from raw Bernstein form')
            error = relative_error(mass(1, 1), radial_weight * abs(jacobian) &
                * raw_norm, expected_norm)
            call require(error < tolerance, 'point mass differs from cylindrical norm')
        end do
        call check_integral(reduced_energy, expected_magnetic, &
            'two-component energy', maximum_error)
        call check_integral(compressible_energy, expected_magnetic + expected_fluid, &
            'compressible energy', maximum_error)
        call check_integral(kinetic_norm, expected_norm, 'physical norm', maximum_error)
        call check_integral(raw_energy, expected_magnetic + expected_fluid, &
            'raw Bernstein integral', maximum_error)
        call check_integral(sum(term_energy(:4)), expected_magnetic, &
            'magnetic term decomposition', maximum_error)
        if (adiabatic_index > 0.0_dp) then
            call check_integral(term_energy(5), expected_fluid, &
                'fluid term', maximum_error)
        else
            call require(abs(term_energy(5)) <= tiny(1.0_dp), &
                'zero gamma retains fluid energy')
        end if
        if (abs(blinear) > 0.0_dp) then
            call require(term_energy(4) < 0.0_dp, 'pinch drive must be destabilizing')
        end if

    end subroutine check_equilibrium

    subroutine check_zero_toroidal_flux()
        integer, parameter :: mode_m(1) = [0], mode_n(1) = [0], parity(1) = [1]
        real(dp), parameter :: sentinel = 2.75_dp
        real(dp) :: fields(1, 1, 13), zero(1, 1), gamma_pressure(1, 1)
        real(dp) :: h1(1, 1), dh1(1, 1), l2(1, 1), full(3, 3)
        real(dp) :: terms(3, 3, 5), mass(3, 3), perpendicular(2, 2), two(2, 2)
        real(dp) :: point_energy, point_mass, relative
        integer :: phase_index, info

        ! Exact constant Cartesian field B=e_theta: T=0, P=1 is valid for
        ! the point kernels and perpendicular form. The three-component
        ! surface variables use nu=mu-(P/T)*J*eta and cannot cover T=0.
        fields = 0.0_dp
        fields(:, :, 2) = 1.0_dp
        fields(:, :, 6:9) = 1.0_dp
        zero = 0.0_dp
        gamma_pressure = 1.0_dp
        h1 = 1.0_dp
        dh1 = 0.3_dp
        l2 = 0.0_dp
        do phase_index = 1, size(phases)
            full = sentinel
            terms = sentinel
            call assemble_compatible_compressible_stiffness_surface(fields, zero, &
                zero, zero, zero, gamma_pressure, mode_m, mode_n, parity, 1, &
                h1, dh1, l2, l2, 1.0_dp, phases(phase_index), full, info, terms)
            call require(info /= 0, 'zero toroidal flux admitted by nu stiffness')
            call require(maxval(abs(full - sentinel)) <= tiny(1.0_dp), &
                'invalid nu stiffness changed its output')
            call require(maxval(abs(terms - sentinel)) <= tiny(1.0_dp), &
                'invalid nu stiffness changed its terms')
            mass = sentinel
            call assemble_compatible_physical_mass_surface(fields, 1.0_dp, mode_m, &
                mode_n, parity, 1, h1, l2, l2, 1.0_dp, phases(phase_index), mass, info)
            call require(info /= 0, 'zero toroidal flux admitted by nu mass')
            call require(maxval(abs(mass - sentinel)) <= tiny(1.0_dp), &
                'invalid nu mass changed its output')
            perpendicular = 0.0_dp
            call assemble_compatible_perpendicular_mass_surface(fields, 1.0_dp, &
                mode_m, mode_n, parity, 1, h1, l2, 1.0_dp, phases(phase_index), &
                perpendicular, info)
            call require(info == 0, 'zero toroidal flux rejected by perpendicular mass')
            call require(all(ieee_is_finite(perpendicular)), &
                'zero toroidal flux produced nonfinite perpendicular mass')
            relative = relative_error(perpendicular(1, 1), 1.0_dp, 1.0_dp)
            call require(relative < tolerance, 'constant-field perpendicular norm is wrong')
        end do
        two = 0.0_dp
        call assemble_compatible_transformed_surface(fields, zero, mode_m, mode_n, &
            parity, 1, h1, dh1, l2, two, info)
        call require(info == 0, 'zero toroidal flux rejected by two-component energy')
        relative = relative_error(two(1, 1), 0.3_dp**2, 1.0_dp)
        call require(relative < tolerance, 'constant-field two-component energy is wrong')
        point_energy = three_component_density_value(0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.3_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.3_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
        relative = relative_error(point_energy, &
            0.3_dp**2 * (1.0_dp / vacuum_permeability + 1.0_dp), abs(point_energy))
        call require(relative < tolerance, 'constant-field point energy is wrong')
        point_mass = 2.0_dp * physical_mass_energy(0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
            [1.0_dp, 0.0_dp, 0.0_dp])
        relative = relative_error(point_mass, 1.0_dp, 1.0_dp)
        call require(relative < tolerance, 'constant-field point mass is wrong')
    end subroutine check_zero_toroidal_flux

    subroutine check_integral(actual, expected, label, maximum_error)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        real(dp), intent(inout) :: maximum_error
        real(dp) :: relative

        relative = relative_error(actual, expected, abs(expected))
        maximum_error = max(maximum_error, relative)
        if (relative >= tolerance) then
            write (error_unit, '(a,1x,3es22.12)') label, actual, expected, relative
            error stop 1
        end if
    end subroutine check_integral

    pure function relative_error(actual, expected, scale) result(error)
        real(dp), intent(in) :: actual, expected, scale
        real(dp) :: error

        error = abs(actual - expected) / scale
        if (.not. ieee_is_finite(error)) error = huge(1.0_dp)
    end function relative_error

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (condition) return
        write (error_unit, '(a)') message
        error stop 1
    end subroutine require
end program test_screw_pinch_analytic_energy
