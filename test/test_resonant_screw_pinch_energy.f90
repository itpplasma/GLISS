program test_resonant_screw_pinch_energy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_family_point_assembly, only: &
        assemble_compatible_transformed_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use compatible_radial_quadrature, only: build_gauss_legendre, &
        compatible_quadrature_ok
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tolerance = 4.0e-12_dp
    real(dp), parameter :: a = 0.7_dp, blinear = 0.8_dp, bz = 1.2_dp
    real(dp), parameter :: p_edge = 17000.0_dp, rho = 1.7e-7_dp
    integer, parameter :: orders(3) = [16, 32, 64]
    integer, parameter :: phases(2) = [phase_assembly_direct, &
        phase_assembly_transformed]
    ! Independent finite-energy instability trial for the exact equilibrium
    ! Btheta=A*r, Bz=Z, p=pe+A^2*(a^2-r^2)/mu0, with physical radius r.
    ! Bernstein et al. (1958), printed p.23, Eq.3.13, no gravity:
    ! https://www.hao.ucar.edu/people/phil-judge/Book/preparations/bernstein58.pdf
    ! 2*dW=int[Q^2/mu0-jhat.(Q cross xi)/mu0+p_r*xi_r*div(xi)
    !         +gamma*p*div(xi)^2]dV, Q=curl(xi cross B), jhat=curl(B).
    ! For xi=exp(i*(m*theta+k*z))*(f, i*g, i*h), h=-A*r*g/Z,
    ! k=-m*A/Z, and D=div(xi), the raw curl gives Q=-B*D. Thus
    ! q=(B^2/mu0+gamma*p)*D^2+2*p_r*f*D. Choosing D=-p_r*f/H,
    ! H=B^2/mu0+gamma*p>0, gives q=-p_r^2*f^2/H<0.
    ! f=a*(r/a)^(m-1)*(1-r^2/a^2)^2 and
    ! g=Z^2*(r*f'+f-r*D)/(m*B^2) realize that divergence. f-g is
    ! r^(m+1) times a smooth function of r^2, f+g is r^(m-1) times
    ! such a function, and h is O(r^m). These are Cartesian regular
    ! (including the m=1 axis tie), and f=f'=g=h=0 at r=a.
    ! For a real nonzero harmonic the angle integral is pi*L, giving
    ! 2*dW=pi*L*int(r*q dr) and N=pi*L*int(r*rho*(f^2+g^2+h^2)dr).
    ! At gamma=1, H=(Z^2+mu0*pe+A^2*a^2)/mu0 is constant, and
    ! 2*dW=-48*pi*L*A^4*a^6/[mu0*(Z^2+mu0*pe+A^2*a^2)
    !                            *(m+1)*(m+2)*(m+3)*(m+4)*(m+5)].
    ! Constants below independently integrate the raw q, its fluid part
    ! gamma*p*D^2, and the physical norm using mpmath 1.3.0 at 80 digits.
    ! Inputs are the exact decimal parameters above, gamma=1 or 5/3,
    ! L=2*pi*abs(Z/A); 22 significant digits are retained in each literal.
    ! The gamma=1 total also agrees with the closed beta-function integral.
    ! This qualifies a prescribed globally resonant negative-energy trial
    ! through production surface assemblies, not an isolated eigenvalue,
    ! continuum growth rate, toroidal primitive import or global FEEC solve.
    real(dp), parameter :: oracle_energy(3, 2) = reshape([ &
        -42646.04906810001839661_dp, &
        -12184.58544802857668475_dp, &
        -4569.21954301071625678_dp, &
        -39059.2601575646109253_dp, &
        -11268.41463231817697414_dp, &
        -4257.800827707277346395_dp], [3, 2])
    real(dp), parameter :: oracle_magnetic(3, 2) = reshape([ &
        -48541.24592092651415739_dp, &
        -13676.71583581543959604_dp, &
        -5072.706730466416884241_dp, &
        -47244.92820305564996218_dp, &
        -13378.47707254683333349_dp, &
        -4980.319315596499051766_dp], [3, 2])
    real(dp), parameter :: oracle_fluid(3, 2) = reshape([ &
        5895.196852826495760777_dp, &
        1492.13038778686291129_dp, &
        503.4871874557006274618_dp, &
        8185.668045491039036879_dp, &
        2110.062440228656359349_dp, &
        722.5184878892217053714_dp], [3, 2])
    real(dp), parameter :: oracle_mass(3, 2) = reshape([ &
        3.162879677164415558192e-7_dp, &
        3.903558731270335941694e-8_dp, &
        9.857497941848811145871e-9_dp, &
        3.154928324905639513341e-7_dp, &
        3.89965516874636840581e-8_dp, &
        9.853299483784979360974e-9_dp], [3, 2])

    real(dp), allocatable :: nodes(:), weights(:)
    real(dp) :: gamma, worst_error, worst_all_error
    integer :: order_index, fixture, m, orientation, phase, parity, status

    worst_error = 0.0_dp
    worst_all_error = 0.0_dp
    do order_index = 1, size(orders)
        call build_gauss_legendre(orders(order_index), nodes, weights, status)
        call require(status == compatible_quadrature_ok, 'quadrature failed')
        nodes = 0.5_dp * (nodes + 1.0_dp)
        weights = 0.5_dp * weights
        do fixture = 1, 2
            gamma = 1.0_dp
            if (fixture == 2) gamma = 5.0_dp / 3.0_dp
            do m = 1, 3
                do orientation = -1, 1, 2
                    do phase = 1, 2
                        do parity = 1, 2
                            call check_trial(orders(order_index), fixture, m, gamma, &
                                orientation, phases(phase), parity, worst_error)
                        end do
                    end do
                end do
            end do
        end do
    end do
    write (*, '(a,es12.4)') 'worst relative error: ', worst_error
    write (*, '(a,es12.4)') 'worst all-check error: ', worst_all_error
    write (*, '(a)') 'PASS'

contains

    subroutine check_trial(order, fixture, m, gamma, handedness, phase, kind, &
            maximum_error)
        integer, intent(in) :: order, fixture, m, handedness, phase, kind
        real(dp), intent(in) :: gamma
        real(dp), intent(inout) :: maximum_error
        integer, parameter :: angular = 8
        real(dp) :: fields(angular, angular, 13), drive(angular, angular)
        real(dp) :: jr(angular, angular), zero(angular, angular)
        real(dp) :: gamma_pressure(angular, angular)
        real(dp) :: h1(1, 1), dh1(1, 1), eta(1, 1), nu(1, 1)
        real(dp) :: two(2, 2), two_terms(2, 2, 4), full(3, 3)
        real(dp) :: terms(3, 3, 5), mass(3, 3)
        real(dp) :: r, u, f, df, g, h, b2, pressure, pr, d, raw_d
        real(dp) :: length, jacobian, flux_norm, sign_phase, kphysical, weight
        real(dp) :: qtheta, qz, raw_mag, raw_fluid, raw_mass
        real(dp) :: native_energy, native_magnetic, native_fluid, native_mass
        real(dp) :: expected_energy, expected_mag, expected_fluid, expected_mass
        real(dp) :: point_scale, closed, beta_factor
        integer :: point, info, mode_m(1), mode_n(1), parity(1), j

        length = 2.0_dp * pi * abs(bz / blinear)
        mode_m = m
        mode_n = handedness * m
        parity = kind
        kphysical = -2.0_dp * pi * real(mode_n(1), dp) &
            / (real(handedness, dp) * length)
        call require(abs(real(m, dp) * blinear + kphysical * bz) &
            < 8.0_dp * epsilon(1.0_dp), 'resonance failed')
        sign_phase = -1.0_dp
        if (kind == 2) sign_phase = 1.0_dp
        zero = 0.0_dp
        native_energy = 0.0_dp
        native_magnetic = 0.0_dp
        native_fluid = 0.0_dp
        native_mass = 0.0_dp
        expected_energy = 0.0_dp
        expected_mag = 0.0_dp
        expected_fluid = 0.0_dp
        expected_mass = 0.0_dp
        do point = 1, size(nodes)
            ! Gauss in area s=(r/a)^2; the surface chart remains physical r.
            r = a * sqrt(nodes(point))
            weight = a**2 * weights(point) / (2.0_dp * r)
            u = r / a
            f = a * u**(m - 1) * (1.0_dp - u**2)**2
            df = -4.0_dp * u**m * (1.0_dp - u**2)
            if (m > 1) df = df + real(m - 1, dp) * u**(m - 2) &
                * (1.0_dp - u**2)**2
            b2 = bz**2 + (blinear * r)**2
            pressure = p_edge + blinear**2 * (a**2 - r**2) &
                / vacuum_permeability
            pr = -2.0_dp * blinear**2 * r / vacuum_permeability
            d = -pr * f / (b2 / vacuum_permeability + gamma * pressure)
            g = bz**2 * (r * df + f - r * d) / (real(m, dp) * b2)
            h = -blinear * r * g / bz
            raw_d = df + f / r - real(m, dp) * g / r - kphysical * h
            call check_error(raw_d, d, 1.0_dp, 'divergence realization')
            call require(abs(blinear * r * g + bz * h) &
                < 8.0_dp * epsilon(1.0_dp), 'perpendicular trial failed')
            jacobian = real(handedness, dp) * 2.0_dp * pi * length * r
            fields(:, :, 1) = 2.0_dp * pi * r * bz
            fields(:, :, 2) = real(handedness, dp) * length * blinear * r
            fields(:, :, 3) = 2.0_dp * pi * bz
            fields(:, :, 4) = real(handedness, dp) * length * blinear
            fields(:, :, 5) = real(handedness, dp) * length * bz
            fields(:, :, 6) = 2.0_dp * pi * blinear * r**2
            fields(:, :, 7) = jacobian
            fields(:, :, 8) = sqrt(b2)
            fields(:, :, 9) = 1.0_dp
            fields(:, :, 10) = 2.0_dp * blinear * bz
            fields(:, :, 11) = vacuum_permeability * pr
            fields(:, :, 12:13) = 0.0_dp
            drive = 4.0_dp * blinear**2
            jr = real(handedness, dp) * 2.0_dp * pi * length
            gamma_pressure = gamma * pressure
            flux_norm = fields(1, 1, 1)**2 + fields(1, 1, 2)**2
            h1 = f
            dh1 = df
            ! eta=T*xi^theta-P*xi^zeta, nu=J*(T^2+P^2)*xi^zeta/T.
            ! Physical tangential components are sign_phase*(g,h) times
            ! the opposite angular parity to the normal component f.
            eta = sign_phase * b2 * g / bz
            nu = sign_phase * flux_norm * h / bz
            two = 0.0_dp
            two_terms = 0.0_dp
            full = 0.0_dp
            terms = 0.0_dp
            mass = 0.0_dp
            call assemble_compatible_transformed_surface(fields, drive, mode_m, &
                mode_n, parity, 1, h1, dh1, eta, two, info, two_terms)
            call require(info == 0, 'two-component assembly failed')
            call assemble_compatible_compressible_stiffness_surface(fields, drive, &
                jr, zero, zero, gamma_pressure, mode_m, mode_n, parity, 1, &
                h1, dh1, eta, nu, weight, phase, full, info, terms)
            call require(info == 0, 'compressible assembly failed')
            call assemble_compatible_physical_mass_surface(fields, rho, mode_m, &
                mode_n, parity, 1, h1, eta, nu, weight, phase, mass, info)
            call require(info == 0, 'mass assembly failed')
            ! Independent raw cylindrical curl of xi cross B, not the
            ! reduced resonant identity or the completed-square drive.
            qtheta = -kphysical * (bz * g - blinear * r * h) &
                - blinear * (f + r * df)
            qz = -bz * (df + f / r) &
                + real(m, dp) * (bz * g - blinear * r * h) / r
            call check_error(qtheta, -blinear * r * raw_d, 1.0_dp, &
                'raw curl theta identity')
            call check_error(qz, -bz * raw_d, 1.0_dp, &
                'raw curl z identity')
            raw_mag = (qtheta**2 + qz**2 + 2.0_dp * blinear * f * qtheta) &
                / vacuum_permeability + pr * f * raw_d
            raw_fluid = gamma * pressure * raw_d**2
            raw_mass = rho * (f**2 + g**2 + h**2)
            point_scale = weight * abs(jacobian) / 2.0_dp
            call check_error(sum(full), point_scale * (raw_mag + raw_fluid), &
                abs(oracle_energy(m, fixture)), 'point raw Bernstein energy')
            call check_error(sum(mass), point_scale * raw_mass, &
                oracle_mass(m, fixture), 'point physical norm')
            call check_error(sum(full), -point_scale * pr**2 * f**2 &
                / (b2 / vacuum_permeability + gamma * pressure), &
                abs(oracle_energy(m, fixture)), 'point optimal negative form')
            native_energy = native_energy + sum(full)
            native_magnetic = native_magnetic + weight * sum(two) &
                / vacuum_permeability
            native_fluid = native_fluid + sum(terms(:, :, 5))
            native_mass = native_mass + sum(mass)
            expected_energy = expected_energy + point_scale * (raw_mag + raw_fluid)
            expected_mag = expected_mag + point_scale * raw_mag
            expected_fluid = expected_fluid + point_scale * raw_fluid
            expected_mass = expected_mass + point_scale * raw_mass
        end do
        call require(native_energy < 0.0_dp, 'strict negative energy missing')
        call check_error(native_energy, expected_energy, &
            abs(oracle_energy(m, fixture)), &
            'integrated raw energy')
        call check_error(native_magnetic, expected_mag, &
            abs(oracle_magnetic(m, fixture)), &
            'integrated raw magnetic')
        call check_error(native_fluid, expected_fluid, oracle_fluid(m, fixture), &
            'integrated raw fluid')
        call check_error(native_mass, expected_mass, oracle_mass(m, fixture), &
            'integrated raw mass')
        call check_error(native_energy, oracle_energy(m, fixture), &
            abs(oracle_energy(m, fixture)), 'independent total oracle')
        call check_error(native_magnetic, oracle_magnetic(m, fixture), &
            abs(oracle_magnetic(m, fixture)), 'independent magnetic oracle')
        call check_error(native_fluid, oracle_fluid(m, fixture), &
            oracle_fluid(m, fixture), 'independent fluid oracle')
        call check_error(native_mass, oracle_mass(m, fixture), &
            oracle_mass(m, fixture), 'independent mass oracle')
        if (fixture == 1) then
            beta_factor = 1.0_dp
            do j = 1, 5
                beta_factor = beta_factor * real(m + j, dp)
            end do
            closed = -48.0_dp * pi * length * blinear**4 * a**6 &
                / (vacuum_permeability * (bz**2 + vacuum_permeability * p_edge &
                + blinear**2 * a**2) * beta_factor)
            call check_error(native_energy, closed, abs(closed), 'closed beta integral')
        end if
        maximum_error = max(maximum_error, abs(native_energy &
            / oracle_energy(m, fixture) - 1.0_dp))
        if (handedness == 1 .and. phase == phase_assembly_transformed .and. kind == 1) &
            write (*, '(3i3,4es25.16)') order, fixture, m, native_energy, &
            native_magnetic, native_fluid, native_mass
    end subroutine check_trial

    subroutine check_error(actual, expected, scale, label)
        real(dp), intent(in) :: actual, expected, scale
        character(*), intent(in) :: label
        real(dp) :: error

        error = abs(actual - expected) / scale
        worst_all_error = max(worst_all_error, error)
        if (.not. ieee_is_finite(error) .or. error >= tolerance) then
            write (error_unit, '(a,4es25.16)') trim(label), actual, expected, &
                scale, error
            error stop 'oracle mismatch'
        end if
    end subroutine check_error

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (.not. condition) then
            write (error_unit, '(a)') trim(label)
            error stop 'resonant pinch check failed'
        end if
    end subroutine require
end program test_resonant_screw_pinch_energy
