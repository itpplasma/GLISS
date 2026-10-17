! Independent, exact zero-pressure vacuum equilibrium.  The meridional
! conformal map is W=R+i Z=D(1+q exp(-2 pi i theta))/(1-q exp(-2 pi i theta)).
! Its circles have centre R0 and radius r, and its positive Cartesian
! Jacobian is 4 pi^2 C R^2.  GLISS uses B^zeta=-Phi'/J, so
! Phi'=-2 pi F C gives B=F/R e_phi.
module vacuum_torus_analytic_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_half
    implicit none
    private

    real(dp), parameter, public :: torus_major_radius = 3.0_dp
    real(dp), parameter, public :: torus_minor_radius = 1.0_dp
    real(dp), parameter, public :: torus_field_constant = 3.0_dp
    real(dp), parameter, public :: torus_c = torus_major_radius &
        - sqrt(torus_major_radius**2 - torus_minor_radius**2)
    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: make_vacuum_torus, vacuum_torus_point

contains

    subroutine make_vacuum_torus(surfaces, harmonics, equilibrium)
        integer, intent(in) :: surfaces, harmonics
        type(gvec_cas3d_equilibrium_t), intent(out) :: equilibrium
        real(dp) :: s, d, radius, q
        real(dp) :: samples(256, 10), values(10), angle, factor
        complex(dp) :: w, ws, wt, wtt
        integer :: radial, mode, point, field

        equilibrium%field_periods = 1
        equilibrium%winding = 1
        equilibrium%radial_grid = radial_grid_half
        equilibrium%stellarator_symmetric = .true.
        equilibrium%has_boozer_position_frame = .true.
        equilibrium%has_chart_metric = .true.
        allocate (equilibrium%poloidal_modes(harmonics + 1))
        do mode = 0, harmonics
            equilibrium%poloidal_modes(mode + 1) = mode
        end do
        equilibrium%toroidal_modes = [0]
        allocate (equilibrium%rho(surfaces), equilibrium%s(surfaces), &
            equilibrium%pressure(surfaces), &
            equilibrium%b_theta_average(surfaces), &
            equilibrium%b_zeta_average(surfaces), &
            equilibrium%toroidal_flux(surfaces), &
            equilibrium%poloidal_flux(surfaces), &
            equilibrium%rotational_transform(surfaces))
        equilibrium%pressure = 0.0_dp
        equilibrium%b_theta_average = 0.0_dp
        equilibrium%b_zeta_average = 2.0_dp * pi * torus_field_constant
        equilibrium%poloidal_flux = 0.0_dp
        equilibrium%rotational_transform = 0.0_dp
        call allocate_pair(equilibrium%xhat, surfaces, harmonics)
        call allocate_pair(equilibrium%yhat, surfaces, harmonics)
        call allocate_pair(equilibrium%zhat, surfaces, harmonics)
        call allocate_pair(equilibrium%mod_b, surfaces, harmonics)
        call allocate_pair(equilibrium%jacobian, surfaces, harmonics)
        call allocate_pair(equilibrium%g_tt, surfaces, harmonics)
        call allocate_pair(equilibrium%g_tz, surfaces, harmonics)
        call allocate_pair(equilibrium%g_zz, surfaces, harmonics)
        call allocate_pair(equilibrium%g_st, surfaces, harmonics)
        call allocate_pair(equilibrium%g_sz, surfaces, harmonics)
        call allocate_pair(equilibrium%second_form_tt, surfaces, harmonics)
        call allocate_pair(equilibrium%second_form_tz, surfaces, harmonics)
        call allocate_pair(equilibrium%second_form_zz, surfaces, harmonics)
        call allocate_pair(equilibrium%b_contravariant_theta, surfaces, harmonics)
        call allocate_pair(equilibrium%b_contravariant_zeta, surfaces, harmonics)

        do radial = 1, surfaces
            equilibrium%rho(radial) = sqrt((real(radial, dp) - 0.5_dp) &
                / real(surfaces, dp))
            s = equilibrium%rho(radial)**2
            equilibrium%s(radial) = s
            equilibrium%toroidal_flux(radial) = &
                -2.0_dp * pi * torus_field_constant * torus_c * s
            d = torus_major_radius - torus_c * s
            radius = sqrt(2.0_dp * torus_major_radius * torus_c * s &
                - torus_c**2 * s**2)
            q = radius / (torus_major_radius + d)
            equilibrium%xhat%cosine(radial, 1, 1) = d
            do mode = 1, harmonics
                equilibrium%xhat%cosine(radial, mode + 1, 1) = &
                    2.0_dp * d * q**mode
                equilibrium%zhat%sine(radial, mode + 1, 1) = &
                    -2.0_dp * d * q**mode
            end do
            do point = 1, size(samples, 1)
                call meridional_point(s, real(point - 1, dp) &
                    / real(size(samples, 1), dp), w, ws, wt, wtt)
                call exported_fields(w, ws, wt, wtt, values)
                samples(point, :) = values
            end do
            do mode = 0, harmonics
                values = 0.0_dp
                factor = 2.0_dp / real(size(samples, 1), dp)
                if (mode == 0) factor = 0.5_dp * factor
                do point = 1, size(samples, 1)
                    angle = 2.0_dp * pi * real(mode * (point - 1), dp) &
                        / real(size(samples, 1), dp)
                    do field = 1, 9
                        values(field) = values(field) &
                            + factor * cos(angle) * samples(point, field)
                    end do
                    values(10) = values(10) &
                        + factor * sin(angle) * samples(point, 10)
                end do
                call store_fields(equilibrium, radial, mode + 1, values)
            end do
        end do
    end subroutine make_vacuum_torus

    subroutine allocate_pair(pair, surfaces, harmonics)
        type(harmonic_pair_t), intent(out) :: pair
        integer, intent(in) :: surfaces, harmonics

        allocate (pair%cosine(surfaces, harmonics + 1, 1), source=0.0_dp)
        allocate (pair%sine(surfaces, harmonics + 1, 1), source=0.0_dp)
    end subroutine allocate_pair

    subroutine meridional_point(s, theta, w, ws, wt, wtt)
        real(dp), intent(in) :: s, theta
        complex(dp), intent(out) :: w, ws, wt, wtt
        real(dp) :: d, q, logarithmic_slope
        complex(dp) :: v, omega

        d = torus_major_radius - torus_c * s
        q = sqrt(torus_c * s / (2.0_dp * torus_major_radius - torus_c * s))
        v = q * exp(cmplx(0.0_dp, -2.0_dp * pi * theta, dp))
        omega = cmplx(0.0_dp, -2.0_dp * pi, dp)
        logarithmic_slope = torus_major_radius &
            / (s * (2.0_dp * torus_major_radius - torus_c * s))
        w = d * (1.0_dp + v) / (1.0_dp - v)
        ws = -torus_c * (1.0_dp + v) / (1.0_dp - v) &
            + 2.0_dp * d * logarithmic_slope * v / (1.0_dp - v)**2
        wt = 2.0_dp * d * omega * v / (1.0_dp - v)**2
        wtt = 2.0_dp * d * omega**2 * v * (1.0_dp + v) &
            / (1.0_dp - v)**3
    end subroutine meridional_point

    subroutine exported_fields(w, ws, wt, wtt, values)
        complex(dp), intent(in) :: w, ws, wt, wtt
        real(dp), intent(out) :: values(10)
        real(dp) :: r, normal_r, normal_z

        r = real(w, dp)
        normal_r = -aimag(wt) / abs(wt)
        normal_z = real(wt, dp) / abs(wt)
        values(1) = torus_field_constant / r
        values(2) = 4.0_dp * pi**2 * torus_c * r**2
        values(3) = abs(wt)**2
        values(4) = 4.0_dp * pi**2 * r**2
        values(5) = normal_r * real(wtt, dp) + normal_z * aimag(wtt)
        values(6) = -4.0_dp * pi**2 * r * normal_r
        values(7) = torus_field_constant / (2.0_dp * pi * r**2)
        values(8:9) = 0.0_dp
        values(10) = real(conjg(ws) * wt, dp)
    end subroutine exported_fields

    subroutine store_fields(equilibrium, radial, mode, values)
        type(gvec_cas3d_equilibrium_t), intent(inout) :: equilibrium
        integer, intent(in) :: radial, mode
        real(dp), intent(in) :: values(10)

        equilibrium%mod_b%cosine(radial, mode, 1) = values(1)
        equilibrium%jacobian%cosine(radial, mode, 1) = values(2)
        equilibrium%g_tt%cosine(radial, mode, 1) = values(3)
        equilibrium%g_zz%cosine(radial, mode, 1) = values(4)
        equilibrium%second_form_tt%cosine(radial, mode, 1) = values(5)
        equilibrium%second_form_zz%cosine(radial, mode, 1) = values(6)
        equilibrium%b_contravariant_zeta%cosine(radial, mode, 1) = values(7)
        equilibrium%g_st%sine(radial, mode, 1) = values(10)
    end subroutine store_fields

    subroutine vacuum_torus_point(s, theta, zeta, position, jacobian, field)
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: position(3), jacobian, field(3)
        complex(dp) :: w, ws, wt, wtt
        real(dp) :: phi, r

        call meridional_point(s, theta, w, ws, wt, wtt)
        r = real(w, dp)
        phi = 2.0_dp * pi * zeta
        position(1) = r * cos(phi)
        position(2) = r * sin(phi)
        position(3) = aimag(w)
        jacobian = 4.0_dp * pi**2 * torus_c * r**2
        field(1) = -torus_field_constant / r * sin(phi)
        field(2) = torus_field_constant / r * cos(phi)
        field(3) = 0.0_dp
    end subroutine vacuum_torus_point

end module vacuum_torus_analytic_fixture
