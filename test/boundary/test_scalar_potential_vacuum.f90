program test_scalar_potential_vacuum
    ! Exact vacuum energies of a toroidal harmonic (toroidal_harmonic_oracle.py):
    ! Phi = sqrt(cosh eta - cos theta) P^1_{3/2}(cosh eta) cos(2 theta) cos(phi)
    ! outside the torus R0 = 3, a = 1, alone and inside the ideal wall of
    ! radius 1.6 of the same toroidal coordinates.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use scalar_potential_vacuum, only: assemble_exterior_vacuum, &
        edge_triangle_centroids, vacuum_bie_degenerate, vacuum_bie_invalid, &
        vacuum_bie_not_nested, vacuum_bie_ok
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major = 3.0_dp, minor = 1.0_dp, wall_minor = 1.6_dp
    ! Cosine series of the normal flux per d theta d phi on the torus.
    real(dp), parameter :: flux_series(0:23) = [ &
        -1.88127356162503134e+00_dp, -1.13197750697218495e+01_dp, &
        -8.54001446054974309e+01_dp, -1.09949051012997412e+01_dp, &
        -1.88322714647422518e+00_dp, -3.36095720435106493e-01_dp, &
        -6.04825207469112461e-02_dp, -1.08869274940653175e-02_dp, &
        -1.95556507011906099e-03_dp, -3.50315685020544272e-04_dp, &
        -6.25819581908054030e-05_dp, -1.11509099025952337e-05_dp, &
        -1.98215988717460242e-06_dp, -3.51586857129386898e-07_dp, &
        -6.22421185514173505e-08_dp, -1.09996047614881340e-08_dp, &
        -1.94081683983883700e-09_dp, -3.41957182483806789e-10_dp, &
        -6.01721683003024624e-11_dp, -1.05756354454090003e-11_dp, &
        -1.85672801316702056e-12_dp, -3.25656740809755853e-13_dp, &
        -5.70658535772567539e-14_dp, -9.99143072484011584e-15_dp]
    real(dp), parameter :: exterior_energy = 1.25761242286838060e+04_dp
    real(dp), parameter :: annulus_energy = 1.81394433228681228e+04_dp

    call test_exterior_convergence()
    call test_wall_convergence()
    call test_symmetry_and_scaling()
    call test_sector_solution()
    call test_rejections()
    write (*, "(a)") "PASS"

contains

    subroutine test_exterior_convergence()
        real(dp) :: coarse, fine

        coarse = relative_error(12, 36, .false.)
        fine = relative_error(16, 48, .false.)
        ! Flat triangles and constant potentials converge at second order.
        call require(fine < 1.5e-2_dp, "exterior energy misses the harmonic")
        call require(coarse / fine > 1.6_dp, &
            "exterior energy does not converge at second order")
    end subroutine test_exterior_convergence

    subroutine test_wall_convergence()
        real(dp) :: coarse, fine

        coarse = relative_error(12, 36, .true.)
        fine = relative_error(16, 48, .true.)
        call require(fine < 1.5e-2_dp, "walled energy misses the harmonic")
        call require(coarse / fine > 1.6_dp, &
            "walled energy does not converge at second order")
        call require(annulus_energy > exterior_energy, &
            "the oracle wall does not stiffen the vacuum")
    end subroutine test_wall_convergence

    real(dp) function relative_error(nu, nv, walled) result(error)
        integer, intent(in) :: nu, nv
        logical, intent(in) :: walled
        real(dp), allocatable :: energy(:, :), flux(:, :), plasma(:, :, :)
        real(dp), allocatable :: wall(:, :, :)
        integer :: info

        call coordinate_torus(minor, nu, nv, plasma)
        call harmonic_flux(nu, nv, 0.0_dp, flux)
        if (walled) then
            call coordinate_torus(wall_minor, nu, nv, wall)
            call assemble_exterior_vacuum(plasma, flux, energy, info, wall)
            call require(info == vacuum_bie_ok, "walled assembly failed")
            error = abs(energy(1, 1) / annulus_energy - 1.0_dp)
        else
            call assemble_exterior_vacuum(plasma, flux, energy, info)
            call require(info == vacuum_bie_ok, "exterior assembly failed")
            error = abs(energy(1, 1) / exterior_energy - 1.0_dp)
        end if
    end function relative_error

    subroutine test_symmetry_and_scaling()
        integer, parameter :: nu = 12, nv = 36
        real(dp), allocatable :: energy(:, :), flux(:, :), plasma(:, :, :)
        real(dp), allocatable :: rotated(:, :), pair(:, :)
        integer :: info

        call coordinate_torus(minor, nu, nv, plasma)
        call harmonic_flux(nu, nv, 0.0_dp, flux)
        call harmonic_flux(nu, nv, 0.25_dp, rotated)
        allocate (pair(size(flux, 1), 3))
        pair(:, 1) = flux(:, 1)
        pair(:, 2) = rotated(:, 1)
        pair(:, 3) = 2.0_dp * flux(:, 1)
        call assemble_exterior_vacuum(plasma, pair, energy, info)
        call require(info == vacuum_bie_ok, "multi-datum assembly failed")
        call require(all(abs(energy - transpose(energy)) <= 1.0e-12_dp &
            * maxval(abs(energy))), "energy form is not symmetric")
        ! A quarter toroidal period turns cos(phi) into sin(phi): the
        ! axisymmetric vacuum gives equal energies and no coupling.
        call require(abs(energy(2, 2) / energy(1, 1) - 1.0_dp) < 1.0e-10_dp, &
            "toroidally rotated datum changed the energy")
        call require(abs(energy(1, 2)) < 1.0e-10_dp * energy(1, 1), &
            "orthogonal toroidal phases couple")
        call require(abs(energy(3, 3) / energy(1, 1) - 4.0_dp) < 1.0e-12_dp, &
            "energy is not quadratic in the datum")
    end subroutine test_symmetry_and_scaling

    ! The block-circulant solve over rotation sectors reproduces the full
    ! system: the axisymmetric torus with and without a wall (one sector per
    ! toroidal node) and a five-period rotating ellipse, whose data mix the
    ! Fourier harmonics over sectors.
    subroutine test_sector_solution()
        integer, parameter :: nu = 10, nv = 20
        real(dp), allocatable :: dense(:, :), flux(:, :), plasma(:, :, :)
        real(dp), allocatable :: sectors(:, :), wall(:, :, :)
        integer :: info

        call coordinate_torus(minor, nu, nv, plasma)
        call harmonic_flux(nu, nv, 0.1_dp, flux)
        call coordinate_torus(wall_minor, nu, nv, wall)
        call assemble_exterior_vacuum(plasma, flux, sectors, info)
        call require(info == vacuum_bie_ok, "sector assembly failed")
        call assemble_exterior_vacuum(plasma, flux, dense, info, &
            symmetric=.false.)
        call require(info == vacuum_bie_ok, "dense assembly failed")
        call require(agree(sectors, dense), &
            "axisymmetric sector solve differs from the full system")
        call assemble_exterior_vacuum(plasma, flux, sectors, info, wall)
        call require(info == vacuum_bie_ok, "walled sector assembly failed")
        call assemble_exterior_vacuum(plasma, flux, dense, info, wall, &
            symmetric=.false.)
        call require(info == vacuum_bie_ok, "walled dense assembly failed")
        call require(agree(sectors, dense), &
            "walled sector solve differs from the full system")

        call rotating_ellipse(0.0_dp, nu, 4 * nv, plasma)
        call rotating_ellipse(0.4_dp, nu, 4 * nv, wall)
        call mixed_flux(nu, 4 * nv, flux)
        call assemble_exterior_vacuum(plasma, flux, sectors, info, wall)
        call require(info == vacuum_bie_ok, "helical sector assembly failed")
        call assemble_exterior_vacuum(plasma, flux, dense, info, wall, &
            symmetric=.false.)
        call require(info == vacuum_bie_ok, "helical dense assembly failed")
        call require(agree(sectors, dense), &
            "five-period sector solve differs from the full system")
    end subroutine test_sector_solution

    logical function agree(first, second)
        real(dp), intent(in) :: first(:, :), second(:, :)

        agree = all(abs(first - second) <= 1.0e-10_dp * maxval(abs(second)))
    end function agree

    ! Ellipse of semi-axes 1 + offset and 0.6 + offset rotating five times
    ! poloidally per toroidal turn about the circle R = 4.
    subroutine rotating_ellipse(offset, nu, nv, surface)
        real(dp), intent(in) :: offset
        integer, intent(in) :: nu, nv
        real(dp), allocatable, intent(out) :: surface(:, :, :)
        real(dp) :: along, across, phi, r, theta, tilt
        integer :: i, k

        allocate (surface(3, nu, nv))
        do k = 1, nv
            phi = 2.0_dp * pi * real(k - 1, dp) / real(nv, dp)
            tilt = 5.0_dp * phi
            do i = 1, nu
                theta = 2.0_dp * pi * real(i - 1, dp) / real(nu, dp)
                along = (1.0_dp + offset) * cos(theta)
                across = (0.6_dp + offset) * sin(theta)
                r = 4.0_dp + along * cos(tilt) - across * sin(tilt)
                surface(1, i, k) = r * cos(phi)
                surface(2, i, k) = r * sin(phi)
                surface(3, i, k) = along * sin(tilt) + across * cos(tilt)
            end do
        end do
    end subroutine rotating_ellipse

    ! Fluxes with poloidal and toroidal numbers that are not multiples of
    ! the five periods, and one datum with zero net flux per sector pair.
    subroutine mixed_flux(nu, nv, flux)
        integer, intent(in) :: nu, nv
        real(dp), allocatable, intent(out) :: flux(:, :)
        real(dp), allocatable :: centroid(:, :)
        real(dp) :: u, v
        integer :: t

        allocate (centroid(2, 2 * nu * nv), flux(2 * nu * nv, 3))
        call edge_triangle_centroids(nu, nv, centroid)
        do t = 1, size(centroid, 2)
            u = 2.0_dp * pi * centroid(1, t)
            v = 2.0_dp * pi * centroid(2, t)
            flux(t, 1) = cos(u - v)
            flux(t, 2) = sin(2.0_dp * u + 3.0_dp * v)
            flux(t, 3) = cos(u) * sin(5.0_dp * v)
        end do
    end subroutine mixed_flux

    subroutine test_rejections()
        integer, parameter :: nu = 8, nv = 24
        real(dp), allocatable :: energy(:, :), flux(:, :), plasma(:, :, :)
        real(dp), allocatable :: wall(:, :, :)
        integer :: info

        call coordinate_torus(minor, nu, nv, plasma)
        call harmonic_flux(nu, nv, 0.0_dp, flux)
        call assemble_exterior_vacuum(plasma, flux(1:10, :), energy, info)
        call require(info == vacuum_bie_invalid, "short datum was accepted")
        call coordinate_torus(0.8_dp, nu, nv, wall)
        call assemble_exterior_vacuum(plasma, flux, energy, info, wall)
        call require(info == vacuum_bie_not_nested, &
            "a wall inside the plasma was accepted")
        wall = plasma
        wall(:, 2, :) = wall(:, 1, :)
        call assemble_exterior_vacuum(wall, flux, energy, info)
        call require(info == vacuum_bie_degenerate, &
            "a collapsed edge was accepted")
    end subroutine test_rejections

    ! Nodes of the toroidal-coordinate torus through radius r about the same
    ! foci as the plasma, indexed by the toroidal-coordinate angle theta.
    subroutine coordinate_torus(radius, nu, nv, surface)
        real(dp), intent(in) :: radius
        integer, intent(in) :: nu, nv
        real(dp), allocatable, intent(out) :: surface(:, :, :)
        real(dp) :: focus, eta, theta, phi, denominator, r
        integer :: i, k

        focus = sqrt(major**2 - minor**2)
        eta = asinh(focus / radius)
        allocate (surface(3, nu, nv))
        do k = 1, nv
            phi = 2.0_dp * pi * real(k - 1, dp) / real(nv, dp)
            do i = 1, nu
                theta = 2.0_dp * pi * real(i - 1, dp) / real(nu, dp)
                denominator = cosh(eta) - cos(theta)
                r = focus * sinh(eta) / denominator
                surface(1, i, k) = r * cos(phi)
                surface(2, i, k) = r * sin(phi)
                surface(3, i, k) = focus * sin(theta) / denominator
            end do
        end do
    end subroutine coordinate_torus

    ! B.n dA / (du dv) of the harmonic at the triangle centroids, with the
    ! toroidal phase shifted by shift periods.
    subroutine harmonic_flux(nu, nv, shift, flux)
        integer, intent(in) :: nu, nv
        real(dp), intent(in) :: shift
        real(dp), allocatable, intent(out) :: flux(:, :)
        real(dp), allocatable :: centroid(:, :)
        real(dp) :: series, theta
        integer :: k, t

        allocate (centroid(2, 2 * nu * nv), flux(2 * nu * nv, 1))
        call edge_triangle_centroids(nu, nv, centroid)
        do t = 1, size(centroid, 2)
            theta = 2.0_dp * pi * centroid(1, t)
            series = 0.0_dp
            do k = 0, size(flux_series) - 1
                series = series + flux_series(k) * cos(real(k, dp) * theta)
            end do
            flux(t, 1) = (2.0_dp * pi)**2 * series &
                * cos(2.0_dp * pi * (centroid(2, t) + shift))
        end do
    end subroutine harmonic_flux

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: "//message
            error stop 1
        end if
    end subroutine require

end program test_scalar_potential_vacuum
