program test_fixed_boundary_symmetry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cylinder_fixture, only: create_cylinder_fixture
    use export_surface_geometry, only: build_angular_grids
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_ok, fixed_boundary_problem_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        evaluate_primitive_equilibrium, primitive_equilibrium_spline_t, &
        primitive_equilibrium_ok
    use primitive_geometry_grid, only: primitive_geometry_grid_t
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium, changed
    type(fixed_boundary_problem_t) :: problem
    integer :: info, unit
    real(dp) :: coupling
    character(len=*), parameter :: filename = "fixed_boundary_symmetry.nc"

    call create_cylinder_fixture(filename, surfaces=5)
    call read_gvec_cas3d_file(filename, equilibrium, info)
    call require(info == reader_ok, "read failed")
    call require(.not. equilibrium%stellarator_symmetric, &
        "control must exercise full sine/cosine storage")
    call build(equilibrium)
    call require(info == fixed_boundary_ok, "symmetric full storage rejected")
    call cross_mass(equilibrium, coupling)
    call require(coupling < 1.0e-12_dp, "symmetric cross-parity mass is nonzero")

    changed = equilibrium
    changed%winding = 0
    changed%has_boozer_position_frame = .false.
    call build(changed)
    call require(info == fixed_boundary_ok, "Cartesian position frame rejected")
    call cross_mass(changed, coupling)
    call require(coupling < 1.0e-12_dp, "Cartesian symmetric cross mass is nonzero")

    ! Opposite toroidal harmonics at m=0 are redundant. Sine contributions
    ! cancel physically; inspecting forbidden coefficients individually fails.
    changed = equilibrium
    changed%xhat%sine(:, 1, 2) = 0.3_dp
    changed%xhat%sine(:, 1, 3) = 0.3_dp
    call build(changed)
    call require(info == fixed_boundary_ok, "redundant cancellation rejected")

    ! A rigid vertical translation leaves every physical operator invariant.
    changed = equilibrium
    changed%zhat%cosine(:, 1, 1) = 0.4_dp
    call build(changed)
    call require(info == fixed_boundary_ok, "translated symmetric geometry rejected")

    changed = equilibrium
    changed%xhat%sine(:, 2, 2) = 0.02_dp
    call cross_mass(changed, coupling)
    write (*, '(a,es16.8)') "odd-harmonic relative cross mass: ", coupling
    call require(coupling > 1.0e-4_dp, "odd-harmonic oracle has no parity coupling")
    call build(changed)
    call require(info /= fixed_boundary_ok, "odd position harmonic accepted")
    changed%stellarator_symmetric = .true.
    call build(changed)
    call require(info /= fixed_boundary_ok, "false symmetry declaration trusted")

    ! The identity chart transform retains parity. A radial-dependent poloidal
    ! shift is physically symmetric but its fixed Fourier classes couple.
    call create_cylinder_fixture(filename, surfaces=5, chart_shift=0.0_dp)
    call read_gvec_cas3d_file(filename, changed, info)
    call require(info == reader_ok, "identity chart read failed")
    call build(changed)
    call require(info == fixed_boundary_ok, "identity chart rejected")
    call create_cylinder_fixture(filename, surfaces=5, chart_shift=0.2_dp)
    call read_gvec_cas3d_file(filename, changed, info)
    call require(info == reader_ok, "shifted chart read failed")
    call cross_mass(changed, coupling)
    write (*, '(a,es16.8)') "shifted-chart relative cross mass: ", coupling
    call require(coupling > 1.0e-4_dp, "shifted chart oracle has no coupling")
    call build(changed)
    call require(info /= fixed_boundary_ok, "parity-coupling chart accepted")

    call check_field_scaling()

    open (newunit=unit, file=filename, status="old")
    close (unit, status="delete")
    write (*, '(a)') "PASS"
contains
    subroutine build(input)
        type(gvec_cas3d_equilibrium_t), intent(in) :: input

        call build_fixed_boundary_problem(input, 5.0_dp / 3.0_dp, &
            2.0_dp, 1.0_dp, [1], [1], 1, problem, info)
    end subroutine build

    subroutine cross_mass(input, coupling)
        type(gvec_cas3d_equilibrium_t), intent(in) :: input
        real(dp), intent(out) :: coupling
        type(primitive_equilibrium_spline_t) :: spline
        type(primitive_geometry_grid_t) :: geometry
        real(dp), allocatable :: theta(:), zeta(:), weight(:, :)
        real(dp) :: pressure, slope, phase, value, normalization
        integer :: status, j, k, m, n

        call fit_primitive_equilibrium(input, spline, status)
        call require(status == primitive_equilibrium_ok, "oracle spline failed")
        call build_angular_grids(96, 96, theta, zeta)
        call evaluate_primitive_equilibrium(spline, 0.5_dp, theta, zeta, &
            geometry, pressure, slope, status)
        call require(status == primitive_equilibrium_ok, "oracle geometry failed")
        ! Independent physical kinetic energy for parallel displacement:
        ! xi = f(theta,zeta) B. Its mass is integral |J| |B|^2 f*g.
        ! The constant and sine trial factors belong to opposite classes.
        ! This oracle bypasses GLISS mass kernels and parity assembly entirely.
        weight = abs(geometry%signed_jacobian) * geometry%mod_b**2
        normalization = sum(weight)
        coupling = 0.0_dp
        do n = -2, 2
            do m = 0, 2
                value = 0.0_dp
                do k = 1, size(zeta)
                    do j = 1, size(theta)
                        phase = 2.0_dp * acos(-1.0_dp) &
                            * (real(m, dp) * theta(j) - real(n, dp) * zeta(k))
                        value = value + weight(j, k) * sin(phase)
                    end do
                end do
                coupling = max(coupling, abs(value) / normalization)
            end do
        end do
    end subroutine cross_mass

    subroutine check_field_scaling()
        type(gvec_cas3d_equilibrium_t) :: sheared, scaled
        real(dp), allocatable :: original(:, :)
        real(dp) :: angle, amplitude, reference, measured, length
        integer :: surface, exponent

        ! A toroidal shear along a purely toroidal field gives even beta and
        ! zero sigma. Its cross-parity normal/parallel mass persists as B->0.
        sheared = equilibrium
        sheared%winding = 0
        sheared%pressure = 100.0_dp
        sheared%poloidal_flux = 0.0_dp
        do surface = 1, size(sheared%s)
            angle = 0.2_dp * sheared%s(surface)
            original = sheared%xhat%cosine(surface, :, :)
            sheared%xhat%cosine(surface, :, :) = cos(angle) * original &
                - sin(angle) * sheared%yhat%cosine(surface, :, :)
            sheared%yhat%cosine(surface, :, :) = sin(angle) * original &
                + cos(angle) * sheared%yhat%cosine(surface, :, :)
            original = sheared%xhat%sine(surface, :, :)
            sheared%xhat%sine(surface, :, :) = cos(angle) * original &
                - sin(angle) * sheared%yhat%sine(surface, :, :)
            sheared%yhat%sine(surface, :, :) = sin(angle) * original &
                + cos(angle) * sheared%yhat%sine(surface, :, :)
        end do
        do exponent = 0, -12, -6
            amplitude = 10.0_dp**exponent
            scaled = equilibrium
            call rescale(scaled, 1.0_dp, amplitude)
            call build(scaled)
            call require(info == fixed_boundary_ok, &
                "field-rescaled symmetric geometry rejected")
            scaled = sheared
            scaled%toroidal_flux = amplitude * sheared%toroidal_flux
            scaled%pressure = amplitude**2 * sheared%pressure
            call parallel_normal_cross_mass(scaled, measured)
            if (exponent == 0) reference = measured
            write (*, '(a,i4,es16.8)') "field scaling cross mass: ", &
                exponent, measured
            call require(measured > 1.0e-3_dp, "shear cross mass vanished")
            call require(abs(measured - reference) < 1.0e-10_dp, &
                "dimensionless shear oracle changed under field scaling")
            call build(scaled)
            call require(info /= fixed_boundary_ok, &
                "field-rescaled parity coupling accepted")
        end do
        do exponent = -3, 3, 3
            length = 10.0_dp**exponent
            amplitude = 10.0_dp**(2 * exponent)
            scaled = equilibrium
            call rescale(scaled, length, amplitude)
            call cross_mass(scaled, measured)
            call require(measured < 1.0e-12_dp, &
                "rescaled symmetric cross mass is nonzero")
            call build(scaled)
            call require(info == fixed_boundary_ok, &
                "rescaled symmetric geometry rejected")
            scaled = sheared
            call rescale(scaled, length, amplitude)
            call parallel_normal_cross_mass(scaled, measured)
            call require(abs(measured - reference) < 1.0e-10_dp, &
                "shear oracle changed under length/field scaling")
            call build(scaled)
            call require(info /= fixed_boundary_ok, &
                "length/field-rescaled parity coupling accepted")
        end do
    end subroutine check_field_scaling

    subroutine rescale(input, length, magnetic)
        type(gvec_cas3d_equilibrium_t), intent(inout) :: input
        real(dp), intent(in) :: length, magnetic

        input%xhat%cosine = length * input%xhat%cosine
        input%xhat%sine = length * input%xhat%sine
        input%yhat%cosine = length * input%yhat%cosine
        input%yhat%sine = length * input%yhat%sine
        input%zhat%cosine = length * input%zhat%cosine
        input%zhat%sine = length * input%zhat%sine
        input%toroidal_flux = magnetic * length**2 * input%toroidal_flux
        input%poloidal_flux = magnetic * length**2 * input%poloidal_flux
        input%pressure = magnetic**2 * input%pressure
    end subroutine rescale

    subroutine parallel_normal_cross_mass(input, coupling)
        type(gvec_cas3d_equilibrium_t), intent(in) :: input
        real(dp), intent(out) :: coupling
        type(primitive_equilibrium_spline_t) :: spline
        type(primitive_geometry_grid_t) :: geometry
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp) :: pressure, slope, cross, normal, parallel, weight
        integer :: status, j, k

        call fit_primitive_equilibrium(input, spline, status)
        call require(status == primitive_equilibrium_ok, "shear spline failed")
        call build_angular_grids(96, 96, theta, zeta)
        call evaluate_primitive_equilibrium(spline, 0.5_dp, theta, zeta, &
            geometry, pressure, slope, status)
        call require(status == primitive_equilibrium_ok, "shear geometry failed")
        ! Coordinate radial displacement e_s and parallel B displacement have
        ! opposite fixed parity when both use cos(theta). Compute the Gram
        ! cross term directly from the physical metric, normalized by norms.
        cross = 0.0_dp
        normal = 0.0_dp
        parallel = 0.0_dp
        do k = 1, size(zeta)
            do j = 1, size(theta)
                weight = abs(geometry%signed_jacobian(j, k)) &
                    * cos(2.0_dp * acos(-1.0_dp) * theta(j))**2
                cross = cross + weight * (geometry%metric(j, k, 1, 2) &
                    * geometry%b_contravariant(j, k, 1) &
                    + geometry%metric(j, k, 1, 3) &
                    * geometry%b_contravariant(j, k, 2))
                normal = normal + weight * geometry%metric(j, k, 1, 1)
                parallel = parallel + weight * geometry%mod_b(j, k)**2
            end do
        end do
        coupling = abs(cross) / sqrt(normal * parallel)
    end subroutine parallel_normal_cross_mass

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (*, '(a)') message
        error stop 1
    end subroutine require
end program test_fixed_boundary_symmetry
