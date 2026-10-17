! Public import -> primitive spline -> fixed-boundary assembly -> certified
! solve on a genuine analytical toroidal equilibrium, independent of GVEC.
! For p=0 and curl B=0, 2 delta W = integral |curl(xi cross B)|^2/mu0 >= 0.
! In particular xi^s=s(1-s), xi^theta=xi^zeta=0 has exact quadratic energy
! Phi_s^2/(3 C mu0).  Eigenfrequencies are not asserted to be closed form.
program test_vacuum_torus_analytic
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cartesian_harmonic_spline, only: cartesian_harmonic_ok, &
        cartesian_jet_grid_t, evaluate_cartesian_harmonic_spline
    use compatible_radial_quadrature, only: build_gauss_legendre
    use field_periodic_cartesian, only: convert_field_periodic_jet, &
        field_periodic_cartesian_ok
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        diagnose_fixed_boundary_energy, fixed_boundary_energy_terms_t, &
        fixed_boundary_full_spectrum_t, fixed_boundary_ok, &
        fixed_boundary_problem_t, fixed_boundary_spectrum_result_t, &
        solve_fixed_boundary_class, solve_fixed_boundary_full_spectrum
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use gvec_cas3d_writer, only: write_gvec_cas3d_file, writer_ok
    use physical_constants, only: vacuum_permeability
    use primitive_equilibrium_spline, only: evaluate_primitive_equilibrium, &
        fit_primitive_equilibrium, primitive_equilibrium_ok, &
        primitive_equilibrium_spline_t
    use primitive_geometry_grid, only: primitive_geometry_grid_t
    use radial_feec_complex, only: build_radial_feec_complex, &
        evaluate_radial_feec_complex, radial_feec_complex_t
    use vacuum_torus_analytic_fixture, only: make_vacuum_torus, torus_c, &
        torus_field_constant, torus_major_radius, vacuum_torus_point
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: density = 1.0_dp / vacuum_permeability
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    real(dp) :: coarse(4), medium(4), fine(4)

    call check_import(9, 4, coarse)
    call check_import(17, 8, medium)
    call check_import(33, 12, fine)
    call require(all(medium < 0.02_dp * coarse), &
        "joint Fourier/radial refinement did not reduce all errors")
    call require(all(fine < 0.05_dp * medium), &
        "fine import has no refinement plateau")
    call require(fine(1) < 2.0e-8_dp, "fine Cartesian position error")
    call require(fine(2) < 2.0e-7_dp, "fine determinant error")
    call require(fine(3) < 2.0e-7_dp, "fine Cartesian field error")
    call require(fine(4) < 2.0e-6_dp, "fine force-free covariant curl error")
    call check_absolute_energy()
    call check_spectra()
    call check_corrupt_controls()
    write (*, '(a)') "vacuum torus analytical public path: PASS"

contains

    subroutine check_import(surfaces, harmonics, errors)
        integer, intent(in) :: surfaces, harmonics
        real(dp), intent(out) :: errors(4)
        type(gvec_cas3d_equilibrium_t) :: source
        type(primitive_equilibrium_spline_t) :: spline
        type(primitive_geometry_grid_t) :: geometry
        type(cartesian_jet_grid_t) :: jet
        real(dp), parameter :: coordinates(7) = &
            [1.0e-8_dp, 1.0e-5_dp, 0.0137_dp, 0.217_dp, &
             0.503_dp, 0.847_dp, 0.997_dp]
        real(dp) :: theta(37), zeta(7), position(3), field(3), actual(3)
        real(dp) :: jacobian, pressure, slope, curl_error, norm
        character(len=80) :: filename
        integer :: status, radial, j, k, unit

        call make_vacuum_torus(surfaces, harmonics, source)
        write (filename, '("vacuum_torus_",i0,"_",i0,".nc")') &
            surfaces, harmonics
        call write_gvec_cas3d_file(filename, source, status)
        call require(status == writer_ok, "analytical writer rejected fixture")
        call read_gvec_cas3d_file(filename, equilibrium, status)
        call require(status == reader_ok, "analytical reader rejected fixture")
        open (newunit=unit, file=filename, status="old")
        close (unit, status="delete")
        call fit_primitive_equilibrium(equilibrium, spline, status)
        call require(status == primitive_equilibrium_ok, "primitive fit failed")
        do j = 1, size(theta)
            theta(j) = real(j - 1, dp) / real(size(theta), dp)
        end do
        do k = 1, size(zeta)
            zeta(k) = real(k - 1, dp) / real(size(zeta), dp)
        end do
        errors = 0.0_dp
        do radial = 1, size(coordinates)
            call evaluate_cartesian_harmonic_spline(spline%radial_grid, &
                spline%position, coordinates(radial), theta, zeta, jet, status)
            call require(status == cartesian_harmonic_ok, "position jet failed")
            call convert_field_periodic_jet(zeta, 1, 1, jet, status)
            call require(status == field_periodic_cartesian_ok, "frame failed")
            call evaluate_primitive_equilibrium(spline, coordinates(radial), &
                theta, zeta, geometry, pressure, slope, status)
            call require(status == primitive_equilibrium_ok, "geometry failed")
            call require(pressure == 0.0_dp .and. slope == 0.0_dp, &
                "vacuum pressure gained a source")
            do k = 1, size(zeta)
                do j = 1, size(theta)
                    call vacuum_torus_point(coordinates(radial), theta(j), &
                        zeta(k), position, jacobian, field)
                    actual = geometry%b_contravariant(j, k, 1) &
                        * jet%poloidal(j, k, :) &
                        + geometry%b_contravariant(j, k, 2) &
                        * jet%toroidal(j, k, :)
                    errors(1) = max(errors(1), &
                        maxval(abs(jet%value(j, k, :) - position)))
                    errors(2) = max(errors(2), &
                        abs(geometry%signed_jacobian(j, k) / jacobian - 1.0_dp))
                    norm = sqrt(dot_product(field, field))
                    errors(3) = max(errors(3), &
                        sqrt(sum((actual - field)**2)) / norm)
                end do
            end do
            ! Exact vacuum covariants are (B_s,B_theta,B_zeta)=(0,0,2 pi F).
            ! Their radial derivative is an independent force-free oracle;
            ! axisymmetry and g_s,zeta=0 make the other curl terms vanish.
            curl_error = maxval(abs(geometry%b_covariant_radial)) &
                / (2.0_dp * pi * torus_field_constant)
            errors(4) = max(errors(4), curl_error)
        end do
        write (*, '(a,i0,a,i0,4(a,es12.4))') "surfaces ", surfaces, &
            " harmonics ", harmonics, " position ", errors(1), &
            " Jacobian ", errors(2), " B ", errors(3), " curl ", errors(4)
    end subroutine check_import

    subroutine build(input, problem, status, coupled, material_density)
        type(gvec_cas3d_equilibrium_t), intent(in) :: input
        type(fixed_boundary_problem_t), intent(out) :: problem
        integer, intent(out) :: status
        logical, optional, intent(in) :: coupled
        real(dp), optional, intent(in) :: material_density
        real(dp) :: rho
        logical :: both

        rho = density
        if (present(material_density)) rho = material_density
        both = .false.
        if (present(coupled)) both = coupled
        call build_fixed_boundary_problem(input, 5.0_dp / 3.0_dp, rho, &
            1.0e-6_dp, [0, 2], [0, 1], 2, problem, status, 96, 8, &
            coupled=both, radial_cells=3)
    end subroutine build

    subroutine check_absolute_energy()
        type(fixed_boundary_problem_t) :: problem
        type(fixed_boundary_energy_terms_t) :: energy
        type(fixed_boundary_spectrum_result_t) :: result
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: h1(:), dh1(:), l2(:), coefficients(:), vector(:)
        real(dp) :: sample(61, 5), rhs(61), s, flux, expected, expected_mass
        integer :: status, row, index

        call build_radial_feec_complex([0.0_dp, 1.0_dp / 3.0_dp, &
            2.0_dp / 3.0_dp, 1.0_dp], 2, .true., .true., complex, status)
        call require(status == 0, "oracle radial complex failed")
        call require(complex%h1_dofs == size(sample, 2), "oracle basis dimension")
        do row = 1, size(rhs)
            s = (real(row, dp) - 0.5_dp) / real(size(rhs), dp)
            call evaluate_radial_feec_complex(complex, s, h1, dh1, l2, status)
            call require(status == 0, "oracle basis evaluation failed")
            sample(row, :) = h1
            rhs(row) = s * (1.0_dp - s)
        end do
        allocate (coefficients(complex%h1_dofs))
        call least_squares(sample, rhs, coefficients)
        call build(equilibrium, problem, status)
        if (status /= fixed_boundary_ok) write (*, '(a,i0)') &
            "analytical public build status ", status
        call require(status == fixed_boundary_ok, "analytical public build failed")
        call solve_fixed_boundary_class(problem, 1, result, status)
        call require(status == fixed_boundary_ok, "analytical certified solve failed")
        allocate (vector(result%unknowns), source=0.0_dp)
        ! Normal trials are ordered by radial function then mode.  Only the
        ! constant (m=n=0) cosine trial is populated; it has no axis tie.
        do index = 1, size(coefficients)
            vector(2 * index - 1) = coefficients(index)
        end do
        call diagnose_fixed_boundary_energy(problem, 1, vector, energy, status)
        call require(status == fixed_boundary_ok, "public energy diagnosis failed")
        flux = 2.0_dp * pi * torus_field_constant * torus_c
        expected = flux**2 / (3.0_dp * torus_c * vacuum_permeability)
        expected_mass = independent_normal_mass()
        write (*, '(a,es12.4,a,es12.4)') "absolute energy error ", &
            abs(energy%potential_energy / expected - 1.0_dp), &
            " physical mass error ", abs(energy%kinetic_energy / expected_mass &
            - 1.0_dp)
        call require(abs(energy%potential_energy / expected - 1.0_dp) &
            < 2.0e-7_dp, "absolute analytical magnetic energy differs")
        call require(abs(energy%kinetic_energy / expected_mass - 1.0_dp) &
            < 2.0e-7_dp, "absolute analytical physical mass differs")
        call require(abs(energy%pressure_drive) < 1.0e-12_dp * expected, &
            "vacuum pressure drive is nonzero")
        call require(abs(energy%plasma_compressibility) < 1.0e-12_dp * expected, &
            "zero pressure has compressive plasma energy")
    end subroutine check_absolute_energy

    real(dp) function independent_normal_mass() result(mass)
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: s, d, radius2, angle, a, r, metric_ss, jacobian, x
        integer :: status, point, theta

        call build_gauss_legendre(32, nodes, weights, status)
        call require(status == 0, "independent integration rule failed")
        mass = 0.0_dp
        do point = 1, size(nodes)
            s = 0.5_dp * (nodes(point) + 1.0_dp)
            d = torus_major_radius - torus_c * s
            radius2 = 2.0_dp * torus_major_radius * torus_c * s &
                - torus_c**2 * s**2
            x = s * (1.0_dp - s)
            do theta = 0, 255
                angle = 2.0_dp * pi * real(theta, dp) / 256.0_dp
                a = torus_major_radius - sqrt(radius2) * cos(angle)
                r = d**2 / a
                metric_ss = torus_c**2 * (d**2 / radius2 &
                    + torus_major_radius**2 * sin(angle)**2 / a**2)
                jacobian = 4.0_dp * pi**2 * torus_c * r**2
                mass = mass + 0.5_dp * weights(point) / 256.0_dp &
                    * density * jacobian * metric_ss * x**2
            end do
        end do
    end function independent_normal_mass

    subroutine check_spectra()
        type(fixed_boundary_problem_t) :: problem
        type(fixed_boundary_full_spectrum_t) :: first, second, coupled, scaled
        type(fixed_boundary_spectrum_result_t) :: certificate
        type(gvec_cas3d_equilibrium_t) :: changed
        real(dp), allocatable :: expected(:)
        real(dp) :: scale
        integer :: status, parity

        call build(equilibrium, problem, status)
        call require(status == fixed_boundary_ok, "spectrum build failed")
        call solve_fixed_boundary_full_spectrum(problem, 1, first, status)
        call require(status == fixed_boundary_ok, "cosine full spectrum failed")
        call solve_fixed_boundary_full_spectrum(problem, 2, second, status)
        call require(status == fixed_boundary_ok, "sine full spectrum failed")
        do parity = 1, 2
            call solve_fixed_boundary_class(problem, parity, certificate, status)
            call require(status == fixed_boundary_ok, "vacuum certificate failed")
            call require(certificate%negative_count == 0, &
                "exact vacuum theorem violated by certified negative inertia")
            call require(certificate%floor_count >= certificate%mu_unknowns, &
                "parallel vacuum neutral space is missing")
        end do
        call require(minval(first%eigenvalues) >= -1.0e-6_dp, &
            "cosine vacuum spectrum has negative energy")
        call require(minval(second%eigenvalues) >= -1.0e-6_dp, &
            "sine vacuum spectrum has negative energy")
        call require(maxval(first%eigenvalues) > 1.0_dp, &
            "vacuum benchmark contains no positive branch")
        call build(equilibrium, problem, status, coupled=.true.)
        call require(status == fixed_boundary_ok, "coupled public build failed")
        call solve_fixed_boundary_full_spectrum(problem, 0, coupled, status)
        call require(status == fixed_boundary_ok, "coupled full spectrum failed")
        allocate (expected(size(first%eigenvalues) + size(second%eigenvalues)))
        expected(:size(first%eigenvalues)) = first%eigenvalues
        expected(size(first%eigenvalues) + 1:) = second%eigenvalues
        call sort_values(expected)
        call require(size(expected) == size(coupled%eigenvalues), &
            "coupled vacuum spectrum dimension differs")
        scale = max(1.0_dp, maxval(abs(expected)))
        call require(maxval(abs(expected - coupled%eigenvalues)) &
            < 1.0e-8_dp * scale, "coupled spectrum differs from parity union")
        changed = equilibrium
        changed%toroidal_flux = 1.7_dp * changed%toroidal_flux
        call build(changed, problem, status, material_density=2.3_dp * density)
        call require(status == fixed_boundary_ok, "scaled public build failed")
        call solve_fixed_boundary_full_spectrum(problem, 1, scaled, status)
        call require(status == fixed_boundary_ok, "scaled full spectrum failed")
        expected = first%eigenvalues * 1.7_dp**2 / 2.3_dp
        scale = max(1.0_dp, maxval(abs(expected)))
        call require(maxval(abs(expected - scaled%eigenvalues)) &
            < 1.0e-8_dp * scale, "physical B squared/density scaling differs")
        write (*, '(a,i0,a,i0,a)') "certified vacuum spectra: ", &
            size(first%eigenvalues), " + ", size(second%eigenvalues), &
            " parity eigenpairs, coupled union and physical scaling PASS"
    end subroutine check_spectra

    subroutine check_corrupt_controls()
        type(gvec_cas3d_equilibrium_t) :: changed
        type(fixed_boundary_problem_t) :: problem
        integer :: status

        changed = equilibrium
        changed%zhat%sine = 0.0_dp
        call build(changed, problem, status)
        call require(status /= fixed_boundary_ok, "collapsed torus was accepted")
        changed = equilibrium
        changed%toroidal_flux = 0.0_dp
        call build(changed, problem, status)
        call require(status /= fixed_boundary_ok, "zero field was accepted")
    end subroutine check_corrupt_controls

    subroutine least_squares(matrix, rhs, solution)
        real(dp), intent(in) :: matrix(:, :), rhs(:)
        real(dp), intent(out) :: solution(:)
        real(dp) :: a(size(matrix, 1), size(matrix, 2)), b(size(rhs), 1)
        real(dp) :: work(4096), reconstructed(size(rhs))
        integer :: status, row, column
        interface
            subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
                import dp
                character, intent(in) :: trans
                integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
                real(dp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
                integer, intent(out) :: info
            end subroutine dgels
        end interface

        a = matrix
        b(:, 1) = rhs
        call dgels('N', size(matrix, 1), size(matrix, 2), 1, a, size(matrix, 1), &
            b, size(rhs), work, size(work), status)
        call require(status == 0, "independent trial fit failed")
        solution = b(:size(solution), 1)
        reconstructed = 0.0_dp
        do column = 1, size(solution)
            do row = 1, size(rhs)
                reconstructed(row) = reconstructed(row) &
                    + matrix(row, column) * solution(column)
            end do
        end do
        call require(maxval(abs(reconstructed - rhs)) < 1.0e-13_dp, &
            "analytical displacement is not represented exactly")
    end subroutine least_squares

    subroutine sort_values(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: key
        integer :: first, previous

        do first = 2, size(values)
            key = values(first)
            previous = first - 1
            do while (previous >= 1)
                if (values(previous) <= key) exit
                values(previous + 1) = values(previous)
                previous = previous - 1
            end do
            values(previous + 1) = key
        end do
    end subroutine sort_values

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, '(a)') "FAIL: " // message
        error stop 1
    end subroutine require

end program test_vacuum_torus_analytic
