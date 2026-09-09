program test_fixed_boundary_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_problem_ok, &
        compatible_two_component_problem_t
    use cylinder_fixture, only: create_cylinder_fixture
    use dense_spectrum_support, only: dense_spectrum_invalid, &
        dense_spectrum_ok, unpermute_dense_vectors
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_invalid, fixed_boundary_ok, &
        fixed_boundary_energy_terms_t, diagnose_fixed_boundary_energy, &
        fixed_boundary_full_spectrum_t, fixed_boundary_problem_t, &
        fixed_boundary_rayleigh_gradient, fixed_boundary_spectrum_result_t, &
        solve_fixed_boundary_class, solve_fixed_boundary_full_spectrum
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_invalid, primitive_kernel_ok
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none

    character(len=*), parameter :: fixture = "fixed_boundary_spectrum.nc"
    real(dp), parameter :: reference_lowest(2) = &
        [-7.9144227183717817e1_dp, -7.9144227377194269e1_dp]
    real(dp), parameter :: reference_certificate_limit = 1.2e-3_dp
    real(dp), parameter :: reference_relative_limit = 1.0e-8_dp
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(fixed_boundary_problem_t) :: problem
    type(fixed_boundary_spectrum_result_t) :: first, second, repeated
    type(fixed_boundary_full_spectrum_t) :: full, full_repeated
    type(fixed_boundary_energy_terms_t) :: energy
    integer :: info

    call create_cylinder_fixture(fixture)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "spectrum fixture read failed")

    call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
        2.0_dp, 1.0_dp, [1, 2], [1, 1], 1, problem, info)
    call require(info == fixed_boundary_ok, "problem construction failed")
    call solve_fixed_boundary_class(problem, 1, first, info)
    call require(info == fixed_boundary_ok, "class-one solve failed")
    call solve_fixed_boundary_class(problem, 2, second, info)
    call require(info == fixed_boundary_ok, "class-two solve failed")
    call check_result(first, 1)
    call check_result(second, 2)
    call diagnose_fixed_boundary_energy(problem, 1, first%eigenvector, &
        energy, info)
    call require(info == fixed_boundary_ok, "energy decomposition failed")
    call check_energy(energy, first%lowest_eigenvalue)
    call check_rayleigh_gradient()
    call check_reference_certificate(first, reference_lowest(1))
    call check_reference_certificate(second, reference_lowest(2))

    call solve_fixed_boundary_full_spectrum(problem, 1, full, info)
    call require(info == fixed_boundary_ok, "full-spectrum solve failed")
    call check_full_spectrum(full, first)
    call solve_fixed_boundary_full_spectrum(problem, 1, full_repeated, info)
    call require(info == fixed_boundary_ok, "repeated full-spectrum solve failed")
    call require(all(full_repeated%eigenvalues == full%eigenvalues), &
        "repeated full-spectrum solve changed the eigenvalues")
    call require(all(full_repeated%eigenvectors == full%eigenvectors), &
        "repeated full-spectrum solve changed the eigenvectors")
    call require(all(full_repeated%rayleigh_quotients &
        == full%rayleigh_quotients), &
        "repeated full-spectrum solve changed the Rayleigh quotients")
    call require(all(full_repeated%residuals == full%residuals), &
        "repeated full-spectrum solve changed the residuals")
    call require(all(full_repeated%resolutions == full%resolutions), &
        "repeated full-spectrum solve changed the resolutions")

    call solve_fixed_boundary_class(problem, 1, repeated, info)
    call require(info == fixed_boundary_ok, "repeated solve failed")
    call require(repeated%lowest_eigenvalue == first%lowest_eigenvalue, &
        "repeated solve changed the eigenvalue")
    call require(all(repeated%eigenvector == first%eigenvector), &
        "repeated solve changed the eigenvector")
    call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
        2.0_dp, 1.0_dp, [1, 2], [1, 1], 1, problem, info)
    call require(info == fixed_boundary_ok, "same-object rebuild failed")
    call solve_fixed_boundary_class(problem, 1, repeated, info)
    call require(info == fixed_boundary_ok, "rebuilt problem solve failed")
    call require(repeated%lowest_eigenvalue == first%lowest_eigenvalue, &
        "same-object rebuild changed the eigenvalue")
    call require(all(repeated%eigenvector == first%eigenvector), &
        "same-object rebuild changed the eigenvector")

    call check_geometry_orientation()
    call check_angular_convergence()
    call check_material_derivative()
    call check_invalid_inputs(equilibrium)
    call check_row_permutation()
    call check_density_scaling()
    call delete_fixture()
    write (*, "(a)") "PASS"

contains

    subroutine check_geometry_orientation()
        type(gvec_cas3d_equilibrium_t) :: mapped
        type(primitive_equilibrium_spline_t) :: spline
        type(fixed_boundary_problem_t) :: folded_problem
        type(compatible_two_component_problem_t) :: folded_marginality
        real(dp), allocatable :: fields(:, :, :), drive(:, :)
        real(dp) :: theta(16), zeta(8), factor
        integer :: i, j, status, orientation

        do i = 1, size(theta)
            theta(i) = real(i - 1, dp) / real(size(theta), dp)
        end do
        do i = 1, size(zeta)
            zeta(i) = real(i - 1, dp) / real(size(zeta), dp)
        end do
        ! The analytic torus map has J=-8*pi^2*r*r'*(R+r*cos(theta)).
        ! r=.5*sqrt(s)*(1-.8*s) is represented exactly: its m=1 quotient
        ! is linear. Its determinant changes sign at s=5/12, while each
        ! individual angular surface retains a coherent orientation.
        mapped = equilibrium
        do i = 1, size(mapped%s)
            factor = 1.0_dp - 0.8_dp * mapped%s(i)
            do j = 1, size(mapped%poloidal_modes)
                if (mapped%poloidal_modes(j) /= 1) cycle
                mapped%xhat%cosine(i, j, :) = &
                    factor * mapped%xhat%cosine(i, j, :)
                mapped%xhat%sine(i, j, :) = &
                    factor * mapped%xhat%sine(i, j, :)
                mapped%yhat%cosine(i, j, :) = &
                    factor * mapped%yhat%cosine(i, j, :)
                mapped%yhat%sine(i, j, :) = &
                    factor * mapped%yhat%sine(i, j, :)
                mapped%zhat%cosine(i, j, :) = &
                    factor * mapped%zhat%cosine(i, j, :)
                mapped%zhat%sine(i, j, :) = &
                    factor * mapped%zhat%sine(i, j, :)
            end do
        end do
        call fit_primitive_equilibrium(mapped, spline, status)
        call require(status == primitive_equilibrium_ok, "polynomial map fit failed")
        orientation = 0
        call evaluate_primitive_kernel_surface(spline, 0.2_dp, theta, zeta, &
            fields, drive, status, orientation=orientation)
        call require(status == primitive_kernel_ok, "negative chart rejected")
        call require(orientation == -1, "analytic negative determinant differs")
        call evaluate_primitive_kernel_surface(spline, 0.8_dp, theta, zeta, &
            fields, drive, status, orientation=orientation)
        call require(status == primitive_kernel_invalid, "radial fold accepted")
        orientation = 0
        call evaluate_primitive_kernel_surface(spline, 0.8_dp, theta, zeta, &
            fields, drive, status, orientation=orientation)
        call require(status == primitive_kernel_ok, "positive chart rejected")
        call require(orientation == 1, "analytic positive determinant differs")
        call evaluate_primitive_kernel_surface(spline, 0.9_dp, theta, zeta, &
            fields, drive, status, orientation=orientation)
        call require(status == primitive_kernel_ok, &
            "consistent positive chart rejected")
        call build_fixed_boundary_problem(mapped, 5.0_dp / 3.0_dp, &
            2.0_dp, 1.0_dp, [1], [1], 1, folded_problem, status, 16, 8)
        call require(status /= fixed_boundary_ok, "production radial fold accepted")
        call build_compatible_two_component_problem(mapped, [1], [1], &
            [0.5_dp], 1, 1, 16, 8, folded_marginality, status)
        call require(status /= compatible_problem_ok, &
            "marginality radial fold accepted")

        ! Moving the major radius R=3+2*s adds 2*cos(theta) to r' in J.
        ! At s=.5 it dominates r', so an angular surface crosses J=0.
        mapped = equilibrium
        do i = 1, size(mapped%s)
            factor = (3.0_dp + 2.0_dp * mapped%s(i)) / 3.0_dp
            mapped%xhat%cosine(i, 1, :) = factor * mapped%xhat%cosine(i, 1, :)
            mapped%yhat%sine(i, 1, :) = factor * mapped%yhat%sine(i, 1, :)
        end do
        call fit_primitive_equilibrium(mapped, spline, status)
        call require(status == primitive_equilibrium_ok, "moving-axis fit failed")
        call evaluate_primitive_kernel_surface(spline, 0.5_dp, theta, zeta, &
            fields, drive, status)
        call require(status == primitive_kernel_invalid, "angular fold accepted")

        mapped = equilibrium
        mapped%zhat%cosine = 0.0_dp
        mapped%zhat%sine = 0.0_dp
        call fit_primitive_equilibrium(mapped, spline, status)
        call require(status == primitive_equilibrium_ok, "degenerate map fit failed")
        call evaluate_primitive_kernel_surface(spline, 0.5_dp, theta, zeta, &
            fields, drive, status)
        call require(status == primitive_kernel_invalid, "zero determinant accepted")
    end subroutine check_geometry_orientation

    subroutine check_material_derivative()
        type(fixed_boundary_problem_t) :: shifted_problem
        type(fixed_boundary_full_spectrum_t) :: shifted
        type(fixed_boundary_energy_terms_t) :: terms
        real(dp), parameter :: gamma = 5.0_dp / 3.0_dp
        real(dp), parameter :: steps(3) = [1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp]
        real(dp) :: gradient, finite_difference, values(2), scale
        integer :: step, side, status, index

        ! Trace of the complete mass-orthonormal spectrum is independent of
        ! eigenvector ordering/crossings. Central differences rebuild the
        ! production pencil; no differentiated code supplies this oracle.
        gradient = 0.0_dp
        do index = 1, size(full%eigenvalues)
            call diagnose_fixed_boundary_energy(problem, 1, &
                full%eigenvectors(:, index), terms, status)
            call require(status == fixed_boundary_ok, "derivative energy failed")
            gradient = gradient + terms%plasma_compressibility &
                / (gamma * terms%kinetic_energy)
        end do
        call require(ieee_is_finite(gradient), "nonfinite gamma gradient")
        call require(gradient > 0.0_dp, "missing compressibility derivative")
        do step = 1, size(steps)
            do side = 1, 2
                call build_fixed_boundary_problem(equilibrium, &
                    gamma + real(2 * side - 3, dp) * steps(step), &
                    2.0_dp, 1.0_dp, [1, 2], [1, 1], 1, shifted_problem, status)
                call require(status == fixed_boundary_ok, "gamma build failed")
                call solve_fixed_boundary_full_spectrum(shifted_problem, &
                    1, shifted, status)
                call require(status == fixed_boundary_ok, "gamma solve failed")
                values(side) = sum(shifted%eigenvalues)
            end do
            finite_difference = (values(2) - values(1)) / (2.0_dp * steps(step))
            call require(ieee_is_finite(finite_difference), "nonfinite gamma FD")
            scale = max(1.0_dp, abs(gradient))
            call require(abs(finite_difference - gradient) < 1.0e-7_dp * scale, &
                "production gamma derivative lacks a step plateau")
        end do
    end subroutine check_material_derivative

    subroutine check_angular_convergence()
        type(fixed_boundary_problem_t) :: refined
        type(fixed_boundary_energy_terms_t) :: refined_energy
        type(fixed_boundary_spectrum_result_t) :: refined_result
        integer :: status

        ! The cylinder is toroidally homogeneous: increasing a resolved
        ! trapezoidal grid preserves its Fourier inner products exactly.
        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            2.0_dp, 1.0_dp, [1, 2], [1, 1], 1, refined, status, 128, 96)
        call require(status == fixed_boundary_ok, "refined grid construction failed")
        call diagnose_fixed_boundary_energy(refined, 1, first%eigenvector, &
            refined_energy, status)
        call require(status == fixed_boundary_ok, "refined energy failed")
        call diagnose_fixed_boundary_energy(problem, 1, first%eigenvector, &
            energy, status)
        call require(status == fixed_boundary_ok, "baseline energy failed")
        call require(abs(refined_energy%kinetic_energy &
            - energy%kinetic_energy) < 1.0e-10_dp, &
            "resolved cylinder Fourier mass changed with angular grid")
        call solve_fixed_boundary_class(refined, 1, refined_result, status)
        call require(status == fixed_boundary_ok, "refined solve failed")
        call require(refined_result%angular_theta == 128, "theta metadata lost")
        call require(refined_result%angular_zeta == 96, "zeta metadata lost")
    end subroutine check_angular_convergence

    subroutine check_density_scaling()
        type(fixed_boundary_problem_t) :: scaled_problem
        type(fixed_boundary_spectrum_result_t) :: scaled
        integer :: status

        ! K is independent of density and M is linear in density, so all
        ! physical omega^2 values scale inversely while inertia is unchanged.
        call build_fixed_boundary_problem(equilibrium, 5.0_dp / 3.0_dp, &
            8.0_dp, 0.25_dp, [1, 2], [1, 1], 1, scaled_problem, status)
        call require(status == fixed_boundary_ok, "scaled-density build failed")
        call solve_fixed_boundary_class(scaled_problem, 1, scaled, status)
        call require(status == fixed_boundary_ok, "scaled-density solve failed")
        call require(abs(4.0_dp * scaled%lowest_eigenvalue &
            - first%lowest_eigenvalue) < 1.0e-7_dp, &
            "physical spectrum violates inverse-density scaling")
        call require(scaled%negative_count == first%negative_count, &
            "positive density rescaling changed the instability count")
    end subroutine check_density_scaling

    subroutine check_reference_certificate(result, reference)
        type(fixed_boundary_spectrum_result_t), intent(in) :: result
        real(dp), intent(in) :: reference

        call require(result%certificate <= reference_certificate_limit, &
            "reference eigenvalue certificate is too wide")
        call require(abs(result%lowest_eigenvalue - reference) &
            <= result%certificate, &
            "reference eigenvalue lies outside the certificate")
        call require(abs(result%lowest_eigenvalue - reference) &
            <= reference_relative_limit * abs(reference), &
            "reference eigenvalue exceeds the portable relative limit")
    end subroutine check_reference_certificate

    subroutine check_rayleigh_gradient()
        type(fixed_boundary_energy_terms_t) :: plus, minus
        real(dp), allocatable :: gradient(:), primal(:), tangent(:), shifted(:)
        real(dp) :: centered, exact, step, scale
        integer :: index, status

        allocate (primal(size(first%eigenvector)), &
            tangent(size(first%eigenvector)), &
            shifted(size(first%eigenvector)))
        do index = 1, size(tangent)
            tangent(index) = sin(real(index, dp))
            primal(index) = first%eigenvector(index) &
                + 0.01_dp * cos(real(index, dp))
        end do
        call fixed_boundary_rayleigh_gradient(problem, 1, primal, &
            gradient, status)
        call require(status == fixed_boundary_ok, &
            "Rayleigh gradient action failed")
        call require(all(ieee_is_finite(gradient)), &
            "Rayleigh gradient contains nonfinite values")
        scale = max(1.0_dp, sqrt(dot_product(gradient, gradient)) &
            * sqrt(dot_product(primal, primal)))
        call require(abs(dot_product(gradient, primal)) &
            < 1.0e-11_dp * scale, &
            "Rayleigh gradient violates scale invariance")
        step = 1.0e-6_dp
        shifted = primal + step * tangent
        call diagnose_fixed_boundary_energy(problem, 1, shifted, plus, status)
        call require(status == fixed_boundary_ok, &
            "positive Rayleigh perturbation failed")
        shifted = primal - step * tangent
        call diagnose_fixed_boundary_energy(problem, 1, shifted, minus, status)
        call require(status == fixed_boundary_ok, &
            "negative Rayleigh perturbation failed")
        centered = (plus%rayleigh_quotient - minus%rayleigh_quotient) &
            / (2.0_dp * step)
        exact = dot_product(gradient, tangent)
        if (abs(centered - exact) >= 2.0e-7_dp &
            * max(1.0_dp, abs(centered), abs(exact))) &
            write (error_unit, "(a,3es24.15)") &
            "Rayleigh centered/exact/difference: ", centered, exact, &
            centered - exact
        call require(abs(centered - exact) &
            < 2.0e-7_dp * max(1.0_dp, abs(centered), abs(exact)), &
            "Rayleigh gradient disagrees with a centered reevaluation")
    end subroutine check_rayleigh_gradient

    subroutine check_energy(terms, eigenvalue)
        type(fixed_boundary_energy_terms_t), intent(in) :: terms
        real(dp), intent(in) :: eigenvalue
        real(dp) :: positive_scale

        positive_scale = max(1.0_dp, abs(terms%potential_energy))
        call require(terms%kinetic_energy > 0.0_dp, &
            "kinetic energy is not positive")
        call require(abs(terms%kinetic_energy - 1.0_dp) < 1.0e-11_dp, &
            "eigenvector is not mass normalized")
        call require(abs(terms%rayleigh_quotient - eigenvalue) &
            <= 1.0e-10_dp * max(1.0_dp, abs(eigenvalue)), &
            "energy quotient differs from the eigenvalue")
        call require(terms%closure_error <= terms%closure_tolerance, &
            "energy terms do not close")
        call require(terms%field_line_bending >= -1.0e-12_dp * positive_scale, &
            "field-line bending energy is negative")
        call require(terms%magnetic_shear >= -1.0e-12_dp * positive_scale, &
            "magnetic shear energy is negative")
        call require(terms%magnetic_compression &
            >= -1.0e-12_dp * positive_scale, &
            "magnetic compression energy is negative")
        call require(terms%plasma_compressibility &
            >= -1.0e-12_dp * positive_scale, &
            "plasma compressibility energy is negative")
    end subroutine check_energy

    subroutine check_result(result, parity_class)
        type(fixed_boundary_spectrum_result_t), intent(in) :: result
        integer, intent(in) :: parity_class

        call require(result%parity_class == parity_class, &
            "result parity class is wrong")
        call require(result%field_periods == 1, &
            "result field-period count is wrong")
        call require(result%mode_count == 2, "result mode count is wrong")
        call require(result%unknowns == 196, "result unknown count is wrong")
        call require(result%normal_unknowns == 64, &
            "normal unknown count is wrong")
        call require(result%eta_unknowns == 66, &
            "eta unknown count is wrong")
        call require(result%mu_unknowns == 66, &
            "mu unknown count is wrong")
        call require(size(result%eigenvector) == result%unknowns, &
            "eigenvector size is wrong")
        call require(all(ieee_is_finite(result%eigenvector)), &
            "eigenvector contains nonfinite values")
        call require(ieee_is_finite(result%lowest_eigenvalue), &
            "eigenvalue is nonfinite")
        call require(result%certificate == result%inertia_interval &
            + result%eigenpair_residual + result%eigenpair_resolution, &
            "certificate components do not close")
        call require(result%negative_count >= 0, "negative count is invalid")
        call require(result%floor_count >= 0, "floor count is invalid")
    end subroutine check_result

    subroutine check_full_spectrum(full, certified)
        type(fixed_boundary_full_spectrum_t), intent(in) :: full
        type(fixed_boundary_spectrum_result_t), intent(in) :: certified
        real(dp) :: overlap
        integer :: certified_index

        call require(size(full%eigenvalues) == certified%unknowns, &
            "full-spectrum eigenvalue count is wrong")
        call require(size(full%eigenvectors, 1) == certified%unknowns &
            .and. size(full%eigenvectors, 2) == certified%unknowns, &
            "full-spectrum eigenvector shape is wrong")
        call require(all(ieee_is_finite(full%eigenvalues)), &
            "full spectrum contains nonfinite eigenvalues")
        call require(all(ieee_is_finite(full%eigenvectors)), &
            "full spectrum contains nonfinite eigenvectors")
        call require(size(full%rayleigh_quotients) == certified%unknowns &
            .and. all(ieee_is_finite(full%rayleigh_quotients)), &
            "full spectrum contains invalid Rayleigh quotients")
        call require(size(full%residuals) == certified%unknowns &
            .and. all(ieee_is_finite(full%residuals)) &
            .and. all(full%residuals >= 0.0_dp), &
            "full spectrum contains invalid residuals")
        call require(size(full%resolutions) == certified%unknowns &
            .and. all(ieee_is_finite(full%resolutions)) &
            .and. all(full%resolutions >= 0.0_dp), &
            "full spectrum contains invalid resolutions")
        call require(all(full%eigenvalues(2:) >= &
            full%eigenvalues(:size(full%eigenvalues) - 1)), &
            "full spectrum is not sorted")
        call require(count(full%eigenvalues < -certified%zero_floor) &
            == certified%negative_count, "full-spectrum negative count differs")
        call require(count(abs(full%eigenvalues) <= certified%zero_floor) &
            == certified%floor_count, "full-spectrum floor count differs")
        certified_index = 1
        if (certified%negative_count == 0) &
            certified_index = certified%floor_count + 1
        call require(abs(full%eigenvalues(certified_index) &
            - certified%lowest_eigenvalue) &
            <= certified%certificate + 1.0e-12_dp &
            * abs(certified%lowest_eigenvalue), &
            "full spectrum disagrees with the certified active eigenvalue")
        overlap = abs(dot_product(full%eigenvectors(:, certified_index), &
            certified%eigenvector)) / sqrt(dot_product( &
            full%eigenvectors(:, certified_index), &
            full%eigenvectors(:, certified_index)) &
            * dot_product(certified%eigenvector, certified%eigenvector))
        call require(overlap > 1.0_dp - 1.0e-10_dp, &
            "full-spectrum vector is not in dynamic component order")
    end subroutine check_full_spectrum

    subroutine check_row_permutation()
        real(dp) :: vectors(3, 2)
        integer :: status

        vectors(1, :) = [1.0_dp, 10.0_dp]
        vectors(2, :) = [2.0_dp, 20.0_dp]
        vectors(3, :) = [3.0_dp, 30.0_dp]
        call unpermute_dense_vectors(vectors, [2, 3, 1], status)
        call require(status == dense_spectrum_ok, &
            "valid row permutation failed")
        call require(all(vectors(1, :) == [3.0_dp, 30.0_dp]) &
            .and. all(vectors(2, :) == [1.0_dp, 10.0_dp]) &
            .and. all(vectors(3, :) == [2.0_dp, 20.0_dp]), &
            "dense eigenvectors were unpermuted incorrectly")
        call unpermute_dense_vectors(vectors, [1, 1, 3], status)
        call require(status == dense_spectrum_invalid, &
            "duplicate row permutation was accepted")
    end subroutine check_row_permutation

    subroutine check_invalid_inputs(local_equilibrium)
        type(gvec_cas3d_equilibrium_t), intent(in) :: local_equilibrium
        type(fixed_boundary_problem_t) :: invalid_problem
        type(fixed_boundary_spectrum_result_t) :: invalid_result
        type(fixed_boundary_energy_terms_t) :: invalid_energy
        real(dp) :: nan
        integer :: status

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call build_fixed_boundary_problem(local_equilibrium, nan, 2.0_dp, &
            1.0_dp, [1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "nonfinite adiabatic index was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            0.0_dp, 1.0_dp, [1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "zero density was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 0.0_dp, [1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "zero floor was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [1], 0, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "degree zero was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [1], 5, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "degree five was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1, 1], [1, 1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "duplicate mode was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [-1], [1], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "negative poloidal mode was accepted")
        ! sin(32*theta) vanishes on the 64-point grid, while modes 31 and
        ! 33 have coincident cosines. Neither table represents its continuum
        ! Fourier inner product on that grid.
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [32], [0], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "poloidal Nyquist mode was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [31, 33], [0, 0], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "aliased poloidal mode pair was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [-32], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "negative toroidal Nyquist mode was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [32], 1, invalid_problem, status)
        call require(status == fixed_boundary_invalid, &
            "positive toroidal Nyquist mode was accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [32], 1, invalid_problem, status, 64, 128)
        call require(status == fixed_boundary_ok, &
            "resolved toroidal mode was rejected")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [1], 1, invalid_problem, status, 0, 64)
        call require(status == fixed_boundary_invalid, "zero angular count accepted")
        call build_fixed_boundary_problem(local_equilibrium, 1.0_dp, &
            2.0_dp, 1.0_dp, [1], [1], 1, invalid_problem, status, huge(1), 2)
        call require(status == fixed_boundary_invalid, &
            "angular product overflow accepted")
        call solve_fixed_boundary_class(problem, 0, invalid_result, status)
        call require(status == fixed_boundary_invalid, &
            "invalid parity class was accepted")
        call solve_fixed_boundary_full_spectrum(problem, 0, full, status)
        call require(status == fixed_boundary_invalid, &
            "full spectrum accepted an invalid parity class")
        call diagnose_fixed_boundary_energy(problem, 0, first%eigenvector, &
            invalid_energy, status)
        call require(status == fixed_boundary_invalid, &
            "energy decomposition accepted an invalid parity class")
        call diagnose_fixed_boundary_energy(problem, 1, &
            first%eigenvector(:size(first%eigenvector) - 1), invalid_energy, &
            status)
        call require(status == fixed_boundary_invalid, &
            "energy decomposition accepted the wrong vector size")
        first%eigenvector(1) = nan
        call diagnose_fixed_boundary_energy(problem, 1, first%eigenvector, &
            invalid_energy, status)
        call require(status == fixed_boundary_invalid, &
            "energy decomposition accepted a nonfinite vector")
    end subroutine check_invalid_inputs

    subroutine delete_fixture()
        integer :: unit, status

        open (newunit=unit, file=fixture, status="old", iostat=status)
        call require(status == 0, "failed to open spectrum fixture")
        close (unit, status="delete", iostat=status)
        call require(status == 0, "failed to delete spectrum fixture")
    end subroutine delete_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_fixed_boundary_spectrum
