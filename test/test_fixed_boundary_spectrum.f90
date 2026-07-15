program test_fixed_boundary_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    use dense_spectrum_support, only: dense_spectrum_invalid, &
        dense_spectrum_ok, unpermute_dense_vectors
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_invalid, fixed_boundary_ok, &
        fixed_boundary_energy_terms_t, diagnose_fixed_boundary_energy, &
        fixed_boundary_full_spectrum_t, fixed_boundary_problem_t, &
        fixed_boundary_rayleigh_gradient, fixed_boundary_spectrum_result_t, &
        solve_fixed_boundary_class, solve_fixed_boundary_full_spectrum
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

    call check_invalid_inputs(equilibrium)
    call check_row_permutation()
    call delete_fixture()
    write (*, "(a)") "PASS"

contains

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
