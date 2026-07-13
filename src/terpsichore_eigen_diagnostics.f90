module terpsichore_eigen_diagnostics
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use variable_block_tridiagonal, only: &
        apply_variable_block_tridiagonal, variable_block_ok, &
        variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: terpsichore_eigen_diagnostics_ok = 0
    integer, parameter, public :: terpsichore_eigen_diagnostics_invalid = -1

    type, public :: terpsichore_eigen_diagnostics_t
        real(dp) :: growth_rate = 0.0_dp
        real(dp) :: reference_quotient = 0.0_dp
        real(dp) :: computed_potential = 0.0_dp
        real(dp) :: computed_kinetic = 0.0_dp
        real(dp) :: reference_residual = 0.0_dp
        real(dp) :: mode_overlap = 0.0_dp
    end type terpsichore_eigen_diagnostics_t

    public :: compute_terpsichore_eigen_diagnostics

contains

    subroutine compute_terpsichore_eigen_diagnostics(stiffness, mass, &
            eigenvalue, eigenvector, reference, reference_potential, &
            reference_kinetic, alfven_normalization, diagnostics, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: eigenvalue, eigenvector(:), reference(:)
        real(dp), intent(in) :: reference_potential, reference_kinetic
        real(dp), intent(in) :: alfven_normalization
        type(terpsichore_eigen_diagnostics_t), intent(out) :: diagnostics
        integer, intent(out) :: info
        real(dp) :: eigen_mass(size(reference)), reference_k(size(reference))
        real(dp) :: reference_m(size(reference)), eigen_norm, scale

        diagnostics = terpsichore_eigen_diagnostics_t()
        info = terpsichore_eigen_diagnostics_invalid
        if (.not. valid_scalars(eigenvalue, reference_potential, &
            reference_kinetic, alfven_normalization)) return
        if (size(eigenvector) /= size(reference)) return
        if (.not. all(ieee_is_finite(eigenvector))) return
        if (.not. all(ieee_is_finite(reference))) return
        call apply_variable_block_tridiagonal(stiffness, reference, &
            reference_k, info)
        if (info /= variable_block_ok) return
        call apply_variable_block_tridiagonal(mass, reference, reference_m, &
            info)
        if (info /= variable_block_ok) return
        call apply_variable_block_tridiagonal(mass, eigenvector, eigen_mass, &
            info)
        if (info /= variable_block_ok) return
        diagnostics%computed_potential = stable_dot(reference, reference_k)
        diagnostics%computed_kinetic = stable_dot(reference, reference_m)
        eigen_norm = stable_dot(eigenvector, eigen_mass)
        if (diagnostics%computed_kinetic <= 0.0_dp .or. eigen_norm <= 0.0_dp) &
            return
        diagnostics%reference_quotient = reference_potential &
            / reference_kinetic
        scale = max(norm2(reference_k), abs(diagnostics%reference_quotient) &
            * norm2(reference_m))
        if (scale <= 0.0_dp) return
        diagnostics%reference_residual = norm2(reference_k &
            - diagnostics%reference_quotient * reference_m) / scale
        diagnostics%mode_overlap = abs(stable_dot(eigenvector, reference_m)) &
            / sqrt(eigen_norm * diagnostics%computed_kinetic)
        diagnostics%growth_rate = -sign(sqrt(abs(eigenvalue) &
            / alfven_normalization), eigenvalue)
        if (.not. diagnostics_are_finite(diagnostics)) return
        info = terpsichore_eigen_diagnostics_ok
    end subroutine compute_terpsichore_eigen_diagnostics

    pure function valid_scalars(eigenvalue, potential, kinetic, normalization) &
            result(valid)
        real(dp), intent(in) :: eigenvalue, potential, kinetic, normalization
        logical :: valid

        valid = ieee_is_finite(eigenvalue) .and. ieee_is_finite(potential) &
            .and. ieee_is_finite(kinetic) &
            .and. ieee_is_finite(normalization)
        if (.not. valid) return
        valid = kinetic > 0.0_dp .and. normalization > 0.0_dp
    end function valid_scalars

    pure function diagnostics_are_finite(diagnostics) result(valid)
        type(terpsichore_eigen_diagnostics_t), intent(in) :: diagnostics
        logical :: valid

        valid = ieee_is_finite(diagnostics%growth_rate) &
            .and. ieee_is_finite(diagnostics%reference_quotient) &
            .and. ieee_is_finite(diagnostics%computed_potential) &
            .and. ieee_is_finite(diagnostics%computed_kinetic) &
            .and. ieee_is_finite(diagnostics%reference_residual) &
            .and. ieee_is_finite(diagnostics%mode_overlap)
    end function diagnostics_are_finite

    pure function stable_dot(first, second) result(product)
        real(dp), intent(in) :: first(:), second(:)
        real(dp) :: product, correction, term, updated
        integer :: i

        product = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(first)
            term = first(i) * second(i)
            updated = product + term
            if (abs(product) >= abs(term)) then
                correction = correction + (product - updated) + term
            else
                correction = correction + (term - updated) + product
            end if
            product = updated
        end do
        product = product + correction
    end function stable_dot

end module terpsichore_eigen_diagnostics
