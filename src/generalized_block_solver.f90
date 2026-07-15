module generalized_block_solver
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use block_tridiagonal, only: apply_block_tridiagonal_into, block_factor_t, &
        block_tridiagonal_t, factorize_shifted, solve_factored
    implicit none
    private

    integer, parameter, public :: generalized_ok = 0
    integer, parameter, public :: generalized_invalid = -1
    integer, parameter, public :: generalized_mass_not_spd = -2
    integer, parameter, public :: generalized_no_convergence = -3

    type :: generalized_factor_t
        type(block_tridiagonal_t) :: shifted
        type(block_factor_t) :: factor
    end type generalized_factor_t

    public :: generalized_inertia
    public :: generalized_eigenpair_diagnostics
    public :: iterate_generalized_eigenvalue

contains

    subroutine generalized_inertia(stiffness, mass, shift, count, info)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        integer, intent(out) :: count, info
        type(generalized_factor_t) :: shifted_factor

        count = -1
        call factorize_generalized_shift(stiffness, mass, shift, &
            shifted_factor, info)
        if (info /= generalized_ok) return
        count = shifted_factor%factor%negative_count
    end subroutine generalized_inertia

    subroutine generalized_eigenpair_diagnostics(stiffness, mass, vector, &
            eigenvalue, quotient, residual, info)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:, :), eigenvalue
        real(dp), intent(out) :: quotient, residual
        integer, intent(out) :: info
        real(dp) :: mass_image(size(vector, 1), size(vector, 2))
        real(dp) :: squared_norm

        call validate_generalized_problem(stiffness, mass, info)
        if (info /= generalized_ok) return
        info = generalized_invalid
        if (.not. ieee_is_finite(eigenvalue)) return
        if (size(vector, 1) /= size(stiffness%diag, 1)) return
        if (size(vector, 2) /= size(stiffness%diag, 3)) return
        if (.not. all(ieee_is_finite(vector))) return
        call apply_block_tridiagonal_into(mass, vector, mass_image)
        squared_norm = sum(vector * mass_image)
        if (.not. ieee_is_finite(squared_norm)) return
        if (squared_norm <= 0.0_dp) return
        quotient = generalized_rayleigh_quotient(stiffness, mass, vector)
        residual = generalized_residual(stiffness, mass, vector, eigenvalue)
        info = generalized_ok
    end subroutine generalized_eigenpair_diagnostics

    subroutine iterate_generalized_eigenvalue(stiffness, mass, shift, &
            eigenvalue, vector, residual, info)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        real(dp), intent(out) :: eigenvalue, residual
        real(dp), allocatable, intent(out) :: vector(:, :)
        integer, intent(out) :: info
        type(generalized_factor_t) :: shifted_factor
        real(dp), allocatable :: iterate(:, :)
        real(dp) :: previous
        integer :: i, j, iteration

        call factorize_generalized_shift(stiffness, mass, shift, &
            shifted_factor, info)
        if (info /= generalized_ok) return
        allocate (vector(size(stiffness%diag, 1), &
            size(stiffness%diag, 3)), iterate(size(stiffness%diag, 1), &
            size(stiffness%diag, 3)))
        do j = 1, size(vector, 2)
            do i = 1, size(vector, 1)
                vector(i, j) = 1.0_dp + 0.1_dp * real(i, dp) &
                    + 0.01_dp * real(j, dp)
            end do
        end do
        call normalize_in_mass(vector, mass, info)
        if (info /= generalized_ok) return
        previous = huge(1.0_dp)
        do iteration = 1, 500
            call apply_block_tridiagonal_into(mass, vector, iterate)
            call solve_factored(shifted_factor%shifted, &
                shifted_factor%factor, iterate, info)
            if (info /= 0) then
                info = generalized_invalid
                return
            end if
            call normalize_in_mass(iterate, mass, info)
            if (info /= generalized_ok) return
            vector = iterate
            eigenvalue = generalized_rayleigh_quotient(stiffness, mass, &
                vector)
            residual = generalized_residual(stiffness, mass, vector, &
                eigenvalue)
            if (abs(eigenvalue - previous) <= 1.0e-13_dp &
                * max(1.0_dp, abs(eigenvalue)) &
                .and. residual <= 1.0e-12_dp) exit
            previous = eigenvalue
        end do
        if (iteration > 500) then
            info = generalized_no_convergence
            return
        end if
        info = generalized_ok
    end subroutine iterate_generalized_eigenvalue

    function generalized_rayleigh_quotient(stiffness, mass, vector) &
            result(quotient)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:, :)
        real(dp) :: quotient
        real(dp) :: stiffness_image(size(vector, 1), size(vector, 2))
        real(dp) :: mass_image(size(vector, 1), size(vector, 2))

        call apply_block_tridiagonal_into(stiffness, vector, stiffness_image)
        call apply_block_tridiagonal_into(mass, vector, mass_image)
        quotient = sum(vector * stiffness_image) / sum(vector * mass_image)
    end function generalized_rayleigh_quotient

    function generalized_residual(stiffness, mass, vector, eigenvalue) &
            result(residual)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:, :), eigenvalue
        real(dp) :: residual
        real(dp) :: stiffness_image(size(vector, 1), size(vector, 2))
        real(dp) :: mass_image(size(vector, 1), size(vector, 2)), scale

        call apply_block_tridiagonal_into(stiffness, vector, stiffness_image)
        call apply_block_tridiagonal_into(mass, vector, mass_image)
        scale = max(1.0_dp, norm2(stiffness_image), &
            abs(eigenvalue) * norm2(mass_image))
        residual = norm2(stiffness_image - eigenvalue * mass_image) / scale
    end function generalized_residual

    subroutine factorize_generalized_shift(stiffness, mass, shift, result, &
            info)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        type(generalized_factor_t), intent(out) :: result
        integer, intent(out) :: info

        info = generalized_invalid
        if (.not. ieee_is_finite(shift)) return
        call validate_generalized_problem(stiffness, mass, info)
        if (info /= generalized_ok) return
        call form_shifted_matrix(stiffness, mass, shift, result%shifted)
        call factorize_shifted(result%shifted, 0.0_dp, result%factor, info)
        if (info /= 0) then
            info = generalized_invalid
            return
        end if
        info = generalized_ok
    end subroutine factorize_generalized_shift

    subroutine validate_generalized_problem(stiffness, mass, info)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        integer, intent(out) :: info
        type(block_factor_t) :: mass_factor

        info = generalized_invalid
        if (.not. valid_block_pair(stiffness, mass)) return
        call factorize_shifted(mass, 0.0_dp, mass_factor, info)
        if (info /= 0) then
            info = generalized_mass_not_spd
            return
        end if
        if (mass_factor%negative_count /= 0) then
            info = generalized_mass_not_spd
            return
        end if
        info = generalized_ok
    end subroutine validate_generalized_problem

    pure subroutine form_shifted_matrix(stiffness, mass, shift, shifted)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        type(block_tridiagonal_t), intent(out) :: shifted

        allocate (shifted%diag, source=stiffness%diag - shift * mass%diag)
        allocate (shifted%off, source=stiffness%off - shift * mass%off)
    end subroutine form_shifted_matrix

    subroutine normalize_in_mass(vector, mass, info)
        real(dp), intent(inout) :: vector(:, :)
        type(block_tridiagonal_t), intent(in) :: mass
        integer, intent(out) :: info
        real(dp) :: mass_image(size(vector, 1), size(vector, 2))
        real(dp) :: squared_norm

        info = generalized_invalid
        call apply_block_tridiagonal_into(mass, vector, mass_image)
        squared_norm = sum(vector * mass_image)
        if (.not. ieee_is_finite(squared_norm)) return
        if (squared_norm <= 0.0_dp) return
        vector = vector / sqrt(squared_norm)
        info = generalized_ok
    end subroutine normalize_in_mass

    pure function valid_block_pair(stiffness, mass) result(valid)
        type(block_tridiagonal_t), intent(in) :: stiffness, mass
        logical :: valid

        valid = valid_block_matrix(stiffness)
        if (.not. valid) return
        valid = valid_block_matrix(mass)
        if (.not. valid) return
        valid = all(shape(stiffness%diag) == shape(mass%diag))
        if (.not. valid) return
        valid = all(shape(stiffness%off) == shape(mass%off))
    end function valid_block_pair

    pure function valid_block_matrix(matrix) result(valid)
        type(block_tridiagonal_t), intent(in) :: matrix
        logical :: valid
        real(dp) :: scale
        integer :: i, blocks, width

        valid = allocated(matrix%diag)
        if (.not. valid) return
        valid = allocated(matrix%off)
        if (.not. valid) return
        width = size(matrix%diag, 1)
        blocks = size(matrix%diag, 3)
        valid = width > 0 .and. blocks > 0
        if (.not. valid) return
        valid = size(matrix%diag, 2) == width
        if (.not. valid) return
        valid = size(matrix%off, 1) == width
        if (.not. valid) return
        valid = size(matrix%off, 2) == width
        if (.not. valid) return
        valid = size(matrix%off, 3) == blocks - 1
        if (.not. valid) return
        valid = all(ieee_is_finite(matrix%diag))
        if (.not. valid) return
        valid = all(ieee_is_finite(matrix%off))
        if (.not. valid) return
        do i = 1, blocks
            scale = max(1.0_dp, maxval(abs(matrix%diag(:, :, i))))
            valid = maxval(abs(matrix%diag(:, :, i) &
                - transpose(matrix%diag(:, :, i)))) &
                <= 128.0_dp * epsilon(1.0_dp) * scale
            if (.not. valid) return
        end do
    end function valid_block_matrix

end module generalized_block_solver
