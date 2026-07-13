module variable_generalized_solver
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use stable_reduction, only: stable_dot_product
    use variable_block_tridiagonal, only: &
        apply_variable_block_tridiagonal, factorize_variable_shifted, &
        solve_variable_factored, validate_variable_blocks, &
        variable_block_factor_t, variable_block_ok, &
        variable_block_tridiagonal_t, variable_matrix_block_t
    implicit none
    private

    integer, parameter, public :: variable_generalized_ok = 0
    integer, parameter, public :: variable_generalized_invalid = -1
    integer, parameter, public :: variable_generalized_mass_not_spd = -2
    integer, parameter, public :: variable_generalized_no_convergence = -3

    public :: iterate_variable_generalized_eigenvalue
    public :: variable_generalized_diagnostics
    public :: variable_generalized_inertia

contains

    subroutine variable_generalized_inertia(stiffness, mass, shift, count, &
            info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        integer, intent(out) :: count, info
        type(variable_block_factor_t) :: factor

        count = -1
        call factorize_generalized_shift(stiffness, mass, shift, factor, info)
        if (info /= variable_generalized_ok) return
        count = factor%negative_count
    end subroutine variable_generalized_inertia

    subroutine variable_generalized_diagnostics(stiffness, mass, vector, &
            eigenvalue, quotient, residual, resolution, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:), eigenvalue
        real(dp), intent(out) :: quotient, residual, resolution
        integer, intent(out) :: info
        real(dp) :: mass_image(size(vector)), squared_norm

        call validate_generalized_problem(stiffness, mass, info)
        if (info /= variable_generalized_ok) return
        info = variable_generalized_invalid
        if (.not. ieee_is_finite(eigenvalue)) return
        if (size(vector) /= sum(stiffness%widths)) return
        if (.not. all(ieee_is_finite(vector))) return
        call apply_variable_block_tridiagonal(mass, vector, mass_image, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        squared_norm = dot_product(vector, mass_image)
        if (.not. ieee_is_finite(squared_norm)) return
        if (squared_norm <= 0.0_dp) return
        call variable_rayleigh_quotient(stiffness, mass, vector, quotient, info)
        if (info /= variable_generalized_ok) return
        call variable_residual(stiffness, mass, vector, eigenvalue, residual, &
            info)
        if (info /= variable_generalized_ok) return
        call variable_resolution(stiffness, mass, vector, eigenvalue, &
            resolution, info)
    end subroutine variable_generalized_diagnostics

    subroutine iterate_variable_generalized_eigenvalue(stiffness, mass, &
            shift, eigenvalue, vector, residual, resolution, info, controls)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        real(dp), intent(out) :: eigenvalue, residual, resolution
        real(dp), allocatable, intent(out) :: vector(:)
        integer, intent(out) :: info
        type(fixed_boundary_solver_controls_t), intent(in), optional :: controls
        type(variable_block_factor_t) :: factor
        type(fixed_boundary_solver_controls_t) :: stopping
        real(dp), allocatable :: iterate(:)
        real(dp) :: previous
        integer :: iteration, n

        stopping = fixed_boundary_solver_controls_t()
        if (present(controls)) stopping = controls
        call factorize_generalized_shift(stiffness, mass, shift, factor, info)
        if (info /= variable_generalized_ok) return
        n = sum(stiffness%widths)
        allocate (vector(n), iterate(n))
        call initialize_variable_iterate(vector)
        call normalize_variable_mass(vector, mass, info)
        if (info /= variable_generalized_ok) return
        previous = huge(1.0_dp)
        do iteration = 1, stopping%inverse_iteration_limit
            call update_variable_iterate(stiffness, mass, factor, vector, &
                iterate, eigenvalue, residual, resolution, info)
            if (info /= variable_generalized_ok) return
            if (iteration_converged(eigenvalue, previous, residual, &
                resolution, stopping)) exit
            previous = eigenvalue
        end do
        if (iteration > stopping%inverse_iteration_limit) then
            info = variable_generalized_no_convergence
            return
        end if
        info = variable_generalized_ok
    end subroutine iterate_variable_generalized_eigenvalue

    subroutine update_variable_iterate(stiffness, mass, factor, vector, &
            iterate, eigenvalue, residual, resolution, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        type(variable_block_factor_t), intent(in) :: factor
        real(dp), intent(inout) :: vector(:), iterate(:)
        real(dp), intent(out) :: eigenvalue, residual, resolution
        integer, intent(out) :: info

        call apply_variable_block_tridiagonal(mass, vector, iterate, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        call solve_variable_factored(factor, iterate, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        call normalize_variable_mass(iterate, mass, info)
        if (info /= variable_generalized_ok) return
        vector = iterate
        call variable_rayleigh_quotient(stiffness, mass, vector, eigenvalue, &
            info)
        if (info /= variable_generalized_ok) return
        call variable_residual(stiffness, mass, vector, eigenvalue, residual, &
            info)
        if (info /= variable_generalized_ok) return
        call variable_resolution(stiffness, mass, vector, eigenvalue, &
            resolution, info)
    end subroutine update_variable_iterate

    pure subroutine initialize_variable_iterate(vector)
        real(dp), intent(out) :: vector(:)
        integer :: i

        do i = 1, size(vector)
            vector(i) = 1.0_dp + 0.1_dp * real(i, dp)
        end do
    end subroutine initialize_variable_iterate

    subroutine factorize_generalized_shift(stiffness, mass, shift, factor, &
            info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        type(variable_block_factor_t), intent(out) :: factor
        integer, intent(out) :: info
        type(variable_block_tridiagonal_t) :: shifted

        info = variable_generalized_invalid
        if (.not. ieee_is_finite(shift)) return
        call validate_generalized_problem(stiffness, mass, info)
        if (info /= variable_generalized_ok) return
        call form_generalized_shift(stiffness, mass, shift, shifted)
        call factorize_variable_shifted(shifted, 0.0_dp, factor, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        info = variable_generalized_ok
    end subroutine factorize_generalized_shift

    subroutine validate_generalized_problem(stiffness, mass, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        integer, intent(out) :: info
        type(variable_block_factor_t) :: mass_factor

        info = variable_generalized_invalid
        if (.not. matching_variable_blocks(stiffness, mass)) return
        call factorize_variable_shifted(mass, 0.0_dp, mass_factor, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_mass_not_spd
            return
        end if
        if (mass_factor%negative_count /= 0) then
            info = variable_generalized_mass_not_spd
            return
        end if
        info = variable_generalized_ok
    end subroutine validate_generalized_problem

    function matching_variable_blocks(stiffness, mass) result(matches)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        logical :: matches
        integer :: info

        call validate_variable_blocks(stiffness, info)
        matches = info == variable_block_ok
        if (.not. matches) return
        call validate_variable_blocks(mass, info)
        matches = info == variable_block_ok
        if (.not. matches) return
        matches = size(stiffness%widths) == size(mass%widths)
        if (.not. matches) return
        matches = all(stiffness%widths == mass%widths)
    end function matching_variable_blocks

    subroutine form_generalized_shift(stiffness, mass, shift, shifted)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        type(variable_block_tridiagonal_t), intent(out) :: shifted
        integer :: block

        allocate (shifted%widths, source=stiffness%widths)
        allocate (shifted%diagonal(size(stiffness%diagonal)))
        allocate (shifted%lower(size(stiffness%lower)))
        do block = 1, size(stiffness%widths)
            allocate (shifted%diagonal(block)%values, &
                source=stiffness%diagonal(block)%values &
                - shift * mass%diagonal(block)%values)
            if (block < size(stiffness%widths)) then
                allocate (shifted%lower(block)%values, &
                    source=stiffness%lower(block)%values &
                    - shift * mass%lower(block)%values)
            end if
        end do
    end subroutine form_generalized_shift

    subroutine normalize_variable_mass(vector, mass, info)
        real(dp), intent(inout) :: vector(:)
        type(variable_block_tridiagonal_t), intent(in) :: mass
        integer, intent(out) :: info
        real(dp) :: image(size(vector)), squared_norm

        call apply_variable_block_tridiagonal(mass, vector, image, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        squared_norm = stable_dot_product(vector, image)
        if (.not. ieee_is_finite(squared_norm)) then
            info = variable_generalized_invalid
            return
        end if
        if (squared_norm <= 0.0_dp) then
            info = variable_generalized_invalid
            return
        end if
        vector = vector / sqrt(squared_norm)
        info = variable_generalized_ok
    end subroutine normalize_variable_mass

    subroutine variable_rayleigh_quotient(stiffness, mass, vector, quotient, &
            info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:)
        real(dp), intent(out) :: quotient
        integer, intent(out) :: info
        real(dp) :: stiffness_image(size(vector)), mass_image(size(vector))
        real(dp) :: denominator, numerator

        call apply_variable_block_tridiagonal(stiffness, vector, &
            stiffness_image, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        call apply_variable_block_tridiagonal(mass, vector, mass_image, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        numerator = stable_dot_product(vector, stiffness_image)
        denominator = stable_dot_product(vector, mass_image)
        if (.not. ieee_is_finite(numerator)) then
            info = variable_generalized_invalid
            return
        end if
        if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_dp) then
            info = variable_generalized_invalid
            return
        end if
        quotient = numerator / denominator
        if (.not. ieee_is_finite(quotient)) then
            info = variable_generalized_invalid
            return
        end if
        info = variable_generalized_ok
    end subroutine variable_rayleigh_quotient

    subroutine variable_residual(stiffness, mass, vector, eigenvalue, &
            residual, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:), eigenvalue
        real(dp), intent(out) :: residual
        integer, intent(out) :: info
        real(dp) :: stiffness_image(size(vector)), mass_image(size(vector))
        real(dp) :: mass_norm

        call apply_variable_block_tridiagonal(stiffness, vector, &
            stiffness_image, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        call apply_variable_block_tridiagonal(mass, vector, mass_image, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        mass_norm = norm2(mass_image)
        if (.not. ieee_is_finite(mass_norm) .or. mass_norm <= 0.0_dp) then
            info = variable_generalized_invalid
            return
        end if
        residual = norm2(stiffness_image - eigenvalue * mass_image) / mass_norm
        if (.not. ieee_is_finite(residual)) then
            info = variable_generalized_invalid
            return
        end if
        info = variable_generalized_ok
    end subroutine variable_residual

    subroutine variable_resolution(stiffness, mass, vector, eigenvalue, &
            resolution, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:), eigenvalue
        real(dp), intent(out) :: resolution
        integer, intent(out) :: info
        real(dp) :: absolute_image(size(vector)), mass_image(size(vector))
        real(dp) :: factor, mass_norm
        integer :: operation_count

        call apply_variable_block_tridiagonal(mass, vector, mass_image, info)
        if (info /= variable_block_ok) then
            info = variable_generalized_invalid
            return
        end if
        mass_norm = norm2(mass_image)
        if (.not. ieee_is_finite(mass_norm) .or. mass_norm <= 0.0_dp) then
            info = variable_generalized_invalid
            return
        end if
        call absolute_shifted_action(stiffness, mass, vector, eigenvalue, &
            absolute_image, operation_count)
        factor = (2.0_dp * real(operation_count, dp) &
            + 4.0_dp * real(size(vector), dp) + 16.0_dp) &
            * epsilon(1.0_dp)
        if (factor >= 1.0_dp) then
            info = variable_generalized_invalid
            return
        end if
        factor = factor / (1.0_dp - factor)
        resolution = factor * norm2(absolute_image) / mass_norm
        if (.not. ieee_is_finite(resolution)) then
            info = variable_generalized_invalid
            return
        end if
        info = variable_generalized_ok
    end subroutine variable_resolution

    pure subroutine absolute_shifted_action(stiffness, mass, vector, &
            eigenvalue, image, maximum_terms)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: vector(:), eigenvalue
        real(dp), intent(out) :: image(:)
        integer, intent(out) :: maximum_terms
        integer :: block, first, last, next_first, next_last, terms

        first = 1
        maximum_terms = 0
        do block = 1, size(stiffness%widths)
            last = first + stiffness%widths(block) - 1
            image(first:last) = matmul(abs(stiffness%diagonal(block)%values) &
                + abs(eigenvalue) * abs(mass%diagonal(block)%values), &
                abs(vector(first:last)))
            terms = stiffness%widths(block)
            if (block > 1) then
                next_first = first - stiffness%widths(block - 1)
                next_last = first - 1
                image(first:last) = image(first:last) + matmul( &
                    abs(stiffness%lower(block - 1)%values) + abs(eigenvalue) &
                    * abs(mass%lower(block - 1)%values), &
                    abs(vector(next_first:next_last)))
                terms = terms + stiffness%widths(block - 1)
            end if
            if (block < size(stiffness%widths)) then
                next_first = last + 1
                next_last = last + stiffness%widths(block + 1)
                image(first:last) = image(first:last) + matmul(transpose( &
                    abs(stiffness%lower(block)%values) + abs(eigenvalue) &
                    * abs(mass%lower(block)%values)), &
                    abs(vector(next_first:next_last)))
                terms = terms + stiffness%widths(block + 1)
            end if
            maximum_terms = max(maximum_terms, terms)
            first = last + 1
        end do
    end subroutine absolute_shifted_action

    pure function iteration_converged(eigenvalue, previous, residual, &
            resolution, controls) &
            result(converged)
        real(dp), intent(in) :: eigenvalue, previous, residual, resolution
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        logical :: converged

        converged = abs(eigenvalue - previous) <= max( &
            controls%eigenvalue_relative &
            * max(1.0_dp, abs(eigenvalue)), resolution)
        if (.not. converged) return
        converged = residual <= max(controls%residual_relative &
            * max(1.0_dp, abs(eigenvalue)), resolution)
    end function iteration_converged

end module variable_generalized_solver
