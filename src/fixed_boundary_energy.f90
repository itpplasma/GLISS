module fixed_boundary_energy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use stable_reduction, only: stable_dot_product
    use variable_block_tridiagonal, only: &
        apply_variable_block_tridiagonal, pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: fixed_boundary_energy_term_count = 5
    integer, parameter, public :: fixed_boundary_energy_ok = 0
    integer, parameter, public :: fixed_boundary_energy_invalid = -1
    integer, parameter, public :: fixed_boundary_energy_allocation = -2
    integer, parameter, public :: fixed_boundary_energy_inconsistent = -3

    type, public :: fixed_boundary_energy_store_t
        type(variable_block_tridiagonal_t) :: terms( &
            fixed_boundary_energy_term_count)
    end type fixed_boundary_energy_store_t

    type, public :: fixed_boundary_energy_terms_t
        real(dp) :: field_line_bending = 0.0_dp
        real(dp) :: magnetic_shear = 0.0_dp
        real(dp) :: magnetic_compression = 0.0_dp
        real(dp) :: pressure_drive = 0.0_dp
        real(dp) :: plasma_compressibility = 0.0_dp
        real(dp) :: potential_energy = 0.0_dp
        real(dp) :: kinetic_energy = 0.0_dp
        real(dp) :: rayleigh_quotient = 0.0_dp
        real(dp) :: closure_error = 0.0_dp
        real(dp) :: closure_tolerance = 0.0_dp
    end type fixed_boundary_energy_terms_t

    public :: diagnose_fixed_boundary_energy_store
    public :: pack_fixed_boundary_energy_store
    public :: rayleigh_gradient_fixed_boundary_store

contains

    subroutine pack_fixed_boundary_energy_store(dense_terms, permutation, &
            widths, store, info)
        real(dp), intent(in) :: dense_terms(:, :, :)
        integer, intent(in) :: permutation(:), widths(:)
        type(fixed_boundary_energy_store_t), intent(out) :: store
        integer, intent(out) :: info
        integer :: term

        info = fixed_boundary_energy_invalid
        if (size(dense_terms, 3) /= fixed_boundary_energy_term_count) return
        do term = 1, fixed_boundary_energy_term_count
            call pack_permuted_variable_blocks(dense_terms(:, :, term), &
                permutation, widths, store%terms(term), info)
            if (info /= variable_block_ok) then
                info = fixed_boundary_energy_invalid
                return
            end if
        end do
        info = fixed_boundary_energy_ok
    end subroutine pack_fixed_boundary_energy_store

    subroutine diagnose_fixed_boundary_energy_store(stiffness, mass, store, &
            permutation, vector, result, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        type(fixed_boundary_energy_store_t), intent(in) :: store
        integer, intent(in) :: permutation(:)
        real(dp), intent(in) :: vector(:)
        type(fixed_boundary_energy_terms_t), intent(out) :: result
        integer, intent(out) :: info
        real(dp), allocatable :: component_image_sum(:), image(:)
        real(dp), allocatable :: permuted(:), potential_image(:)
        real(dp) :: components(fixed_boundary_energy_term_count)
        real(dp) :: matrix_closure, ones(fixed_boundary_energy_term_count)
        real(dp) :: roundoff_tolerance, scale
        integer :: allocation_status, index

        info = fixed_boundary_energy_invalid
        if (size(vector) < 1 .or. size(vector) /= size(permutation)) return
        if (.not. all(ieee_is_finite(vector))) return
        allocate (permuted(size(vector)), image(size(vector)), &
            potential_image(size(vector)), component_image_sum(size(vector)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = fixed_boundary_energy_allocation
            return
        end if
        do index = 1, size(vector)
            if (permutation(index) < 1 .or. &
                permutation(index) > size(vector)) return
            permuted(index) = vector(permutation(index))
        end do
        call quadratic_form(stiffness, permuted, potential_image, &
            result%potential_energy, info)
        if (info /= fixed_boundary_energy_ok) return
        call quadratic_form(mass, permuted, image, result%kinetic_energy, info)
        if (info /= fixed_boundary_energy_ok) return
        if (result%kinetic_energy <= 0.0_dp) then
            info = fixed_boundary_energy_invalid
            return
        end if
        component_image_sum = 0.0_dp
        do index = 1, fixed_boundary_energy_term_count
            call quadratic_form(store%terms(index), permuted, image, &
                components(index), info)
            if (info /= fixed_boundary_energy_ok) return
            component_image_sum = component_image_sum + image
        end do
        call assign_components(components, result)
        result%rayleigh_quotient = result%potential_energy &
            / result%kinetic_energy
        ones = 1.0_dp
        result%closure_error = abs(result%potential_energy &
            - stable_dot_product(ones, components))
        scale = max(1.0_dp, abs(result%potential_energy), &
            sum(abs(components)))
        matrix_closure = stable_dot_product(abs(permuted), &
            abs(potential_image - component_image_sum))
        roundoff_tolerance = 128.0_dp * epsilon(1.0_dp) &
            * real(size(vector), dp) * scale
        if (matrix_closure > 64.0_dp * roundoff_tolerance) then
            info = fixed_boundary_energy_inconsistent
            return
        end if
        result%closure_tolerance = matrix_closure + roundoff_tolerance
        if (.not. valid_result(result)) then
            info = fixed_boundary_energy_inconsistent
            return
        end if
        info = fixed_boundary_energy_ok
    end subroutine diagnose_fixed_boundary_energy_store

    subroutine rayleigh_gradient_fixed_boundary_store(stiffness, mass, &
            permutation, vector, gradient, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        integer, intent(in) :: permutation(:)
        real(dp), intent(in) :: vector(:)
        real(dp), allocatable, intent(out) :: gradient(:)
        integer, intent(out) :: info
        real(dp), allocatable :: permuted(:), stiffness_image(:)
        real(dp), allocatable :: mass_image(:), permuted_gradient(:)
        real(dp) :: kinetic_energy, potential_energy, quotient
        integer :: allocation_status, index

        info = fixed_boundary_energy_invalid
        if (size(vector) < 1 .or. size(vector) /= size(permutation)) return
        if (.not. all(ieee_is_finite(vector))) return
        allocate (permuted(size(vector)), stiffness_image(size(vector)), &
            mass_image(size(vector)), permuted_gradient(size(vector)), &
            gradient(size(vector)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = fixed_boundary_energy_allocation
            return
        end if
        do index = 1, size(vector)
            if (permutation(index) < 1 .or. &
                permutation(index) > size(vector)) return
            permuted(index) = vector(permutation(index))
        end do
        call quadratic_form(stiffness, permuted, stiffness_image, &
            potential_energy, info)
        if (info /= fixed_boundary_energy_ok) return
        call quadratic_form(mass, permuted, mass_image, kinetic_energy, info)
        if (info /= fixed_boundary_energy_ok) return
        if (kinetic_energy <= 0.0_dp) then
            info = fixed_boundary_energy_invalid
            return
        end if
        quotient = potential_energy / kinetic_energy
        permuted_gradient = 2.0_dp * (stiffness_image &
            - quotient * mass_image) / kinetic_energy
        gradient = 0.0_dp
        do index = 1, size(vector)
            gradient(permutation(index)) = permuted_gradient(index)
        end do
        if (.not. all(ieee_is_finite(gradient))) then
            info = fixed_boundary_energy_inconsistent
            return
        end if
        info = fixed_boundary_energy_ok
    end subroutine rayleigh_gradient_fixed_boundary_store

    subroutine quadratic_form(matrix, vector, image, value, info)
        type(variable_block_tridiagonal_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:)
        real(dp), intent(out) :: image(:), value
        integer, intent(out) :: info

        call apply_variable_block_tridiagonal(matrix, vector, image, info)
        if (info /= variable_block_ok) then
            info = fixed_boundary_energy_invalid
            return
        end if
        value = stable_dot_product(vector, image)
        if (.not. ieee_is_finite(value)) then
            info = fixed_boundary_energy_invalid
            return
        end if
        info = fixed_boundary_energy_ok
    end subroutine quadratic_form

    pure subroutine assign_components(components, result)
        real(dp), intent(in) :: components(fixed_boundary_energy_term_count)
        type(fixed_boundary_energy_terms_t), intent(inout) :: result

        result%field_line_bending = components(1)
        result%magnetic_shear = components(2)
        result%magnetic_compression = components(3)
        result%pressure_drive = components(4)
        result%plasma_compressibility = components(5)
    end subroutine assign_components

    pure function valid_result(result) result(valid)
        type(fixed_boundary_energy_terms_t), intent(in) :: result
        real(dp) :: values(10)
        logical :: valid

        values = [result%field_line_bending, result%magnetic_shear, &
            result%magnetic_compression, result%pressure_drive, &
            result%plasma_compressibility, result%potential_energy, &
            result%kinetic_energy, result%rayleigh_quotient, &
            result%closure_error, result%closure_tolerance]
        valid = all(ieee_is_finite(values))
        if (.not. valid) return
        valid = result%closure_error <= result%closure_tolerance
    end function valid_result

end module fixed_boundary_energy
