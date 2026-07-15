module gliss_energy_capi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
        c_f_pointer, c_int, c_ptr, c_size_t, c_sizeof
    use fixed_boundary_energy, only: fixed_boundary_energy_terms_t
    use fixed_boundary_spectrum, only: diagnose_fixed_boundary_energy, &
        fixed_boundary_allocation_error, fixed_boundary_invalid, &
        fixed_boundary_ok, fixed_boundary_rayleigh_gradient, &
        fixed_boundary_unknown_count
    use gliss_c_abi_support, only: error_buffer_status, status_allocation_error, &
        status_compute_error, status_invalid_argument, status_ok, write_error
    use gliss_c_contexts, only: stability_problem_context_t
    implicit none
    private

    type, bind(c) :: energy_terms_c
        integer(c_size_t) :: struct_size
        real(c_double) :: field_line_bending
        real(c_double) :: magnetic_shear
        real(c_double) :: magnetic_compression
        real(c_double) :: pressure_drive
        real(c_double) :: plasma_compressibility
        real(c_double) :: potential_energy
        real(c_double) :: kinetic_energy
        real(c_double) :: rayleigh_quotient
        real(c_double) :: closure_error
        real(c_double) :: closure_tolerance
    end type energy_terms_c

    public :: gliss_stability_problem_energy_c
    public :: gliss_stability_problem_rayleigh_vjp_c

contains

    function gliss_stability_problem_energy_c(handle, parity_class, &
            vector_count, vector_pointer, terms_pointer, error_pointer, &
            error_capacity) bind(c, name="gliss_stability_problem_energy") &
            result(status)
        type(c_ptr), value, intent(in) :: handle, vector_pointer, terms_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: vector_count, error_capacity
        integer(c_int) :: status
        type(stability_problem_context_t), pointer :: context
        type(energy_terms_c), pointer :: terms
        type(fixed_boundary_energy_terms_t) :: result
        real(c_double), pointer :: vector(:)
        integer :: count, info, pointer_shape(1)

        status = prepare_energy_output(terms_pointer, error_pointer, &
            error_capacity, terms)
        if (status /= status_ok) return
        if (.not. c_associated(handle)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "stability problem handle is null")
            return
        end if
        call c_f_pointer(handle, context)
        call fixed_boundary_unknown_count(context%problem, int(parity_class), &
            count, info)
        if (info /= fixed_boundary_ok) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "parity_class must be 1 or 2")
            return
        end if
        if (vector_count /= int(count, c_size_t)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "energy vector_count must equal the unknown count")
            return
        end if
        if (.not. c_associated(vector_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "energy vector pointer is null")
            return
        end if
        pointer_shape(1) = count
        call c_f_pointer(vector_pointer, vector, pointer_shape)
        call diagnose_fixed_boundary_energy(context%problem, int(parity_class), &
            vector, result, info)
        if (info /= fixed_boundary_ok) then
            call report_energy_error(info, status, error_pointer, error_capacity)
            return
        end if
        call fill_energy_terms(result, terms)
        status = status_ok
    end function gliss_stability_problem_energy_c

    function gliss_stability_problem_rayleigh_vjp_c(handle, parity_class, &
            vector_count, vector_pointer, cotangent, gradient_capacity, &
            gradient_pointer, error_pointer, error_capacity) &
            bind(c, name="gliss_stability_problem_rayleigh_vjp") result(status)
        type(c_ptr), value, intent(in) :: handle, vector_pointer
        type(c_ptr), value, intent(in) :: gradient_pointer, error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: vector_count
        integer(c_size_t), value, intent(in) :: gradient_capacity
        integer(c_size_t), value, intent(in) :: error_capacity
        real(c_double), value, intent(in) :: cotangent
        integer(c_int) :: status
        type(stability_problem_context_t), pointer :: context
        real(c_double), pointer :: vector(:), gradient(:)
        real(c_double), allocatable :: native_gradient(:)
        integer :: count, info, pointer_shape(1)

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(handle)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "stability problem handle is null")
            return
        end if
        call c_f_pointer(handle, context)
        call fixed_boundary_unknown_count(context%problem, int(parity_class), &
            count, info)
        if (info /= fixed_boundary_ok) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "parity_class must be 1 or 2")
            return
        end if
        if (vector_count /= int(count, c_size_t) .or. &
            gradient_capacity /= int(count, c_size_t)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "Rayleigh vector and gradient counts must equal the unknown count")
            return
        end if
        if (.not. c_associated(vector_pointer) .or. &
            .not. c_associated(gradient_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "Rayleigh vector and gradient pointers must not be null")
            return
        end if
        if (.not. ieee_is_finite(cotangent)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "Rayleigh cotangent must be finite")
            return
        end if
        pointer_shape(1) = count
        call c_f_pointer(vector_pointer, vector, pointer_shape)
        call fixed_boundary_rayleigh_gradient(context%problem, &
            int(parity_class), vector, native_gradient, info)
        if (info /= fixed_boundary_ok) then
            call report_rayleigh_error(info, status, error_pointer, &
                error_capacity)
            return
        end if
        native_gradient = cotangent * native_gradient
        if (.not. all(ieee_is_finite(native_gradient))) then
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "Rayleigh cotangent scaling produced a nonfinite gradient")
            return
        end if
        call c_f_pointer(gradient_pointer, gradient, pointer_shape)
        gradient = native_gradient
        status = status_ok
    end function gliss_stability_problem_rayleigh_vjp_c

    function prepare_energy_output(terms_pointer, error_pointer, &
            error_capacity, terms) result(status)
        type(c_ptr), value, intent(in) :: terms_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        type(energy_terms_c), pointer, intent(out) :: terms
        integer(c_int) :: status

        nullify (terms)
        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(terms_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "energy terms pointer is null")
            return
        end if
        call c_f_pointer(terms_pointer, terms)
        if (terms%struct_size /= c_sizeof(terms)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "energy terms struct_size is incompatible")
            return
        end if
        status = status_ok
    end function prepare_energy_output

    subroutine fill_energy_terms(result, terms)
        type(fixed_boundary_energy_terms_t), intent(in) :: result
        type(energy_terms_c), intent(out) :: terms

        terms%struct_size = c_sizeof(terms)
        terms%field_line_bending = result%field_line_bending
        terms%magnetic_shear = result%magnetic_shear
        terms%magnetic_compression = result%magnetic_compression
        terms%pressure_drive = result%pressure_drive
        terms%plasma_compressibility = result%plasma_compressibility
        terms%potential_energy = result%potential_energy
        terms%kinetic_energy = result%kinetic_energy
        terms%rayleigh_quotient = result%rayleigh_quotient
        terms%closure_error = result%closure_error
        terms%closure_tolerance = result%closure_tolerance
    end subroutine fill_energy_terms

    subroutine report_energy_error(info, status, error_pointer, error_capacity)
        integer, intent(in) :: info
        integer(c_int), intent(out) :: status
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity

        if (info == fixed_boundary_invalid) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "energy vector must be finite with positive kinetic norm")
        else if (info == fixed_boundary_allocation_error) then
            status = status_allocation_error
            call write_error(error_pointer, error_capacity, &
                "failed to allocate energy diagnostic storage")
        else
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "fixed-boundary energy decomposition failed")
        end if
    end subroutine report_energy_error

    subroutine report_rayleigh_error(info, status, error_pointer, &
            error_capacity)
        integer, intent(in) :: info
        integer(c_int), intent(out) :: status
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity

        if (info == fixed_boundary_invalid) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "Rayleigh vector must be finite with positive kinetic norm")
        else if (info == fixed_boundary_allocation_error) then
            status = status_allocation_error
            call write_error(error_pointer, error_capacity, &
                "failed to allocate Rayleigh derivative storage")
        else
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "fixed-boundary Rayleigh derivative failed")
        end if
    end subroutine report_rayleigh_error

end module gliss_energy_capi
