module gliss_full_spectrum_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
        c_f_pointer, c_int, c_ptr, c_size_t
    use fixed_boundary_spectrum, only: fixed_boundary_allocation_error, &
        fixed_boundary_full_spectrum_t, fixed_boundary_ok, &
        fixed_boundary_unknown_count, solve_fixed_boundary_full_spectrum
    use gliss_c_abi_support, only: error_buffer_status, status_allocation_error, &
        status_capacity, status_compute_error, status_invalid_argument, &
        status_ok, write_error
    use gliss_c_contexts, only: stability_problem_context_t
    implicit none
    private

    public :: gliss_stability_problem_full_spectrum_c

contains

    function gliss_stability_problem_full_spectrum_c(handle, parity_class, &
            eigenvalue_capacity, eigenvalue_pointer, residual_pointer, &
            resolution_pointer, rayleigh_pointer, eigenvector_capacity, &
            eigenvector_pointer, eigenvalues_written_pointer, &
            eigenvectors_written_pointer, error_pointer, error_capacity) &
            bind(c, name="gliss_stability_problem_full_spectrum") result(status)
        type(c_ptr), value, intent(in) :: handle, eigenvalue_pointer
        type(c_ptr), value, intent(in) :: residual_pointer, resolution_pointer
        type(c_ptr), value, intent(in) :: rayleigh_pointer, eigenvector_pointer
        type(c_ptr), value, intent(in) :: eigenvalues_written_pointer
        type(c_ptr), value, intent(in) :: eigenvectors_written_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: eigenvalue_capacity
        integer(c_size_t), value, intent(in) :: eigenvector_capacity
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        integer(c_size_t), pointer :: eigenvalues_written, eigenvectors_written
        type(stability_problem_context_t), pointer :: context
        type(fixed_boundary_full_spectrum_t) :: result
        real(c_double), pointer :: eigenvalues(:), residuals(:), resolutions(:)
        real(c_double), pointer :: rayleigh_quotients(:), eigenvectors(:, :)
        integer(c_size_t) :: dimension, required_vectors
        integer :: count, info, matrix_shape(2), vector_shape(1)

        status = prepare_full_spectrum_outputs(eigenvalues_written_pointer, &
            eigenvectors_written_pointer, error_pointer, error_capacity, &
            eigenvalues_written, eigenvectors_written)
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
        dimension = int(count, c_size_t)
        if (dimension > huge(required_vectors) / dimension) then
            status = status_capacity
            call write_error(error_pointer, error_capacity, &
                "full spectrum exceeds the ABI size limit")
            return
        end if
        required_vectors = dimension * dimension
        eigenvalues_written = dimension
        eigenvectors_written = required_vectors
        if (eigenvalue_capacity < dimension &
            .or. eigenvector_capacity < required_vectors) then
            status = status_capacity
            call write_error(error_pointer, error_capacity, &
                "full-spectrum output capacity is too small")
            return
        end if
        if (.not. valid_output_pointers(eigenvalue_pointer, residual_pointer, &
            resolution_pointer, rayleigh_pointer, eigenvector_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "full-spectrum output pointer is null")
            return
        end if
        call solve_fixed_boundary_full_spectrum(context%problem, &
            int(parity_class), result, info)
        if (info /= fixed_boundary_ok) then
            call report_full_spectrum_error(info, status, error_pointer, &
                error_capacity)
            return
        end if
        vector_shape(1) = count
        matrix_shape(1) = count
        matrix_shape(2) = count
        call c_f_pointer(eigenvalue_pointer, eigenvalues, vector_shape)
        call c_f_pointer(residual_pointer, residuals, vector_shape)
        call c_f_pointer(resolution_pointer, resolutions, vector_shape)
        call c_f_pointer(rayleigh_pointer, rayleigh_quotients, vector_shape)
        call c_f_pointer(eigenvector_pointer, eigenvectors, matrix_shape)
        eigenvalues = result%eigenvalues
        residuals = result%residuals
        resolutions = result%resolutions
        rayleigh_quotients = result%rayleigh_quotients
        eigenvectors = result%eigenvectors
        status = status_ok
    end function gliss_stability_problem_full_spectrum_c

    function prepare_full_spectrum_outputs(eigenvalues_written_pointer, &
            eigenvectors_written_pointer, error_pointer, error_capacity, &
            eigenvalues_written, eigenvectors_written) result(status)
        type(c_ptr), value, intent(in) :: eigenvalues_written_pointer
        type(c_ptr), value, intent(in) :: eigenvectors_written_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_size_t), pointer, intent(out) :: eigenvalues_written
        integer(c_size_t), pointer, intent(out) :: eigenvectors_written
        integer(c_int) :: status

        nullify (eigenvalues_written, eigenvectors_written)
        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(eigenvalues_written_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "eigenvalues-written output pointer is null")
            return
        end if
        call c_f_pointer(eigenvalues_written_pointer, eigenvalues_written)
        eigenvalues_written = 0_c_size_t
        if (.not. c_associated(eigenvectors_written_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "eigenvectors-written output pointer is null")
            return
        end if
        call c_f_pointer(eigenvectors_written_pointer, eigenvectors_written)
        eigenvectors_written = 0_c_size_t
        status = status_ok
    end function prepare_full_spectrum_outputs

    function valid_output_pointers(eigenvalues, residuals, resolutions, &
            rayleigh_quotients, eigenvectors) result(valid)
        type(c_ptr), value, intent(in) :: eigenvalues, residuals, resolutions
        type(c_ptr), value, intent(in) :: rayleigh_quotients, eigenvectors
        logical :: valid

        valid = c_associated(eigenvalues) .and. c_associated(residuals) &
            .and. c_associated(resolutions) &
            .and. c_associated(rayleigh_quotients) &
            .and. c_associated(eigenvectors)
    end function valid_output_pointers

    subroutine report_full_spectrum_error(info, status, error_pointer, &
            error_capacity)
        integer, intent(in) :: info
        integer(c_int), intent(out) :: status
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity

        if (info == fixed_boundary_allocation_error) then
            status = status_allocation_error
            call write_error(error_pointer, error_capacity, &
                "failed to allocate fixed-boundary spectrum storage")
        else
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "fixed-boundary full-spectrum computation failed")
        end if
    end subroutine report_full_spectrum_error

end module gliss_full_spectrum_capi
