module gliss_pressure_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
        c_f_pointer, c_int, c_ptr, c_size_t
    use, intrinsic :: iso_fortran_env, only: int64
    use fixed_boundary_spectrum, only: fixed_boundary_allocation_error, &
        fixed_boundary_invalid, fixed_boundary_is_free, fixed_boundary_ok, &
        fixed_boundary_pressure_trace_jvp, fixed_boundary_pressure_trace_vjp, &
        fixed_boundary_unknown_count
    use gliss_c_abi_support, only: error_buffer_status, status_allocation_error, &
        status_capacity, status_compute_error, status_invalid_argument, &
        status_ok, write_error
    use gliss_c_contexts, only: stability_problem_context_t
    implicit none
    private

    public :: gliss_stability_problem_pressure_samples_c
    public :: gliss_stability_problem_pressure_trace_jvp_c
    public :: gliss_stability_problem_pressure_trace_vjp_c

contains

    function gliss_stability_problem_pressure_samples_c(handle, capacity, &
            nodes_pointer, pressure_pointer, written_pointer, error_pointer, &
            error_capacity) bind(c, name="gliss_stability_problem_pressure_samples") &
            result(status)
        type(c_ptr), value, intent(in) :: handle, nodes_pointer, pressure_pointer
        type(c_ptr), value, intent(in) :: written_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: capacity, error_capacity
        integer(c_int) :: status
        type(stability_problem_context_t), pointer :: context
        integer(c_size_t), pointer :: written
        real(c_double), pointer, contiguous :: nodes(:), pressure(:)
        integer :: count, shape(1), sample

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(written_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, "written pointer is null")
            return
        end if
        call c_f_pointer(written_pointer, written)
        written = 0_c_size_t
        status = pressure_context(handle, context, error_pointer, error_capacity)
        if (status /= status_ok) return
        count = size(context%equilibrium%equilibrium%s)
        written = int(count, c_size_t)
        if (capacity < written) then
            status = status_capacity
            call write_error(error_pointer, error_capacity, &
                "pressure buffer is too small")
            return
        end if
        if (.not. c_associated(nodes_pointer) &
            .or. .not. c_associated(pressure_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "pressure output pointer is null")
            return
        end if
        shape(1) = count
        call c_f_pointer(nodes_pointer, nodes, shape)
        call c_f_pointer(pressure_pointer, pressure, shape)
        do sample = 1, count
            nodes(sample) = context%equilibrium%equilibrium%s(sample)
            pressure(sample) = context%equilibrium%equilibrium%pressure(sample)
        end do
        status = status_ok
    end function gliss_stability_problem_pressure_samples_c

    function gliss_stability_problem_pressure_trace_jvp_c(handle, parity_class, &
            vector_count, basis_count, vectors_pointer, pressure_count, &
            direction_pointer, derivative_pointer, error_pointer, error_capacity) &
            bind(c, name="gliss_stability_problem_pressure_trace_jvp") result(status)
        type(c_ptr), value, intent(in) :: handle, vectors_pointer, direction_pointer
        type(c_ptr), value, intent(in) :: derivative_pointer, error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: vector_count, basis_count
        integer(c_size_t), value, intent(in) :: pressure_count
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(stability_problem_context_t), pointer :: context
        real(c_double), pointer, contiguous :: vectors(:, :), direction(:)
        real(c_double), pointer :: output
        real(c_double) :: derivative
        integer :: info, shape(1)

        status = prepare_basis(handle, parity_class, vector_count, basis_count, &
            vectors_pointer, context, vectors, error_pointer, error_capacity)
        if (status /= status_ok) return
        status = status_invalid_argument
        if (pressure_count /= &
            int(size(context%equilibrium%equilibrium%s), c_size_t)) then
            call write_error(error_pointer, error_capacity, &
                "pressure direction has the wrong sample count")
            return
        end if
        if (.not. c_associated(direction_pointer) &
            .or. .not. c_associated(derivative_pointer)) then
            call write_error(error_pointer, error_capacity, &
                "pressure direction or derivative pointer is null")
            return
        end if
        shape(1) = int(pressure_count)
        call c_f_pointer(direction_pointer, direction, shape)
        call fixed_boundary_pressure_trace_jvp(context%problem, &
            context%equilibrium%equilibrium, int(parity_class), vectors, direction, &
            derivative, info)
        if (info /= fixed_boundary_ok) then
            call report_pressure_error(info, status, error_pointer, error_capacity)
            return
        end if
        call c_f_pointer(derivative_pointer, output)
        output = derivative
        status = status_ok
    end function gliss_stability_problem_pressure_trace_jvp_c

    function gliss_stability_problem_pressure_trace_vjp_c(handle, parity_class, &
            vector_count, basis_count, vectors_pointer, cotangent, gradient_capacity, &
            gradient_pointer, error_pointer, error_capacity) &
            bind(c, name="gliss_stability_problem_pressure_trace_vjp") result(status)
        type(c_ptr), value, intent(in) :: handle, vectors_pointer, gradient_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: vector_count, basis_count
        integer(c_size_t), value, intent(in) :: gradient_capacity
        integer(c_size_t), value, intent(in) :: error_capacity
        real(c_double), value, intent(in) :: cotangent
        integer(c_int) :: status
        type(stability_problem_context_t), pointer :: context
        real(c_double), pointer, contiguous :: vectors(:, :), output(:)
        real(c_double), allocatable :: gradient(:)
        integer :: info, shape(1)

        status = prepare_basis(handle, parity_class, vector_count, basis_count, &
            vectors_pointer, context, vectors, error_pointer, error_capacity)
        if (status /= status_ok) return
        shape(1) = size(context%equilibrium%equilibrium%s)
        if (gradient_capacity < int(shape(1), c_size_t)) then
            status = status_capacity
            call write_error(error_pointer, error_capacity, &
                "pressure gradient buffer is too small")
            return
        end if
        if (.not. c_associated(gradient_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "pressure gradient pointer is null")
            return
        end if
        call fixed_boundary_pressure_trace_vjp(context%problem, &
            context%equilibrium%equilibrium, int(parity_class), vectors, cotangent, &
            gradient, info)
        if (info /= fixed_boundary_ok) then
            call report_pressure_error(info, status, error_pointer, error_capacity)
            return
        end if
        call c_f_pointer(gradient_pointer, output, shape)
        output = gradient
        status = status_ok
    end function gliss_stability_problem_pressure_trace_vjp_c

    function pressure_context(handle, context, error_pointer, error_capacity) &
            result(status)
        type(c_ptr), value, intent(in) :: handle, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        type(stability_problem_context_t), pointer, intent(out) :: context
        integer(c_int) :: status

        status = status_invalid_argument
        nullify (context)
        if (.not. c_associated(handle)) then
            call write_error(error_pointer, error_capacity, &
                "stability problem handle is null")
            return
        end if
        call c_f_pointer(handle, context)
        if (.not. associated(context%equilibrium)) then
            call write_error(error_pointer, error_capacity, &
                "equilibrium data is unavailable")
            return
        end if
        if (.not. allocated(context%equilibrium%equilibrium%s)) return
        if (.not. allocated(context%equilibrium%equilibrium%pressure)) return
        if (size(context%equilibrium%equilibrium%s) /= &
            size(context%equilibrium%equilibrium%pressure)) return
        status = status_ok
    end function pressure_context

    function prepare_basis(handle, parity_class, vector_count, basis_count, &
            vectors_pointer, context, vectors, error_pointer, error_capacity) &
            result(status)
        type(c_ptr), value, intent(in) :: handle, vectors_pointer, error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: vector_count, basis_count
        integer(c_size_t), value, intent(in) :: error_capacity
        type(stability_problem_context_t), pointer, intent(out) :: context
        real(c_double), pointer, contiguous, intent(out) :: vectors(:, :)
        integer(c_int) :: status
        integer :: count, info, shape(2)

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        status = pressure_context(handle, context, error_pointer, error_capacity)
        if (status /= status_ok) return
        if (fixed_boundary_is_free(context%problem)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "pressure derivatives require a fixed-boundary problem")
            return
        end if
        call fixed_boundary_unknown_count(context%problem, int(parity_class), &
            count, info)
        if (info /= fixed_boundary_ok) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "invalid pressure parity class")
            return
        end if
        status = status_invalid_argument
        if (vector_count /= int(count, c_size_t)) return
        if (basis_count < 1_c_size_t .or. basis_count > vector_count) return
        if (int(count, int64) * int(basis_count, int64) > int(huge(1), int64)) return
        if (.not. c_associated(vectors_pointer)) return
        shape(1) = count
        shape(2) = int(basis_count)
        call c_f_pointer(vectors_pointer, vectors, shape)
        status = status_ok
    end function prepare_basis

    subroutine report_pressure_error(info, status, error_pointer, error_capacity)
        integer, intent(in) :: info
        integer(c_int), intent(out) :: status
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity

        select case (info)
        case (fixed_boundary_allocation_error)
            status = status_allocation_error
        case (fixed_boundary_invalid)
            status = status_invalid_argument
        case default
            status = status_compute_error
        end select
        call write_error(error_pointer, error_capacity, &
            "pressure derivative failed its pressure, resonance or basis contract")
    end subroutine report_pressure_error

end module gliss_pressure_capi
