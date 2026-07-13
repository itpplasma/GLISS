module gliss_solver_controls_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
        c_f_pointer, c_int, c_ptr, c_size_t, c_sizeof
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use fixed_boundary_spectrum, only: fixed_boundary_ok, &
        set_fixed_boundary_solver_controls
    use gliss_c_abi_support, only: error_buffer_status, &
        status_invalid_argument, status_ok, write_error
    use gliss_c_contexts, only: stability_problem_context_t
    implicit none
    private

    type, bind(c) :: solver_tolerances_c
        integer(c_size_t) :: struct_size
        real(c_double) :: eigenvalue_relative
        real(c_double) :: residual_relative
        real(c_double) :: negative_bracket_relative
        real(c_double) :: negative_bracket_floor
        integer(c_int) :: inverse_iteration_limit
        integer(c_int) :: bracket_iteration_limit
    end type solver_tolerances_c

    public :: gliss_stability_problem_set_solver_tolerances_c

contains

    function gliss_stability_problem_set_solver_tolerances_c(handle, &
            tolerances_pointer, error_pointer, error_capacity) bind(c, &
            name="gliss_stability_problem_set_solver_tolerances") result(status)
        type(c_ptr), value, intent(in) :: handle, tolerances_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(stability_problem_context_t), pointer :: context
        type(solver_tolerances_c), pointer :: supplied
        type(fixed_boundary_solver_controls_t) :: controls
        integer :: info

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(handle)) then
            call reject("stability problem handle is null")
            return
        end if
        if (.not. c_associated(tolerances_pointer)) then
            call reject("solver tolerances pointer is null")
            return
        end if
        call c_f_pointer(handle, context)
        call c_f_pointer(tolerances_pointer, supplied)
        if (supplied%struct_size /= c_sizeof(supplied)) then
            call reject("solver tolerances struct_size is incompatible")
            return
        end if
        controls = fixed_boundary_solver_controls_t( &
            supplied%eigenvalue_relative, supplied%residual_relative, &
            supplied%negative_bracket_relative, &
            supplied%negative_bracket_floor, &
            int(supplied%inverse_iteration_limit), &
            int(supplied%bracket_iteration_limit))
        call set_fixed_boundary_solver_controls(context%problem, controls, info)
        if (info /= fixed_boundary_ok) then
            call reject("solver tolerances must be finite and positive with iteration limits at least 1")
            return
        end if
        status = status_ok

    contains

        subroutine reject(message)
            character(len=*), intent(in) :: message

            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, message)
        end subroutine reject

    end function gliss_stability_problem_set_solver_tolerances_c

end module gliss_solver_controls_capi
