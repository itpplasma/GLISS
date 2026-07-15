module gliss_spectrum_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
        c_f_pointer, c_int, c_loc, c_null_ptr, c_ptr, c_size_t, c_sizeof
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_allocation_error, &
        fixed_boundary_invalid, fixed_boundary_ok, fixed_boundary_problem_t, &
        fixed_boundary_spectrum_result_t, fixed_boundary_unknown_count, &
        solve_fixed_boundary_class
    use gliss_c_abi_support, only: error_buffer_status, status_allocation_error, &
        status_capacity, status_compute_error, status_invalid_argument, &
        status_internal_error, status_ok, write_error
    use gliss_c_contexts, only: equilibrium_context_t, &
        stability_problem_context_t
    implicit none
    private

    type, bind(c) :: spectrum_summary_c
        integer(c_size_t) :: struct_size
        integer(c_int) :: has_chart_metric
        integer(c_int) :: has_eigenvector
        integer(c_int) :: field_periods
        integer(c_int) :: parity_class
        integer(c_int) :: radial_quadrature
        integer(c_int) :: angular_theta
        integer(c_int) :: angular_zeta
        integer(c_size_t) :: mode_count
        integer(c_size_t) :: unknowns
        integer(c_size_t) :: normal_unknowns
        integer(c_size_t) :: eta_unknowns
        integer(c_size_t) :: mu_unknowns
        integer(c_size_t) :: negative_count
        integer(c_size_t) :: floor_count
        real(c_double) :: adiabatic_index
        real(c_double) :: density_kg_m3
        real(c_double) :: zero_floor
        real(c_double) :: lowest_eigenvalue
        real(c_double) :: certificate
        real(c_double) :: eigenpair_residual
        real(c_double) :: eigenpair_resolution
        real(c_double) :: inertia_interval
    end type spectrum_summary_c

    public :: gliss_stability_problem_create_c
    public :: gliss_stability_problem_destroy_c
    public :: gliss_stability_problem_unknown_count_c
    public :: gliss_stability_problem_solve_class_c

contains

    function gliss_stability_problem_create_c(equilibrium_handle, &
            adiabatic_index, density_kg_m3, zero_floor, mode_count, &
            mode_m_pointer, mode_n_pointer, radial_quadrature, handle_pointer, &
            error_pointer, error_capacity) &
            bind(c, name="gliss_stability_problem_create") result(status)
        type(c_ptr), value, intent(in) :: equilibrium_handle
        real(c_double), value, intent(in) :: adiabatic_index, density_kg_m3
        real(c_double), value, intent(in) :: zero_floor
        integer(c_size_t), value, intent(in) :: mode_count
        type(c_ptr), value, intent(in) :: mode_m_pointer, mode_n_pointer
        integer(c_int), value, intent(in) :: radial_quadrature
        type(c_ptr), value, intent(in) :: handle_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(equilibrium_context_t), pointer :: equilibrium
        type(stability_problem_context_t), pointer :: context
        type(c_ptr), pointer :: handle
        integer, allocatable :: mode_m(:), mode_n(:)
        integer :: allocation_status, info

        status = prepare_output_handle(handle_pointer, error_pointer, &
            error_capacity, handle)
        if (status /= status_ok) return
        status = equilibrium_from_handle(equilibrium_handle, equilibrium)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        status = decode_modes(mode_count, mode_m_pointer, mode_n_pointer, &
            mode_m, mode_n)
        if (status /= status_ok) then
            if (status == status_allocation_error) then
                call write_error(error_pointer, error_capacity, &
                    "failed to allocate mode storage")
            else
                call write_error(error_pointer, error_capacity, &
                    "modes must be nonempty valid int32 arrays")
            end if
            return
        end if
        allocate (context, stat=allocation_status)
        if (allocation_status /= 0) then
            status = status_allocation_error
            call write_error(error_pointer, error_capacity, &
                "failed to allocate stability problem")
            return
        end if
        call build_fixed_boundary_problem(equilibrium%equilibrium, &
            adiabatic_index, density_kg_m3, zero_floor, mode_m, mode_n, &
            int(radial_quadrature), context%problem, info)
        if (info /= fixed_boundary_ok) then
            deallocate (context)
            call report_problem_error(info, status, error_pointer, &
                error_capacity)
            return
        end if
        handle = c_loc(context)
        status = status_ok
    end function gliss_stability_problem_create_c

    function prepare_output_handle(handle_pointer, error_pointer, &
            error_capacity, handle) result(status)
        type(c_ptr), value, intent(in) :: handle_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        type(c_ptr), pointer, intent(out) :: handle
        integer(c_int) :: status

        nullify (handle)
        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(handle_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "stability problem output pointer is null")
            return
        end if
        call c_f_pointer(handle_pointer, handle)
        handle = c_null_ptr
        status = status_ok
    end function prepare_output_handle

    function decode_modes(mode_count, mode_m_pointer, mode_n_pointer, mode_m, &
            mode_n) result(status)
        integer(c_size_t), value, intent(in) :: mode_count
        type(c_ptr), value, intent(in) :: mode_m_pointer, mode_n_pointer
        integer, allocatable, intent(out) :: mode_m(:), mode_n(:)
        integer(c_int), pointer :: c_mode_m(:), c_mode_n(:)
        integer(c_int) :: status
        integer :: allocation_status, count, index, pointer_shape(1)

        status = status_invalid_argument
        if (mode_count < 1_c_size_t) return
        if (mode_count > int(huge(count), c_size_t)) return
        if (.not. c_associated(mode_m_pointer)) return
        if (.not. c_associated(mode_n_pointer)) return
        count = int(mode_count)
        pointer_shape(1) = count
        call c_f_pointer(mode_m_pointer, c_mode_m, pointer_shape)
        call c_f_pointer(mode_n_pointer, c_mode_n, pointer_shape)
        allocate (mode_m(count), mode_n(count), stat=allocation_status)
        if (allocation_status /= 0) then
            status = status_allocation_error
            return
        end if
        do index = 1, count
            mode_m(index) = int(c_mode_m(index))
            mode_n(index) = int(c_mode_n(index))
        end do
        status = status_ok
    end function decode_modes

    function gliss_stability_problem_destroy_c(handle_pointer, error_pointer, &
            error_capacity) bind(c, name="gliss_stability_problem_destroy") &
            result(status)
        type(c_ptr), value, intent(in) :: handle_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(c_ptr), pointer :: handle
        type(stability_problem_context_t), pointer :: context
        integer :: allocation_status

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(handle_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "stability problem handle pointer is null")
            return
        end if
        call c_f_pointer(handle_pointer, handle)
        if (.not. c_associated(handle)) then
            status = status_ok
            return
        end if
        call c_f_pointer(handle, context)
        deallocate (context, stat=allocation_status)
        if (allocation_status /= 0) then
            status = status_internal_error
            call write_error(error_pointer, error_capacity, &
                "failed to destroy stability problem")
            return
        end if
        handle = c_null_ptr
        status = status_ok
    end function gliss_stability_problem_destroy_c

    function gliss_stability_problem_unknown_count_c(handle, parity_class, &
            unknown_pointer, error_pointer, error_capacity) &
            bind(c, name="gliss_stability_problem_unknown_count") &
            result(status)
        type(c_ptr), value, intent(in) :: handle, unknown_pointer, error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        integer(c_size_t), pointer :: unknowns
        type(stability_problem_context_t), pointer :: context
        integer :: count, info

        status = prepare_required_output(unknown_pointer, error_pointer, &
            error_capacity, "unknown-count output pointer is null")
        if (status /= status_ok) return
        call c_f_pointer(unknown_pointer, unknowns)
        unknowns = 0_c_size_t
        status = problem_from_handle(handle, context)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, &
                "stability problem handle is null")
            return
        end if
        call fixed_boundary_unknown_count(context%problem, int(parity_class), &
            count, info)
        if (info /= fixed_boundary_ok) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "parity_class must be 1 or 2")
            return
        end if
        unknowns = int(count, c_size_t)
        status = status_ok
    end function gliss_stability_problem_unknown_count_c

    function prepare_required_output(output_pointer, error_pointer, &
            error_capacity, message) result(status)
        type(c_ptr), value, intent(in) :: output_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        character(len=*), intent(in) :: message
        integer(c_int) :: status

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(output_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, message)
            return
        end if
        status = status_ok
    end function prepare_required_output

    function gliss_stability_problem_solve_class_c(handle, parity_class, &
            capacity, vector_pointer, written_pointer, summary_pointer, &
            error_pointer, error_capacity) &
            bind(c, name="gliss_stability_problem_solve_class") result(status)
        type(c_ptr), value, intent(in) :: handle, vector_pointer
        type(c_ptr), value, intent(in) :: written_pointer, summary_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_int), value, intent(in) :: parity_class
        integer(c_size_t), value, intent(in) :: capacity, error_capacity
        integer(c_int) :: status
        integer(c_size_t), pointer :: written
        type(spectrum_summary_c), pointer :: summary
        type(stability_problem_context_t), pointer :: context
        type(fixed_boundary_spectrum_result_t) :: result
        real(c_double), pointer :: vector(:)
        integer :: info, pointer_shape(1), required

        status = prepare_solve_outputs(written_pointer, summary_pointer, &
            error_pointer, error_capacity, written, summary)
        if (status /= status_ok) return
        status = problem_from_handle(handle, context)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, &
                "stability problem handle is null")
            return
        end if
        call solve_fixed_boundary_class(context%problem, int(parity_class), &
            result, info)
        if (info /= fixed_boundary_ok) then
            call report_problem_error(info, status, error_pointer, &
                error_capacity)
            return
        end if
        call fill_summary(result, summary)
        required = size(result%eigenvector)
        written = int(required, c_size_t)
        if (capacity < written) then
            status = status_capacity
            call write_error(error_pointer, error_capacity, &
                "eigenvector capacity is smaller than the unknown count")
            return
        end if
        if (required == 0) then
            status = status_ok
            return
        end if
        if (.not. c_associated(vector_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "eigenvector output pointer is null")
            return
        end if
        pointer_shape(1) = required
        call c_f_pointer(vector_pointer, vector, pointer_shape)
        vector = result%eigenvector
        status = status_ok
    end function gliss_stability_problem_solve_class_c


    function prepare_solve_outputs(written_pointer, summary_pointer, &
            error_pointer, error_capacity, written, summary) result(status)
        type(c_ptr), value, intent(in) :: written_pointer, summary_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_size_t), pointer, intent(out) :: written
        type(spectrum_summary_c), pointer, intent(out) :: summary
        integer(c_int) :: status

        nullify (written, summary)
        status = prepare_required_output(written_pointer, error_pointer, &
            error_capacity, "written output pointer is null")
        if (status /= status_ok) return
        call c_f_pointer(written_pointer, written)
        written = 0_c_size_t
        if (.not. c_associated(summary_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "spectrum summary pointer is null")
            return
        end if
        call c_f_pointer(summary_pointer, summary)
        if (summary%struct_size /= c_sizeof(summary)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "spectrum summary struct_size is incompatible")
            return
        end if
        status = status_ok
    end function prepare_solve_outputs

    subroutine fill_summary(result, summary)
        type(fixed_boundary_spectrum_result_t), intent(in) :: result
        type(spectrum_summary_c), intent(out) :: summary

        summary%struct_size = c_sizeof(summary)
        summary%has_chart_metric = merge(1_c_int, 0_c_int, &
            result%has_chart_metric)
        summary%has_eigenvector = merge(1_c_int, 0_c_int, &
            result%has_eigenvector)
        summary%field_periods = int(result%field_periods, c_int)
        summary%parity_class = int(result%parity_class, c_int)
        summary%radial_quadrature = int(result%radial_quadrature, c_int)
        summary%angular_theta = int(result%angular_theta, c_int)
        summary%angular_zeta = int(result%angular_zeta, c_int)
        summary%mode_count = int(result%mode_count, c_size_t)
        summary%unknowns = int(result%unknowns, c_size_t)
        summary%normal_unknowns = int(result%normal_unknowns, c_size_t)
        summary%eta_unknowns = int(result%eta_unknowns, c_size_t)
        summary%mu_unknowns = int(result%mu_unknowns, c_size_t)
        summary%negative_count = int(result%negative_count, c_size_t)
        summary%floor_count = int(result%floor_count, c_size_t)
        summary%adiabatic_index = result%adiabatic_index
        summary%density_kg_m3 = result%density_kg_m3
        summary%zero_floor = result%zero_floor
        summary%lowest_eigenvalue = result%lowest_eigenvalue
        summary%certificate = result%certificate
        summary%eigenpair_residual = result%eigenpair_residual
        summary%eigenpair_resolution = result%eigenpair_resolution
        summary%inertia_interval = result%inertia_interval
    end subroutine fill_summary

    function equilibrium_from_handle(handle, context) result(status)
        type(c_ptr), value, intent(in) :: handle
        type(equilibrium_context_t), pointer, intent(out) :: context
        integer(c_int) :: status

        nullify (context)
        if (.not. c_associated(handle)) then
            status = status_invalid_argument
            return
        end if
        call c_f_pointer(handle, context)
        status = status_ok
    end function equilibrium_from_handle

    function problem_from_handle(handle, context) result(status)
        type(c_ptr), value, intent(in) :: handle
        type(stability_problem_context_t), pointer, intent(out) :: context
        integer(c_int) :: status

        nullify (context)
        if (.not. c_associated(handle)) then
            status = status_invalid_argument
            return
        end if
        call c_f_pointer(handle, context)
        status = status_ok
    end function problem_from_handle

    subroutine report_problem_error(info, status, error_pointer, error_capacity)
        integer, intent(in) :: info
        integer(c_int), intent(out) :: status
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity

        if (info == fixed_boundary_invalid) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "invalid fixed-boundary stability configuration")
        else if (info == fixed_boundary_allocation_error) then
            status = status_allocation_error
            call write_error(error_pointer, error_capacity, &
                "failed to allocate fixed-boundary spectrum storage")
        else
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "fixed-boundary stability computation failed")
        end if
    end subroutine report_problem_error

end module gliss_spectrum_capi
