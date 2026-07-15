module gliss_axisymmetric_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
        c_f_pointer, c_int, c_ptr, c_size_t, c_sizeof
    use axisymmetric_spectrum, only: axisymmetric_spectrum_compute_error, &
        axisymmetric_spectrum_invalid_input, axisymmetric_spectrum_ok, &
        axisymmetric_spectrum_result_t, compute_axisymmetric_spectrum
    use gliss_c_abi_support, only: error_buffer_status, status_compute_error, &
        status_invalid_argument, status_ok, write_error
    use gliss_c_contexts, only: equilibrium_context_t
    implicit none
    private

    type, bind(c) :: axisymmetric_spectrum_result_c
        integer(c_size_t) :: struct_size
        integer(c_int) :: has_eigenpair
        integer(c_int) :: field_periods
        integer(c_int) :: toroidal_mode
        integer(c_int) :: poloidal_max
        integer(c_size_t) :: mode_count
        integer(c_size_t) :: radial_surfaces
        integer(c_int) :: parity_class
        integer(c_int) :: degree
        integer(c_size_t) :: negative_count
        real(c_double) :: lowest_eigenvalue
        real(c_double) :: certificate
        real(c_double) :: eigenpair_residual
        real(c_double) :: force_balance_residual
    end type axisymmetric_spectrum_result_c

    public :: gliss_axisymmetric_spectrum_c

contains

    function gliss_axisymmetric_spectrum_c(equilibrium_handle, &
            toroidal_mode, poloidal_max, degree, solve_eigenpair, &
            result_pointer, error_pointer, error_capacity) bind(c, &
            name="gliss_axisymmetric_spectrum") result(status)
        type(c_ptr), value, intent(in) :: equilibrium_handle, result_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_int), value, intent(in) :: toroidal_mode, poloidal_max
        integer(c_int), value, intent(in) :: degree, solve_eigenpair
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(axisymmetric_spectrum_result_c), pointer :: result
        type(axisymmetric_spectrum_result_t) :: native
        type(equilibrium_context_t), pointer :: equilibrium
        character(len=128) :: message
        integer :: info

        status = prepare_call(equilibrium_handle, result_pointer, &
            error_pointer, error_capacity, equilibrium, result)
        if (status /= status_ok) return
        if (solve_eigenpair /= 0_c_int .and. solve_eigenpair /= 1_c_int) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "solve_eigenpair must be 0 or 1")
            return
        end if
        call compute_axisymmetric_spectrum(equilibrium%equilibrium, &
            int(toroidal_mode), int(poloidal_max), int(degree), &
            solve_eigenpair == 1_c_int, native, info, message)
        select case (info)
        case (axisymmetric_spectrum_ok)
            call fill_result(native, result)
            status = status_ok
        case (axisymmetric_spectrum_invalid_input)
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, trim(message))
        case (axisymmetric_spectrum_compute_error)
            status = status_compute_error
            call write_error(error_pointer, error_capacity, trim(message))
        case default
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "axisymmetric solve returned an unknown status")
        end select
    end function gliss_axisymmetric_spectrum_c

    function prepare_call(equilibrium_handle, result_pointer, error_pointer, &
            error_capacity, equilibrium, result) result(status)
        type(c_ptr), value, intent(in) :: equilibrium_handle, result_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        type(equilibrium_context_t), pointer, intent(out) :: equilibrium
        type(axisymmetric_spectrum_result_c), pointer, intent(out) :: result
        integer(c_int) :: status

        nullify (equilibrium, result)
        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(equilibrium_handle)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        if (.not. c_associated(result_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "axisymmetric result pointer is null")
            return
        end if
        call c_f_pointer(equilibrium_handle, equilibrium)
        call c_f_pointer(result_pointer, result)
        if (.not. associated(equilibrium)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        if (result%struct_size /= c_sizeof(result)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "axisymmetric result struct_size is incompatible")
            return
        end if
        status = status_ok
    end function prepare_call

    subroutine fill_result(native, result)
        type(axisymmetric_spectrum_result_t), intent(in) :: native
        type(axisymmetric_spectrum_result_c), intent(out) :: result

        result%struct_size = c_sizeof(result)
        result%has_eigenpair = merge(1_c_int, 0_c_int, native%has_eigenpair)
        result%field_periods = int(native%field_periods, c_int)
        result%toroidal_mode = int(native%toroidal_mode, c_int)
        result%poloidal_max = int(native%poloidal_max, c_int)
        result%mode_count = int(native%mode_count, c_size_t)
        result%radial_surfaces = int(native%radial_surfaces, c_size_t)
        result%parity_class = int(native%parity_class, c_int)
        result%degree = int(native%degree, c_int)
        result%negative_count = int(native%negative_count, c_size_t)
        result%lowest_eigenvalue = native%lowest_eigenvalue
        result%certificate = native%certificate
        result%eigenpair_residual = native%eigenpair_residual
        result%force_balance_residual = native%force_balance_residual
    end subroutine fill_result

end module gliss_axisymmetric_capi
