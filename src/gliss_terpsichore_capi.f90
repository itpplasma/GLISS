module gliss_terpsichore_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_f_pointer, &
        c_int, c_ptr, c_size_t, c_sizeof
    use gliss_c_abi_support, only: decode_path, error_buffer_status, &
        status_compute_error, status_invalid_argument, status_ok, &
        status_read_error, write_error
    use terpsichore_fixed_boundary_spectrum, only: &
        solve_terpsichore_fixed_boundary_file, &
        terpsichore_fixed_boundary_result_t, &
        terpsichore_fixed_spectrum_compute_error, &
        terpsichore_fixed_spectrum_ok, terpsichore_fixed_spectrum_read_error
    implicit none
    private

    type, bind(c) :: terpsichore_fixed_boundary_result_c
        integer(c_size_t) :: struct_size
        integer(c_size_t) :: unknowns
        integer(c_size_t) :: negative_count
        real(c_double) :: eigenvalue
        real(c_double) :: certificate
        real(c_double) :: residual
        real(c_double) :: resolution
    end type terpsichore_fixed_boundary_result_c

    public :: gliss_terpsichore_fixed_boundary_c

contains

    function gliss_terpsichore_fixed_boundary_c(path_pointer, path_length, &
            result_pointer, error_pointer, error_capacity) bind(c, &
            name="gliss_terpsichore_fixed_boundary") result(status)
        type(c_ptr), value, intent(in) :: path_pointer, result_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: path_length, error_capacity
        integer(c_int) :: status
        type(terpsichore_fixed_boundary_result_c), pointer :: result
        type(terpsichore_fixed_boundary_result_t) :: native
        character(len=:), allocatable :: filename
        character(len=128) :: message
        integer :: info

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(result_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "TERPSICHORE result pointer is null")
            return
        end if
        call c_f_pointer(result_pointer, result)
        if (result%struct_size /= c_sizeof(result)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "TERPSICHORE result struct_size is incompatible")
            return
        end if
        status = decode_path(path_pointer, path_length, filename, message)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, trim(message))
            return
        end if
        call solve_terpsichore_fixed_boundary_file(filename, native, info, &
            message)
        select case (info)
        case (terpsichore_fixed_spectrum_ok)
            call fill_result(native, result)
            status = status_ok
        case (terpsichore_fixed_spectrum_read_error)
            status = status_read_error
            call write_error(error_pointer, error_capacity, trim(message))
        case (terpsichore_fixed_spectrum_compute_error)
            status = status_compute_error
            call write_error(error_pointer, error_capacity, trim(message))
        case default
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "TERPSICHORE solve returned an unknown status")
        end select
    end function gliss_terpsichore_fixed_boundary_c

    subroutine fill_result(native, result)
        type(terpsichore_fixed_boundary_result_t), intent(in) :: native
        type(terpsichore_fixed_boundary_result_c), intent(out) :: result

        result%struct_size = c_sizeof(result)
        result%unknowns = int(native%unknowns, c_size_t)
        result%negative_count = int(native%negative_count, c_size_t)
        result%eigenvalue = native%eigenvalue
        result%certificate = native%certificate
        result%residual = native%residual
        result%resolution = native%resolution
    end subroutine fill_result

end module gliss_terpsichore_capi
