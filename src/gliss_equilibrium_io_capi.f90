module gliss_equilibrium_io_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_int, &
        c_ptr, c_size_t
    use gliss_c_abi_support, only: decode_path, error_buffer_status, &
        status_internal_error, status_invalid_argument, status_ok, &
        status_read_error, write_error
    use gliss_c_contexts, only: equilibrium_context_t
    use gvec_cas3d_writer, only: write_gvec_cas3d_file, writer_invalid, &
        writer_netcdf_error, writer_ok, writer_open_error
    implicit none
    private

    public :: gliss_equilibrium_schema_version_c
    public :: gliss_equilibrium_write_c

contains

    function gliss_equilibrium_schema_version_c(handle, version_pointer, &
            error_pointer, error_capacity) bind(c, &
            name="gliss_equilibrium_schema_version") result(status)
        type(c_ptr), value, intent(in) :: handle, version_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        integer(c_int), pointer :: version
        type(equilibrium_context_t), pointer :: context

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(version_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "schema_version output pointer is null")
            return
        end if
        call c_f_pointer(version_pointer, version)
        version = 0_c_int
        status = context_from_handle(handle, context)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        version = int(context%equilibrium%schema_version, c_int)
    end function gliss_equilibrium_schema_version_c

    function gliss_equilibrium_write_c(handle, path_pointer, path_length, &
            error_pointer, error_capacity) bind(c, &
            name="gliss_equilibrium_write") result(status)
        type(c_ptr), value, intent(in) :: handle, path_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: path_length, error_capacity
        integer(c_int) :: status
        type(equilibrium_context_t), pointer :: context
        character(len=:), allocatable :: filename
        character(len=128) :: message
        integer :: info

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        status = context_from_handle(handle, context)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        status = decode_path(path_pointer, path_length, filename, message)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, trim(message))
            return
        end if
        call write_gvec_cas3d_file(filename, context%equilibrium, info)
        select case (info)
        case (writer_ok)
            status = status_ok
        case (writer_invalid)
            status = status_internal_error
            call write_error(error_pointer, error_capacity, &
                "equilibrium context is not serializable")
        case (writer_open_error)
            status = status_read_error
            call write_error(error_pointer, error_capacity, &
                "failed to create equilibrium export; output may exist")
        case (writer_netcdf_error)
            status = status_read_error
            call write_error(error_pointer, error_capacity, &
                "failed to write equilibrium export")
        case default
            status = status_internal_error
            call write_error(error_pointer, error_capacity, &
                "equilibrium export returned an unknown status")
        end select
    end function gliss_equilibrium_write_c

    function context_from_handle(handle, context) result(status)
        type(c_ptr), value, intent(in) :: handle
        type(equilibrium_context_t), pointer, intent(out) :: context
        integer(c_int) :: status

        nullify (context)
        if (.not. c_associated(handle)) then
            status = status_invalid_argument
            return
        end if
        call c_f_pointer(handle, context)
        if (.not. associated(context)) then
            status = status_invalid_argument
            return
        end if
        status = status_ok
    end function context_from_handle

end module gliss_equilibrium_io_capi
