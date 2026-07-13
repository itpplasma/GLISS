module gliss_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
        c_f_pointer, c_int, c_loc, c_null_char, c_null_ptr, c_ptr, c_size_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: compute_mercier, mercier_ok, &
        mercier_result_t
    implicit none
    private

    character(len=*), parameter :: version_string = "0.0.1"
    integer(c_int), parameter :: abi_version_number = 1
    integer(c_int), parameter :: status_ok = 0
    integer(c_int), parameter :: status_read_error = 1
    integer(c_int), parameter :: status_compute_error = 2
    integer(c_int), parameter :: status_capacity = 3
    integer(c_int), parameter :: status_invalid_argument = 4
    integer(c_int), parameter :: status_allocation_error = 5
    integer(c_int), parameter :: status_internal_error = 6

    type :: equilibrium_context_t
        type(gvec_cas3d_equilibrium_t) :: equilibrium
    end type equilibrium_context_t

    public :: gliss_version_c
    public :: gliss_abi_version_c
    public :: gliss_mercier_profile_c
    public :: gliss_equilibrium_create_c
    public :: gliss_equilibrium_destroy_c
    public :: gliss_equilibrium_surface_count_c
    public :: gliss_mercier_profile_context_c

contains

    subroutine gliss_version_c(buffer, length) bind(c, name="gliss_version")
        ! Contract: buffer is length bytes owned by the caller; the
        ! written string is truncated to fit and always null-terminated.
        character(c_char), intent(out) :: buffer(*)
        integer(c_int), value, intent(in) :: length
        integer :: index, copy_length

        if (length < 1) return
        copy_length = min(length - 1, len(version_string))
        do index = 1, copy_length
            buffer(index) = version_string(index:index)
        end do
        buffer(copy_length + 1) = c_null_char
    end subroutine gliss_version_c

    function gliss_abi_version_c() bind(c, name="gliss_abi_version") &
            result(version)
        integer(c_int) :: version

        version = abi_version_number
    end function gliss_abi_version_c

    subroutine gliss_mercier_profile_c(path, path_length, n_theta, n_zeta, &
            capacity, surfaces, s_values, d_mercier, status) &
            bind(c, name="gliss_mercier_profile")
        character(kind=c_char), intent(in) :: path(*)
        integer(c_int), value, intent(in) :: path_length, n_theta, n_zeta
        integer(c_int), value, intent(in) :: capacity
        integer(c_int), intent(out) :: surfaces, status
        real(c_double), intent(out) :: s_values(*), d_mercier(*)
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        type(mercier_result_t) :: mercier
        character(len=:), allocatable :: filename
        integer :: info, i

        surfaces = 0
        status = status_invalid_argument
        if (path_length < 1) return
        if (n_theta < 1 .or. n_zeta < 1) return
        if (capacity < 0) return
        allocate (character(len=path_length) :: filename)
        do i = 1, path_length
            filename(i:i) = path(i)
        end do
        call read_gvec_cas3d_file(filename, equilibrium, info)
        if (info /= reader_ok) then
            status = status_read_error
            return
        end if
        call compute_mercier(equilibrium, int(n_theta), int(n_zeta), &
            mercier, info)
        if (info /= mercier_ok) then
            status = status_compute_error
            return
        end if
        surfaces = int(size(mercier%s), c_int)
        if (surfaces > capacity) then
            status = status_capacity
            return
        end if
        do i = 1, size(mercier%s)
            s_values(i) = mercier%s(i)
            d_mercier(i) = mercier%d_mercier(i)
        end do
        status = status_ok
    end subroutine gliss_mercier_profile_c

    function gliss_equilibrium_create_c(path_pointer, path_length, &
            handle_pointer, &
            error_pointer, error_capacity) &
            bind(c, name="gliss_equilibrium_create") result(status)
        type(c_ptr), value, intent(in) :: path_pointer, handle_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: path_length, error_capacity
        integer(c_int) :: status
        type(c_ptr), pointer :: handle
        type(equilibrium_context_t), pointer :: context
        character(len=:), allocatable :: filename
        character(len=128) :: message
        integer :: allocation_status, info

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(handle_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "equilibrium output pointer is null")
            return
        end if
        call c_f_pointer(handle_pointer, handle)
        handle = c_null_ptr
        status = decode_path(path_pointer, path_length, filename, message)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, trim(message))
            return
        end if
        allocate (context, stat=allocation_status)
        if (allocation_status /= 0) then
            status = status_allocation_error
            call write_error(error_pointer, error_capacity, &
                "failed to allocate equilibrium context")
            return
        end if
        call read_gvec_cas3d_file(filename, context%equilibrium, info)
        if (info /= reader_ok) then
            deallocate (context)
            status = status_read_error
            call write_error(error_pointer, error_capacity, &
                "failed to read equilibrium export")
            return
        end if
        handle = c_loc(context)
        status = status_ok
    end function gliss_equilibrium_create_c

    function gliss_equilibrium_destroy_c(handle_pointer, error_pointer, &
            error_capacity) bind(c, name="gliss_equilibrium_destroy") &
            result(status)
        type(c_ptr), value, intent(in) :: handle_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(c_ptr), pointer :: handle
        type(equilibrium_context_t), pointer :: context
        integer :: allocation_status

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(handle_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle pointer is null")
            return
        end if
        call c_f_pointer(handle_pointer, handle)
        if (.not. c_associated(handle)) then
            status = status_ok
            return
        end if
        call c_f_pointer(handle, context)
        if (.not. associated(context)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is invalid")
            return
        end if
        deallocate (context, stat=allocation_status)
        if (allocation_status /= 0) then
            status = status_internal_error
            call write_error(error_pointer, error_capacity, &
                "failed to destroy equilibrium context")
            return
        end if
        handle = c_null_ptr
        status = status_ok
    end function gliss_equilibrium_destroy_c

    function gliss_equilibrium_surface_count_c(handle, surfaces_pointer, &
            error_pointer, error_capacity) &
            bind(c, name="gliss_equilibrium_surface_count") result(status)
        type(c_ptr), value, intent(in) :: handle, surfaces_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        integer(c_size_t), pointer :: surfaces
        type(equilibrium_context_t), pointer :: context

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        if (.not. c_associated(surfaces_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "surface_count output pointer is null")
            return
        end if
        call c_f_pointer(surfaces_pointer, surfaces)
        surfaces = 0_c_size_t
        call write_error(error_pointer, error_capacity, "")
        status = context_from_handle(handle, context)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        surfaces = int(size(context%equilibrium%s), c_size_t)
    end function gliss_equilibrium_surface_count_c

    function gliss_mercier_profile_context_c(handle, n_theta, n_zeta, &
            capacity, s_pointer, d_pointer, written, error_pointer, &
            error_capacity) bind(c, name="gliss_mercier_profile_context") &
            result(status)
        type(c_ptr), value, intent(in) :: handle, s_pointer, d_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_int), value, intent(in) :: n_theta, n_zeta
        integer(c_size_t), value, intent(in) :: capacity, error_capacity
        type(c_ptr), value, intent(in) :: written
        integer(c_int) :: status
        type(equilibrium_context_t), pointer :: context
        type(mercier_result_t) :: mercier
        real(c_double), pointer :: s_values(:), d_mercier(:)
        integer(c_size_t), pointer :: written_value
        integer :: info, result_size

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        if (.not. c_associated(written)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "written output pointer is null")
            return
        end if
        call c_f_pointer(written, written_value)
        written_value = 0_c_size_t
        call write_error(error_pointer, error_capacity, "")
        status = context_from_handle(handle, context)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        if (n_theta < 1 .or. n_zeta < 1) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "n_theta and n_zeta must be positive")
            return
        end if
        call compute_mercier(context%equilibrium, int(n_theta), int(n_zeta), &
            mercier, info)
        if (info /= mercier_ok) then
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "Mercier stability computation failed")
            return
        end if
        result_size = size(mercier%s)
        written_value = int(result_size, c_size_t)
        if (capacity < written_value) then
            status = status_capacity
            call write_error(error_pointer, error_capacity, &
                "output capacity is smaller than the surface count")
            return
        end if
        if (result_size == 0) then
            status = status_ok
            return
        end if
        if (.not. c_associated(s_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "s_values output pointer is null")
            return
        end if
        if (.not. c_associated(d_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "d_mercier output pointer is null")
            return
        end if
        call c_f_pointer(s_pointer, s_values, [result_size])
        call c_f_pointer(d_pointer, d_mercier, [result_size])
        s_values = mercier%s
        d_mercier = mercier%d_mercier
        status = status_ok
    end function gliss_mercier_profile_context_c

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

    function decode_path(path_pointer, path_length, filename, message) &
            result(status)
        type(c_ptr), value, intent(in) :: path_pointer
        integer(c_size_t), value, intent(in) :: path_length
        character(len=:), allocatable, intent(out) :: filename
        character(len=*), intent(out) :: message
        integer(c_int) :: status
        character(c_char), pointer :: path(:)
        integer :: allocation_status, i, length

        message = ""
        if (.not. c_associated(path_pointer)) then
            status = status_invalid_argument
            message = "path pointer is null"
            return
        end if
        if (path_length < 1_c_size_t) then
            status = status_invalid_argument
            message = "path must not be empty"
            return
        end if
        if (path_length > int(huge(length), c_size_t)) then
            status = status_invalid_argument
            message = "path length exceeds the Fortran string limit"
            return
        end if
        length = int(path_length)
        call c_f_pointer(path_pointer, path, [length])
        allocate (character(len=length) :: filename, stat=allocation_status)
        if (allocation_status /= 0) then
            status = status_allocation_error
            message = "failed to allocate path storage"
            return
        end if
        do i = 1, length
            if (path(i) == c_null_char) then
                deallocate (filename)
                status = status_invalid_argument
                message = "path contains an embedded null byte"
                return
            end if
            filename(i:i) = path(i)
        end do
        status = status_ok
    end function decode_path

    subroutine write_error(buffer_pointer, capacity, message)
        type(c_ptr), value, intent(in) :: buffer_pointer
        integer(c_size_t), value, intent(in) :: capacity
        character(len=*), intent(in) :: message
        character(c_char), pointer :: buffer(:)
        integer :: copy_length, i, mapped_length

        if (.not. c_associated(buffer_pointer)) return
        if (capacity < 1_c_size_t) return
        mapped_length = int(min(capacity, &
            int(len_trim(message) + 1, c_size_t)))
        call c_f_pointer(buffer_pointer, buffer, [mapped_length])
        copy_length = mapped_length - 1
        do i = 1, copy_length
            buffer(i) = message(i:i)
        end do
        buffer(mapped_length) = c_null_char
    end subroutine write_error

    function error_buffer_status(buffer_pointer, capacity) result(status)
        type(c_ptr), value, intent(in) :: buffer_pointer
        integer(c_size_t), value, intent(in) :: capacity
        integer(c_int) :: status

        if (capacity < 1_c_size_t) then
            status = status_ok
            return
        end if
        if (.not. c_associated(buffer_pointer)) then
            status = status_invalid_argument
            return
        end if
        status = status_ok
    end function error_buffer_status

end module gliss_capi
