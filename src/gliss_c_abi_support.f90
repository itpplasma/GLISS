module gliss_c_abi_support
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_pointer, &
        c_int, c_null_char, c_ptr, c_size_t
    implicit none
    private

    integer(c_int), parameter, public :: status_ok = 0
    integer(c_int), parameter, public :: status_read_error = 1
    integer(c_int), parameter, public :: status_compute_error = 2
    integer(c_int), parameter, public :: status_capacity = 3
    integer(c_int), parameter, public :: status_invalid_argument = 4
    integer(c_int), parameter, public :: status_allocation_error = 5
    integer(c_int), parameter, public :: status_internal_error = 6

    public :: error_buffer_status
    public :: decode_path
    public :: write_error

contains

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

end module gliss_c_abi_support
