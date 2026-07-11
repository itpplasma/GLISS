module gliss_capi
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    implicit none
    private

    character(len=*), parameter :: version_string = "0.1.0"
    integer(c_int), parameter :: abi_version_number = 1

    public :: gliss_version_c
    public :: gliss_abi_version_c

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

end module gliss_capi
