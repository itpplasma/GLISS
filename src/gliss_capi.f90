module gliss_capi
    use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, &
        c_null_char
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: compute_mercier, mercier_ok, &
        mercier_result_t
    implicit none
    private

    character(len=*), parameter :: version_string = "0.1.0"
    integer(c_int), parameter :: abi_version_number = 1

    public :: gliss_version_c
    public :: gliss_abi_version_c
    public :: gliss_mercier_profile_c

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
        status = 4
        if (path_length < 1) return
        if (n_theta < 1 .or. n_zeta < 1) return
        if (capacity < 0) return
        allocate (character(len=path_length) :: filename)
        do i = 1, path_length
            filename(i:i) = path(i)
        end do
        call read_gvec_cas3d_file(filename, equilibrium, info)
        if (info /= reader_ok) then
            status = 1
            return
        end if
        call compute_mercier(equilibrium, int(n_theta), int(n_zeta), &
            mercier, info)
        if (info /= mercier_ok) then
            status = 2
            return
        end if
        surfaces = int(size(mercier%s), c_int)
        if (surfaces > capacity) then
            status = 3
            return
        end if
        do i = 1, size(mercier%s)
            s_values(i) = mercier%s(i)
            d_mercier(i) = mercier%d_mercier(i)
        end do
        status = 0
    end subroutine gliss_mercier_profile_c

end module gliss_capi
