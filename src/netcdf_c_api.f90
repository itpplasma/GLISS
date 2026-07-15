module netcdf_c_api
    use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_loc, &
        c_null_char, c_null_ptr, c_size_t
    use netcdf_c_bindings
    implicit none
    private

    public :: nc_noerr
    public :: nc_enotatt
    public :: nc_enotvar
    public :: nc_global
    public :: nc_double
    public :: nc_int64

    public :: nc_close_file
    public :: nc_create_netcdf4
    public :: nc_create_netcdf4_exclusive
    public :: nc_def_dimension
    public :: nc_def_scalar
    public :: nc_def_variable
    public :: nc_end_definitions
    public :: nc_get_global_text
    public :: nc_get_integer
    public :: nc_get_real
    public :: nc_inquire_dimension_id
    public :: nc_inquire_dimension_length
    public :: nc_inquire_dimension_name
    public :: nc_inquire_variable_dimensions
    public :: nc_inquire_variable_id
    public :: nc_inquire_variable_rank
    public :: nc_open_read
    public :: nc_open_write
    public :: nc_put_global_text
    public :: nc_put_integer
    public :: nc_put_real

    interface nc_get_integer
        module procedure nc_get_integer_scalar
        module procedure nc_get_integer_vector
    end interface

    interface nc_get_real
        module procedure nc_get_real_scalar
        module procedure nc_get_real_vector
        module procedure nc_get_real_tensor
    end interface

    interface nc_put_integer
        module procedure nc_put_integer_scalar
        module procedure nc_put_integer_vector
    end interface

    interface nc_put_real
        module procedure nc_put_real_scalar
        module procedure nc_put_real_vector
        module procedure nc_put_real_tensor
    end interface

contains

    integer function nc_open_read(path, ncid) result(status)
        character(len=*), intent(in) :: path
        integer, intent(out) :: ncid
        character(c_char), allocatable :: c_path(:)
        integer(c_int) :: c_ncid

        call build_c_string(path, c_path)
        status = c_nc_open(c_path, 0_c_int, c_ncid)
        ncid = int(c_ncid)
    end function nc_open_read

    integer function nc_open_write(path, ncid) result(status)
        character(len=*), intent(in) :: path
        integer, intent(out) :: ncid
        character(c_char), allocatable :: c_path(:)
        integer(c_int) :: c_ncid

        call build_c_string(path, c_path)
        status = c_nc_open(c_path, 1_c_int, c_ncid)
        ncid = int(c_ncid)
    end function nc_open_write

    integer function nc_create_netcdf4(path, ncid) result(status)
        character(len=*), intent(in) :: path
        integer, intent(out) :: ncid
        character(c_char), allocatable :: c_path(:)
        integer(c_int) :: c_ncid

        call build_c_string(path, c_path)
        status = c_nc_create(c_path, int(z'1000', c_int), c_ncid)
        ncid = int(c_ncid)
    end function nc_create_netcdf4

    integer function nc_create_netcdf4_exclusive(path, ncid) result(status)
        character(len=*), intent(in) :: path
        integer, intent(out) :: ncid
        character(c_char), allocatable :: c_path(:)
        integer(c_int) :: c_ncid

        call build_c_string(path, c_path)
        status = c_nc_create(c_path, int(z'1004', c_int), c_ncid)
        ncid = int(c_ncid)
    end function nc_create_netcdf4_exclusive

    integer function nc_close_file(ncid) result(status)
        integer, intent(in) :: ncid

        status = c_nc_close(int(ncid, c_int))
    end function nc_close_file

    integer function nc_inquire_dimension_id(ncid, name, dimid) result(status)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: dimid
        character(c_char), allocatable :: c_name(:)
        integer(c_int) :: c_dimid

        call build_c_string(name, c_name)
        status = c_nc_inq_dimid(int(ncid, c_int), c_name, c_dimid)
        dimid = int(c_dimid)
    end function nc_inquire_dimension_id

    integer function nc_inquire_dimension_length(ncid, dimid, length) &
            result(status)
        integer, intent(in) :: ncid, dimid
        integer, intent(out) :: length
        integer(c_size_t) :: c_length

        status = c_nc_inq_dimlen(int(ncid, c_int), int(dimid, c_int), c_length)
        length = int(c_length)
    end function nc_inquire_dimension_length

    integer function nc_inquire_dimension_name(ncid, dimid, name) result(status)
        integer, intent(in) :: ncid, dimid
        character(len=*), intent(out) :: name
        character(c_char), target :: c_name(257)

        c_name = c_null_char
        status = c_nc_inq_dimname(int(ncid, c_int), int(dimid, c_int), c_name)
        if (status == nc_noerr) call copy_c_string(c_name, name)
    end function nc_inquire_dimension_name

    integer function nc_inquire_variable_id(ncid, name, varid) result(status)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: varid
        character(c_char), allocatable :: c_name(:)
        integer(c_int) :: c_varid

        call build_c_string(name, c_name)
        status = c_nc_inq_varid(int(ncid, c_int), c_name, c_varid)
        varid = int(c_varid)
    end function nc_inquire_variable_id

    integer function nc_inquire_variable_rank(ncid, varid, rank) result(status)
        integer, intent(in) :: ncid, varid
        integer, intent(out) :: rank
        integer(c_int) :: c_rank

        status = c_nc_inq_varndims(int(ncid, c_int), int(varid, c_int), c_rank)
        rank = int(c_rank)
    end function nc_inquire_variable_rank

    integer function nc_inquire_variable_dimensions(ncid, varid, dimids) &
            result(status)
        integer, intent(in) :: ncid, varid
        integer, intent(out) :: dimids(:)
        integer(c_int), allocatable, target :: c_dimids(:)

        allocate (c_dimids(size(dimids)))
        status = c_nc_inq_vardimid(int(ncid, c_int), int(varid, c_int), &
            c_loc(c_dimids(1)))
        if (status == nc_noerr) dimids = int(c_dimids)
    end function nc_inquire_variable_dimensions

    integer function nc_get_integer_scalar(ncid, varid, value) result(status)
        integer, intent(in) :: ncid, varid
        integer, intent(out) :: value
        integer(c_int), target :: c_value

        status = c_nc_get_var_int(int(ncid, c_int), int(varid, c_int), &
            c_loc(c_value))
        if (status == nc_noerr) value = int(c_value)
    end function nc_get_integer_scalar

    integer function nc_get_integer_vector(ncid, varid, values) result(status)
        integer, intent(in) :: ncid, varid
        integer, intent(out) :: values(:)
        integer(c_int), allocatable, target :: c_values(:)

        allocate (c_values(size(values)))
        status = c_nc_get_var_int(int(ncid, c_int), int(varid, c_int), &
            c_loc(c_values(1)))
        if (status == nc_noerr) values = int(c_values)
    end function nc_get_integer_vector

    integer function nc_get_real_scalar(ncid, varid, value) result(status)
        integer, intent(in) :: ncid, varid
        real(c_double), intent(out), target :: value

        status = c_nc_get_var_double(int(ncid, c_int), int(varid, c_int), &
            c_loc(value))
    end function nc_get_real_scalar

    integer function nc_get_real_vector(ncid, varid, values) result(status)
        integer, intent(in) :: ncid, varid
        real(c_double), intent(out), target, contiguous :: values(:)

        status = c_nc_get_var_double(int(ncid, c_int), int(varid, c_int), &
            c_loc(values(1)))
    end function nc_get_real_vector

    integer function nc_get_real_tensor(ncid, varid, values) result(status)
        integer, intent(in) :: ncid, varid
        real(c_double), intent(out), target, contiguous :: values(:, :, :)

        status = c_nc_get_var_double(int(ncid, c_int), int(varid, c_int), &
            c_loc(values(1, 1, 1)))
    end function nc_get_real_tensor

    integer function nc_get_global_text(ncid, name, value) result(status)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        character(c_char), allocatable, target :: buffer(:)
        character(c_char), allocatable :: c_name(:)
        integer(c_size_t) :: length
        integer :: index

        value = ""
        call build_c_string(name, c_name)
        status = c_nc_inq_attlen(int(ncid, c_int), int(nc_global, c_int), &
            c_name, length)
        if (status /= nc_noerr) return
        if (length > len(value)) then
            status = -1
            return
        end if
        if (length == 0) return
        allocate (buffer(int(length)))
        status = c_nc_get_att_text(int(ncid, c_int), int(nc_global, c_int), &
            c_name, c_loc(buffer(1)))
        if (status /= nc_noerr) return
        do index = 1, int(length)
            value(index:index) = buffer(index)
        end do
    end function nc_get_global_text

    integer function nc_def_dimension(ncid, name, length, dimid) result(status)
        integer, intent(in) :: ncid, length
        character(len=*), intent(in) :: name
        integer, intent(out) :: dimid
        character(c_char), allocatable :: c_name(:)
        integer(c_int) :: c_dimid

        call build_c_string(name, c_name)
        status = c_nc_def_dim(int(ncid, c_int), c_name, &
            int(length, c_size_t), c_dimid)
        dimid = int(c_dimid)
    end function nc_def_dimension

    integer function nc_def_scalar(ncid, name, xtype, varid) result(status)
        integer, intent(in) :: ncid, xtype
        character(len=*), intent(in) :: name
        integer, intent(out) :: varid
        character(c_char), allocatable :: c_name(:)
        integer(c_int) :: c_varid

        call build_c_string(name, c_name)
        status = c_nc_def_var(int(ncid, c_int), c_name, int(xtype, c_int), &
            0_c_int, c_null_ptr, c_varid)
        varid = int(c_varid)
    end function nc_def_scalar

    integer function nc_def_variable(ncid, name, xtype, dimids, varid) &
            result(status)
        integer, intent(in) :: ncid, xtype, dimids(:)
        character(len=*), intent(in) :: name
        integer, intent(out) :: varid
        character(c_char), allocatable :: c_name(:)
        integer(c_int), allocatable, target :: c_dimids(:)
        integer(c_int) :: c_varid

        call build_c_string(name, c_name)
        c_dimids = int(dimids, c_int)
        status = c_nc_def_var(int(ncid, c_int), c_name, int(xtype, c_int), &
            int(size(dimids), c_int), c_loc(c_dimids(1)), c_varid)
        varid = int(c_varid)
    end function nc_def_variable

    integer function nc_put_global_text(ncid, name, value) result(status)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name, value
        character(c_char), allocatable, target :: c_value(:)
        character(c_char), allocatable :: c_name(:)
        integer :: index, length

        length = len_trim(value)
        call build_c_string(name, c_name)
        if (length == 0) then
            status = c_nc_put_att_text(int(ncid, c_int), &
                int(nc_global, c_int), c_name, 0_c_size_t, c_null_ptr)
            return
        end if
        allocate (c_value(length))
        do index = 1, length
            c_value(index) = value(index:index)
        end do
        status = c_nc_put_att_text(int(ncid, c_int), int(nc_global, c_int), &
            c_name, int(length, c_size_t), c_loc(c_value(1)))
    end function nc_put_global_text

    integer function nc_end_definitions(ncid) result(status)
        integer, intent(in) :: ncid

        status = c_nc_enddef(int(ncid, c_int))
    end function nc_end_definitions

    integer function nc_put_integer_scalar(ncid, varid, value) result(status)
        integer, intent(in) :: ncid, varid
        integer, intent(in) :: value
        integer(c_int), target :: c_value

        c_value = int(value, c_int)
        status = c_nc_put_var_int(int(ncid, c_int), int(varid, c_int), &
            c_loc(c_value))
    end function nc_put_integer_scalar

    integer function nc_put_integer_vector(ncid, varid, values) result(status)
        integer, intent(in) :: ncid, varid
        integer, intent(in) :: values(:)
        integer(c_int), allocatable, target :: c_values(:)

        c_values = int(values, c_int)
        status = c_nc_put_var_int(int(ncid, c_int), int(varid, c_int), &
            c_loc(c_values(1)))
    end function nc_put_integer_vector

    integer function nc_put_real_scalar(ncid, varid, value) result(status)
        integer, intent(in) :: ncid, varid
        real(c_double), intent(in), target :: value

        status = c_nc_put_var_double(int(ncid, c_int), int(varid, c_int), &
            c_loc(value))
    end function nc_put_real_scalar

    integer function nc_put_real_vector(ncid, varid, values) result(status)
        integer, intent(in) :: ncid, varid
        real(c_double), intent(in), target, contiguous :: values(:)

        status = c_nc_put_var_double(int(ncid, c_int), int(varid, c_int), &
            c_loc(values(1)))
    end function nc_put_real_vector

    integer function nc_put_real_tensor(ncid, varid, values) result(status)
        integer, intent(in) :: ncid, varid
        real(c_double), intent(in), target, contiguous :: values(:, :, :)

        status = c_nc_put_var_double(int(ncid, c_int), int(varid, c_int), &
            c_loc(values(1, 1, 1)))
    end function nc_put_real_tensor

    subroutine build_c_string(value, c_value)
        character(len=*), intent(in) :: value
        character(c_char), allocatable, intent(out) :: c_value(:)
        integer :: index, length

        length = len_trim(value)
        allocate (c_value(length + 1))
        do index = 1, length
            c_value(index) = value(index:index)
        end do
        c_value(length + 1) = c_null_char
    end subroutine build_c_string

    subroutine copy_c_string(c_value, value)
        character(c_char), intent(in) :: c_value(:)
        character(len=*), intent(out) :: value
        integer :: index

        value = ""
        do index = 1, min(size(c_value), len(value))
            if (c_value(index) == c_null_char) return
            value(index:index) = c_value(index)
        end do
    end subroutine copy_c_string

end module netcdf_c_api
