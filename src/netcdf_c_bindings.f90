module netcdf_c_bindings
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr, c_size_t
    implicit none
    private

    integer, parameter, public :: nc_noerr = 0
    integer, parameter, public :: nc_enotatt = -43
    integer, parameter, public :: nc_enotvar = -49
    integer, parameter, public :: nc_global = -1
    integer, parameter, public :: nc_double = 6
    integer, parameter, public :: nc_int64 = 10

    public :: c_nc_close
    public :: c_nc_create
    public :: c_nc_def_dim
    public :: c_nc_def_var
    public :: c_nc_enddef
    public :: c_nc_get_att_text
    public :: c_nc_get_var_double
    public :: c_nc_get_var_int
    public :: c_nc_inq_attlen
    public :: c_nc_inq_dimid
    public :: c_nc_inq_dimlen
    public :: c_nc_inq_dimname
    public :: c_nc_inq_vardimid
    public :: c_nc_inq_varid
    public :: c_nc_inq_varndims
    public :: c_nc_open
    public :: c_nc_put_att_text
    public :: c_nc_put_var_double
    public :: c_nc_put_var_int

    interface
        function c_nc_open(path, mode, ncid) bind(c, name="nc_open") result(status)
            import :: c_char, c_int
            character(c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int), intent(out) :: ncid
            integer(c_int) :: status
        end function c_nc_open

        function c_nc_create(path, mode, ncid) bind(c, name="nc_create") &
                result(status)
            import :: c_char, c_int
            character(c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int), intent(out) :: ncid
            integer(c_int) :: status
        end function c_nc_create

        function c_nc_close(ncid) bind(c, name="nc_close") result(status)
            import :: c_int
            integer(c_int), value :: ncid
            integer(c_int) :: status
        end function c_nc_close

        function c_nc_inq_dimid(ncid, name, dimid) bind(c, &
                name="nc_inq_dimid") result(status)
            import :: c_char, c_int
            integer(c_int), value :: ncid
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: dimid
            integer(c_int) :: status
        end function c_nc_inq_dimid

        function c_nc_inq_dimlen(ncid, dimid, length) bind(c, &
                name="nc_inq_dimlen") result(status)
            import :: c_int, c_size_t
            integer(c_int), value :: ncid, dimid
            integer(c_size_t), intent(out) :: length
            integer(c_int) :: status
        end function c_nc_inq_dimlen

        function c_nc_inq_dimname(ncid, dimid, name) bind(c, &
                name="nc_inq_dimname") result(status)
            import :: c_char, c_int
            integer(c_int), value :: ncid, dimid
            character(c_char), intent(out) :: name(*)
            integer(c_int) :: status
        end function c_nc_inq_dimname

        function c_nc_inq_varid(ncid, name, varid) bind(c, &
                name="nc_inq_varid") result(status)
            import :: c_char, c_int
            integer(c_int), value :: ncid
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: varid
            integer(c_int) :: status
        end function c_nc_inq_varid

        function c_nc_inq_varndims(ncid, varid, rank) bind(c, &
                name="nc_inq_varndims") result(status)
            import :: c_int
            integer(c_int), value :: ncid, varid
            integer(c_int), intent(out) :: rank
            integer(c_int) :: status
        end function c_nc_inq_varndims

        function c_nc_inq_vardimid(ncid, varid, dimids) bind(c, &
                name="nc_inq_vardimid") result(status)
            import :: c_int, c_ptr
            integer(c_int), value :: ncid, varid
            type(c_ptr), value :: dimids
            integer(c_int) :: status
        end function c_nc_inq_vardimid

        function c_nc_get_var_int(ncid, varid, values) bind(c, &
                name="nc_get_var_int") result(status)
            import :: c_int, c_ptr
            integer(c_int), value :: ncid, varid
            type(c_ptr), value :: values
            integer(c_int) :: status
        end function c_nc_get_var_int

        function c_nc_get_var_double(ncid, varid, values) bind(c, &
                name="nc_get_var_double") result(status)
            import :: c_int, c_ptr
            integer(c_int), value :: ncid, varid
            type(c_ptr), value :: values
            integer(c_int) :: status
        end function c_nc_get_var_double

        function c_nc_inq_attlen(ncid, varid, name, length) bind(c, &
                name="nc_inq_attlen") result(status)
            import :: c_char, c_int, c_size_t
            integer(c_int), value :: ncid, varid
            character(c_char), intent(in) :: name(*)
            integer(c_size_t), intent(out) :: length
            integer(c_int) :: status
        end function c_nc_inq_attlen

        function c_nc_get_att_text(ncid, varid, name, value) bind(c, &
                name="nc_get_att_text") result(status)
            import :: c_char, c_int, c_ptr
            integer(c_int), value :: ncid, varid
            character(c_char), intent(in) :: name(*)
            type(c_ptr), value :: value
            integer(c_int) :: status
        end function c_nc_get_att_text

        function c_nc_def_dim(ncid, name, length, dimid) bind(c, &
                name="nc_def_dim") result(status)
            import :: c_char, c_int, c_size_t
            integer(c_int), value :: ncid
            character(c_char), intent(in) :: name(*)
            integer(c_size_t), value :: length
            integer(c_int), intent(out) :: dimid
            integer(c_int) :: status
        end function c_nc_def_dim

        function c_nc_def_var(ncid, name, xtype, rank, dimids, varid) bind(c, &
                name="nc_def_var") result(status)
            import :: c_char, c_int, c_ptr
            integer(c_int), value :: ncid, xtype, rank
            character(c_char), intent(in) :: name(*)
            type(c_ptr), value :: dimids
            integer(c_int), intent(out) :: varid
            integer(c_int) :: status
        end function c_nc_def_var

        function c_nc_put_att_text(ncid, varid, name, length, value) bind(c, &
                name="nc_put_att_text") result(status)
            import :: c_char, c_int, c_size_t, c_ptr
            integer(c_int), value :: ncid, varid
            character(c_char), intent(in) :: name(*)
            integer(c_size_t), value :: length
            type(c_ptr), value :: value
            integer(c_int) :: status
        end function c_nc_put_att_text

        function c_nc_enddef(ncid) bind(c, name="nc_enddef") result(status)
            import :: c_int
            integer(c_int), value :: ncid
            integer(c_int) :: status
        end function c_nc_enddef

        function c_nc_put_var_int(ncid, varid, values) bind(c, &
                name="nc_put_var_int") result(status)
            import :: c_int, c_ptr
            integer(c_int), value :: ncid, varid
            type(c_ptr), value :: values
            integer(c_int) :: status
        end function c_nc_put_var_int

        function c_nc_put_var_double(ncid, varid, values) bind(c, &
                name="nc_put_var_double") result(status)
            import :: c_int, c_ptr
            integer(c_int), value :: ncid, varid
            type(c_ptr), value :: values
            integer(c_int) :: status
        end function c_nc_put_var_double
    end interface

end module netcdf_c_bindings
