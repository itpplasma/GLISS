module gvec_cas3d_netcdf
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf_c_api, only: nc_enotvar, nc_get_integer, nc_get_real, &
        nc_inquire_dimension_id, nc_inquire_dimension_length, &
        nc_inquire_dimension_name, nc_inquire_variable_dimensions, &
        nc_inquire_variable_id, nc_inquire_variable_rank, nc_noerr
    implicit none
    private

    integer, parameter, public :: reader_ok = 0
    integer, parameter, public :: reader_open_error = 1
    integer, parameter, public :: reader_schema_error = 2
    integer, parameter, public :: reader_data_error = 3
    integer, parameter, public :: reader_coordinate_error = 4

    public :: read_dimension
    public :: read_harmonic_component
    public :: read_integer_scalar
    public :: read_integer_vector
    public :: read_optional_real_vector
    public :: read_real_scalar
    public :: read_real_tensor
    public :: read_real_vector

contains

    subroutine read_harmonic_component(ncid, name, ns, nm, nn, values, &
            present, info)
        integer, intent(in) :: ncid, ns, nm, nn
        character(len=*), intent(in) :: name
        real(dp), allocatable, intent(out) :: values(:, :, :)
        logical, intent(out) :: present
        integer, intent(out) :: info
        integer :: nc_status, varid

        present = .false.
        nc_status = nc_inquire_variable_id(ncid, name, varid)
        if (nc_status == nc_enotvar) then
            allocate (values(ns, nm, nn))
            values = 0.0_dp
            info = reader_ok
            return
        end if
        if (nc_status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        call read_real_tensor(ncid, name, ["s", "m", "n"], [ns, nm, nn], &
            values, info)
        if (info /= reader_ok) return
        present = .true.
    end subroutine read_harmonic_component

    subroutine require_variable_dimensions(ncid, varid, names, lengths, info)
        integer, intent(in) :: ncid, varid
        character(len=*), intent(in) :: names(:)
        integer, intent(in) :: lengths(:)
        integer, intent(out) :: info
        integer, allocatable :: dimids(:)
        integer :: ndims, index, status

        status = nc_inquire_variable_rank(ncid, varid, ndims)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        if (ndims /= size(names)) then
            info = reader_schema_error
            return
        end if
        allocate (dimids(ndims))
        if (ndims > 0) then
            status = nc_inquire_variable_dimensions(ncid, varid, dimids)
            if (status /= nc_noerr) then
                info = reader_data_error
                return
            end if
        end if
        do index = 1, ndims
            call require_dimension_id(ncid, dimids(index), names(index), &
                lengths(index), info)
            if (info /= reader_ok) return
        end do
        info = reader_ok
    end subroutine require_variable_dimensions

    subroutine read_dimension(ncid, name, length, info)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: length, info
        integer :: dimid, status

        status = nc_inquire_dimension_id(ncid, name, dimid)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        status = nc_inquire_dimension_length(ncid, dimid, length)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        if (length < 1) then
            info = reader_schema_error
            return
        end if
        info = reader_ok
    end subroutine read_dimension

    subroutine require_dimension_id(ncid, dimid, expected_name, &
            expected_length, info)
        integer, intent(in) :: ncid, dimid, expected_length
        character(len=*), intent(in) :: expected_name
        integer, intent(out) :: info
        character(len=256) :: actual_name
        integer :: actual_length, status

        status = nc_inquire_dimension_name(ncid, dimid, actual_name)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        status = nc_inquire_dimension_length(ncid, dimid, actual_length)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        if (trim(actual_name) /= trim(expected_name)) then
            info = reader_schema_error
            return
        end if
        if (actual_length /= expected_length) then
            info = reader_schema_error
            return
        end if
        info = reader_ok
    end subroutine require_dimension_id

    subroutine read_integer_scalar(ncid, name, value, info)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: value, info
        integer :: varid, status

        status = nc_inquire_variable_id(ncid, name, varid)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        call require_variable_dimensions(ncid, varid, [character(len=1) ::], &
            [integer ::], info)
        if (info /= reader_ok) return
        status = nc_get_integer(ncid, varid, value)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        info = reader_ok
    end subroutine read_integer_scalar

    subroutine read_real_scalar(ncid, name, value, info)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: value
        integer, intent(out) :: info
        integer :: varid, status

        status = nc_inquire_variable_id(ncid, name, varid)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        call require_variable_dimensions(ncid, varid, [character(len=1) ::], &
            [integer ::], info)
        if (info /= reader_ok) return
        status = nc_get_real(ncid, varid, value)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        if (.not. ieee_is_finite(value)) then
            info = reader_data_error
            return
        end if
        info = reader_ok
    end subroutine read_real_scalar

    subroutine read_integer_vector(ncid, name, dimension_name, count, &
            values, info)
        integer, intent(in) :: ncid, count
        character(len=*), intent(in) :: name, dimension_name
        integer, allocatable, intent(out) :: values(:)
        integer, intent(out) :: info
        integer :: varid, status

        status = nc_inquire_variable_id(ncid, name, varid)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        call require_variable_dimensions(ncid, varid, [dimension_name], &
            [count], info)
        if (info /= reader_ok) return
        allocate (values(count))
        status = nc_get_integer(ncid, varid, values)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        info = reader_ok
    end subroutine read_integer_vector

    subroutine read_real_vector(ncid, name, dimension_name, count, values, info)
        integer, intent(in) :: ncid, count
        character(len=*), intent(in) :: name, dimension_name
        real(dp), allocatable, intent(out) :: values(:)
        integer, intent(out) :: info
        integer :: varid, status

        status = nc_inquire_variable_id(ncid, name, varid)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        call require_variable_dimensions(ncid, varid, [dimension_name], &
            [count], info)
        if (info /= reader_ok) return
        allocate (values(count))
        status = nc_get_real(ncid, varid, values)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        if (.not. all(ieee_is_finite(values))) then
            info = reader_data_error
            return
        end if
        info = reader_ok
    end subroutine read_real_vector

    subroutine read_real_tensor(ncid, name, dimension_names, counts, values, &
            info)
        integer, intent(in) :: ncid, counts(3)
        character(len=*), intent(in) :: name, dimension_names(3)
        real(dp), allocatable, intent(out) :: values(:, :, :)
        integer, intent(out) :: info
        real(dp), allocatable :: file_values(:, :, :)
        integer :: first, second, third, status, varid

        status = nc_inquire_variable_id(ncid, name, varid)
        if (status /= nc_noerr) then
            info = reader_schema_error
            return
        end if
        call require_variable_dimensions(ncid, varid, dimension_names, &
            counts, info)
        if (info /= reader_ok) return
        allocate (file_values(counts(3), counts(2), counts(1)))
        status = nc_get_real(ncid, varid, file_values)
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        if (.not. all(ieee_is_finite(file_values))) then
            info = reader_data_error
            return
        end if
        allocate (values(counts(1), counts(2), counts(3)))
        do third = 1, counts(3)
            do second = 1, counts(2)
                do first = 1, counts(1)
                    values(first, second, third) = &
                        file_values(third, second, first)
                end do
            end do
        end do
        info = reader_ok
    end subroutine read_real_tensor

    subroutine read_optional_real_vector(ncid, name, dimension_name, count, &
            values, present, info)
        integer, intent(in) :: ncid, count
        character(len=*), intent(in) :: name, dimension_name
        real(dp), allocatable, intent(out) :: values(:)
        logical, intent(out) :: present
        integer, intent(out) :: info
        integer :: varid, status

        status = nc_inquire_variable_id(ncid, name, varid)
        if (status == nc_enotvar) then
            allocate (values(0))
            present = .false.
            info = reader_ok
            return
        end if
        if (status /= nc_noerr) then
            info = reader_data_error
            return
        end if
        call read_real_vector(ncid, name, dimension_name, count, values, info)
        present = info == reader_ok
    end subroutine read_optional_real_vector

end module gvec_cas3d_netcdf
