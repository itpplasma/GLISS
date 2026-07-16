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
    use terpsichore_pseudoplasma_spectrum, only: &
        solve_terpsichore_pseudoplasma_files, &
        terpsichore_pseudoplasma_result_t, &
        terpsichore_pseudoplasma_spectrum_compute_error, &
        terpsichore_pseudoplasma_spectrum_ok, &
        terpsichore_pseudoplasma_spectrum_read_error
    implicit none
    private

    type, bind(c) :: terpsichore_fixed_boundary_result_legacy_c
        integer(c_size_t) :: struct_size
        integer(c_size_t) :: unknowns
        integer(c_size_t) :: negative_count
        real(c_double) :: eigenvalue
        real(c_double) :: certificate
        real(c_double) :: residual
        real(c_double) :: resolution
    end type terpsichore_fixed_boundary_result_legacy_c

    type, bind(c) :: terpsichore_fixed_boundary_result_c
        integer(c_size_t) :: struct_size
        integer(c_size_t) :: unknowns
        integer(c_size_t) :: negative_count
        real(c_double) :: eigenvalue
        real(c_double) :: certificate
        real(c_double) :: residual
        real(c_double) :: resolution
        real(c_double) :: reference_eigenvalue
        real(c_double) :: reference_potential
        real(c_double) :: computed_potential
        real(c_double) :: reference_kinetic
        real(c_double) :: computed_kinetic
        real(c_double) :: reference_residual
        real(c_double) :: mode_overlap
    end type terpsichore_fixed_boundary_result_c

    type, bind(c) :: result_size_prefix_c
        integer(c_size_t) :: struct_size
    end type result_size_prefix_c

    type, bind(c) :: terpsichore_pseudoplasma_result_c
        integer(c_size_t) :: struct_size
        integer(c_size_t) :: unknowns
        integer(c_size_t) :: negative_count
        real(c_double) :: eigenvalue
        real(c_double) :: certificate
        real(c_double) :: residual
        real(c_double) :: resolution
        real(c_double) :: growth_rate
        real(c_double) :: reference_eigenvalue
        real(c_double) :: reference_potential
        real(c_double) :: computed_potential
        real(c_double) :: reference_kinetic
        real(c_double) :: computed_kinetic
        real(c_double) :: reference_residual
        real(c_double) :: mode_overlap
    end type terpsichore_pseudoplasma_result_c

    public :: gliss_terpsichore_fixed_boundary_c
    public :: gliss_terpsichore_pseudoplasma_c

contains

    function gliss_terpsichore_fixed_boundary_c(path_pointer, path_length, &
            result_pointer, error_pointer, error_capacity) bind(c, &
            name="gliss_terpsichore_fixed_boundary") result(status)
        type(c_ptr), value, intent(in) :: path_pointer, result_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: path_length, error_capacity
        integer(c_int) :: status
        type(terpsichore_fixed_boundary_result_legacy_c), pointer :: legacy_result
        type(terpsichore_fixed_boundary_result_c), pointer :: result
        type(result_size_prefix_c), pointer :: prefix
        type(terpsichore_fixed_boundary_result_legacy_c) :: legacy_probe
        type(terpsichore_fixed_boundary_result_c) :: result_probe
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
        call c_f_pointer(result_pointer, prefix)
        if (prefix%struct_size /= c_sizeof(legacy_probe) .and. &
            prefix%struct_size /= c_sizeof(result_probe)) then
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
            if (prefix%struct_size == c_sizeof(legacy_probe)) then
                call c_f_pointer(result_pointer, legacy_result)
                call fill_legacy_result(native, legacy_result)
            else
                call c_f_pointer(result_pointer, result)
                call fill_result(native, result)
            end if
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

    function gliss_terpsichore_pseudoplasma_c(matrix_path_pointer, &
            matrix_path_length, vacuum_intervals, vacuum_path_pointer, &
            vacuum_path_length, result_pointer, error_pointer, error_capacity) &
            bind(c, name="gliss_terpsichore_pseudoplasma") result(status)
        type(c_ptr), value, intent(in) :: matrix_path_pointer
        type(c_ptr), value, intent(in) :: vacuum_path_pointer, result_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: matrix_path_length
        integer(c_size_t), value, intent(in) :: vacuum_path_length
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int), value, intent(in) :: vacuum_intervals
        integer(c_int) :: status
        type(terpsichore_pseudoplasma_result_c), pointer :: result
        type(terpsichore_pseudoplasma_result_t) :: native
        character(len=:), allocatable :: matrix_filename, vacuum_filename
        character(len=128) :: message
        integer :: info

        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(result_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "TERPSICHORE pseudoplasma result pointer is null")
            return
        end if
        call c_f_pointer(result_pointer, result)
        if (result%struct_size /= c_sizeof(result)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "TERPSICHORE pseudoplasma result struct_size is incompatible")
            return
        end if
        if (vacuum_intervals <= 0_c_int) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "TERPSICHORE IVAC must be positive")
            return
        end if
        status = decode_path(matrix_path_pointer, matrix_path_length, &
            matrix_filename, message)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, trim(message))
            return
        end if
        status = decode_path(vacuum_path_pointer, vacuum_path_length, &
            vacuum_filename, message)
        if (status /= status_ok) then
            call write_error(error_pointer, error_capacity, trim(message))
            return
        end if
        call solve_terpsichore_pseudoplasma_files(matrix_filename, &
            int(vacuum_intervals), vacuum_filename, native, info, message)
        select case (info)
        case (terpsichore_pseudoplasma_spectrum_ok)
            call fill_pseudoplasma_result(native, result)
            status = status_ok
        case (terpsichore_pseudoplasma_spectrum_read_error)
            status = status_read_error
            call write_error(error_pointer, error_capacity, trim(message))
        case (terpsichore_pseudoplasma_spectrum_compute_error)
            status = status_compute_error
            call write_error(error_pointer, error_capacity, trim(message))
        case default
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "TERPSICHORE pseudoplasma solve returned an unknown status")
        end select
    end function gliss_terpsichore_pseudoplasma_c

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
        result%reference_eigenvalue = native%reference_eigenvalue
        result%reference_potential = native%reference_potential
        result%computed_potential = native%computed_potential
        result%reference_kinetic = native%reference_kinetic
        result%computed_kinetic = native%computed_kinetic
        result%reference_residual = native%reference_residual
        result%mode_overlap = native%mode_overlap
    end subroutine fill_result

    subroutine fill_legacy_result(native, result)
        type(terpsichore_fixed_boundary_result_t), intent(in) :: native
        type(terpsichore_fixed_boundary_result_legacy_c), intent(out) :: result

        result%struct_size = c_sizeof(result)
        result%unknowns = int(native%unknowns, c_size_t)
        result%negative_count = int(native%negative_count, c_size_t)
        result%eigenvalue = native%eigenvalue
        result%certificate = native%certificate
        result%residual = native%residual
        result%resolution = native%resolution
    end subroutine fill_legacy_result

    subroutine fill_pseudoplasma_result(native, result)
        type(terpsichore_pseudoplasma_result_t), intent(in) :: native
        type(terpsichore_pseudoplasma_result_c), intent(out) :: result

        result%struct_size = c_sizeof(result)
        result%unknowns = int(native%unknowns, c_size_t)
        result%negative_count = int(native%negative_count, c_size_t)
        result%eigenvalue = native%eigenvalue
        result%certificate = native%certificate
        result%residual = native%residual
        result%resolution = native%resolution
        result%growth_rate = native%growth_rate
        result%reference_eigenvalue = native%reference_eigenvalue
        result%reference_potential = native%reference_potential
        result%computed_potential = native%computed_potential
        result%reference_kinetic = native%reference_kinetic
        result%computed_kinetic = native%computed_kinetic
        result%reference_residual = native%reference_residual
        result%mode_overlap = native%mode_overlap
    end subroutine fill_pseudoplasma_result

end module gliss_terpsichore_capi
