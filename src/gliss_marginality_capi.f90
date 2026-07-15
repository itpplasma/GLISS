module gliss_marginality_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
        c_f_pointer, c_int, c_ptr, c_size_t, c_sizeof
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gliss_c_abi_support, only: error_buffer_status, &
        status_allocation_error, status_compute_error, &
        status_invalid_argument, status_ok, write_error
    use gliss_c_contexts, only: equilibrium_context_t
    use two_component_spectrum, only: compute_phase_envelope_spectrum, &
        compute_two_component_spectrum, two_component_spectrum_compute_error, &
        two_component_spectrum_invalid, two_component_spectrum_ok, &
        two_component_spectrum_result_t
    implicit none
    private

    type, bind(c) :: marginality_result_c
        integer(c_size_t) :: struct_size
        integer(c_int) :: has_eigenpair
        integer(c_int) :: field_periods
        integer(c_size_t) :: mode_count
        integer(c_size_t) :: radial_surfaces
        integer(c_int) :: parity_class
        integer(c_int) :: radial_quadrature
        integer(c_int) :: angular_theta
        integer(c_int) :: angular_zeta
        integer(c_size_t) :: negative_count
        real(c_double) :: lowest_eigenvalue
        real(c_double) :: certificate
        real(c_double) :: eigenpair_residual
        real(c_double) :: force_balance_residual
    end type marginality_result_c

    public :: gliss_cas3d_marginality_c
    public :: gliss_cas3d_phase_envelope_c

contains

    function gliss_cas3d_marginality_c(equilibrium_handle, mode_count, &
            mode_m_pointer, mode_n_pointer, parity_class, radial_quadrature, &
            angular_theta, angular_zeta, solve_eigenpair, result_pointer, &
            error_pointer, error_capacity) bind(c, &
            name="gliss_cas3d_marginality") result(status)
        type(c_ptr), value, intent(in) :: equilibrium_handle
        integer(c_size_t), value, intent(in) :: mode_count
        type(c_ptr), value, intent(in) :: mode_m_pointer, mode_n_pointer
        integer(c_int), value, intent(in) :: parity_class, radial_quadrature
        integer(c_int), value, intent(in) :: angular_theta, angular_zeta
        integer(c_int), value, intent(in) :: solve_eigenpair
        type(c_ptr), value, intent(in) :: result_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(equilibrium_context_t), pointer :: equilibrium
        type(marginality_result_c), pointer :: result
        type(two_component_spectrum_result_t) :: native
        integer, allocatable :: mode_m(:), mode_n(:)
        real(dp), allocatable :: stored_power(:)
        character(len=128) :: message
        integer :: info

        status = prepare_call(equilibrium_handle, result_pointer, &
            error_pointer, error_capacity, equilibrium, result)
        if (status /= status_ok) return
        if (solve_eigenpair /= 0_c_int .and. solve_eigenpair /= 1_c_int) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "solve_eigenpair must be 0 or 1")
            return
        end if
        status = decode_mode_table(mode_count, mode_m_pointer, mode_n_pointer, &
            mode_m, mode_n, stored_power)
        if (status /= status_ok) then
            if (status == status_allocation_error) then
                call write_error(error_pointer, error_capacity, &
                    "failed to allocate mode storage")
            else
                call write_error(error_pointer, error_capacity, &
                    "modes must be nonempty valid int32 arrays")
            end if
            return
        end if
        call compute_two_component_spectrum(equilibrium%equilibrium, mode_m, &
            mode_n, stored_power, int(parity_class), int(radial_quadrature), &
            int(angular_theta), int(angular_zeta), &
            solve_eigenpair == 1_c_int, native, info, message)
        select case (info)
        case (two_component_spectrum_ok)
            call fill_result(native, int(angular_theta), int(angular_zeta), &
                result)
            status = status_ok
        case (two_component_spectrum_invalid)
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, trim(message))
        case (two_component_spectrum_compute_error)
            status = status_compute_error
            call write_error(error_pointer, error_capacity, trim(message))
        case default
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "marginality solve returned an unknown status")
        end select
    end function gliss_cas3d_marginality_c

    function gliss_cas3d_phase_envelope_c(equilibrium_handle, base_m, &
            base_n, envelope_count, envelope_m_pointer, envelope_n_pointer, &
            parity_class, radial_quadrature, angular_theta, angular_zeta, &
            solve_eigenpair, result_pointer, error_pointer, error_capacity) &
            bind(c, name="gliss_cas3d_phase_envelope") result(status)
        type(c_ptr), value, intent(in) :: equilibrium_handle
        integer(c_int), value, intent(in) :: base_m, base_n
        integer(c_size_t), value, intent(in) :: envelope_count
        type(c_ptr), value, intent(in) :: envelope_m_pointer
        type(c_ptr), value, intent(in) :: envelope_n_pointer
        integer(c_int), value, intent(in) :: parity_class, radial_quadrature
        integer(c_int), value, intent(in) :: angular_theta, angular_zeta
        integer(c_int), value, intent(in) :: solve_eigenpair
        type(c_ptr), value, intent(in) :: result_pointer, error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        integer(c_int) :: status
        type(equilibrium_context_t), pointer :: equilibrium
        type(marginality_result_c), pointer :: result
        type(two_component_spectrum_result_t) :: native
        integer, allocatable :: envelope_m(:), envelope_n(:)
        real(dp), allocatable :: unused_power(:)
        character(len=128) :: message
        integer :: info

        status = prepare_call(equilibrium_handle, result_pointer, &
            error_pointer, error_capacity, equilibrium, result)
        if (status /= status_ok) return
        if (solve_eigenpair /= 0_c_int .and. solve_eigenpair /= 1_c_int) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "solve_eigenpair must be 0 or 1")
            return
        end if
        status = decode_mode_table(envelope_count, envelope_m_pointer, &
            envelope_n_pointer, envelope_m, envelope_n, unused_power)
        if (status /= status_ok) then
            if (status == status_allocation_error) then
                call write_error(error_pointer, error_capacity, &
                    "failed to allocate phase-envelope storage")
            else
                call write_error(error_pointer, error_capacity, &
                    "phase-envelope modes must be nonempty int32 arrays")
            end if
            return
        end if
        call compute_phase_envelope_spectrum(equilibrium%equilibrium, &
            int(base_m), int(base_n), envelope_m, envelope_n, &
            int(parity_class), int(radial_quadrature), int(angular_theta), &
            int(angular_zeta), solve_eigenpair == 1_c_int, native, info, &
            message)
        select case (info)
        case (two_component_spectrum_ok)
            call fill_result(native, int(angular_theta), int(angular_zeta), &
                result)
            status = status_ok
        case (two_component_spectrum_invalid)
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, trim(message))
        case (two_component_spectrum_compute_error)
            status = status_compute_error
            call write_error(error_pointer, error_capacity, trim(message))
        case default
            status = status_compute_error
            call write_error(error_pointer, error_capacity, &
                "phase-envelope solve returned an unknown status")
        end select
    end function gliss_cas3d_phase_envelope_c

    function prepare_call(equilibrium_handle, result_pointer, error_pointer, &
            error_capacity, equilibrium, result) result(status)
        type(c_ptr), value, intent(in) :: equilibrium_handle, result_pointer
        type(c_ptr), value, intent(in) :: error_pointer
        integer(c_size_t), value, intent(in) :: error_capacity
        type(equilibrium_context_t), pointer, intent(out) :: equilibrium
        type(marginality_result_c), pointer, intent(out) :: result
        integer(c_int) :: status

        nullify (equilibrium, result)
        status = error_buffer_status(error_pointer, error_capacity)
        if (status /= status_ok) return
        call write_error(error_pointer, error_capacity, "")
        if (.not. c_associated(equilibrium_handle)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        if (.not. c_associated(result_pointer)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "marginality result pointer is null")
            return
        end if
        call c_f_pointer(equilibrium_handle, equilibrium)
        call c_f_pointer(result_pointer, result)
        if (.not. associated(equilibrium)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "equilibrium handle is null")
            return
        end if
        if (result%struct_size /= c_sizeof(result)) then
            status = status_invalid_argument
            call write_error(error_pointer, error_capacity, &
                "marginality result struct_size is incompatible")
            return
        end if
        status = status_ok
    end function prepare_call

    function decode_mode_table(mode_count, mode_m_pointer, mode_n_pointer, &
            mode_m, mode_n, stored_power) result(status)
        integer(c_size_t), value, intent(in) :: mode_count
        type(c_ptr), value, intent(in) :: mode_m_pointer, mode_n_pointer
        integer, allocatable, intent(out) :: mode_m(:), mode_n(:)
        real(dp), allocatable, intent(out) :: stored_power(:)
        integer(c_int), pointer :: c_mode_m(:), c_mode_n(:)
        integer(c_int) :: status
        integer :: allocation_status, count, index, pointer_shape(1)

        status = status_invalid_argument
        if (mode_count < 1_c_size_t) return
        if (mode_count > int(huge(count), c_size_t)) return
        if (.not. c_associated(mode_m_pointer)) return
        if (.not. c_associated(mode_n_pointer)) return
        count = int(mode_count)
        pointer_shape(1) = count
        call c_f_pointer(mode_m_pointer, c_mode_m, pointer_shape)
        call c_f_pointer(mode_n_pointer, c_mode_n, pointer_shape)
        allocate (mode_m(count), mode_n(count), stored_power(count), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            status = status_allocation_error
            return
        end if
        do index = 1, count
            mode_m(index) = int(c_mode_m(index))
            mode_n(index) = int(c_mode_n(index))
            stored_power(index) = 0.0_dp
            if (mode_m(index) > 0) stored_power(index) = &
                1.0_dp - 0.5_dp * real(mode_m(index), dp)
        end do
        status = status_ok
    end function decode_mode_table

    subroutine fill_result(native, angular_theta, angular_zeta, result)
        type(two_component_spectrum_result_t), intent(in) :: native
        integer, intent(in) :: angular_theta, angular_zeta
        type(marginality_result_c), intent(out) :: result

        result%struct_size = c_sizeof(result)
        result%has_eigenpair = merge(1_c_int, 0_c_int, native%has_eigenpair)
        result%field_periods = int(native%field_periods, c_int)
        result%mode_count = int(native%mode_count, c_size_t)
        result%radial_surfaces = int(native%radial_surfaces, c_size_t)
        result%parity_class = int(native%parity_class, c_int)
        result%radial_quadrature = int(native%radial_quadrature, c_int)
        result%angular_theta = int(angular_theta, c_int)
        result%angular_zeta = int(angular_zeta, c_int)
        result%negative_count = int(native%negative_count, c_size_t)
        result%lowest_eigenvalue = native%lowest_eigenvalue
        result%certificate = native%certificate
        result%eigenpair_residual = native%eigenpair_residual
        result%force_balance_residual = native%force_balance_residual
    end subroutine fill_result

end module gliss_marginality_capi
