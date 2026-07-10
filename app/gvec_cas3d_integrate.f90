program gvec_cas3d_integrate
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_integrals, only: integrate_half_mesh_volume, integration_ok
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    character(len=1024) :: filename
    real(dp) :: full_device_volume, signed_period_volume
    integer :: argument_status, info

    call get_command_argument(1, filename, status=argument_status)
    call require(argument_status == 0 .and. len_trim(filename) > 0, &
        "usage: gvec_cas3d_integrate FOURIER")
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    call require(info == reader_ok, "Fourier export could not be read")
    call integrate_half_mesh_volume(equilibrium, signed_period_volume, &
        full_device_volume, info)
    call require(info == integration_ok, "volume integration failed")
    write (*, "(a)") "quantity,value"
    write (*, "(a,',',es24.16e3)") &
        "signed_one_period_volume", signed_period_volume
    write (*, "(a,',',es24.16e3)") &
        "full_device_volume", full_device_volume

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program gvec_cas3d_integrate
