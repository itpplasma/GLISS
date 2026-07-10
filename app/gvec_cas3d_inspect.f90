program gvec_cas3d_inspect
    use, intrinsic :: iso_fortran_env, only: error_unit
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    character(len=1024) :: filename
    integer :: argument_status, info

    call get_command_argument(1, filename, status=argument_status)
    if (argument_status /= 0) then
        write (error_unit, "(a)") "usage: gvec_cas3d_inspect FILE"
        error stop 2
    end if
    if (len_trim(filename) == 0) then
        write (error_unit, "(a)") "usage: gvec_cas3d_inspect FILE"
        error stop 2
    end if
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) then
        write (error_unit, "(a,i0)") "GVEC CAS3D reader status: ", info
        error stop 1
    end if

    write (*, "(a,i0)") "field periods: ", equilibrium%field_periods
    write (*, "(a,i0)") "radial surfaces: ", size(equilibrium%s)
    write (*, "(a,i0)") "poloidal modes: ", size(equilibrium%poloidal_modes)
    write (*, "(a,i0)") "toroidal modes: ", size(equilibrium%toroidal_modes)
    write (*, "(a,l1)") "stellarator symmetric: ", &
        equilibrium%stellarator_symmetric
end program gvec_cas3d_inspect
