program gliss_mercier
    use, intrinsic :: iso_fortran_env, only: error_unit
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: compute_mercier, mercier_ok, &
        mercier_result_t
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(mercier_result_t) :: result
    character(len=1024) :: filename
    integer :: info, i

    if (command_argument_count() /= 1) then
        write (error_unit, "(a)") "usage: gliss_mercier EXPORT_FILE"
        error stop 2
    end if
    call get_command_argument(1, filename)
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) then
        write (error_unit, "(a, i0)") "reader error ", info
        error stop 1
    end if
    call compute_mercier(equilibrium, 64, 64, result, info)
    if (info /= mercier_ok) then
        write (error_unit, "(a, i0)") "mercier error ", info
        error stop 1
    end if

    write (*, "(a)") "s,D_shear,D_current,D_well,D_geodesic,D_Mercier," // &
        "iota_deviation,boozer_deviation,force_balance_residual," // &
        "jacobian_identity_deviation"
    do i = 1, size(result%s)
        write (*, "(es24.16, 9(',', es24.16))") result%s(i), &
            result%d_shear(i), result%d_current(i), result%d_well(i), &
            result%d_geodesic(i), result%d_mercier(i), &
            result%iota_deviation(i), result%boozer_deviation(i), &
            result%force_balance_residual(i), &
            result%jacobian_identity_deviation(i)
    end do
end program gliss_mercier
