program gliss_geometry_admission
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use export_surface_geometry, only: build_angular_grids
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use primitive_equilibrium_spline, only: evaluate_primitive_equilibrium, &
        fit_primitive_equilibrium, primitive_equilibrium_ok, &
        primitive_equilibrium_spline_t
    use primitive_geometry_grid, only: primitive_geometry_grid_t
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(primitive_equilibrium_spline_t) :: spline
    type(primitive_geometry_grid_t) :: geometry
    character(len=1024) :: filename
    real(dp), allocatable :: theta(:), zeta(:)
    real(dp) :: s, pressure, slope
    integer :: info, i, count_s

    if (command_argument_count() /= 1) &
        error stop 'usage: gliss_geometry_admission EXPORT'
    call get_command_argument(1, filename)
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) error stop 'export could not be read'
    call fit_primitive_equilibrium(equilibrium, spline, info)
    if (info /= primitive_equilibrium_ok) error stop 'primitive fit failed'
    call build_angular_grids(256, 8, theta, zeta)
    count_s = 4 * size(equilibrium%s)
    write (*, '(a)') 's,status,min_signed_jacobian,max_signed_jacobian,' // &
        'max_abs_jacobian_s,min_mod_b,max_mod_b'
    do i = 1, count_s
        s = (real(i, dp) - 0.5_dp) / real(count_s, dp)
        call evaluate_primitive_equilibrium(spline, s, theta, zeta, &
            geometry, pressure, slope, info)
        if (info /= primitive_equilibrium_ok) then
            write (*, '(es24.16,a,i0)') s, ',', info
            cycle
        end if
        write (*, '(es24.16,a,i0,5(a,es24.16))') s, ',', info, ',', &
            minval(geometry%signed_jacobian), ',', &
            maxval(geometry%signed_jacobian), ',', &
            maxval(abs(geometry%jacobian_s)), ',', &
            minval(geometry%mod_b), ',', maxval(geometry%mod_b)
    end do
end program gliss_geometry_admission
