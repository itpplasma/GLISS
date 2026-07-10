program gliss_family
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use eigenvalue_tracking, only: certified_lowest_eigenvalue
    use family_assembly, only: family_negative_count, &
        surface_geometry_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none

    integer, parameter :: n_theta = 64, n_zeta = 64
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(surface_geometry_t), allocatable :: geometry(:)
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    integer, allocatable :: mode_m(:), mode_n(:)
    character(len=1024) :: filename, token
    real(dp) :: step, lowest, width
    integer :: info, i, ns, count, arguments, comma, selector

    arguments = command_argument_count()
    if (arguments < 2) then
        write (error_unit, "(a)") &
            "usage: gliss_family EXPORT_FILE m,n [m,n ...]"
        error stop 1
    end if
    call get_command_argument(1, filename)
    allocate (mode_m(arguments - 1), mode_n(arguments - 1))
    do i = 2, arguments
        call get_command_argument(i, token)
        comma = index(token, ",")
        if (comma <= 1) then
            write (error_unit, "(a)") "modes must be given as m,n"
            error stop 1
        end if
        read (token(1:comma - 1), *) mode_m(i - 1)
        read (token(comma + 1:), *) mode_n(i - 1)
    end do

    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) then
        write (error_unit, "(a, i0)") "reader error ", info
        error stop 1
    end if
    call build_kernel_geometry(equilibrium, n_theta, n_zeta, fields, &
        drive, info)
    if (info /= mercier_ok) then
        write (error_unit, "(a, i0)") "geometry error ", info
        error stop 1
    end if
    ns = size(equilibrium%s)
    allocate (geometry(ns))
    do i = 1, ns
        geometry(i)%fields = fields(:, :, :, i)
        geometry(i)%drive = drive(:, :, i)
    end do
    step = 1.0_dp / real(ns, dp)

    write (*, "(a)") "chart_metric,modes,parity_class," // &
        "lowest_eigenvalue,negative_count"
    do selector = 0, 2
        call certified_lowest_eigenvalue(geometry, mode_m, mode_n, &
            step, lowest, width, info, selector)
        if (info /= 0) then
            write (error_unit, "(a, i0)") "eigensolver error ", info
            error stop 1
        end if
        write (error_unit, "(a, i0, a, es9.2)") &
            "certified class ", selector, " window ", width
        call family_negative_count(geometry, mode_m, mode_n, step, &
            0.0_dp, count, info, selector)
        if (info /= 0) then
            write (error_unit, "(a, i0)") "inertia error ", info
            error stop 1
        end if
        write (*, "(l1, a, i0, a, i0, a, es24.16, a, i0)") &
            equilibrium%has_chart_metric, ",", size(mode_m), ",", &
            selector, ",", lowest, ",", count
    end do
end program gliss_family
