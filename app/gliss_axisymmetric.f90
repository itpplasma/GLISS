program gliss_axisymmetric
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: error_unit
    use axisymmetric_spectrum, only: axisymmetric_spectrum_invalid_input, &
        axisymmetric_spectrum_ok, axisymmetric_spectrum_result_t, &
        compute_axisymmetric_spectrum
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(axisymmetric_spectrum_result_t) :: result
    character(len=1024) :: filename, message, token
    integer :: arguments, degree, info, poloidal_max, toroidal_mode
    logical :: count_only

    interface
        subroutine terminate_process(status) bind(C, name="exit")
            import c_int
            integer(c_int), value :: status
        end subroutine terminate_process
    end interface

    arguments = command_argument_count()
    if (arguments < 3 .or. arguments > 6) &
        call fail_usage("expected three arguments plus an optional rule")
    call read_argument(1, "EXPORT_FILE", filename)
    call read_integer_argument(2, "N", toroidal_mode)
    call read_integer_argument(3, "MMAX", poloidal_max)
    if (toroidal_mode <= 0) call fail_usage("N must be positive")
    if (poloidal_max < 1) call fail_usage("MMAX must be positive")
    degree = 2
    count_only = .false.
    if (arguments >= 5) then
        call read_argument(4, "option", token)
        if (trim(token) /= "--degree") &
            call fail_usage("the optional argument must be --degree")
        call read_integer_argument(5, "DEGREE", degree)
        if (degree < 1 .or. degree > 4) &
            call fail_usage("DEGREE must be between 1 and 4")
        if (arguments == 6) then
            call read_argument(6, "option", token)
            if (trim(token) /= "--count-only") &
                call fail_usage("the final option must be --count-only")
            count_only = .true.
        end if
    else if (arguments == 4) then
        call read_argument(4, "option", token)
        if (trim(token) /= "--count-only") &
            call fail_usage("the optional argument must be --count-only")
        count_only = .true.
    end if

    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) call fail("equilibrium export could not be read")
    call compute_axisymmetric_spectrum(equilibrium, toroidal_mode, &
        poloidal_max, degree, .not. count_only, result, info, message)
    if (info == axisymmetric_spectrum_invalid_input .and. &
        index(message, "aliases") > 0) &
        call fail_usage("MMAX aliases the fixed angular quadrature")
    if (info /= axisymmetric_spectrum_ok) call fail(trim(message))

    write (*, "(a)") "chart_metric,field_periods,toroidal_mode," // &
        "poloidal_max,modes,radial_surfaces,parity_class," // &
        "lowest_eigenvalue,inertia_certificate,eigenpair_residual," // &
        "negative_count,force_balance_residual,degree"
    write (*, "(l1, 6(a, i0), 3(a, es24.16), a, i0, a, es24.16, a, i0)") &
        equilibrium%has_chart_metric, ",", equilibrium%field_periods, ",", &
        toroidal_mode, ",", poloidal_max, ",", result%mode_count, ",", &
        result%radial_surfaces, ",", result%parity_class, ",", &
        result%lowest_eigenvalue, ",", result%certificate, ",", &
        result%eigenpair_residual, ",", result%negative_count, ",", &
        result%force_balance_residual, ",", result%degree

contains

    subroutine fail_usage(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_axisymmetric: " // trim(message)
        write (error_unit, "(a)") &
            "usage: gliss_axisymmetric EXPORT_FILE N MMAX " // &
            "[--degree DEGREE] [--count-only]"
        call terminate_process(2_c_int)
    end subroutine fail_usage

    subroutine fail(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_axisymmetric: " // trim(message)
        call terminate_process(1_c_int)
    end subroutine fail

    subroutine read_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        integer :: status

        call get_command_argument(position, value, status=status)
        if (status /= 0) call fail_usage(trim(name) // " is too long")
        if (len_trim(value) == 0) call fail_usage(trim(name) // " is empty")
    end subroutine read_argument

    subroutine read_integer_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        integer, intent(out) :: value
        integer :: status

        call read_argument(position, name, token)
        read (token, *, iostat=status) value
        if (status /= 0) call fail_usage(trim(name) // " must be an integer")
    end subroutine read_integer_argument

end program gliss_axisymmetric
