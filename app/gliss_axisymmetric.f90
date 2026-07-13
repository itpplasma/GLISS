program gliss_axisymmetric
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use eigenvalue_tracking, only: certified_lowest_eigenvalue
    use family_assembly, only: family_assembly_options_t, &
        family_negative_count, surface_geometry_t
    use field_profile_identities, only: compute_field_profile_identities, &
        field_profile_identities_ok, field_profile_identity_result_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: equilibrium_is_axisymmetric, &
        gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none

    integer, parameter :: n_theta = 64, n_zeta = 8
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(family_assembly_options_t) :: options
    type(field_profile_identity_result_t) :: identities
    type(surface_geometry_t), allocatable :: geometry(:)
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    real(dp), allocatable :: normal_stored_power(:)
    integer, allocatable :: mode_m(:), mode_n(:)
    character(len=1024) :: filename, token
    real(dp) :: certificate, force_balance, lowest, residual, step
    integer :: arguments, count, info, m, mode, ns, poloidal_max
    integer :: toroidal_mode
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
    options%radial_space%quadrature_points = 1
    count_only = .false.
    if (arguments >= 5) then
        call read_argument(4, "option", token)
        if (trim(token) /= "--quadrature") &
            call fail_usage("the optional argument must be --quadrature")
        call read_argument(5, "quadrature rule", token)
        select case (trim(token))
        case ("midpoint")
            options%radial_space%quadrature_points = 1
        case ("gauss2")
            options%radial_space%quadrature_points = 2
        case default
            call fail_usage("quadrature RULE must be midpoint or gauss2")
        end select
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
    if (.not. equilibrium%has_chart_metric) &
        call fail("equilibrium export lacks g_st/g_sz chart metrics")
    if (equilibrium%field_periods /= 1) &
        call fail("axisymmetric comparison requires N_FP=1")
    if (.not. equilibrium_is_axisymmetric(equilibrium)) &
        call fail("equilibrium contains nonaxisymmetric harmonics")
    if (2 * poloidal_max + maxval(abs(equilibrium%poloidal_modes)) &
        >= n_theta) &
        call fail_usage("MMAX aliases the fixed angular quadrature")

    allocate (mode_m(2 * poloidal_max + 1))
    allocate (mode_n(2 * poloidal_max + 1))
    mode_m(1) = 0
    mode_n(1) = toroidal_mode
    mode = 1
    do m = 1, poloidal_max
        mode = mode + 1
        mode_m(mode) = m
        mode_n(mode) = -toroidal_mode
        mode = mode + 1
        mode_m(mode) = m
        mode_n(mode) = toroidal_mode
    end do
    allocate (normal_stored_power(size(mode_m)), source=0.0_dp)
    where (mode_m > 0)
        normal_stored_power = 1.0_dp - 0.5_dp * real(mode_m, dp)
    end where

    call build_kernel_geometry(equilibrium, n_theta, n_zeta, fields, &
        drive, info)
    if (info /= mercier_ok) call fail("kernel geometry could not be built")
    ns = size(equilibrium%s)
    allocate (geometry(ns))
    do m = 1, ns
        geometry(m)%fields = fields(:, :, :, m)
        geometry(m)%drive = drive(:, :, m)
    end do
    step = 1.0_dp / real(ns, dp)
    options%field_periods = 1
    options%parity_class = 1
    call family_negative_count(geometry, mode_m, mode_n, step, 0.0_dp, &
        count, info, options, normal_stored_power)
    if (info /= 0) call fail("zero-shift inertia failed")
    if (count_only) then
        lowest = ieee_value(lowest, ieee_quiet_nan)
        certificate = ieee_value(certificate, ieee_quiet_nan)
        residual = ieee_value(residual, ieee_quiet_nan)
    else
        call certified_lowest_eigenvalue(geometry, mode_m, mode_n, step, &
            lowest, certificate, info, options, normal_stored_power, residual)
        if (info /= 0) call fail("certified eigensolve failed")
    end if
    call compute_field_profile_identities(equilibrium, n_theta, n_zeta, &
        identities, info)
    if (info /= field_profile_identities_ok) &
        call fail("force-balance diagnostic failed")
    force_balance = maxval(identities%general_force_balance_deviation)

    write (*, "(a)") "chart_metric,field_periods,toroidal_mode," // &
        "poloidal_max,modes,radial_surfaces,parity_class," // &
        "lowest_eigenvalue,inertia_certificate,eigenpair_residual," // &
        "negative_count,force_balance_residual,radial_quadrature_points"
    write (*, "(l1, 6(a, i0), 3(a, es24.16), a, i0, a, es24.16, a, i0)") &
        equilibrium%has_chart_metric, ",", equilibrium%field_periods, ",", &
        toroidal_mode, ",", poloidal_max, ",", size(mode_m), ",", ns, &
        ",", options%parity_class, ",", lowest, ",", certificate, ",", &
        residual, ",", count, ",", force_balance, ",", &
        options%radial_space%quadrature_points

contains

    subroutine fail_usage(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_axisymmetric: " // trim(message)
        write (error_unit, "(a)") &
            "usage: gliss_axisymmetric EXPORT_FILE N MMAX " // &
            "[--quadrature midpoint|gauss2] [--count-only]"
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
