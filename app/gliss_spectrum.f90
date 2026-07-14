program gliss_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        fixed_boundary_ok, fixed_boundary_problem_t, &
        fixed_boundary_spectrum_result_t, solve_fixed_boundary_class
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mode_topology, only: build_mode_family, mode_family_t
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(fixed_boundary_problem_t) :: problem
    type(mode_family_t) :: family
    character(len=1024) :: filename, token
    integer, allocatable :: mode_m(:), mode_n(:)
    real(dp) :: adiabatic_index, density_kg_m3, zero_floor
    integer :: info, i, j, arguments, comma, selector_position, mode_index
    integer :: family_index, poloidal_max, toroidal_max, radial_quadrature
    logical :: generated_family

    interface
        subroutine terminate_process(status) bind(C, name="exit")
            import c_int
            integer(c_int), value :: status
        end subroutine terminate_process
    end interface

    arguments = command_argument_count()
    if (arguments < 5) call fail_usage("missing required arguments")
    call read_argument(1, "EXPORT_FILE", filename)
    call read_real_argument(2, "GAMMA", adiabatic_index)
    call read_real_argument(3, "DENSITY", density_kg_m3)
    call read_real_argument(4, "FLOOR", zero_floor)
    if (adiabatic_index < 0.0_dp) call fail_usage("GAMMA must be nonnegative")
    if (density_kg_m3 <= 0.0_dp) call fail_usage("DENSITY must be positive")
    if (zero_floor <= 0.0_dp) call fail_usage("FLOOR must be positive")
    radial_quadrature = 1
    selector_position = 5
    call read_argument(selector_position, "mode or option", token)
    if (trim(token) == "--quadrature") then
        call parse_quadrature(arguments, selector_position, radial_quadrature)
    end if
    generated_family = trim(token) == "--family"
    if (generated_family) then
        call parse_family(arguments, selector_position, family_index, &
            poloidal_max, toroidal_max)
    else
        call parse_modes(arguments, selector_position, mode_m, mode_n)
    end if

    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) call fail_solver("reader", info)
    if (generated_family) then
        call build_mode_family(equilibrium%field_periods, family_index, &
            poloidal_max, toroidal_max, family, info)
        if (info /= 0) call fail_solver("mode-family configuration", info)
        mode_m = family%poloidal
        mode_n = family%toroidal
    end if
    call build_fixed_boundary_problem(equilibrium, adiabatic_index, &
        density_kg_m3, zero_floor, mode_m, mode_n, radial_quadrature, &
        problem, info)
    if (info /= fixed_boundary_ok) call fail_solver("spectrum problem", info)

    write (*, "(a)") "chart_metric,field_periods,modes,parity_class," // &
        "adiabatic_index,density_kg_m3,unknowns,negative_count," // &
        "floor_count,lowest_eigenvalue,certificate,eigenpair_residual," // &
        "eigenpair_resolution,inertia_interval,radial_quadrature_points"
    do i = 1, 2
        call report_class(i)
    end do

contains

    subroutine parse_quadrature(count, selector, quadrature)
        integer, intent(in) :: count
        integer, intent(inout) :: selector, quadrature

        if (count < 7) &
            call fail_usage("--quadrature requires RULE and a mode selector")
        call read_argument(6, "quadrature rule", token)
        select case (trim(token))
        case ("midpoint")
            quadrature = 1
        case default
            call fail_usage("quadrature RULE must be midpoint")
        end select
        selector = 7
        call read_argument(selector, "mode or --family", token)
    end subroutine parse_quadrature

    subroutine parse_family(count, selector, index_value, m_max, n_max)
        integer, intent(in) :: count, selector
        integer, intent(out) :: index_value, m_max, n_max

        if (count /= selector + 3) &
            call fail_usage("--family requires exactly INDEX MMAX NMAX")
        call read_integer_argument(selector + 1, "INDEX", index_value)
        call read_integer_argument(selector + 2, "MMAX", m_max)
        call read_integer_argument(selector + 3, "NMAX", n_max)
        if (index_value < 0) call fail_usage("INDEX must be nonnegative")
        if (m_max < 0) call fail_usage("MMAX must be nonnegative")
        if (n_max < 0) call fail_usage("NMAX must be nonnegative")
    end subroutine parse_family

    subroutine parse_modes(count, selector, poloidal, toroidal)
        integer, intent(in) :: count, selector
        integer, allocatable, intent(out) :: poloidal(:), toroidal(:)

        allocate (poloidal(count - selector + 1), &
            toroidal(count - selector + 1))
        do i = selector, count
            mode_index = i - selector + 1
            call read_argument(i, "mode", token)
            comma = index(token, ",")
            if (comma <= 1 .or. comma == len_trim(token)) &
                call fail_usage("modes must be given as m,n")
            if (index(token(comma + 1:), ",") > 0) &
                call fail_usage("modes must contain exactly one comma")
            call parse_integer(token(:comma - 1), "poloidal mode", &
                poloidal(mode_index))
            call parse_integer(token(comma + 1:), "toroidal mode", &
                toroidal(mode_index))
            if (poloidal(mode_index) < 0) &
                call fail_usage("poloidal mode m must be nonnegative")
            if (poloidal(mode_index) == 0 .and. toroidal(mode_index) < 0) &
                call fail_usage("axis modes require nonnegative n")
            do j = 1, mode_index - 1
                if (poloidal(j) == poloidal(mode_index) .and. &
                    toroidal(j) == toroidal(mode_index)) &
                    call fail_usage("duplicate mode")
            end do
        end do
    end subroutine parse_modes

    subroutine report_class(parity_class)
        integer, intent(in) :: parity_class
        type(fixed_boundary_spectrum_result_t) :: result

        call solve_fixed_boundary_class(problem, parity_class, result, info)
        if (info /= fixed_boundary_ok) call fail_solver("spectrum solve", info)
        write (*, "(l1, a, i0, a, i0, a, i0, a, es9.2, a, es9.2, a, i0, " // &
            "a, i0, a, i0, a, es24.16, 4(a, es9.2), a, i0)") &
            result%has_chart_metric, ",", result%field_periods, ",", &
            result%mode_count, ",", result%parity_class, ",", &
            result%adiabatic_index, ",", result%density_kg_m3, ",", &
            result%unknowns, ",", result%negative_count, ",", &
            result%floor_count, ",", result%lowest_eigenvalue, ",", &
            result%certificate, ",", result%eigenpair_residual, ",", &
            result%eigenpair_resolution, ",", result%inertia_interval, ",", &
            result%radial_quadrature
    end subroutine report_class

    subroutine fail_usage(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_spectrum: " // trim(message)
        write (error_unit, "(a)") &
            "usage: gliss_spectrum EXPORT_FILE GAMMA DENSITY FLOOR " // &
            "[--quadrature midpoint] " // &
            "m,n [m,n ...] | --family INDEX MMAX NMAX"
        call terminate_process(2_c_int)
    end subroutine fail_usage

    subroutine fail_solver(operation, status)
        character(len=*), intent(in) :: operation
        integer, intent(in) :: status

        write (error_unit, "(a, i0)") trim(operation) // " error ", status
        error stop 1
    end subroutine fail_solver

    subroutine read_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        integer :: status

        call get_command_argument(position, value, status=status)
        if (status /= 0) call fail_usage(trim(name) // " is too long")
        if (len_trim(value) == 0) call fail_usage(trim(name) // " is empty")
    end subroutine read_argument

    subroutine read_real_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: value

        call read_argument(position, name, token)
        call parse_real(token, name, value)
    end subroutine read_real_argument

    subroutine read_integer_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        integer, intent(out) :: value

        call read_argument(position, name, token)
        call parse_integer(token, name, value)
    end subroutine read_integer_argument

    subroutine parse_real(text, name, value)
        character(len=*), intent(in) :: text, name
        real(dp), intent(out) :: value
        integer :: status

        read (text, *, iostat=status) value
        if (status /= 0) call fail_usage(trim(name) // " must be a number")
        if (.not. ieee_is_finite(value)) &
            call fail_usage(trim(name) // " must be finite")
    end subroutine parse_real

    subroutine parse_integer(text, name, value)
        character(len=*), intent(in) :: text, name
        integer, intent(out) :: value
        integer :: status

        read (text, *, iostat=status) value
        if (status /= 0) call fail_usage(trim(name) // " must be an integer")
    end subroutine parse_integer

end program gliss_spectrum
