program gliss_compatible_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_three_component_problem, only: &
        build_compatible_three_component_problem, &
        compatible_three_component_ok, compatible_three_component_problem_t
    use compatible_problem_assembly_support, only: &
        evaluate_generalized_eigenpair, quadratic_form
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated, &
        symmetric_eigensolver_ok
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(compatible_three_component_problem_t) :: problem
    character(len=1024) :: filename, token
    integer, allocatable :: mode_m(:), mode_n(:)
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :), stored_power(:)
    real(dp), allocatable :: stiffness_copy(:, :), mass_copy(:, :)
    real(dp) :: adiabatic_index, density, floor, kinetic, potential, residual
    real(dp) :: term_energy(5)
    integer :: allocation_status, arguments, comma, degree, info, mode
    integer :: negative_count
    integer :: n_theta, n_zeta, parity, trial

    interface
        subroutine terminate_process(status) bind(C, name="exit")
            import c_int
            integer(c_int), value :: status
        end subroutine terminate_process
    end interface

    arguments = command_argument_count()
    if (arguments < 9) call fail_usage("missing required arguments")
    call read_argument(1, "EXPORT_FILE", filename)
    call read_integer_argument(2, "DEGREE", degree)
    call read_integer_argument(3, "NTHETA", n_theta)
    call read_integer_argument(4, "NZETA", n_zeta)
    call read_real_argument(5, "GAMMA", adiabatic_index)
    call read_real_argument(6, "DENSITY", density)
    call read_real_argument(7, "FLOOR", floor)
    call read_integer_argument(8, "PARITY", parity)
    if (degree < 1 .or. degree > 4) &
        call fail_usage("DEGREE must be between 1 and 4")
    if (n_theta < 8 .or. n_zeta < 8) &
        call fail_usage("NTHETA and NZETA must be at least 8")
    if (adiabatic_index <= 0.0_dp) call fail_usage("GAMMA must be positive")
    if (density <= 0.0_dp) call fail_usage("DENSITY must be positive")
    if (floor <= 0.0_dp) call fail_usage("FLOOR must be positive")
    if (parity < 1 .or. parity > 2) call fail_usage("PARITY must be 1 or 2")
    allocate (mode_m(arguments - 8), mode_n(arguments - 8), &
        stored_power(arguments - 8), stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("mode allocation", -1)
    do trial = 9, arguments
        mode = trial - 8
        call read_argument(trial, "mode", token)
        comma = index(token, ",")
        if (comma <= 1 .or. comma == len_trim(token)) &
            call fail_usage("modes must be given as m,n")
        if (index(token(comma + 1:), ",") > 0) &
            call fail_usage("modes must contain exactly one comma")
        call parse_integer(token(:comma - 1), "poloidal mode", mode_m(mode))
        call parse_integer(token(comma + 1:), "toroidal mode", mode_n(mode))
        if (mode_m(mode) < 0) &
            call fail_usage("poloidal mode m must be nonnegative")
        if (mode_m(mode) == 0 .and. mode_n(mode) < 0) &
            call fail_usage("axis modes require nonnegative n")
        do info = 1, mode - 1
            if (mode_m(info) == mode_m(mode) &
                .and. mode_n(info) == mode_n(mode)) &
                call fail_usage("duplicate mode")
        end do
        stored_power(mode) = 0.0_dp
        if (mode_m(mode) > 0) stored_power(mode) = &
            1.0_dp - 0.5_dp * real(mode_m(mode), dp)
    end do
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) call fail_solver("reader", info)
    call build_compatible_three_component_problem(equilibrium, &
        adiabatic_index, density, mode_m, mode_n, stored_power, parity, &
        degree, n_theta, n_zeta, problem, info)
    if (info /= compatible_three_component_ok) &
        call fail_solver("compatible spectrum problem", info)
    allocate (stiffness_copy, source=problem%stiffness, stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("stiffness allocation", -1)
    allocate (mass_copy, source=problem%mass, stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("mass allocation", -1)
    call solve_symmetric_generalized_allocated(stiffness_copy, mass_copy, &
        eigenvalues, eigenvectors, info, .true.)
    if (info /= symmetric_eigensolver_ok) &
        call fail_solver("compatible spectrum solve", info)
    negative_count = count(eigenvalues < -floor)
    call evaluate_generalized_eigenpair(problem%stiffness, problem%mass, &
        eigenvectors(:, 1), eigenvalues(1), kinetic, potential, residual, info)
    if (info /= 0) call fail_solver("compatible spectrum diagnostics", info)
    do info = 1, 5
        call quadratic_form(problem%stiffness_terms(:, :, info), &
            eigenvectors(:, 1), term_energy(info), trial)
        if (trial /= 0) call fail_solver("compatible term energy", trial)
    end do
    write (*, "(a)") "degree,parity,unknowns,normal_unknowns," // &
        "eta_unknowns,mu_unknowns,negative_count,lowest_eigenvalue," // &
        "residual,kinetic,potential,bending,shear,compression,drive,fluid"
    write (*, "(i0, 6(a, i0), 9(a, es24.16))") degree, ",", parity, ",", &
        size(eigenvalues), ",", problem%normal_unknowns, ",", &
        problem%eta_unknowns, ",", problem%mu_unknowns, ",", negative_count, &
        ",", eigenvalues(1), ",", residual, ",", kinetic, ",", potential, &
        ",", term_energy(1), ",", term_energy(2), ",", term_energy(3), &
        ",", term_energy(4), ",", term_energy(5)

contains

    subroutine fail_usage(message)
        character(len=*), intent(in) :: message
        write (error_unit, "(a)") "gliss_compatible_spectrum: " // trim(message)
        write (error_unit, "(a)") "usage: gliss_compatible_spectrum " // &
            "EXPORT_FILE DEGREE NTHETA NZETA GAMMA DENSITY FLOOR PARITY " // &
            "m,n [m,n ...]"
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

end program gliss_compatible_spectrum
