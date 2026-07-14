program gliss_compatible_marginality
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_problem_ok, &
        compatible_two_component_problem_t
    use compatible_problem_assembly_support, only: &
        evaluate_generalized_eigenpair, quadratic_form, sum_tensor
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use symmetric_pivot_inertia, only: pivot_negative_count
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated, &
        symmetric_eigensolver_ok
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(compatible_two_component_problem_t) :: problem
    character(len=1024) :: filename, token
    integer, allocatable :: mode_m(:), mode_n(:)
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
    real(dp), allocatable :: stored_power(:), stiffness(:, :), mass(:, :)
    real(dp), allocatable :: stiffness_operator(:, :), work(:)
    integer, allocatable :: pivots(:)
    real(dp) :: density, floor, kinetic, m0_stored_power, potential, residual
    real(dp) :: term_energy(4)
    integer :: allocation_status, arguments, comma, degree, first_mode, info
    integer :: inertia_negative_count, mode, negative_count, work_size
    integer :: n_theta, n_zeta, parity, trial
    logical :: has_m0_power, inertia_only, physical_mass

    interface
        subroutine terminate_process(status) bind(C, name="exit")
            import c_int
            integer(c_int), value :: status
        end subroutine terminate_process
        subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
            import dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine dsytrf
    end interface

    arguments = command_argument_count()
    if (arguments < 7) call fail_usage("missing required arguments")
    call read_argument(1, "EXPORT_FILE", filename)
    call read_integer_argument(2, "DEGREE", degree)
    call read_integer_argument(3, "NTHETA", n_theta)
    call read_integer_argument(4, "NZETA", n_zeta)
    call read_real_argument(5, "FLOOR", floor)
    call read_integer_argument(6, "PARITY", parity)
    if (degree < 1 .or. degree > 4) &
        call fail_usage("DEGREE must be between 1 and 4")
    if (n_theta < 8 .or. n_zeta < 8) &
        call fail_usage("NTHETA and NZETA must be at least 8")
    if (floor <= 0.0_dp) call fail_usage("FLOOR must be positive")
    if (parity < 1 .or. parity > 2) call fail_usage("PARITY must be 1 or 2")
    first_mode = 7
    inertia_only = .false.
    physical_mass = .false.
    has_m0_power = .false.
    density = 0.0_dp
    m0_stored_power = 0.0_dp
    do while (first_mode <= arguments)
        call read_argument(first_mode, "mode", token)
        if (index(token, "--") /= 1) exit
        if (trim(token) == "--inertia-only") then
            if (inertia_only) call fail_usage("duplicate --inertia-only")
            inertia_only = .true.
        else if (index(token, "--physical-density=") == 1) then
            if (physical_mass) call fail_usage("duplicate --physical-density")
            if (len_trim(token) == len("--physical-density=")) &
                call fail_usage("--physical-density requires a value")
            call parse_real(token(len("--physical-density=") + 1:), &
                "physical density", density)
            if (density <= 0.0_dp) &
                call fail_usage("physical density must be positive")
            physical_mass = .true.
        else if (index(token, "--m0-stored-power=") == 1) then
            if (has_m0_power) call fail_usage("duplicate --m0-stored-power")
            if (len_trim(token) == len("--m0-stored-power=")) &
                call fail_usage("--m0-stored-power requires a value")
            call parse_real(token(len("--m0-stored-power=") + 1:), &
                "m=0 stored power", m0_stored_power)
            if (m0_stored_power < 0.0_dp .or. m0_stored_power > 1.0_dp) &
                call fail_usage("m=0 stored power must lie in [0,1]")
            has_m0_power = .true.
        else
            call fail_usage("unknown option " // trim(token))
        end if
        first_mode = first_mode + 1
    end do
    if (arguments < first_mode) call fail_usage("at least one mode is required")
    allocate (mode_m(arguments - first_mode + 1), &
        mode_n(arguments - first_mode + 1), &
        stored_power(arguments - first_mode + 1), stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("mode allocation", -1)
    do trial = first_mode, arguments
        mode = trial - first_mode + 1
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
        if (mode_m(mode) == 0 .and. has_m0_power) &
            stored_power(mode) = m0_stored_power
    end do
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) call fail_solver("reader", info)
    if (physical_mass) then
        call build_compatible_two_component_problem(equilibrium, mode_m, &
            mode_n, stored_power, parity, degree, n_theta, n_zeta, problem, &
            info, density_kg_m3=density)
    else
        call build_compatible_two_component_problem(equilibrium, mode_m, &
            mode_n, stored_power, parity, degree, n_theta, n_zeta, problem, &
            info)
    end if
    if (info /= compatible_problem_ok) &
        call fail_solver("compatible marginality problem", info)
    allocate (stiffness_operator(size(problem%stiffness, 1), &
        size(problem%stiffness, 2)), stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("operator allocation", -1)
    call sum_tensor(problem%stiffness_terms, stiffness_operator)
    allocate (stiffness(size(stiffness_operator, 1), &
        size(stiffness_operator, 2)), stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("stiffness allocation", -1)
    call add_mass_shift(stiffness_operator, problem%mass, floor, stiffness)
    if (size(stiffness, 1) > huge(work_size) / 64) &
        call fail_solver("workspace size", -1)
    work_size = 64 * size(stiffness, 1)
    allocate (pivots(size(stiffness, 1)), work(work_size), &
        stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("workspace allocation", -1)
    call dsytrf("U", size(stiffness, 1), stiffness, size(stiffness, 1), &
        pivots, work, size(work), info)
    if (info /= 0) call fail_solver("shifted stiffness inertia", info)
    inertia_negative_count = pivot_negative_count(stiffness, pivots)
    if (inertia_only) then
        write (*, "(a)") "degree,parity,unknowns,normal_unknowns," // &
            "eta_unknowns,inertia_negative_count,physical_mass," // &
            "m0_stored_power"
        write (*, "(i0,6(a,i0),a,es24.16)") degree, ",", parity, ",", &
            size(stiffness, 1), ",", problem%normal_unknowns, ",", &
            problem%eta_unknowns, ",", inertia_negative_count, ",", &
            merge(1, 0, physical_mass), ",", m0_stored_power
        call terminate_process(0_c_int)
    end if
    deallocate (stiffness)
    allocate (stiffness, source=stiffness_operator, stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("stiffness allocation", -1)
    allocate (mass, source=problem%mass, stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("mass allocation", -1)
    call solve_symmetric_generalized_allocated(stiffness, mass, eigenvalues, &
        eigenvectors, info, .true.)
    if (info /= symmetric_eigensolver_ok) &
        call fail_solver("compatible marginality solve", info)
    negative_count = count(eigenvalues < -floor)
    call evaluate_generalized_eigenpair(stiffness_operator, problem%mass, &
        eigenvectors(:, 1), eigenvalues(1), kinetic, potential, residual, info)
    if (info /= 0) call fail_solver("compatible marginality diagnostics", info)
    do info = 1, 4
        call quadratic_form(problem%stiffness_terms(:, :, info), &
            eigenvectors(:, 1), term_energy(info), trial)
        if (trial /= 0) call fail_solver("compatible term energy", trial)
    end do
    write (*, "(a)") "degree,parity,unknowns,normal_unknowns," // &
        "eta_unknowns,negative_count,inertia_negative_count,physical_mass," // &
        "m0_stored_power,lowest_eigenvalue,residual,kinetic," // &
        "potential,bending,shear,compression,drive"
    write (*, "(i0, 7(a, i0), 9(a, es24.16))") degree, ",", parity, ",", &
        size(eigenvalues), ",", problem%normal_unknowns, ",", &
        problem%eta_unknowns, ",", negative_count, ",", &
        inertia_negative_count, ",", merge(1, 0, physical_mass), &
        ",", m0_stored_power, ",", eigenvalues(1), &
        ",", residual, ",", kinetic, ",", potential, ",", term_energy(1), &
        ",", term_energy(2), ",", term_energy(3), ",", term_energy(4)
    write (*, "(a)") "eigenvalue_index,eigenvalue"
    do info = 1, min(20, size(eigenvalues))
        write (*, "(i0,',',es24.16)") info, eigenvalues(info)
    end do

contains

    subroutine add_mass_shift(operator, mass_matrix, shift, shifted)
        real(dp), intent(in) :: operator(:, :), mass_matrix(:, :), shift
        real(dp), intent(out) :: shifted(:, :)
        integer :: column, row

        do column = 1, size(operator, 2)
            do row = 1, size(operator, 1)
                shifted(row, column) = operator(row, column) &
                    + shift * mass_matrix(row, column)
            end do
        end do
    end subroutine add_mass_shift

    subroutine fail_usage(message)
        character(len=*), intent(in) :: message
        write (error_unit, "(a)") &
            "gliss_compatible_marginality: " // trim(message)
        write (error_unit, "(a)") "usage: gliss_compatible_marginality " // &
            "EXPORT_FILE DEGREE NTHETA NZETA FLOOR PARITY " // &
            "[--inertia-only] [--physical-density=VALUE] " // &
            "[--m0-stored-power=VALUE] m,n [m,n ...]"
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

end program gliss_compatible_marginality
