program gliss_compatible_marginality
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, &
        int64, iostat_end
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_problem_ok, &
        compatible_two_component_problem_t, &
        evaluate_compatible_two_component_vector
    use compatible_problem_assembly_support, only: &
        evaluate_generalized_eigenpair, quadratic_form, &
        quadratic_form_with_absolute_sum, sum_tensor
    use dense_generalized_inverse_iteration, only: dense_inverse_ok, &
        solve_dense_generalized_near_shift, &
        solve_dense_generalized_subspace_near_shift
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use symmetric_pivot_inertia, only: pivot_negative_count
    use symmetric_eigensolver, only: solve_symmetric_generalized, &
        solve_symmetric_generalized_allocated, symmetric_eigensolver_ok
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(compatible_two_component_problem_t) :: problem
    character(len=1024) :: count_shift_token, eta_power_token, filename
    character(len=1024) :: external_vector_file
    character(len=1024) :: external_subspace_file
    character(len=1024) :: stored_power_token, token
    integer, allocatable :: mode_m(:), mode_n(:)
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
    real(dp), allocatable :: eigenvector(:)
    real(dp), allocatable :: count_shifts(:)
    real(dp), allocatable :: eta_stored_power(:), stored_power(:)
    real(dp), allocatable :: profile_coordinates(:), profile_eta(:, :)
    real(dp), allocatable :: profile_normal(:, :)
    real(dp), allocatable :: stiffness(:, :), mass(:, :)
    real(dp), allocatable :: stiffness_operator(:, :), work(:)
    integer, allocatable :: pivots(:)
    real(dp) :: bracket_lower, bracket_tolerance, bracket_upper
    real(dp) :: density, floor, kinetic, m0_stored_power, potential, residual
    real(dp) :: subspace_shift
    real(dp) :: term_energy(4)
    integer :: allocation_status, arguments, comma, degree, first_mode, info
    integer :: bracket_kind, inertia_negative_count, mode, negative_count
    integer :: eigenprofile_index, profile_points
    integer :: subspace_iterations
    integer :: work_size
    integer :: n_theta, n_zeta, parity, trial
    integer(int64) :: requested_work_size
    logical :: has_count_shifts, has_eigenprofile, has_eta_powers, has_m0_power
    logical :: has_external_vector
    logical :: has_external_subspace
    logical :: has_subspace_iterations, has_subspace_shift
    logical :: has_profile_points
    logical :: has_stored_powers, inertia_only
    logical :: minimize_compression, minimize_potential
    logical :: physical_mass

    interface
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
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
        subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dsytrs
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
    minimize_compression = .false.
    minimize_potential = .false.
    bracket_kind = 0
    has_m0_power = .false.
    has_stored_powers = .false.
    has_eta_powers = .false.
    has_count_shifts = .false.
    has_eigenprofile = .false.
    has_external_vector = .false.
    has_external_subspace = .false.
    has_subspace_iterations = .false.
    has_subspace_shift = .false.
    has_profile_points = .false.
    eigenprofile_index = 0
    profile_points = 201
    density = 0.0_dp
    m0_stored_power = 0.0_dp
    subspace_iterations = 0
    subspace_shift = 0.0_dp
    do while (first_mode <= arguments)
        call read_argument(first_mode, "mode", token)
        if (index(token, "--") /= 1) exit
        if (trim(token) == "--inertia-only") then
            if (inertia_only) call fail_usage("duplicate --inertia-only")
            inertia_only = .true.
        else if (trim(token) == "--minimize-compression") then
            if (minimize_compression) &
                call fail_usage("duplicate --minimize-compression")
            minimize_compression = .true.
        else if (trim(token) == "--minimize-potential") then
            if (minimize_potential) &
                call fail_usage("duplicate --minimize-potential")
            minimize_potential = .true.
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
        else if (index(token, "--stored-powers=") == 1) then
            if (has_stored_powers) call fail_usage("duplicate --stored-powers")
            if (len_trim(token) == len("--stored-powers=")) &
                call fail_usage("--stored-powers requires a comma-separated list")
            stored_power_token = token(len("--stored-powers=") + 1:)
            has_stored_powers = .true.
        else if (index(token, "--eta-stored-powers=") == 1) then
            if (has_eta_powers) &
                call fail_usage("duplicate --eta-stored-powers")
            if (len_trim(token) == len("--eta-stored-powers=")) &
                call fail_usage( &
                "--eta-stored-powers requires a comma-separated list")
            eta_power_token = token(len("--eta-stored-powers=") + 1:)
            has_eta_powers = .true.
        else if (index(token, "--negative-bracket=") == 1) then
            if (bracket_kind == 1) &
                call fail_usage("duplicate --negative-bracket")
            if (bracket_kind /= 0 .or. has_count_shifts) &
                call fail_usage("conflicting eigenvalue bracket options")
            if (len_trim(token) == len("--negative-bracket=")) &
                call fail_usage("--negative-bracket requires lower,upper,tolerance")
            call parse_bracket(token(len("--negative-bracket=") + 1:), &
                bracket_lower, bracket_upper, bracket_tolerance)
            bracket_kind = 1
        else if (index(token, "--eigenvalue-bracket=") == 1) then
            if (bracket_kind == 2) &
                call fail_usage("duplicate --eigenvalue-bracket")
            if (bracket_kind /= 0 .or. has_count_shifts) &
                call fail_usage("conflicting eigenvalue bracket options")
            if (len_trim(token) == len("--eigenvalue-bracket=")) &
                call fail_usage( &
                "--eigenvalue-bracket requires lower,upper,tolerance")
            call parse_eigenvalue_bracket( &
                token(len("--eigenvalue-bracket=") + 1:), bracket_lower, &
                bracket_upper, bracket_tolerance)
            bracket_kind = 2
        else if (index(token, "--count-shifts=") == 1) then
            if (has_count_shifts) call fail_usage("duplicate --count-shifts")
            if (bracket_kind /= 0) &
                call fail_usage("conflicting eigenvalue count options")
            if (len_trim(token) == len("--count-shifts=")) &
                call fail_usage( &
                "--count-shifts requires a comma-separated list")
            count_shift_token = token(len("--count-shifts=") + 1:)
            has_count_shifts = .true.
        else if (index(token, "--eigenprofile-index=") == 1) then
            if (has_eigenprofile) &
                call fail_usage("duplicate --eigenprofile-index")
            if (len_trim(token) == len("--eigenprofile-index=")) &
                call fail_usage("--eigenprofile-index requires an integer")
            call parse_integer(token(len("--eigenprofile-index=") + 1:), &
                "eigenprofile index", eigenprofile_index)
            if (eigenprofile_index < 1) &
                call fail_usage("eigenprofile index must be positive")
            has_eigenprofile = .true.
        else if (index(token, "--profile-points=") == 1) then
            if (has_profile_points) call fail_usage("duplicate --profile-points")
            if (len_trim(token) == len("--profile-points=")) &
                call fail_usage("--profile-points requires an integer")
            call parse_integer(token(len("--profile-points=") + 1:), &
                "profile points", profile_points)
            if (profile_points < 8 .or. profile_points > 10000) &
                call fail_usage("profile points must lie in [8,10000]")
            has_profile_points = .true.
        else if (index(token, "--evaluate-vector=") == 1) then
            if (has_external_vector) &
                call fail_usage("duplicate --evaluate-vector")
            if (len_trim(token) == len("--evaluate-vector=")) &
                call fail_usage("--evaluate-vector requires a file")
            external_vector_file = token(len("--evaluate-vector=") + 1:)
            has_external_vector = .true.
        else if (index(token, "--evaluate-subspace=") == 1) then
            if (has_external_subspace) &
                call fail_usage("duplicate --evaluate-subspace")
            if (len_trim(token) == len("--evaluate-subspace=")) &
                call fail_usage("--evaluate-subspace requires a path list")
            external_subspace_file = token(len("--evaluate-subspace=") + 1:)
            has_external_subspace = .true.
        else if (index(token, "--subspace-shift=") == 1) then
            if (has_subspace_shift) &
                call fail_usage("duplicate --subspace-shift")
            if (len_trim(token) == len("--subspace-shift=")) &
                call fail_usage("--subspace-shift requires a value")
            call parse_real(token(len("--subspace-shift=") + 1:), &
                "subspace shift", subspace_shift)
            has_subspace_shift = .true.
        else if (index(token, "--subspace-iterations=") == 1) then
            if (has_subspace_iterations) &
                call fail_usage("duplicate --subspace-iterations")
            if (len_trim(token) == len("--subspace-iterations=")) &
                call fail_usage("--subspace-iterations requires an integer")
            call parse_integer(token(len("--subspace-iterations=") + 1:), &
                "subspace iterations", subspace_iterations)
            if (subspace_iterations < 1 .or. subspace_iterations > 1000) &
                call fail_usage("subspace iterations must lie in [1,1000]")
            has_subspace_iterations = .true.
        else
            call fail_usage("unknown option " // trim(token))
        end if
        first_mode = first_mode + 1
    end do
    if (arguments < first_mode) call fail_usage("at least one mode is required")
    if (has_m0_power .and. has_stored_powers) &
        call fail_usage("--m0-stored-power conflicts with --stored-powers")
    if (has_count_shifts .and. inertia_only) &
        call fail_usage("--count-shifts conflicts with --inertia-only")
    if (has_profile_points .and. .not. has_eigenprofile &
        .and. .not. has_external_vector) &
        call fail_usage( &
        "--profile-points requires --eigenprofile-index or --evaluate-vector")
    if (has_eigenprofile .and. has_count_shifts) &
        call fail_usage("--eigenprofile-index conflicts with --count-shifts")
    if (has_eigenprofile .and. inertia_only) &
        call fail_usage("--eigenprofile-index conflicts with --inertia-only")
    if (has_eigenprofile .and. bracket_kind == 0) &
        call fail_usage( &
        "--eigenprofile-index requires --eigenvalue-bracket")
    if (has_eigenprofile .and. bracket_kind == 1) &
        call fail_usage( &
        "--eigenprofile-index conflicts with --negative-bracket")
    if (has_external_vector .and. (has_eigenprofile .or. inertia_only &
        .or. has_count_shifts .or. bracket_kind /= 0)) &
        call fail_usage("--evaluate-vector conflicts with solver options")
    if (has_external_subspace .and. (has_external_vector .or. has_eigenprofile &
        .or. inertia_only .or. has_count_shifts .or. bracket_kind /= 0)) &
        call fail_usage("--evaluate-subspace conflicts with solver options")
    if (has_external_subspace .and. has_profile_points) &
        call fail_usage("--profile-points conflicts with --evaluate-subspace")
    if (has_subspace_shift .neqv. has_subspace_iterations) &
        call fail_usage( &
        "--subspace-shift and --subspace-iterations must be given together")
    if (has_subspace_shift .and. .not. has_external_subspace) &
        call fail_usage("subspace solver options require --evaluate-subspace")
    if (minimize_compression .and. .not. has_external_vector &
        .and. .not. has_external_subspace) &
        call fail_usage("eta minimization requires an external vector or subspace")
    if (minimize_potential .and. .not. has_external_vector &
        .and. .not. has_external_subspace) &
        call fail_usage("eta minimization requires an external vector or subspace")
    if (minimize_compression .and. minimize_potential) &
        call fail_usage("eta minimization options conflict")
    if (has_external_subspace) &
        call preflight_subspace_paths(trim(external_subspace_file))
    allocate (mode_m(arguments - first_mode + 1), &
        mode_n(arguments - first_mode + 1), &
        stored_power(arguments - first_mode + 1), &
        eta_stored_power(arguments - first_mode + 1), stat=allocation_status)
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
    eta_stored_power = 0.0_dp
    if (has_stored_powers) call parse_real_list(stored_power_token, &
        "stored powers", stored_power)
    if (has_eta_powers) call parse_real_list(eta_power_token, &
        "eta stored powers", eta_stored_power)
    if (has_count_shifts) &
        call parse_count_shift_list(count_shift_token, count_shifts)
    if (has_stored_powers) then
        do mode = 1, size(mode_m)
            if (mode_m(mode) == 0) m0_stored_power = stored_power(mode)
        end do
    end if
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) call fail_solver("reader", info)
    if (physical_mass) then
        call build_compatible_two_component_problem(equilibrium, mode_m, &
            mode_n, stored_power, parity, degree, n_theta, n_zeta, problem, &
            info, density_kg_m3=density, &
            eta_stored_power=eta_stored_power)
    else
        call build_compatible_two_component_problem(equilibrium, mode_m, &
            mode_n, stored_power, parity, degree, n_theta, n_zeta, problem, &
            info, eta_stored_power=eta_stored_power)
    end if
    if (info /= compatible_problem_ok) &
        call fail_solver("compatible marginality problem", info)
    allocate (stiffness_operator(size(problem%stiffness, 1), &
        size(problem%stiffness, 2)), stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("operator allocation", -1)
    call sum_tensor(problem%stiffness_terms, stiffness_operator)
    call write_problem_metadata
    if (has_external_vector) then
        call write_external_vector_energy
        call terminate_process(0_c_int)
    end if
    if (has_external_subspace) then
        call write_external_subspace_energy
        call terminate_process(0_c_int)
    end if
    allocate (stiffness(size(stiffness_operator, 1), &
        size(stiffness_operator, 2)), stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("stiffness allocation", -1)
    requested_work_size = 64_int64 * int(size(stiffness, 1), int64)
    if (requested_work_size > int(huge(work_size), int64)) &
        call fail_solver("workspace size", -1)
    work_size = int(requested_work_size)
    allocate (pivots(size(stiffness, 1)), work(work_size), &
        stat=allocation_status)
    if (allocation_status /= 0) call fail_solver("workspace allocation", -1)
    if (has_eigenprofile) then
        call write_eigenprofile
        call terminate_process(0_c_int)
    end if
    if (bracket_kind /= 0) then
        if (bracket_kind == 1) then
            call bracket_negative_eigenvalue
        else
            call bracket_eigenvalue
        end if
        call terminate_process(0_c_int)
    end if
    if (has_count_shifts) then
        call write_shift_counts
        call terminate_process(0_c_int)
    end if
    call shifted_negative_count(floor, inertia_negative_count)
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

    subroutine write_problem_metadata
        write (*, "(a)") "stored_power_index,m,n,power"
        do mode = 1, size(mode_m)
            write (*, "(i0,2(a,i0),a,es24.16)") mode, ",", mode_m(mode), &
                ",", mode_n(mode), ",", stored_power(mode)
        end do
        write (*, "(a)") "eta_stored_power_index,m,n,power"
        do mode = 1, size(mode_m)
            write (*, "(i0,2(a,i0),a,es24.16)") mode, ",", mode_m(mode), &
                ",", mode_n(mode), ",", eta_stored_power(mode)
        end do
        write (*, "(a)") "radial_surfaces,degree,parity,n_theta,n_zeta," // &
            "physical_mass,density_kg_m3,floor,bracket_requested"
        write (*, "(i0,5(a,i0),2(a,es24.16),a,i0)") size(equilibrium%s), &
            ",", degree, ",", parity, ",", n_theta, ",", n_zeta, ",", &
            merge(1, 0, physical_mass), ",", density, ",", floor, ",", &
            merge(1, 0, bracket_kind /= 0)
    end subroutine write_problem_metadata

    subroutine write_external_vector_energy
        real(dp) :: closure, term_sum
        integer :: point, term

        allocate (eigenvector(size(problem%mass, 1)), stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("external vector allocation", -1)
        call read_external_vector(trim(external_vector_file), eigenvector)
        if (minimize_compression) call project_minimizing_eta( &
            problem%stiffness_terms(:, :, 3))
        if (minimize_potential) call project_minimizing_eta(stiffness_operator)
        call quadratic_form(problem%mass, eigenvector, kinetic, info)
        if (info /= 0 .or. kinetic <= 0.0_dp) &
            call fail_external_vector("vector has nonpositive mass norm")
        call quadratic_form(stiffness_operator, eigenvector, potential, info)
        if (info /= 0) call fail_solver("external vector potential", info)
        term_sum = 0.0_dp
        do term = 1, 4
            call quadratic_form(problem%stiffness_terms(:, :, term), &
                eigenvector, term_energy(term), info)
            if (info /= 0) call fail_solver("external vector term", info)
            term_sum = term_sum + term_energy(term)
        end do
        closure = abs(potential - term_sum) / kinetic
        write (*, "(a)") "unknowns,eta_relaxation,kinetic,potential," // &
            "rayleigh_quotient,bending,shear,compression,drive," // &
            "energy_closure_absolute"
        write (*, "(i0,a,i0,8(a,es24.16))") size(eigenvector), ",", &
            merge(1, merge(2, 0, minimize_potential), &
            minimize_compression), ",", kinetic, ",", &
            potential, ",", potential / kinetic, ",", &
            term_energy(1) / kinetic, ",", term_energy(2) / kinetic, ",", &
            term_energy(3) / kinetic, ",", term_energy(4) / kinetic, ",", &
            closure
        allocate (profile_coordinates(profile_points), stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("external profile coordinate allocation", -1)
        do point = 1, profile_points
            profile_coordinates(point) = &
                (real(point, dp) - 0.5_dp) / real(profile_points, dp)
        end do
        call evaluate_compatible_two_component_vector(size(equilibrium%s), &
            mode_m, mode_n, stored_power, eta_stored_power, parity, degree, &
            profile_coordinates, eigenvector, profile_normal, profile_eta, &
            info)
        if (info /= compatible_problem_ok) &
            call fail_solver("external profile vector evaluation", info)
        write (*, "(a)") "profile_index,s,r,mode_index,m,n,normal,eta"
        do point = 1, profile_points
            do mode = 1, size(mode_m)
                write (*, "(i0,2(a,es24.16),3(a,i0),2(a,es24.16))") &
                    point, ",", profile_coordinates(point), ",", &
                    sqrt(profile_coordinates(point)), ",", mode, ",", &
                    mode_m(mode), ",", mode_n(mode), ",", &
                    profile_normal(point, mode), ",", profile_eta(point, mode)
            end do
        end do
    end subroutine write_external_vector_energy

    subroutine write_external_subspace_energy
        character(len=1024), allocatable :: paths(:)
        real(dp), allocatable :: block_eigenvalues(:), block_residuals(:)
        real(dp), allocatable :: converged_mass(:, :)
        real(dp), allocatable :: converged_stiffness(:, :)
        real(dp), allocatable :: converged_vectors(:, :)
        real(dp), allocatable :: initial_final_overlap(:, :)
        real(dp), allocatable :: reduced_mass(:, :), reduced_stiffness(:, :)
        real(dp), allocatable :: reduced_term(:, :, :), reduced_vectors(:, :)
        real(dp), allocatable :: ritz_values(:), ritz_vector(:)
        real(dp), allocatable :: image(:, :), vectors(:, :)
        real(dp) :: closure, value
        integer :: block_iterations, column, dimension, first, row, second, term

        call read_subspace_paths(trim(external_subspace_file), paths)
        dimension = size(paths)
        allocate (vectors(size(problem%mass, 1), dimension), &
            reduced_mass(dimension, dimension), &
            reduced_stiffness(dimension, dimension), &
            reduced_term(dimension, dimension, 4), &
            image(size(problem%mass, 1), dimension), &
            ritz_vector(size(problem%mass, 1)), stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("external subspace allocation", -1)
        do column = 1, dimension
            call read_external_vector(trim(paths(column)), vectors(:, column))
        end do
        if (minimize_compression) call project_minimizing_eta_block( &
            problem%stiffness_terms(:, :, 3), vectors)
        if (minimize_potential) call project_minimizing_eta_block( &
            stiffness_operator, vectors)
        call contract_subspace(problem%mass, vectors, reduced_mass, image)
        call contract_subspace(stiffness_operator, vectors, &
            reduced_stiffness, image)
        do term = 1, 4
            call contract_subspace(problem%stiffness_terms(:, :, term), &
                vectors, reduced_term(:, :, term), image)
        end do
        call solve_symmetric_generalized(reduced_stiffness, reduced_mass, &
            ritz_values, reduced_vectors, info)
        if (info /= symmetric_eigensolver_ok) &
            call fail_external_subspace("reduced mass is not positive definite")
        write (*, "(a)") "subspace_dimension,unknowns,eta_relaxation"
        write (*, "(i0,a,i0,a,i0)") dimension, ",", size(vectors, 1), ",", &
            merge(1, merge(2, 0, minimize_potential), minimize_compression)
        write (*, "(a)") "reduced_matrix,kind,row,column,value"
        do second = 1, dimension
            do first = 1, dimension
                write (*, "(a,2(a,i0),a,es24.16)") "mass", ",", first, &
                    ",", second, ",", reduced_mass(first, second)
                write (*, "(a,2(a,i0),a,es24.16)") "stiffness", ",", &
                    first, ",", second, ",", &
                    reduced_stiffness(first, second)
                do term = 1, 4
                    write (*, "(a,i0,2(a,i0),a,es24.16)") "term", term, &
                        ",", first, ",", second, ",", &
                        reduced_term(first, second, term)
                end do
            end do
        end do
        write (*, "(a)") "ritz_index,eigenvalue,residual,kinetic,potential," // &
            "bending,shear,compression,drive,energy_closure_absolute"
        do column = 1, dimension
            do row = 1, size(ritz_vector)
                ritz_vector(row) = 0.0_dp
                do first = 1, dimension
                    ritz_vector(row) = ritz_vector(row) &
                        + vectors(row, first) * reduced_vectors(first, column)
                end do
            end do
            call evaluate_generalized_eigenpair(stiffness_operator, &
                problem%mass, ritz_vector, ritz_values(column), kinetic, &
                potential, residual, info)
            if (info /= 0) call fail_solver("subspace Ritz diagnostics", info)
            closure = potential
            do term = 1, 4
                call quadratic_form(problem%stiffness_terms(:, :, term), &
                    ritz_vector, value, info)
                if (info /= 0) call fail_solver("subspace Ritz term", info)
                term_energy(term) = value
                closure = closure - value
            end do
            write (*, "(i0,9(a,es24.16))") column, ",", &
                ritz_values(column), ",", residual, ",", kinetic, ",", &
                potential, ",", term_energy(1), ",", term_energy(2), ",", &
                term_energy(3), ",", term_energy(4), ",", abs(closure)
        end do
        if (.not. has_subspace_shift) return
        call solve_dense_generalized_subspace_near_shift( &
            stiffness_operator, problem%mass, subspace_shift, vectors, &
            subspace_iterations, block_eigenvalues, converged_vectors, &
            block_residuals, block_iterations, info)
        if (info /= dense_inverse_ok) &
            call fail_solver("subspace block inverse iteration", info)
        allocate (converged_mass(dimension, dimension), &
            converged_stiffness(dimension, dimension), &
            initial_final_overlap(dimension, dimension), &
            stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("converged subspace allocation", -1)
        call contract_subspace(problem%mass, converged_vectors, &
            converged_mass, image)
        call contract_subspace(stiffness_operator, converged_vectors, &
            converged_stiffness, image)
        call contract_cross_subspace(problem%mass, vectors, converged_vectors, &
            initial_final_overlap, image)
        write (*, "(a)") &
            "subspace_solver,shift,iteration_limit,iterations"
        write (*, "(a,a,es24.16,2(a,i0))") "block_inverse", ",", &
            subspace_shift, ",", subspace_iterations, ",", block_iterations
        write (*, "(a)") "converged_matrix,kind,row,column,value"
        do second = 1, dimension
            do first = 1, dimension
                write (*, "(a,2(a,i0),a,es24.16)") "mass", ",", first, &
                    ",", second, ",", converged_mass(first, second)
                write (*, "(a,2(a,i0),a,es24.16)") "stiffness", ",", &
                    first, ",", second, ",", &
                    converged_stiffness(first, second)
            end do
        end do
        write (*, "(a)") "initial_final_mass_overlap,row,column,value"
        do second = 1, dimension
            do first = 1, dimension
                write (*, "(2(i0,a),es24.16)") first, ",", second, ",", &
                    initial_final_overlap(first, second)
            end do
        end do
        write (*, "(a)") "converged_ritz_index,eigenvalue," // &
            "block_residual,diagnostic_residual,kinetic,potential," // &
            "bending,shear,compression,drive,energy_closure_absolute"
        do column = 1, dimension
            call evaluate_generalized_eigenpair(stiffness_operator, &
                problem%mass, converged_vectors(:, column), &
                block_eigenvalues(column), kinetic, potential, residual, info)
            if (info /= 0) &
                call fail_solver("converged subspace diagnostics", info)
            closure = potential
            do term = 1, 4
                call quadratic_form(problem%stiffness_terms(:, :, term), &
                    converged_vectors(:, column), value, info)
                if (info /= 0) &
                    call fail_solver("converged subspace term", info)
                term_energy(term) = value
                closure = closure - value
            end do
            write (*, "(i0,10(a,es24.16))") column, ",", &
                block_eigenvalues(column), ",", block_residuals(column), &
                ",", residual, ",", kinetic, ",", potential, ",", &
                term_energy(1), ",", term_energy(2), ",", term_energy(3), &
                ",", term_energy(4), ",", abs(closure)
        end do
    end subroutine write_external_subspace_energy

    subroutine read_subspace_paths(path, paths)
        character(len=*), intent(in) :: path
        character(len=1024), allocatable, intent(out) :: paths(:)
        character(len=1024) :: line
        integer :: column, count, io_status, item, unit

        open (newunit=unit, file=path, status="old", action="read", &
            form="formatted", iostat=io_status)
        if (io_status /= 0) call fail_external_subspace("cannot open path list")
        read (unit, "(a)", iostat=io_status) line
        if (io_status /= 0 .or. trim(line) /= "path") &
            call fail_external_subspace("first path-list row must be path")
        count = 0
        do
            read (unit, "(a)", iostat=io_status) line
            if (io_status == iostat_end) exit
            if (io_status /= 0) call fail_external_subspace("cannot read path list")
            if (len_trim(line) == 0) &
                call fail_external_subspace("path list contains an empty row")
            if (len_trim(line) == len(line)) &
                call fail_external_subspace("path-list row is too long")
            count = count + 1
            if (count > 64) &
                call fail_external_subspace("path list exceeds 64 vectors")
        end do
        if (count < 1) call fail_external_subspace("path list contains no vectors")
        rewind (unit)
        read (unit, "(a)", iostat=io_status) line
        allocate (paths(count), stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("external subspace path allocation", -1)
        do item = 1, count
            read (unit, "(a)", iostat=io_status) line
            if (io_status /= 0) call fail_external_subspace("path list changed")
            paths(item) = trim(line)
            do column = 1, item - 1
                if (paths(column) == paths(item)) &
                    call fail_external_subspace("path list contains duplicates")
            end do
        end do
        close (unit)
    end subroutine read_subspace_paths

    subroutine preflight_subspace_paths(path)
        character(len=*), intent(in) :: path
        character(len=1024), allocatable :: paths(:)
        character(len=4096) :: line
        integer :: io_status, item, unit

        call read_subspace_paths(path, paths)
        do item = 1, size(paths)
            open (newunit=unit, file=trim(paths(item)), status="old", &
                action="read", form="formatted", iostat=io_status)
            if (io_status /= 0) &
                call fail_external_subspace("cannot open a listed vector")
            read (unit, "(a)", iostat=io_status) line
            if (io_status /= 0 .or. trim(line) /= "index,value") &
                call fail_external_subspace( &
                "listed vector does not begin with index,value")
            close (unit)
        end do
    end subroutine preflight_subspace_paths

    subroutine read_external_vector(path, vector)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: vector(:)
        character(len=4096) :: line
        integer :: file_index, io_status, point, unit

        open (newunit=unit, file=path, status="old", action="read", &
            form="formatted", iostat=io_status)
        if (io_status /= 0) call fail_external_vector("cannot open file")
        read (unit, "(a)", iostat=io_status) line
        if (io_status /= 0 .or. trim(line) /= "index,value") &
            call fail_external_vector("first row must be index,value")
        do point = 1, size(vector)
            read (unit, "(a)", iostat=io_status) line
            if (io_status /= 0) call fail_external_vector( &
                "file ends before the expected vector size")
            if (len_trim(line) == len(line)) &
                call fail_external_vector("input row is too long")
            call parse_external_vector_row(line, file_index, vector(point))
            if (file_index /= point) call fail_external_vector( &
                "indices must be consecutive from one")
        end do
        read (unit, "(a)", iostat=io_status) line
        if (io_status /= iostat_end) call fail_external_vector( &
            "file contains rows after the expected vector")
        close (unit)
    end subroutine read_external_vector

    subroutine contract_subspace(matrix, vectors, reduced, image)
        real(dp), contiguous, intent(in) :: matrix(:, :), vectors(:, :)
        real(dp), contiguous, intent(out) :: reduced(:, :), image(:, :)
        integer :: first, row, second

        call dgemm("N", "N", size(matrix, 1), size(vectors, 2), &
            size(matrix, 2), 1.0_dp, matrix, size(matrix, 1), vectors, &
            size(vectors, 1), 0.0_dp, image, size(image, 1))
        do second = 1, size(vectors, 2)
            do first = 1, second
                reduced(first, second) = 0.0_dp
                do row = 1, size(vectors, 1)
                    reduced(first, second) = reduced(first, second) &
                        + vectors(row, first) * image(row, second)
                end do
                reduced(second, first) = reduced(first, second)
            end do
        end do
    end subroutine contract_subspace

    subroutine contract_cross_subspace(matrix, left, right, reduced, image)
        real(dp), contiguous, intent(in) :: matrix(:, :), left(:, :)
        real(dp), contiguous, intent(in) :: right(:, :)
        real(dp), contiguous, intent(out) :: reduced(:, :), image(:, :)

        call dgemm("N", "N", size(matrix, 1), size(right, 2), &
            size(matrix, 2), 1.0_dp, matrix, size(matrix, 1), right, &
            size(right, 1), 0.0_dp, image, size(image, 1))
        call dgemm("T", "N", size(left, 2), size(right, 2), &
            size(left, 1), 1.0_dp, left, size(left, 1), image, &
            size(image, 1), 0.0_dp, reduced, size(reduced, 1))
    end subroutine contract_cross_subspace

    subroutine project_minimizing_eta(operator)
        real(dp), intent(in) :: operator(:, :)
        real(dp), allocatable :: eta_matrix(:, :), rhs(:, :), solve_work(:)
        integer, allocatable :: eta_pivots(:)
        integer :: column, eta_unknowns, normal_unknowns, row, solve_work_size

        normal_unknowns = problem%normal_unknowns
        eta_unknowns = problem%eta_unknowns
        solve_work_size = 64 * eta_unknowns
        allocate (eta_matrix(eta_unknowns, eta_unknowns), &
            rhs(eta_unknowns, 1), eta_pivots(eta_unknowns), &
            solve_work(solve_work_size), stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("eta projection allocation", -1)
        do column = 1, eta_unknowns
            do row = 1, eta_unknowns
                eta_matrix(row, column) = operator( &
                    normal_unknowns + row, normal_unknowns + column)
            end do
        end do
        rhs = 0.0_dp
        do column = 1, normal_unknowns
            do row = 1, eta_unknowns
                rhs(row, 1) = rhs(row, 1) - operator( &
                    normal_unknowns + row, column) * eigenvector(column)
            end do
        end do
        call dsytrf("U", eta_unknowns, eta_matrix, eta_unknowns, eta_pivots, &
            solve_work, solve_work_size, info)
        if (info /= 0) call fail_solver("eta projection factor", info)
        call dsytrs("U", eta_unknowns, 1, eta_matrix, eta_unknowns, &
            eta_pivots, rhs, eta_unknowns, info)
        if (info /= 0) call fail_solver("eta projection solve", info)
        do row = 1, eta_unknowns
            eigenvector(normal_unknowns + row) = rhs(row, 1)
        end do
    end subroutine project_minimizing_eta

    subroutine project_minimizing_eta_block(operator, vectors)
        real(dp), intent(in) :: operator(:, :)
        real(dp), intent(inout) :: vectors(:, :)
        real(dp), allocatable :: eta_matrix(:, :), rhs(:, :), solve_work(:)
        integer, allocatable :: eta_pivots(:)
        integer :: column, eta_unknowns, normal_unknowns, row, solve_work_size
        integer :: vector

        normal_unknowns = problem%normal_unknowns
        eta_unknowns = problem%eta_unknowns
        solve_work_size = 64 * eta_unknowns
        allocate (eta_matrix(eta_unknowns, eta_unknowns), &
            rhs(eta_unknowns, size(vectors, 2)), eta_pivots(eta_unknowns), &
            solve_work(solve_work_size), stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("eta subspace projection allocation", -1)
        do column = 1, eta_unknowns
            do row = 1, eta_unknowns
                eta_matrix(row, column) = operator( &
                    normal_unknowns + row, normal_unknowns + column)
            end do
        end do
        rhs = 0.0_dp
        do vector = 1, size(vectors, 2)
            do column = 1, normal_unknowns
                do row = 1, eta_unknowns
                    rhs(row, vector) = rhs(row, vector) - operator( &
                        normal_unknowns + row, column) * vectors(column, vector)
                end do
            end do
        end do
        call dsytrf("U", eta_unknowns, eta_matrix, eta_unknowns, eta_pivots, &
            solve_work, solve_work_size, info)
        if (info /= 0) call fail_solver("eta subspace projection factor", info)
        call dsytrs("U", eta_unknowns, size(vectors, 2), eta_matrix, &
            eta_unknowns, eta_pivots, rhs, eta_unknowns, info)
        if (info /= 0) call fail_solver("eta subspace projection solve", info)
        do vector = 1, size(vectors, 2)
            do row = 1, eta_unknowns
                vectors(normal_unknowns + row, vector) = rhs(row, vector)
            end do
        end do
    end subroutine project_minimizing_eta_block

    subroutine parse_external_vector_row(line, index_value, value)
        character(len=*), intent(in) :: line
        integer, intent(out) :: index_value
        real(dp), intent(out) :: value
        integer :: comma

        comma = index(line, ",")
        if (comma <= 1 .or. comma == len_trim(line) &
            .or. index(line(comma + 1:), ",") > 0) &
            call fail_external_vector("rows must contain exactly index,value")
        call parse_integer(line(:comma - 1), "external vector index", &
            index_value)
        call parse_real(line(comma + 1:), "external vector value", value)
    end subroutine parse_external_vector_row

    subroutine fail_external_vector(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_compatible_marginality: " // &
            "invalid external vector: " // trim(message)
        call terminate_process(2_c_int)
    end subroutine fail_external_vector

    subroutine fail_external_subspace(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_compatible_marginality: " // &
            "invalid external subspace: " // trim(message)
        call terminate_process(2_c_int)
    end subroutine fail_external_subspace

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
            "[--m0-stored-power=VALUE] " // &
            "[--stored-powers=P0,P1,...] " // &
            "[--eta-stored-powers=P0,P1,...] " // &
            "[--negative-bracket=LOWER,UPPER,TOLERANCE] " // &
            "[--eigenvalue-bracket=LOWER,UPPER,TOLERANCE] " // &
            "[--count-shifts=S0,S1,...] " // &
            "[--eigenprofile-index=INDEX] [--profile-points=COUNT] " // &
            "[--evaluate-vector=INDEXED_CSV] " // &
            "[--evaluate-subspace=PATH_LIST] " // &
            "[--subspace-shift=VALUE --subspace-iterations=COUNT] " // &
            "[--minimize-compression] [--minimize-potential] " // &
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

    subroutine parse_real_list(text, name, values)
        character(len=*), intent(in) :: text, name
        real(dp), intent(out) :: values(:)
        integer :: comma, first, item

        first = 1
        do item = 1, size(values)
            comma = index(text(first:), ",")
            if (item < size(values)) then
                if (comma <= 1) &
                    call fail_usage(trim(name) // " requires one value per mode")
                call parse_real(text(first:first + comma - 2), name, values(item))
                first = first + comma
            else
                if (comma > 0 .or. first > len_trim(text)) &
                    call fail_usage(trim(name) // " requires one value per mode")
                call parse_real(text(first:), name, values(item))
            end if
            if (abs(values(item)) > 4.0_dp) &
                call fail_usage(trim(name) // " values must lie in [-4,4]")
        end do
    end subroutine parse_real_list

    subroutine parse_count_shift_list(text, values)
        character(len=*), intent(in) :: text
        real(dp), allocatable, intent(out) :: values(:)
        integer :: allocation_status, comma, first, item, items, position

        items = 1
        do position = 1, len_trim(text)
            if (text(position:position) == ",") items = items + 1
        end do
        allocate (values(items), stat=allocation_status)
        if (allocation_status /= 0) call fail_solver("count shift allocation", -1)
        first = 1
        do item = 1, items
            comma = index(text(first:), ",")
            if (item < items) then
                if (comma <= 1) call fail_usage( &
                    "count shifts require nonempty comma-separated values")
                call parse_real(text(first:first + comma - 2), &
                    "count shift", values(item))
                first = first + comma
            else
                if (comma > 0 .or. first > len_trim(text)) call fail_usage( &
                    "count shifts require nonempty comma-separated values")
                call parse_real(text(first:), "count shift", values(item))
            end if
            if (values(item) < 0.0_dp) &
                call fail_usage("count shifts must be nonnegative")
            if (item > 1) then
                if (values(item) <= values(item - 1)) &
                    call fail_usage("count shifts must be strictly increasing")
            end if
        end do
    end subroutine parse_count_shift_list

    subroutine parse_bracket(text, lower, upper, tolerance)
        character(len=*), intent(in) :: text
        real(dp), intent(out) :: lower, upper, tolerance
        integer :: first_comma, second_comma

        first_comma = index(text, ",")
        if (first_comma <= 1 .or. first_comma == len_trim(text)) &
            call fail_usage("negative bracket requires lower,upper,tolerance")
        second_comma = index(text(first_comma + 1:), ",")
        if (second_comma <= 1) &
            call fail_usage("negative bracket requires lower,upper,tolerance")
        second_comma = first_comma + second_comma
        if (second_comma == len_trim(text)) &
            call fail_usage("negative bracket requires lower,upper,tolerance")
        if (index(text(second_comma + 1:), ",") > 0) &
            call fail_usage("negative bracket requires lower,upper,tolerance")
        call parse_real(text(:first_comma - 1), "negative bracket lower", lower)
        call parse_real(text(first_comma + 1:second_comma - 1), &
            "negative bracket upper", upper)
        call parse_real(text(second_comma + 1:), &
            "negative bracket tolerance", tolerance)
        if (lower <= 0.0_dp .or. upper <= lower) &
            call fail_usage("negative bracket must obey 0 < lower < upper")
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
            call fail_usage("negative bracket tolerance must lie in (0,1)")
    end subroutine parse_bracket

    subroutine parse_eigenvalue_bracket(text, lower, upper, tolerance)
        character(len=*), intent(in) :: text
        real(dp), intent(out) :: lower, upper, tolerance
        integer :: first_comma, second_comma

        first_comma = index(text, ",")
        if (first_comma <= 1 .or. first_comma == len_trim(text)) &
            call fail_usage( &
            "eigenvalue bracket requires lower,upper,tolerance")
        second_comma = index(text(first_comma + 1:), ",")
        if (second_comma <= 1) &
            call fail_usage( &
            "eigenvalue bracket requires lower,upper,tolerance")
        second_comma = first_comma + second_comma
        if (second_comma == len_trim(text)) &
            call fail_usage( &
            "eigenvalue bracket requires lower,upper,tolerance")
        if (index(text(second_comma + 1:), ",") > 0) &
            call fail_usage( &
            "eigenvalue bracket requires lower,upper,tolerance")
        call parse_real(text(:first_comma - 1), &
            "eigenvalue bracket lower", lower)
        call parse_real(text(first_comma + 1:second_comma - 1), &
            "eigenvalue bracket upper", upper)
        call parse_real(text(second_comma + 1:), &
            "eigenvalue bracket tolerance", tolerance)
        if (upper <= lower) &
            call fail_usage("eigenvalue bracket must obey lower < upper")
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
            call fail_usage( &
            "eigenvalue bracket tolerance must lie in (0,1)")
    end subroutine parse_eigenvalue_bracket

    subroutine shifted_negative_count(shift, count)
        real(dp), intent(in) :: shift
        integer, intent(out) :: count

        call add_mass_shift(stiffness_operator, problem%mass, shift, stiffness)
        call dsytrf("U", size(stiffness, 1), stiffness, size(stiffness, 1), &
            pivots, work, size(work), info)
        if (info /= 0) call fail_solver("shifted stiffness inertia", info)
        count = pivot_negative_count(stiffness, pivots)
    end subroutine shifted_negative_count

    subroutine write_shift_counts
        integer :: count, item

        write (*, "(a)") "shift,eigenvalues_below_shift"
        do item = 1, size(count_shifts)
            call shifted_negative_count(-count_shifts(item), count)
            write (*, "(es24.16,',',i0)") count_shifts(item), count
        end do
    end subroutine write_shift_counts

    subroutine write_eigenprofile
        real(dp), parameter :: profile_angle_limit = 1.0e-3_dp
        real(dp) :: action_residual, angle_bound, backward_error
        real(dp) :: certified_eigenvalue, eigenvalue
        real(dp) :: bound_evaluation_factor
        real(dp) :: closure_absolute, closure_bound, closure_bound_relative
        real(dp) :: closure_ratio, closure_scale, closure_upper
        real(dp) :: equilibrated_action_residual, inverse_shift
        real(dp) :: matrix_reduction_bound, matrix_reduction_factor
        real(dp) :: mass_reciprocal_condition, neighbor_gap
        real(dp) :: outer_lower, outer_upper, standard_residual_norm
        real(dp) :: potential_absolute_sum, potential_error_bound
        real(dp) :: quadratic_reduction_factor
        real(dp) :: term_absolute_sum(4), term_error_bound(4)
        real(dp) :: term_energy_absolute_sum, term_energy_sum
        real(dp) :: term_reduction_bound, term_reduction_factor
        real(dp) :: term_weighted_absolute_sum
        real(dp) :: vector_residual_bound
        integer :: iterations, lower_count, point, term, upper_count

        outer_lower = bracket_lower
        outer_upper = bracket_upper
        inverse_shift = outer_lower + 0.5_dp * (outer_upper - outer_lower)
        call refine_eigenprofile_bracket(lower_count, upper_count)
        call solve_dense_generalized_near_shift(stiffness_operator, &
            problem%mass, inverse_shift, eigenvalue, eigenvector, &
            action_residual, equilibrated_action_residual, backward_error, &
            standard_residual_norm, mass_reciprocal_condition, iterations, &
            info)
        if (info /= dense_inverse_ok) then
            write (error_unit, "(a,es24.16,5(a,es24.16),a,i0)") &
                "eigenprofile last eigenvalue=", eigenvalue, &
                " action_residual=", action_residual, &
                " equilibrated_action_residual=", &
                equilibrated_action_residual, &
                " backward_error=", backward_error, &
                " standard_residual_norm=", standard_residual_norm, &
                " mass_reciprocal_condition=", mass_reciprocal_condition, &
                " iterations=", iterations
            call fail_solver("eigenprofile inverse iteration", info)
        end if
        certified_eigenvalue = bracket_lower &
            + 0.5_dp * (bracket_upper - bracket_lower)
        vector_residual_bound = standard_residual_norm &
            + max(abs(eigenvalue - bracket_lower), &
            abs(eigenvalue - bracket_upper))
        neighbor_gap = outer_upper - bracket_upper
        if (eigenprofile_index > 1) neighbor_gap = min(neighbor_gap, &
            bracket_lower - outer_lower)
        if (neighbor_gap <= 0.0_dp) &
            call fail_solver("eigenprofile neighbor separation", -1)
        angle_bound = vector_residual_bound / neighbor_gap
        if (.not. ieee_is_finite(angle_bound) &
            .or. angle_bound > profile_angle_limit) then
            write (error_unit, "(2(a,es24.16))") &
                "eigenprofile angle_bound=", angle_bound, &
                " limit=", profile_angle_limit
            call fail_solver("eigenprofile forward error", -1)
        end if
        call evaluate_generalized_eigenpair(stiffness_operator, problem%mass, &
            eigenvector, eigenvalue, kinetic, potential, residual, info)
        if (info /= 0) call fail_solver("eigenprofile diagnostics", info)
        call quadratic_form_with_absolute_sum(stiffness_operator, eigenvector, &
            potential, potential_absolute_sum, potential_error_bound, trial)
        if (trial /= 0) call fail_solver("eigenprofile potential energy", trial)
        term_energy_sum = 0.0_dp
        term_energy_absolute_sum = 0.0_dp
        term_weighted_absolute_sum = 0.0_dp
        do term = 1, 4
            call quadratic_form_with_absolute_sum( &
                problem%stiffness_terms(:, :, term), eigenvector, &
                term_energy(term), term_absolute_sum(term), &
                term_error_bound(term), trial)
            if (trial /= 0) call fail_solver("eigenprofile term energy", trial)
            term_energy_sum = term_energy_sum + term_energy(term)
            term_energy_absolute_sum = term_energy_absolute_sum &
                + abs(term_energy(term))
            term_weighted_absolute_sum = term_weighted_absolute_sum &
                + term_absolute_sum(term)
        end do
        quadratic_reduction_factor = roundoff_factor( &
            2 * size(eigenvector) + 2)
        matrix_reduction_factor = roundoff_factor(4)
        term_reduction_factor = roundoff_factor(4)
        matrix_reduction_bound = matrix_reduction_factor &
            * term_weighted_absolute_sum &
            / (1.0_dp - quadratic_reduction_factor)
        term_reduction_bound = term_reduction_factor &
            * term_energy_absolute_sum / (1.0_dp - term_reduction_factor)
        closure_bound = potential_error_bound + matrix_reduction_bound &
            + term_reduction_bound
        do term = 1, 4
            closure_bound = closure_bound + term_error_bound(term)
        end do
        ! The independently reduced total and terms can differ through the
        ! four-term matrix reduction, five quadratic-form evaluations and the
        ! reported four-term scalar reduction.  Inflate the positive bound for
        ! the finite-precision arithmetic used to evaluate the bound itself.
        bound_evaluation_factor = roundoff_factor(32)
        closure_bound = closure_bound / (1.0_dp - bound_evaluation_factor)
        closure_absolute = abs(potential - term_energy_sum)
        closure_scale = max(abs(potential), term_energy_absolute_sum, &
            tiny(1.0_dp))
        closure_upper = (closure_absolute + spacing(closure_scale)) &
            / (1.0_dp - matrix_reduction_factor)
        closure_ratio = closure_upper / closure_bound
        closure_bound_relative = closure_bound / closure_scale
        if (.not. ieee_is_finite(closure_ratio) &
            .or. closure_ratio > 1.0_dp) then
            write (error_unit, "(3(a,es24.16))") &
                "energy closure=", closure_absolute, &
                " forward bound=", closure_bound, &
                " ratio=", closure_ratio
            call fail_solver("eigenprofile energy closure", -1)
        end if
        allocate (profile_coordinates(profile_points), stat=allocation_status)
        if (allocation_status /= 0) &
            call fail_solver("eigenprofile coordinate allocation", -1)
        do point = 1, profile_points
            profile_coordinates(point) = &
                (real(point, dp) - 0.5_dp) / real(profile_points, dp)
        end do
        call evaluate_compatible_two_component_vector(size(equilibrium%s), &
            mode_m, mode_n, stored_power, eta_stored_power, parity, degree, &
            profile_coordinates, eigenvector, profile_normal, profile_eta, &
            info)
        if (info /= compatible_problem_ok) &
            call fail_solver("eigenprofile vector evaluation", info)
        write (*, "(a)") "selected_eigenpair_index,inverse_shift," // &
            "outer_lower,outer_upper,lambda_lower,lambda_upper,iterations," // &
            "certified_eigenvalue,rayleigh_quotient,action_residual," // &
            "equilibrated_action_residual,backward_error," // &
            "standard_residual_norm,angle_bound," // &
            "mass_reciprocal_condition,kinetic,potential,bending,shear," // &
            "compression,drive,energy_closure_absolute," // &
            "energy_closure_bound,energy_closure_ratio," // &
            "energy_closure_bound_relative"
        write (*, "(i0,5(a,es24.16),a,i0,18(a,es24.16))") &
            eigenprofile_index, ",", inverse_shift, ",", outer_lower, ",", &
            outer_upper, ",", bracket_lower, ",", bracket_upper, ",", &
            iterations, ",", certified_eigenvalue, ",", eigenvalue, ",", &
            action_residual, ",", equilibrated_action_residual, ",", &
            backward_error, ",", standard_residual_norm, ",", angle_bound, &
            ",", mass_reciprocal_condition, ",", kinetic, ",", potential, &
            ",", term_energy(1), ",", term_energy(2), ",", term_energy(3), &
            ",", term_energy(4), ",", closure_absolute, ",", closure_bound, &
            ",", closure_ratio, ",", closure_bound_relative
        write (*, "(a)") "profile_index,s,r,mode_index,m,n,normal,eta"
        do point = 1, profile_points
            do mode = 1, size(mode_m)
                write (*, "(i0,2(a,es24.16),3(a,i0),2(a,es24.16))") &
                    point, ",", profile_coordinates(point), ",", &
                    sqrt(profile_coordinates(point)), ",", mode, ",", &
                    mode_m(mode), ",", mode_n(mode), ",", &
                    profile_normal(point, mode), ",", profile_eta(point, mode)
            end do
        end do
    end subroutine write_eigenprofile

    pure function roundoff_factor(operations) result(gamma)
        integer, intent(in) :: operations
        real(dp) :: gamma, product

        product = real(operations, dp) * epsilon(1.0_dp)
        gamma = product / (1.0_dp - product)
    end function roundoff_factor

    subroutine refine_eigenprofile_bracket(lower_count, upper_count)
        character(len=192) :: message
        real(dp) :: midpoint, scale, target_tolerance
        integer, intent(out) :: lower_count, upper_count
        integer :: iterations, midpoint_count

        call shifted_negative_count(-bracket_lower, lower_count)
        call shifted_negative_count(-bracket_upper, upper_count)
        if (lower_count /= eigenprofile_index - 1 &
            .or. upper_count /= eigenprofile_index) then
            write (message, "(a,i0,a,i0,a,i0,a,i0,a)") &
                "eigenprofile bracket counts must be ", &
                eigenprofile_index - 1, " and ", eigenprofile_index, &
                " (found ", lower_count, " and ", upper_count, ")"
            call fail_usage(trim(message))
        end if
        iterations = 0
        scale = max(abs(bracket_lower), abs(bracket_upper), tiny(1.0_dp))
        target_tolerance = min(bracket_tolerance, 1.0e-10_dp)
        do while (bracket_upper - bracket_lower > target_tolerance * scale)
            if (iterations >= 80) &
                call fail_solver("eigenprofile bracket iteration", -1)
            midpoint = bracket_lower + 0.5_dp * &
                (bracket_upper - bracket_lower)
            call shifted_negative_count(-midpoint, midpoint_count)
            if (midpoint_count == lower_count) then
                bracket_lower = midpoint
            else if (midpoint_count == upper_count) then
                bracket_upper = midpoint
            else
                call fail_solver( &
                    "eigenprofile bracket changed by multiple levels", -1)
            end if
            iterations = iterations + 1
        end do
    end subroutine refine_eigenprofile_bracket

    subroutine bracket_negative_eigenvalue
        character(len=192) :: message
        real(dp) :: midpoint
        integer :: iterations, lower_count, midpoint_count, upper_count

        call shifted_negative_count(bracket_lower, lower_count)
        call shifted_negative_count(bracket_upper, upper_count)
        if (lower_count /= upper_count + 1) then
            write (message, "(a,i0,a,i0,a,i0,a)") &
                "negative bracket must contain exactly one eigenvalue (contains ", &
                lower_count - upper_count, "; lower count ", lower_count, &
                ", upper count ", upper_count, ")"
            call fail_usage(trim(message))
        end if
        iterations = 0
        do while (bracket_upper - bracket_lower > &
                bracket_tolerance * bracket_lower)
            if (iterations >= 80) call fail_solver("negative bracket iteration", -1)
            midpoint = bracket_lower + 0.5_dp * &
                (bracket_upper - bracket_lower)
            call shifted_negative_count(midpoint, midpoint_count)
            if (midpoint_count == lower_count) then
                bracket_lower = midpoint
            else if (midpoint_count == upper_count) then
                bracket_upper = midpoint
            else
                call fail_solver("negative bracket changed by multiple levels", -1)
            end if
            iterations = iterations + 1
        end do
        write (*, "(a)") "degree,parity,unknowns,normal_unknowns," // &
            "eta_unknowns,lower_count,upper_count,iterations," // &
            "lambda_lower,lambda_upper,relative_tolerance,m0_stored_power"
        write (*, "(i0,7(a,i0),4(a,es24.16))") degree, ",", parity, ",", &
            size(stiffness, 1), ",", problem%normal_unknowns, ",", &
            problem%eta_unknowns, ",", lower_count, ",", upper_count, ",", &
            iterations, ",", -bracket_upper, ",", -bracket_lower, ",", &
            bracket_tolerance, ",", m0_stored_power
    end subroutine bracket_negative_eigenvalue

    subroutine bracket_eigenvalue
        character(len=192) :: message
        real(dp) :: midpoint, scale
        integer :: iterations, lower_count, midpoint_count, upper_count

        call shifted_negative_count(-bracket_lower, lower_count)
        call shifted_negative_count(-bracket_upper, upper_count)
        if (upper_count /= lower_count + 1) then
            write (message, "(a,i0,a,i0,a,i0,a)") &
                "eigenvalue bracket must contain exactly one eigenvalue (contains ", &
                upper_count - lower_count, "; lower count ", lower_count, &
                ", upper count ", upper_count, ")"
            call fail_usage(trim(message))
        end if
        iterations = 0
        scale = max(abs(bracket_lower), abs(bracket_upper), tiny(1.0_dp))
        do while (bracket_upper - bracket_lower > bracket_tolerance * scale)
            if (iterations >= 80) &
                call fail_solver("eigenvalue bracket iteration", -1)
            midpoint = bracket_lower + 0.5_dp * &
                (bracket_upper - bracket_lower)
            call shifted_negative_count(-midpoint, midpoint_count)
            if (midpoint_count == lower_count) then
                bracket_lower = midpoint
            else if (midpoint_count == upper_count) then
                bracket_upper = midpoint
            else
                call fail_solver( &
                    "eigenvalue bracket changed by multiple levels", -1)
            end if
            iterations = iterations + 1
        end do
        write (*, "(a)") "degree,parity,unknowns,normal_unknowns," // &
            "eta_unknowns,lower_count,upper_count,iterations," // &
            "lambda_lower,lambda_upper,relative_tolerance,m0_stored_power"
        write (*, "(i0,7(a,i0),4(a,es24.16))") degree, ",", parity, ",", &
            size(stiffness, 1), ",", problem%normal_unknowns, ",", &
            problem%eta_unknowns, ",", lower_count, ",", upper_count, ",", &
            iterations, ",", bracket_lower, ",", bracket_upper, ",", &
            bracket_tolerance, ",", m0_stored_power
    end subroutine bracket_eigenvalue

end program gliss_compatible_marginality
