program gliss_compatible_operator_trace
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_operator_geometry_trace, only: &
        operator_geometry_trace_ok, write_compatible_operator_geometry
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_cell_trace_t, &
        compatible_problem_ok, compatible_two_component_problem_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        evaluate_radial_cubic_spline_field, fit_radial_cubic_spline_field, &
        radial_cubic_spline_field_t, radial_cubic_spline_grid_t, &
        radial_cubic_spline_ok
    implicit none

    character(len=16), parameter :: field_names(14) = [character(len=16) :: &
        "flux_t_slope", "flux_p_slope", "flux_t_curve", "flux_p_curve", &
        "current_i", "current_j", "signed_sqrtg", "mod_b", "grad_s2", &
        "j_dot_b", "pressure_slope", "sigma_tilde", "beta_tilde", "drive"]
    integer, parameter :: selected_cell_count = 8
    integer, parameter :: q1_selection_rank = 4
    real(dp), parameter :: target_poloidal_radii(7) = [ &
        0.10_dp, 0.25_dp, 0.40_dp, 0.55_dp, 0.70_dp, 0.85_dp, 0.95_dp]
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(compatible_two_component_problem_t) :: problem
    type(compatible_cell_trace_t), allocatable :: traces(:)
    type(radial_cubic_spline_grid_t) :: profile_grid
    type(radial_cubic_spline_field_t) :: profile_spline
    character(len=1024) :: filename
    integer, allocatable :: mode_m(:), mode_n(:), selected_cells(:)
    real(dp), allocatable :: profile_values(:, :), stored_power(:)
    real(dp) :: edge_profiles(3), edge_seconds(3), edge_slopes(3), q1_distance
    integer :: allocation_status, arguments, degree, info, n_theta, n_zeta
    integer :: parity, profile_index, q1_index

    interface
        subroutine terminate_process(status) bind(C, name="exit")
            import c_int
            integer(c_int), value :: status
        end subroutine terminate_process
    end interface

    call read_arguments(filename, degree, n_theta, n_zeta, parity, mode_m, &
        mode_n, stored_power)
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) call fail("reader", info)
    q1_index = 1
    q1_distance = abs(equilibrium%rotational_transform(1) - 1.0_dp)
    do profile_index = 2, size(equilibrium%rotational_transform)
        if (abs(equilibrium%rotational_transform(profile_index) - 1.0_dp) < &
            q1_distance) then
            q1_index = profile_index
            q1_distance = abs( &
                equilibrium%rotational_transform(profile_index) - 1.0_dp)
        end if
    end do
    allocate (profile_values(size(equilibrium%s), 3), stat=allocation_status)
    if (allocation_status /= 0) call fail("profile allocation", -1)
    profile_values(:, 1) = equilibrium%toroidal_flux
    profile_values(:, 2) = equilibrium%poloidal_flux
    profile_values(:, 3) = equilibrium%pressure
    call build_radial_cubic_spline_grid(equilibrium%s, 0.0_dp, 1.0_dp, &
        profile_grid, info)
    if (info /= radial_cubic_spline_ok) call fail("profile grid", info)
    call fit_radial_cubic_spline_field(profile_grid, profile_values, &
        profile_spline, info)
    if (info /= radial_cubic_spline_ok) call fail("profile spline", info)
    call evaluate_radial_cubic_spline_field(profile_grid, profile_spline, &
        1.0_dp, edge_profiles, edge_slopes, edge_seconds, info)
    if (info /= radial_cubic_spline_ok) call fail("edge profiles", info)
    if (edge_profiles(2) == 0.0_dp) call fail("zero edge poloidal flux", -1)
    call select_cells(size(equilibrium%s), q1_index, selected_cells)
    call build_compatible_two_component_problem(equilibrium, mode_m, mode_n, &
        stored_power, parity, degree, n_theta, n_zeta, problem, info, &
        selected_cells, traces)
    if (info /= compatible_problem_ok) call fail("operator assembly", info)
    call write_trace(q1_index)
    call write_compatible_operator_geometry(equilibrium, traces, n_theta, &
        n_zeta, info)
    if (info /= operator_geometry_trace_ok) call fail("geometry trace", info)

contains

    subroutine write_trace(local_q1_index)
        integer, intent(in) :: local_q1_index
        integer :: field, mode, surface, trace_index

        write (*, "(a)") "TRACE_VERSION,1"
        write (*, "(a,9(i0,:,','))") "DIMENSIONS,", &
            size(equilibrium%s), degree, n_theta, n_zeta, size(mode_m), &
            size(problem%stiffness, 1), problem%normal_unknowns, &
            problem%eta_unknowns, size(traces)
        do mode = 1, size(mode_m)
            write (*, "(a,i0,3(',',i0),',',es24.16e3)") "MODE,", mode, &
                mode_m(mode), mode_n(mode), parity, stored_power(mode)
        end do
        do field = 1, size(field_names)
            write (*, "(a,i0,',',a)") "FIELD_NAME,", field, &
                trim(field_names(field))
        end do
        do surface = 1, size(equilibrium%s)
            call write_profile(surface)
        end do
        write (*, "(a,2(i0,','),4(es24.16e3,:,','))") "Q1_SELECTION,", &
            local_q1_index, selected_position(local_q1_index), &
            equilibrium%s(local_q1_index), &
            poloidal_radius(equilibrium%poloidal_flux(local_q1_index)), &
            1.0_dp / equilibrium%rotational_transform(local_q1_index), &
            equilibrium%rotational_transform(local_q1_index)
        call write_angles
        do trace_index = 1, size(traces)
            call write_cell(trace_index)
        end do
    end subroutine write_trace

    subroutine write_profile(surface)
        integer, intent(in) :: surface
        write (*, "(a,i0,9(',',es24.16e3))") "PROFILE,", surface, &
            equilibrium%s(surface), &
            poloidal_radius(equilibrium%poloidal_flux(surface)), &
            1.0_dp / equilibrium%rotational_transform(surface), &
            equilibrium%pressure(surface), equilibrium%toroidal_flux(surface), &
            equilibrium%poloidal_flux(surface), &
            equilibrium%rotational_transform(surface), &
            equilibrium%b_theta_average(surface), &
            equilibrium%b_zeta_average(surface)
    end subroutine write_profile

    subroutine write_angles
        real(dp) :: theta, zeta
        integer :: j, k
        do k = 1, n_zeta
            zeta = real(k - 1, dp) / real(n_zeta, dp)
            do j = 1, n_theta
                theta = real(j - 1, dp) / real(n_theta, dp)
                write (*, "(a,2(i0,','),2(es24.16e3,:,','))") "ANGLE,", &
                    j, k, theta, zeta
            end do
        end do
    end subroutine write_angles

    subroutine write_cell(trace_index)
        integer, intent(in) :: trace_index
        real(dp), allocatable :: cell_mass(:, :), cell_terms(:, :, :)
        real(dp) :: lower, upper
        integer :: allocation_status, cell, point

        cell = traces(trace_index)%cell
        lower = real(cell - 1, dp) / real(size(equilibrium%s), dp)
        upper = real(cell, dp) / real(size(equilibrium%s), dp)
        write (*, "(a,2(i0,','),2(es24.16e3,:,','))") &
            "CELL_SELECTION,", trace_index, cell, lower, upper
        allocate (cell_mass(size(traces(trace_index)%points(1)%map), &
            size(traces(trace_index)%points(1)%map)), source=0.0_dp, &
            stat=allocation_status)
        if (allocation_status /= 0) call fail("cell mass allocation", -1)
        allocate (cell_terms(size(cell_mass, 1), size(cell_mass, 2), 4), &
            source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) call fail("cell term allocation", -1)
        do point = 1, size(traces(trace_index)%points)
            call write_point(trace_index, point, cell_mass, cell_terms)
        end do
        call write_matrix("CELL", cell, "B", cell_mass)
        call write_term_matrices("CELL", cell, cell_terms)
        call write_global_block(trace_index)
    end subroutine write_cell

    subroutine write_point(trace_index, point, cell_mass, cell_terms)
        integer, intent(in) :: trace_index, point
        real(dp), intent(inout) :: cell_mass(:, :), cell_terms(:, :, :)
        real(dp) :: profiles(3), seconds(3), slopes(3)
        character(len=10) :: kind
        integer :: info_local, term

        associate (sample => traces(trace_index)%points(point))
            call evaluate_radial_cubic_spline_field(profile_grid, &
                profile_spline, sample%coordinate, profiles, slopes, seconds, &
                info_local)
            if (info_local /= radial_cubic_spline_ok) &
                call fail("point profiles", info_local)
            kind = merge("accurate  ", "constraint", sample%assembles_mass)
            write (*, "(a,2(i0,','),a,',',11(es24.16e3,:,','))") "POINT,", &
                traces(trace_index)%cell, point, trim(kind), &
                sample%coordinate, poloidal_radius(profiles(2)), &
                slopes(1) / slopes(2), profiles(3), profiles(1), profiles(2), &
                slopes(1), slopes(2), seconds(1), seconds(2), sample%weight
            call write_basis(trace_index, point)
            call write_point_fields(trace_index, point)
            do term = 1, 4
                call write_matrix("POINT", point, term_name(term), &
                    sample%stiffness_terms(:, :, term), &
                    traces(trace_index)%cell)
                if (sample%term_mask(term)) cell_terms(:, :, term) = &
                    cell_terms(:, :, term) &
                    + sample%weight * sample%stiffness_terms(:, :, term)
            end do
            call write_matrix("POINT", point, "B", sample%mass, &
                traces(trace_index)%cell)
            if (sample%assembles_mass) cell_mass = cell_mass &
                + sample%weight * sample%mass
        end associate
    end subroutine write_point

    subroutine write_basis(trace_index, point)
        integer, intent(in) :: trace_index, point
        integer :: basis, column, trial, trials

        associate (sample => traces(trace_index)%points(point))
            trials = size(mode_m)
            do basis = 1, size(sample%h1, 1)
                do trial = 1, trials
                    column = (basis - 1) * trials + trial
                    write (*, "(a,2(i0,','),a,3(',',i0),3(',',es24.16e3))") &
                        "BASIS,", traces(trace_index)%cell, point, "H1", &
                        basis, trial, sample%map(column), sample%h1(basis, trial), &
                        sample%dh1(basis, trial), sample%weight
                end do
            end do
            do basis = 1, size(sample%l2, 1)
                do trial = 1, trials
                    column = size(sample%h1, 1) * trials &
                        + (basis - 1) * trials + trial
                    write (*, "(a,2(i0,','),a,3(',',i0),3(',',es24.16e3))") &
                        "BASIS,", traces(trace_index)%cell, point, "L2", &
                        basis, trial, sample%map(column), sample%l2(basis, trial), &
                        0.0_dp, sample%weight
                end do
            end do
        end associate
    end subroutine write_basis

    subroutine write_point_fields(trace_index, point)
        integer, intent(in) :: trace_index, point
        integer :: field, j, k

        associate (sample => traces(trace_index)%points(point))
            do k = 1, n_zeta
                do j = 1, n_theta
                    do field = 1, 13
                        write (*, "(a,5(i0,','),es24.16e3)") "FIELD,", &
                            traces(trace_index)%cell, point, j, k, field, &
                            sample%fields(j, k, field)
                    end do
                    write (*, "(a,5(i0,','),es24.16e3)") "FIELD,", &
                        traces(trace_index)%cell, point, j, k, 14, &
                        sample%drive(j, k)
                end do
            end do
        end associate
    end subroutine write_point_fields

    subroutine write_global_block(trace_index)
        integer, intent(in) :: trace_index
        integer, allocatable :: indices(:)
        integer :: a, allocation_status, b, cell, count, source, term

        cell = traces(trace_index)%cell
        count = 0
        do source = 1, size(traces(trace_index)%points(1)%map)
            if (traces(trace_index)%points(1)%map(source) > 0) count = count + 1
        end do
        allocate (indices(count), stat=allocation_status)
        if (allocation_status /= 0) call fail("global index allocation", -1)
        count = 0
        do source = 1, size(traces(trace_index)%points(1)%map)
            if (traces(trace_index)%points(1)%map(source) <= 0) cycle
            count = count + 1
            indices(count) = traces(trace_index)%points(1)%map(source)
        end do
        do a = 1, size(indices)
            write (*, "(a,3(i0,:,','))") "GLOBAL_INDEX,", cell, a, indices(a)
        end do
        do b = 1, size(indices)
            do a = 1, size(indices)
                write (*, "(a,i0,',',a,2(',',i0),',',es24.16e3)") &
                    "GLOBAL_MATRIX,", cell, "A", a, b, &
                    problem%stiffness(indices(a), indices(b))
                write (*, "(a,i0,',',a,2(',',i0),',',es24.16e3)") &
                    "GLOBAL_MATRIX,", cell, "B", a, b, &
                    problem%mass(indices(a), indices(b))
                do term = 1, 4
                    write (*, "(a,i0,',',a,2(',',i0),',',es24.16e3)") &
                        "GLOBAL_MATRIX,", cell, term_name(term), a, b, &
                        problem%stiffness_terms(indices(a), indices(b), term)
                end do
            end do
        end do
    end subroutine write_global_block

    subroutine write_term_matrices(scope, index, matrices)
        character(len=*), intent(in) :: scope
        integer, intent(in) :: index
        real(dp), intent(in) :: matrices(:, :, :)
        integer :: term
        do term = 1, size(matrices, 3)
            call write_matrix(scope, index, term_name(term), &
                matrices(:, :, term))
        end do
    end subroutine write_term_matrices

    subroutine write_matrix(scope, index, name, matrix, cell)
        character(len=*), intent(in) :: scope, name
        integer, intent(in) :: index
        real(dp), intent(in) :: matrix(:, :)
        integer, optional, intent(in) :: cell
        integer :: a, b, local_cell

        local_cell = index
        if (present(cell)) local_cell = cell
        do b = 1, size(matrix, 2)
            do a = 1, size(matrix, 1)
                write (*, "(a,a,',',i0,',',i0,',',a,2(',',i0),',',es24.16e3)") &
                    "MATRIX,", trim(scope), local_cell, index, trim(name), &
                    a, b, matrix(a, b)
            end do
        end do
    end subroutine write_matrix

    pure function term_name(term) result(name)
        integer, intent(in) :: term
        character(len=2) :: name
        select case (term)
        case (1)
            name = "K1"
        case (2)
            name = "K2"
        case (3)
            name = "K3"
        case default
            name = "K4"
        end select
    end function term_name

    pure function poloidal_radius(poloidal_flux) result(radius)
        real(dp), intent(in) :: poloidal_flux
        real(dp) :: radius
        radius = sqrt(poloidal_flux / edge_profiles(2))
    end function poloidal_radius

    pure function selected_position(cell) result(position)
        integer, intent(in) :: cell
        integer :: position
        position = findloc(selected_cells, cell, dim=1)
    end function selected_position

    subroutine select_cells(intervals, q1_cell, cells)
        integer, intent(in) :: intervals, q1_cell
        integer, allocatable, intent(out) :: cells(:)
        integer :: allocation_status, cell, previous, rank, target
        real(dp) :: distance, nearest_distance
        logical :: used

        if (intervals < selected_cell_count) &
            call fail("operator trace requires at least eight radial cells", intervals)
        if (q1_cell < 1 .or. q1_cell > intervals) &
            call fail("q=1 trace cell is out of bounds", q1_cell)
        allocate (cells(selected_cell_count), stat=allocation_status)
        if (allocation_status /= 0) call fail("cell selection allocation", -1)
        cells = 0
        cells(q1_selection_rank) = q1_cell
        target = 0
        do rank = 1, selected_cell_count
            if (rank == q1_selection_rank) cycle
            target = target + 1
            nearest_distance = huge(0.0_dp)
            do cell = 1, intervals
                used = .false.
                do previous = 1, selected_cell_count
                    if (cells(previous) == cell) used = .true.
                end do
                if (used) cycle
                distance = abs(poloidal_radius( &
                    equilibrium%poloidal_flux(cell)) - &
                    target_poloidal_radii(target))
                if (distance < nearest_distance) then
                    nearest_distance = distance
                    cells(rank) = cell
                end if
            end do
            if (cells(rank) == 0) call fail("radial trace-cell selection", rank)
        end do
        do rank = 2, selected_cell_count
            if (cells(rank) <= cells(rank - 1)) call fail( &
                "q=1 radius must lie between trace targets 0.40 and 0.55", &
                q1_cell)
        end do
    end subroutine select_cells

    subroutine read_arguments(path, local_degree, poloidal_points, &
            toroidal_points, local_parity, poloidal_modes, toroidal_modes, &
            powers)
        character(len=*), intent(out) :: path
        integer, intent(out) :: local_degree, poloidal_points, toroidal_points
        integer, intent(out) :: local_parity
        integer, allocatable, intent(out) :: poloidal_modes(:), toroidal_modes(:)
        real(dp), allocatable, intent(out) :: powers(:)
        character(len=64) :: token
        integer :: allocation_status, argument, comma, local_info, mode

        arguments = command_argument_count()
        if (arguments < 6) call usage("missing required arguments")
        call read_argument(1, "EXPORT_FILE", path)
        call read_integer(2, "DEGREE", local_degree)
        call read_integer(3, "NTHETA", poloidal_points)
        call read_integer(4, "NZETA", toroidal_points)
        call read_integer(5, "PARITY", local_parity)
        if (local_degree < 1 .or. local_degree > 4) &
            call usage("DEGREE must be between 1 and 4")
        if (poloidal_points < 8 .or. toroidal_points < 8) &
            call usage("NTHETA and NZETA must be at least 8")
        if (local_parity < 1 .or. local_parity > 2) &
            call usage("PARITY must be 1 or 2")
        allocate (poloidal_modes(arguments - 5), &
            toroidal_modes(arguments - 5), powers(arguments - 5), &
            stat=allocation_status)
        if (allocation_status /= 0) call fail("mode allocation", -1)
        do argument = 6, arguments
            mode = argument - 5
            call read_argument(argument, "mode", token)
            comma = index(token, ",")
            if (comma <= 1) call usage("modes must be given once as m,n")
            if (comma == len_trim(token)) &
                call usage("modes must be given once as m,n")
            if (index(token(comma + 1:), ",") > 0) &
                call usage("modes must be given once as m,n")
            call parse_integer(token(:comma - 1), "poloidal mode", &
                poloidal_modes(mode), local_info)
            call parse_integer(token(comma + 1:), "toroidal mode", &
                toroidal_modes(mode), local_info)
            if (poloidal_modes(mode) < 0) &
                call usage("poloidal mode must be nonnegative")
            if (poloidal_modes(mode) == 0 .and. toroidal_modes(mode) < 0) &
                call usage("axis modes require nonnegative n")
            if (mode > 1) then
                if (any(poloidal_modes(:mode - 1) == poloidal_modes(mode) &
                    .and. toroidal_modes(:mode - 1) == toroidal_modes(mode))) &
                    call usage("duplicate mode")
            end if
            powers(mode) = 0.0_dp
            if (poloidal_modes(mode) > 0) powers(mode) = &
                1.0_dp - 0.5_dp * real(poloidal_modes(mode), dp)
        end do
    end subroutine read_arguments

    subroutine read_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        integer :: status
        call get_command_argument(position, value, status=status)
        if (status /= 0) call usage(trim(name) // " is too long")
        if (len_trim(value) == 0) call usage(trim(name) // " is empty")
    end subroutine read_argument

    subroutine read_integer(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        integer, intent(out) :: value
        character(len=64) :: token
        integer :: status
        call read_argument(position, name, token)
        call parse_integer(token, name, value, status)
    end subroutine read_integer

    subroutine parse_integer(text, name, value, status)
        character(len=*), intent(in) :: text, name
        integer, intent(out) :: value, status
        read (text, *, iostat=status) value
        if (status /= 0) call usage(trim(name) // " must be an integer")
    end subroutine parse_integer

    subroutine usage(message)
        character(len=*), intent(in) :: message
        write (error_unit, "(a)") "gliss_compatible_operator_trace: " &
            // trim(message)
        write (error_unit, "(a)") "usage: gliss_compatible_operator_trace " &
            // "EXPORT_FILE DEGREE NTHETA NZETA PARITY m,n [m,n ...]"
        call terminate_process(2_c_int)
    end subroutine usage

    subroutine fail(operation, status)
        character(len=*), intent(in) :: operation
        integer, intent(in) :: status
        write (error_unit, "(a,i0)") &
            "gliss_compatible_operator_trace: " // trim(operation) &
            // " error ", status
        error stop 1
    end subroutine fail

end program gliss_compatible_operator_trace
