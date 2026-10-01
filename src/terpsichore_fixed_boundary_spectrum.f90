module terpsichore_fixed_boundary_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_block_scatter, only: dynamic_block_map_t
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        dynamic_family_layout_t, dynamic_layout_ok
    use terpsichore_eigen_diagnostics, only: &
        compute_terpsichore_eigen_diagnostics, &
        terpsichore_eigen_diagnostics_ok, terpsichore_eigen_diagnostics_t
    use terpsichore_matrix_fixture, only: &
        read_terpsichore_fixed_boundary_potential_fixture, &
        terpsichore_matrix_fixture_ok, terpsichore_matrix_fixture_t
    use terpsichore_noninteracting_stiffness, only: &
        assemble_terpsichore_noninteracting_fixed_boundary_blocks, &
        terpsichore_noninteracting_ok
    use terpsichore_reduced_mass_adapter, only: &
        assemble_terpsichore_fixture_reduced_mass_blocks, &
        terpsichore_reduced_adapter_ok
    use terpsichore_solution_fixture, only: &
        build_terpsichore_plasma_solution, &
        read_terpsichore_solution_fixture, terpsichore_solution_ok, &
        terpsichore_solution_fixture_t
    use variable_block_tridiagonal, only: pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, pencil_roundoff, &
        variable_generalized_inertia, variable_generalized_ok
    implicit none
    private

    integer, parameter, public :: terpsichore_fixed_spectrum_ok = 0
    integer, parameter, public :: terpsichore_fixed_spectrum_read_error = 1
    integer, parameter, public :: terpsichore_fixed_spectrum_compute_error = 2

    type, public :: terpsichore_fixed_boundary_result_t
        integer :: unknowns = 0
        integer :: negative_count = 0
        real(dp) :: eigenvalue = 0.0_dp
        real(dp) :: certificate = 0.0_dp
        real(dp) :: residual = 0.0_dp
        real(dp) :: resolution = 0.0_dp
        real(dp) :: reference_eigenvalue = 0.0_dp
        real(dp) :: reference_potential = 0.0_dp
        real(dp) :: computed_potential = 0.0_dp
        real(dp) :: reference_kinetic = 0.0_dp
        real(dp) :: computed_kinetic = 0.0_dp
        real(dp) :: reference_residual = 0.0_dp
        real(dp) :: mode_overlap = 0.0_dp
    end type terpsichore_fixed_boundary_result_t

    public :: pack_terpsichore_problem
    public :: read_terpsichore_reference
    public :: solve_terpsichore_fixed_boundary_file
    public :: solve_terpsichore_lowest
    public :: terpsichore_layouts_match

contains

    subroutine solve_terpsichore_fixed_boundary_file(path, result, info, &
            message)
        character(len=*), intent(in) :: path
        type(terpsichore_fixed_boundary_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(terpsichore_matrix_fixture_t) :: fixture
        type(dynamic_family_layout_t) :: layout
        type(dynamic_block_map_t) :: map
        type(variable_block_tridiagonal_t) :: stiffness_blocks, mass_blocks
        type(terpsichore_eigen_diagnostics_t) :: diagnostics
        type(terpsichore_solution_fixture_t) :: solution
        real(dp), allocatable :: vector(:), reference(:)

        result = terpsichore_fixed_boundary_result_t()
        call read_fixed_fixture(path, fixture, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        call assemble_fixed_problem(fixture, stiffness_blocks, mass_blocks, &
            layout, map, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        result%unknowns = layout%total_unknowns
        call solve_terpsichore_lowest(stiffness_blocks, mass_blocks, &
            result%eigenvalue, vector, result%residual, result%resolution, &
            result%certificate, result%negative_count, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        call read_terpsichore_reference(path, 0, fixture, layout, &
            map%permutation, solution, reference, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        call compute_terpsichore_eigen_diagnostics(stiffness_blocks, mass_blocks, &
            result%eigenvalue, vector, reference, solution%potential_energy, &
            solution%kinetic_energy, 1.0_dp, diagnostics, info)
        if (info /= terpsichore_eigen_diagnostics_ok) then
            call solve_failure("TERPSICHORE eigen diagnostics failed", info, &
                message)
            return
        end if
        result%reference_eigenvalue = diagnostics%reference_quotient
        result%reference_potential = solution%potential_energy
        result%computed_potential = diagnostics%computed_potential
        result%reference_kinetic = solution%kinetic_energy
        result%computed_kinetic = diagnostics%computed_kinetic
        result%reference_residual = diagnostics%reference_residual
        result%mode_overlap = diagnostics%mode_overlap
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine solve_terpsichore_fixed_boundary_file

    subroutine read_terpsichore_reference(path, vacuum_intervals, matrix, &
            layout, permutation, fixture, reference, info, message)
        character(len=*), intent(in) :: path
        integer, intent(in) :: vacuum_intervals, permutation(:)
        type(terpsichore_matrix_fixture_t), intent(in) :: matrix
        type(dynamic_family_layout_t), intent(in) :: layout
        type(terpsichore_solution_fixture_t), intent(out) :: fixture
        real(dp), allocatable, intent(out) :: reference(:)
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        real(dp), allocatable :: original(:)
        integer :: allocation_status, close_status, i, read_status, unit

        info = terpsichore_fixed_spectrum_read_error
        message = "cannot open TERPSICHORE FORT.23 reference"
        open (newunit=unit, file=trim(path), status="old", action="read", &
            access="sequential", form="unformatted", iostat=read_status)
        if (read_status /= 0) return
        call read_terpsichore_solution_fixture(unit, vacuum_intervals, fixture, &
            read_status)
        close (unit, iostat=close_status)
        if (read_status /= terpsichore_solution_ok) then
            message = "invalid TERPSICHORE FORT.23 reference"
            return
        end if
        if (close_status /= 0) then
            message = "cannot close TERPSICHORE FORT.23 reference"
            return
        end if
        if (.not. solution_modes_match(matrix, fixture)) then
            message = "TERPSICHORE matrix and solution mode tables differ"
            return
        end if
        call build_terpsichore_plasma_solution(fixture, layout, original, &
            read_status)
        if (read_status /= terpsichore_solution_ok) then
            info = terpsichore_fixed_spectrum_compute_error
            message = "TERPSICHORE solution mapping failed"
            return
        end if
        allocate (reference(size(original)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = terpsichore_fixed_spectrum_compute_error
            message = "TERPSICHORE reference allocation failed"
            return
        end if
        do i = 1, size(original)
            reference(i) = original(permutation(i))
        end do
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine read_terpsichore_reference

    pure function solution_modes_match(matrix, solution) result(matches)
        type(terpsichore_matrix_fixture_t), intent(in) :: matrix
        type(terpsichore_solution_fixture_t), intent(in) :: solution
        logical :: matches

        matches = matrix%intervals == solution%plasma_intervals
        if (.not. matches) return
        matches = matrix%modes == solution%modes
        if (.not. matches) return
        matches = allocated(matrix%mode_m)
        if (.not. matches) return
        matches = allocated(matrix%mode_n)
        if (.not. matches) return
        matches = allocated(solution%mode_m)
        if (.not. matches) return
        matches = allocated(solution%mode_n)
        if (.not. matches) return
        matches = all(matrix%mode_m == solution%mode_m)
        if (.not. matches) return
        matches = all(matrix%mode_n == solution%mode_n)
    end function solution_modes_match

    subroutine read_fixed_fixture(path, fixture, info, message)
        character(len=*), intent(in) :: path
        type(terpsichore_matrix_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer :: io_status, unit

        info = terpsichore_fixed_spectrum_read_error
        message = "cannot open TERPSICHORE FORT.23"
        open (newunit=unit, file=trim(path), status="old", action="read", &
            access="sequential", form="unformatted", iostat=io_status)
        if (io_status /= 0) return
        if (has_vacuum_records(unit)) then
            close (unit, iostat=io_status)
            info = terpsichore_fixed_spectrum_read_error
            message = "TERPSICHORE FORT.23 has vacuum intervals (IVAC>0); " &
                // "use the pseudoplasma solver"
            return
        end if
        rewind (unit)
        call read_terpsichore_fixed_boundary_potential_fixture(unit, 0, &
            fixture, io_status)
        close (unit, iostat=info)
        if (io_status /= terpsichore_matrix_fixture_ok) then
            info = terpsichore_fixed_spectrum_read_error
            message = "invalid TERPSICHORE FORT.23"
            return
        end if
        if (info /= 0) then
            info = terpsichore_fixed_spectrum_read_error
            message = "cannot close TERPSICHORE FORT.23"
            return
        end if
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine read_fixed_fixture

    function has_vacuum_records(unit) result(vacuum)
        ! The radial record of an IVAC=0 file holds exactly the plasma grid
        ! and profiles. With IVAC>0 the grid extends into the vacuum, so the
        ! record is longer and the fixed-boundary layout would misread every
        ! later field (the parity first).
        integer, intent(in) :: unit
        logical :: vacuum
        integer :: intervals, poloidal, toroidal, periods, field, modes
        integer :: io_status
        real(dp), allocatable :: probe(:)

        vacuum = .false.
        read (unit, iostat=io_status) intervals, poloidal, toroidal, periods, &
            field, modes
        if (io_status /= 0 .or. intervals < 1) return
        allocate (probe(7 * intervals + 5), stat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) probe
        vacuum = io_status == 0
    end function has_vacuum_records

    ! Stiffness and mass in block-tridiagonal storage, assembled interval
    ! by interval: no dense matrix of the problem order is formed.
    subroutine assemble_fixed_problem(fixture, stiffness, mass, layout, map, &
            info, message)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        type(variable_block_tridiagonal_t), intent(out) :: stiffness, mass
        type(dynamic_family_layout_t), intent(out) :: layout
        type(dynamic_block_map_t), intent(out) :: map
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(dynamic_family_layout_t) :: mass_layout
        type(dynamic_block_map_t) :: mass_map

        info = terpsichore_fixed_spectrum_compute_error
        if (fixture%legacy_modelk /= 0) then
            message = "TERPSICHORE fixed-boundary solve requires MODELK=0"
            return
        end if
        if (fixture%parity /= 0.0_dp) then
            message = "TERPSICHORE fixed-boundary solve requires sine parity"
            return
        end if
        call assemble_terpsichore_noninteracting_fixed_boundary_blocks( &
            fixture, stiffness, layout, map, info)
        if (info /= terpsichore_noninteracting_ok) then
            message = "TERPSICHORE stiffness assembly failed"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        call assemble_terpsichore_fixture_reduced_mass_blocks(fixture, mass, &
            mass_layout, mass_map, info)
        if (info /= terpsichore_reduced_adapter_ok) then
            message = "TERPSICHORE reduced mass assembly failed"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        if (.not. terpsichore_layouts_match(layout, mass_layout) &
            .or. any(map%permutation /= mass_map%permutation)) then
            message = "TERPSICHORE stiffness and mass layouts differ"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine assemble_fixed_problem

    subroutine pack_terpsichore_problem(layout, stiffness, mass, &
            stiffness_blocks, mass_blocks, widths, permutation, info, message)
        type(dynamic_family_layout_t), intent(in) :: layout
        real(dp), allocatable, intent(inout) :: stiffness(:, :), mass(:, :)
        type(variable_block_tridiagonal_t), intent(out) :: stiffness_blocks
        type(variable_block_tridiagonal_t), intent(out) :: mass_blocks
        integer, allocatable, intent(out) :: widths(:), permutation(:)
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        call build_dynamic_block_permutation(layout, widths, permutation, info)
        if (info /= dynamic_layout_ok) then
            message = "TERPSICHORE block permutation failed"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        call pack_permuted_variable_blocks(stiffness, permutation, widths, &
            stiffness_blocks, info)
        if (info /= variable_block_ok) then
            message = "TERPSICHORE stiffness packing failed"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        deallocate (stiffness)
        call pack_permuted_variable_blocks(mass, permutation, widths, &
            mass_blocks, info)
        if (info /= variable_block_ok) then
            message = "TERPSICHORE mass packing failed"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        deallocate (mass)
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine pack_terpsichore_problem

    subroutine solve_terpsichore_lowest(stiffness, mass, eigenvalue, &
            eigenvector, residual, resolution, certificate, negative_count, &
            info, message)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(out) :: eigenvalue, residual, resolution, certificate
        real(dp), allocatable, intent(out) :: eigenvector(:)
        integer, intent(out) :: negative_count, info
        character(len=*), intent(out) :: message
        real(dp) :: zero_floor, lower, upper, middle
        integer :: below, iteration

        ! Null directions within the pencil roundoff are not negative.
        zero_floor = 64.0_dp * pencil_roundoff(stiffness, mass)
        ! This first factorization validates the pencil for every later
        ! probe of solve_terpsichore_lowest and its brackets.
        call variable_generalized_inertia(stiffness, mass, -zero_floor, &
            negative_count, info)
        if (info /= variable_generalized_ok) then
            call solve_failure("TERPSICHORE inertia failed", info, message)
            return
        end if
        ! A stable pencil has no negative eigenvalue; its certified lowest
        ! eigenpair is then the lowest one at or above -zero_floor.
        if (negative_count > 0) then
            call bracket_lowest(stiffness, mass, zero_floor, lower, upper, &
                info, message)
        else
            call bracket_lowest_stable(stiffness, mass, zero_floor, lower, &
                upper, info, message)
        end if
        if (info /= terpsichore_fixed_spectrum_ok) return
        do iteration = 1, 200
            middle = lower + 0.5_dp * (upper - lower)
            if (upper - lower <= 5.0e-5_dp * abs(middle) + zero_floor) exit
            call variable_generalized_inertia(stiffness, mass, middle, below, &
                info, validated=.true.)
            if (info /= variable_generalized_ok) then
                middle = nearest(middle, 1.0_dp)
                call variable_generalized_inertia(stiffness, mass, middle, &
                    below, info, validated=.true.)
            end if
            if (info /= variable_generalized_ok) then
                call solve_failure("TERPSICHORE bisection inertia failed", &
                    info, message)
                return
            end if
            if (below == 0) then
                lower = middle
            else
                upper = middle
            end if
        end do
        if (iteration > 200) then
            call solve_failure("TERPSICHORE eigenvalue bisection failed", &
                info, message)
            return
        end if
        certificate = upper - lower
        call iterate_variable_generalized_eigenvalue(stiffness, mass, lower, &
            eigenvalue, eigenvector, residual, resolution, info, &
            validated=.true.)
        if (info /= variable_generalized_ok) then
            call solve_failure("TERPSICHORE inverse iteration failed", info, &
                message)
            return
        end if
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine solve_terpsichore_lowest

    subroutine bracket_lowest(stiffness, mass, zero_floor, lower, upper, &
            info, message)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: zero_floor
        real(dp), intent(out) :: lower, upper
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer :: below, iteration

        lower = -2.0_dp * zero_floor
        do iteration = 1, 200
            call variable_generalized_inertia(stiffness, mass, lower, below, &
                info, validated=.true.)
            if (info == variable_generalized_ok .and. below == 0) exit
            if (info /= variable_generalized_ok) &
                lower = nearest(lower, -1.0_dp)
            lower = 2.0_dp * lower
        end do
        if (iteration > 200) then
            call solve_failure("TERPSICHORE cannot bracket lowest eigenvalue", &
                info, message)
            return
        end if
        upper = -zero_floor
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine bracket_lowest

    subroutine bracket_lowest_stable(stiffness, mass, zero_floor, lower, &
            upper, info, message)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: zero_floor
        real(dp), intent(out) :: lower, upper
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer :: below, iteration

        lower = -zero_floor
        upper = 2.0_dp * zero_floor
        do iteration = 1, 1100
            call variable_generalized_inertia(stiffness, mass, upper, below, &
                info, validated=.true.)
            if (info == variable_generalized_ok .and. below > 0) exit
            if (info /= variable_generalized_ok) &
                upper = nearest(upper, 1.0_dp)
            upper = 2.0_dp * upper
            if (.not. ieee_is_finite(upper)) exit
        end do
        if (iteration > 1100 .or. .not. ieee_is_finite(upper)) then
            call solve_failure("TERPSICHORE cannot bracket lowest eigenvalue", &
                info, message)
            return
        end if
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine bracket_lowest_stable

    subroutine solve_failure(text, info, message)
        character(len=*), intent(in) :: text
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        info = terpsichore_fixed_spectrum_compute_error
        message = text
    end subroutine solve_failure

    pure function terpsichore_layouts_match(first, second) result(matches)
        type(dynamic_family_layout_t), intent(in) :: first, second
        logical :: matches

        matches = first%trials == second%trials &
            .and. first%intervals == second%intervals &
            .and. first%total_unknowns == second%total_unknowns &
            .and. first%outer_normal_retained .eqv. &
            second%outer_normal_retained
        if (.not. matches) return
        matches = allocated(first%active) .and. allocated(second%active)
        if (.not. matches) return
        matches = all(shape(first%active) == shape(second%active))
        if (.not. matches) return
        matches = all(first%active .eqv. second%active)
    end function terpsichore_layouts_match

end module terpsichore_fixed_boundary_spectrum
