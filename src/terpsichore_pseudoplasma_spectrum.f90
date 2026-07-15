module terpsichore_pseudoplasma_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: dynamic_family_layout_t
    use terpsichore_eigen_diagnostics, only: &
        compute_terpsichore_eigen_diagnostics, &
        terpsichore_eigen_diagnostics_ok, terpsichore_eigen_diagnostics_t
    use terpsichore_fixed_boundary_spectrum, only: &
        pack_terpsichore_problem, solve_terpsichore_lowest_negative, &
        terpsichore_fixed_spectrum_compute_error, &
        terpsichore_fixed_spectrum_ok, terpsichore_fixed_spectrum_read_error, &
        terpsichore_layouts_match
    use terpsichore_matrix_fixture, only: &
        read_terpsichore_potential_fixture, terpsichore_matrix_fixture_ok, &
        terpsichore_matrix_fixture_t
    use terpsichore_noninteracting_stiffness, only: &
        assemble_terpsichore_noninteracting_free_boundary_stiffness, &
        terpsichore_noninteracting_ok
    use terpsichore_pseudoplasma_coupling, only: &
        add_terpsichore_pseudoplasma_schur, pseudoplasma_coupling_ok
    use terpsichore_pseudoplasma_fixture, only: &
        pseudoplasma_fixture_ok, read_terpsichore_pseudoplasma_fixture, &
        terpsichore_pseudoplasma_fixture_t
    use terpsichore_reduced_mass_adapter, only: &
        assemble_terpsichore_fixture_reduced_mass_free_boundary, &
        terpsichore_reduced_adapter_ok
    use terpsichore_solution_fixture, only: &
        build_terpsichore_plasma_solution, &
        read_terpsichore_solution_fixture, terpsichore_solution_ok, &
        terpsichore_solution_fixture_t
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: terpsichore_pseudoplasma_spectrum_ok = &
        terpsichore_fixed_spectrum_ok
    integer, parameter, public :: terpsichore_pseudoplasma_spectrum_read_error = &
        terpsichore_fixed_spectrum_read_error
    integer, parameter, public :: &
        terpsichore_pseudoplasma_spectrum_compute_error = &
        terpsichore_fixed_spectrum_compute_error

    type, public :: terpsichore_pseudoplasma_result_t
        integer :: unknowns = 0
        integer :: negative_count = 0
        real(dp) :: eigenvalue = 0.0_dp
        real(dp) :: certificate = 0.0_dp
        real(dp) :: residual = 0.0_dp
        real(dp) :: resolution = 0.0_dp
        real(dp) :: growth_rate = 0.0_dp
        real(dp) :: reference_eigenvalue = 0.0_dp
        real(dp) :: reference_potential = 0.0_dp
        real(dp) :: computed_potential = 0.0_dp
        real(dp) :: reference_kinetic = 0.0_dp
        real(dp) :: computed_kinetic = 0.0_dp
        real(dp) :: reference_residual = 0.0_dp
        real(dp) :: mode_overlap = 0.0_dp
    end type terpsichore_pseudoplasma_result_t

    public :: solve_terpsichore_pseudoplasma_files

contains

    subroutine solve_terpsichore_pseudoplasma_files(matrix_path, &
            vacuum_intervals, vacuum_path, result, info, message)
        character(len=*), intent(in) :: matrix_path, vacuum_path
        integer, intent(in) :: vacuum_intervals
        type(terpsichore_pseudoplasma_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(terpsichore_matrix_fixture_t) :: fixture
        type(terpsichore_pseudoplasma_fixture_t) :: vacuum
        type(terpsichore_solution_fixture_t) :: solution
        type(dynamic_family_layout_t) :: stiffness_layout, mass_layout
        type(variable_block_tridiagonal_t) :: stiffness_blocks, mass_blocks
        type(terpsichore_eigen_diagnostics_t) :: diagnostics
        real(dp), allocatable :: stiffness(:, :), mass(:, :), vector(:)
        real(dp), allocatable :: reference(:)
        integer, allocatable :: permutation(:), widths(:)

        result = terpsichore_pseudoplasma_result_t()
        info = terpsichore_pseudoplasma_spectrum_compute_error
        message = "TERPSICHORE IVAC must be positive"
        if (vacuum_intervals <= 0) return
        call read_matrix(matrix_path, vacuum_intervals, fixture, info, message)
        if (info /= terpsichore_pseudoplasma_spectrum_ok) return
        call assemble_problem(fixture, vacuum_intervals, vacuum_path, &
            stiffness, mass, stiffness_layout, mass_layout, vacuum, info, &
            message)
        if (info /= terpsichore_pseudoplasma_spectrum_ok) return
        result%unknowns = stiffness_layout%total_unknowns
        call pack_terpsichore_problem(stiffness_layout, stiffness, mass, &
            stiffness_blocks, mass_blocks, widths, permutation, info, message)
        if (info /= terpsichore_pseudoplasma_spectrum_ok) return
        call solve_terpsichore_lowest_negative(stiffness_blocks, mass_blocks, &
            result%eigenvalue, vector, result%residual, result%resolution, &
            result%certificate, result%negative_count, info, message)
        if (info /= terpsichore_pseudoplasma_spectrum_ok) return
        call read_reference(matrix_path, vacuum_intervals, fixture, &
            stiffness_layout, permutation, solution, reference, info, message)
        if (info /= terpsichore_pseudoplasma_spectrum_ok) return
        call compute_terpsichore_eigen_diagnostics(stiffness_blocks, &
            mass_blocks, result%eigenvalue, vector, reference, &
            solution%potential_energy, solution%kinetic_energy, &
            vacuum%alfven_normalization, diagnostics, info)
        if (info /= terpsichore_eigen_diagnostics_ok) then
            info = terpsichore_pseudoplasma_spectrum_compute_error
            message = "TERPSICHORE eigen diagnostics failed"
            return
        end if
        result%growth_rate = diagnostics%growth_rate
        result%reference_eigenvalue = diagnostics%reference_quotient
        result%reference_potential = solution%potential_energy
        result%computed_potential = diagnostics%computed_potential
        result%reference_kinetic = solution%kinetic_energy
        result%computed_kinetic = diagnostics%computed_kinetic
        result%reference_residual = diagnostics%reference_residual
        result%mode_overlap = diagnostics%mode_overlap
        info = terpsichore_pseudoplasma_spectrum_ok
        message = ""
    end subroutine solve_terpsichore_pseudoplasma_files

    subroutine read_matrix(path, vacuum_intervals, fixture, info, message)
        character(len=*), intent(in) :: path
        integer, intent(in) :: vacuum_intervals
        type(terpsichore_matrix_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer :: close_status, read_status, unit

        info = terpsichore_pseudoplasma_spectrum_read_error
        message = "cannot open TERPSICHORE FORT.23"
        open (newunit=unit, file=trim(path), status="old", action="read", &
            access="sequential", form="unformatted", iostat=read_status)
        if (read_status /= 0) return
        call read_terpsichore_potential_fixture(unit, vacuum_intervals, &
            fixture, read_status)
        close (unit, iostat=close_status)
        if (read_status /= terpsichore_matrix_fixture_ok) then
            message = "invalid TERPSICHORE FORT.23"
            return
        end if
        if (close_status /= 0) then
            message = "cannot close TERPSICHORE FORT.23"
            return
        end if
        info = terpsichore_pseudoplasma_spectrum_ok
        message = ""
    end subroutine read_matrix

    subroutine assemble_problem(fixture, vacuum_intervals, vacuum_path, &
            stiffness, mass, stiffness_layout, mass_layout, vacuum, info, &
            message)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        integer, intent(in) :: vacuum_intervals
        character(len=*), intent(in) :: vacuum_path
        real(dp), allocatable, intent(out) :: stiffness(:, :), mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: stiffness_layout
        type(dynamic_family_layout_t), intent(out) :: mass_layout
        type(terpsichore_pseudoplasma_fixture_t), intent(out) :: vacuum
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        real(dp), allocatable :: effective(:, :), response(:, :)

        info = terpsichore_pseudoplasma_spectrum_compute_error
        if (fixture%legacy_modelk /= 0) then
            message = "TERPSICHORE pseudoplasma solve requires MODELK=0"
            return
        end if
        if (fixture%parity /= 0.0_dp) then
            message = "TERPSICHORE pseudoplasma solve requires sine parity"
            return
        end if
        call assemble_terpsichore_noninteracting_free_boundary_stiffness( &
            fixture, stiffness, stiffness_layout, info)
        if (info /= terpsichore_noninteracting_ok) then
            info = terpsichore_pseudoplasma_spectrum_compute_error
            message = "TERPSICHORE pseudoplasma stiffness assembly failed"
            return
        end if
        call read_vacuum(vacuum_path, vacuum, info, message)
        if (info /= terpsichore_pseudoplasma_spectrum_ok) return
        if (vacuum%vacuum_intervals /= vacuum_intervals) then
            info = terpsichore_pseudoplasma_spectrum_read_error
            message = "TERPSICHORE FORT.24 IVAC does not match the request"
            return
        end if
        call add_terpsichore_pseudoplasma_schur(fixture, vacuum, &
            stiffness_layout, stiffness, effective, response, info)
        if (info /= pseudoplasma_coupling_ok) then
            info = terpsichore_pseudoplasma_spectrum_compute_error
            message = "TERPSICHORE pseudoplasma Schur coupling failed"
            return
        end if
        call assemble_terpsichore_fixture_reduced_mass_free_boundary( &
            fixture, mass, mass_layout, info)
        if (info /= terpsichore_reduced_adapter_ok) then
            info = terpsichore_pseudoplasma_spectrum_compute_error
            message = "TERPSICHORE free-boundary reduced mass assembly failed"
            return
        end if
        if (.not. terpsichore_layouts_match(stiffness_layout, mass_layout)) then
            info = terpsichore_pseudoplasma_spectrum_compute_error
            message = "TERPSICHORE stiffness and mass layouts differ"
            return
        end if
        info = terpsichore_pseudoplasma_spectrum_ok
        message = ""
    end subroutine assemble_problem

    subroutine read_vacuum(path, fixture, info, message)
        character(len=*), intent(in) :: path
        type(terpsichore_pseudoplasma_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer :: close_status, read_status, unit

        info = terpsichore_pseudoplasma_spectrum_read_error
        message = "cannot open TERPSICHORE FORT.24"
        open (newunit=unit, file=trim(path), status="old", action="read", &
            access="sequential", form="unformatted", iostat=read_status)
        if (read_status /= 0) return
        call read_terpsichore_pseudoplasma_fixture(unit, fixture, read_status)
        close (unit, iostat=close_status)
        if (read_status /= pseudoplasma_fixture_ok) then
            message = "invalid TERPSICHORE FORT.24"
            return
        end if
        if (close_status /= 0) then
            message = "cannot close TERPSICHORE FORT.24"
            return
        end if
        info = terpsichore_pseudoplasma_spectrum_ok
        message = ""
    end subroutine read_vacuum

    subroutine read_reference(path, vacuum_intervals, matrix, layout, &
            permutation, fixture, reference, info, message)
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

        info = terpsichore_pseudoplasma_spectrum_read_error
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
            info = terpsichore_pseudoplasma_spectrum_compute_error
            message = "TERPSICHORE solution mapping failed"
            return
        end if
        allocate (reference(size(original)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = terpsichore_pseudoplasma_spectrum_compute_error
            message = "TERPSICHORE reference allocation failed"
            return
        end if
        do i = 1, size(original)
            reference(i) = original(permutation(i))
        end do
        info = terpsichore_pseudoplasma_spectrum_ok
        message = ""
    end subroutine read_reference

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

end module terpsichore_pseudoplasma_spectrum
