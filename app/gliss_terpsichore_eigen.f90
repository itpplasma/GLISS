program gliss_terpsichore_eigen
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: dynamic_family_layout_t
    use terpsichore_eigen_diagnostics, only: &
        compute_terpsichore_eigen_diagnostics, &
        terpsichore_eigen_diagnostics_ok, terpsichore_eigen_diagnostics_t
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
    use terpsichore_fixed_boundary_spectrum, only: &
        pack_terpsichore_problem, solve_terpsichore_fixed_boundary_file, &
        solve_terpsichore_lowest_negative, &
        terpsichore_fixed_boundary_result_t, &
        terpsichore_fixed_spectrum_ok, terpsichore_layouts_match
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    implicit none

    type(terpsichore_matrix_fixture_t) :: fixture
    type(terpsichore_pseudoplasma_fixture_t) :: vacuum
    type(terpsichore_solution_fixture_t) :: solution
    type(dynamic_family_layout_t) :: layout, mass_layout
    type(variable_block_tridiagonal_t) :: k_blocks, m_blocks
    type(terpsichore_eigen_diagnostics_t) :: diagnostics
    type(terpsichore_fixed_boundary_result_t) :: fixed_result
    real(dp), allocatable :: stiffness(:, :), mass(:, :), vector(:)
    real(dp), allocatable :: reference(:)
    integer, allocatable :: permutation(:), widths(:)
    character(len=1024) :: matrix_file, vacuum_file
    character(len=128) :: message
    real(dp) :: eigenvalue, residual, resolution, certificate
    integer :: info, vacuum_intervals, negatives

    call parse_arguments(matrix_file, vacuum_intervals, vacuum_file)
    if (vacuum_intervals == 0) then
        call solve_terpsichore_fixed_boundary_file(matrix_file, fixed_result, &
            info, message)
        if (info /= terpsichore_fixed_spectrum_ok) call fail(message, 1)
        write (error_unit, "(a,i0)") "negative count at zero shift: ", &
            fixed_result%negative_count
        call write_fixed_result(fixed_result%unknowns, &
            fixed_result%negative_count, fixed_result%eigenvalue, &
            fixed_result%certificate, fixed_result%residual, &
            fixed_result%resolution)
    else
        call read_matrix(matrix_file, vacuum_intervals, fixture)
        call assemble_problem(fixture, vacuum_intervals, vacuum_file, &
            stiffness, mass, layout, mass_layout, vacuum)
        if (.not. terpsichore_layouts_match(layout, mass_layout)) &
            call fail("stiffness and mass layouts differ", 1)
        call pack_terpsichore_problem(layout, stiffness, mass, k_blocks, &
            m_blocks, widths, permutation, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) call fail(message, 1)
        call solve_terpsichore_lowest_negative(k_blocks, m_blocks, eigenvalue, &
            vector, residual, resolution, certificate, negatives, info, &
            message)
        if (info /= terpsichore_fixed_spectrum_ok) call fail(message, 1)
        write (error_unit, "(a,i0)") &
            "negative count at zero shift: ", negatives
        call read_reference(matrix_file, vacuum_intervals, fixture, layout, &
            permutation, solution, reference)
        call compute_terpsichore_eigen_diagnostics(k_blocks, m_blocks, &
            eigenvalue, vector, reference, solution%potential_energy, &
            solution%kinetic_energy, vacuum%alfven_normalization, &
            diagnostics, info)
        if (info /= terpsichore_eigen_diagnostics_ok) &
            call fail_status("eigen diagnostics", info)
        call write_free_result(layout%total_unknowns, negatives, eigenvalue, &
            certificate, residual, resolution, solution, diagnostics)
    end if

contains

    subroutine parse_arguments(matrix_path, vacuum_count, vacuum_path)
        character(len=*), intent(out) :: matrix_path, vacuum_path
        integer, intent(out) :: vacuum_count
        character(len=1024) :: argument
        integer :: count, status

        matrix_path = ""
        vacuum_path = ""
        count = command_argument_count()
        if (count < 2 .or. count > 3) call usage()
        call get_command_argument(1, matrix_path, status=status)
        if (status /= 0 .or. len_trim(matrix_path) == 0) call usage()
        call get_command_argument(2, argument, status=status)
        if (status /= 0 .or. .not. decimal_integer(argument)) &
            call fail("IVAC must be a nonnegative decimal integer", 2)
        read (argument, *, iostat=status) vacuum_count
        if (status /= 0 .or. vacuum_count < 0) &
            call fail("IVAC is outside the supported integer range", 2)
        if (vacuum_count == 0 .and. count /= 2) &
            call fail("IVAC=0 does not accept FORT.24", 2)
        if (vacuum_count > 0 .and. count /= 3) &
            call fail("IVAC>0 requires FORT.24", 2)
        if (vacuum_count == 0) return
        call get_command_argument(3, vacuum_path, status=status)
        if (status /= 0 .or. len_trim(vacuum_path) == 0) &
            call fail("IVAC>0 requires FORT.24", 2)
    end subroutine parse_arguments

    pure function decimal_integer(argument) result(valid)
        character(len=*), intent(in) :: argument
        logical :: valid
        integer :: i, digit

        valid = len_trim(argument) > 0
        if (.not. valid) return
        do i = 1, len_trim(argument)
            digit = iachar(argument(i:i))
            valid = digit >= iachar("0") .and. digit <= iachar("9")
            if (.not. valid) return
        end do
    end function decimal_integer

    subroutine read_matrix(path, vacuum_count, value)
        character(len=*), intent(in) :: path
        integer, intent(in) :: vacuum_count
        type(terpsichore_matrix_fixture_t), intent(out) :: value
        integer :: info, io_status, unit

        call open_fixture(path, unit)
        call read_terpsichore_potential_fixture(unit, vacuum_count, value, &
            info)
        close (unit, iostat=io_status)
        if (info /= terpsichore_matrix_fixture_ok) &
            call fail_status("TERPSICHORE matrix fixture", info)
        if (io_status /= 0) call fail("cannot close FORT.23", 1)
    end subroutine read_matrix

    subroutine assemble_problem(value, vacuum_count, vacuum_path, k, m, &
            k_layout, m_layout, vacuum_value)
        type(terpsichore_matrix_fixture_t), intent(in) :: value
        integer, intent(in) :: vacuum_count
        character(len=*), intent(in) :: vacuum_path
        real(dp), allocatable, intent(out) :: k(:, :), m(:, :)
        type(dynamic_family_layout_t), intent(out) :: k_layout, m_layout
        type(terpsichore_pseudoplasma_fixture_t), intent(out) :: vacuum_value
        real(dp), allocatable :: effective(:, :), response(:, :)
        integer :: info

        call assemble_terpsichore_noninteracting_free_boundary_stiffness( &
            value, k, k_layout, info)
        if (info /= terpsichore_noninteracting_ok) &
            call fail_status("stiffness assembly", info)
        call read_vacuum(vacuum_path, vacuum_value)
        if (vacuum_value%vacuum_intervals /= vacuum_count) &
            call fail("FORT.24 IVAC does not match the command line", 1)
        call add_terpsichore_pseudoplasma_schur(value, vacuum_value, &
            k_layout, k, effective, response, info)
        if (info /= pseudoplasma_coupling_ok) &
            call fail_status("pseudoplasma Schur coupling", info)
        call assemble_terpsichore_fixture_reduced_mass_free_boundary( &
            value, m, m_layout, info)
        if (info /= terpsichore_reduced_adapter_ok) &
            call fail_status("mass assembly", info)
    end subroutine assemble_problem

    subroutine read_vacuum(path, value)
        character(len=*), intent(in) :: path
        type(terpsichore_pseudoplasma_fixture_t), intent(out) :: value
        integer :: info, io_status, unit

        call open_fixture(path, unit)
        call read_terpsichore_pseudoplasma_fixture(unit, value, info)
        close (unit, iostat=io_status)
        if (info /= pseudoplasma_fixture_ok) &
            call fail_status("TERPSICHORE pseudoplasma fixture", info)
        if (io_status /= 0) call fail("cannot close FORT.24", 1)
    end subroutine read_vacuum

    subroutine read_reference(path, vacuum_count, matrix, value_layout, perm, &
            value, permuted)
        character(len=*), intent(in) :: path
        integer, intent(in) :: vacuum_count, perm(:)
        type(terpsichore_matrix_fixture_t), intent(in) :: matrix
        type(dynamic_family_layout_t), intent(in) :: value_layout
        type(terpsichore_solution_fixture_t), intent(out) :: value
        real(dp), allocatable, intent(out) :: permuted(:)
        real(dp), allocatable :: original(:)
        integer :: i, info, io_status, unit

        call open_fixture(path, unit)
        call read_terpsichore_solution_fixture(unit, vacuum_count, value, info)
        close (unit, iostat=io_status)
        if (info /= terpsichore_solution_ok) &
            call fail_status("TERPSICHORE solution fixture", info)
        if (io_status /= 0) call fail("cannot close FORT.23", 1)
        if (.not. solution_modes_match(matrix, value)) &
            call fail("matrix and solution mode tables differ", 1)
        call build_terpsichore_plasma_solution(value, value_layout, original, &
            info)
        if (info /= terpsichore_solution_ok) &
            call fail_status("TERPSICHORE solution mapping", info)
        allocate (permuted(size(original)))
        do i = 1, size(original)
            permuted(i) = original(perm(i))
        end do
    end subroutine read_reference

    pure function solution_modes_match(matrix, solution_value) result(matches)
        type(terpsichore_matrix_fixture_t), intent(in) :: matrix
        type(terpsichore_solution_fixture_t), intent(in) :: solution_value
        logical :: matches

        matches = matrix%intervals == solution_value%plasma_intervals &
            .and. matrix%modes == solution_value%modes
        if (.not. matches) return
        matches = allocated(matrix%mode_m) .and. allocated(matrix%mode_n)
        if (.not. matches) return
        matches = allocated(solution_value%mode_m) &
            .and. allocated(solution_value%mode_n)
        if (.not. matches) return
        matches = all(matrix%mode_m == solution_value%mode_m) &
            .and. all(matrix%mode_n == solution_value%mode_n)
    end function solution_modes_match

    subroutine open_fixture(path, unit)
        character(len=*), intent(in) :: path
        integer, intent(out) :: unit
        integer :: io_status

        open (newunit=unit, file=trim(path), status="old", action="read", &
            access="sequential", form="unformatted", iostat=io_status)
        if (io_status /= 0) call fail("cannot open " // trim(path), 1)
    end subroutine open_fixture

    subroutine write_fixed_result(unknowns, negative_count, value, width, &
            value_residual, value_resolution)
        integer, intent(in) :: unknowns, negative_count
        real(dp), intent(in) :: value, width, value_residual, value_resolution

        write (*, "(a)") "unknowns,negative_count,eigenvalue,certificate," // &
            "residual,resolution"
        write (*, '(i0,",",i0,4(",",es24.16))') unknowns, &
            negative_count, value, width, value_residual, value_resolution
    end subroutine write_fixed_result

    subroutine write_free_result(unknowns, negative_count, value, width, &
            value_residual, value_resolution, reference_value, diagnostic)
        integer, intent(in) :: unknowns, negative_count
        real(dp), intent(in) :: value, width, value_residual, value_resolution
        type(terpsichore_solution_fixture_t), intent(in) :: reference_value
        type(terpsichore_eigen_diagnostics_t), intent(in) :: diagnostic

        write (*, "(a)") "unknowns,negative_count,eigenvalue,certificate," // &
            "residual,resolution,growth_rate,reference_eigenvalue," // &
            "reference_potential,computed_potential,reference_kinetic," // &
            "computed_kinetic,reference_residual,mode_overlap"
        write (*, '(i0,",",i0,12(",",es24.16))') unknowns, &
            negative_count, value, width, value_residual, value_resolution, &
            diagnostic%growth_rate, diagnostic%reference_quotient, &
            reference_value%potential_energy, diagnostic%computed_potential, &
            reference_value%kinetic_energy, diagnostic%computed_kinetic, &
            diagnostic%reference_residual, diagnostic%mode_overlap
    end subroutine write_free_result

    subroutine usage()
        write (error_unit, "(a)") &
            "usage: gliss_terpsichore_eigen FORT.23 IVAC [FORT.24]"
        flush (error_unit)
        stop 2
    end subroutine usage

    subroutine fail_status(operation, status)
        character(len=*), intent(in) :: operation
        integer, intent(in) :: status

        write (error_unit, "(a,a,i0)") trim(operation), " status: ", status
        flush (error_unit)
        stop 1
    end subroutine fail_status

    subroutine fail(message, status)
        character(len=*), intent(in) :: message
        integer, intent(in) :: status

        write (error_unit, "(a)") trim(message)
        flush (error_unit)
        if (status == 2) stop 2
        stop 1
    end subroutine fail

end program gliss_terpsichore_eigen
