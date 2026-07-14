module terpsichore_fixed_boundary_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        dynamic_family_layout_t, dynamic_layout_ok
    use terpsichore_matrix_fixture, only: &
        read_terpsichore_fixed_boundary_potential_fixture, &
        terpsichore_matrix_fixture_ok, terpsichore_matrix_fixture_t
    use terpsichore_noninteracting_stiffness, only: &
        assemble_terpsichore_noninteracting_fixed_boundary_stiffness, &
        terpsichore_noninteracting_ok
    use terpsichore_reduced_mass_adapter, only: &
        assemble_terpsichore_fixture_reduced_mass, &
        terpsichore_reduced_adapter_ok
    use variable_block_tridiagonal, only: pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_generalized_inertia, &
        variable_generalized_ok
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
    end type terpsichore_fixed_boundary_result_t

    public :: pack_terpsichore_problem
    public :: solve_terpsichore_fixed_boundary_file
    public :: solve_terpsichore_lowest_negative
    public :: terpsichore_layouts_match

contains

    subroutine solve_terpsichore_fixed_boundary_file(path, result, info, &
            message)
        character(len=*), intent(in) :: path
        type(terpsichore_fixed_boundary_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(terpsichore_matrix_fixture_t) :: fixture
        type(dynamic_family_layout_t) :: stiffness_layout, mass_layout
        type(variable_block_tridiagonal_t) :: stiffness_blocks, mass_blocks
        real(dp), allocatable :: stiffness(:, :), mass(:, :), vector(:)
        integer, allocatable :: permutation(:), widths(:)

        result = terpsichore_fixed_boundary_result_t()
        call read_fixed_fixture(path, fixture, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        call assemble_fixed_problem(fixture, stiffness, mass, &
            stiffness_layout, mass_layout, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        result%unknowns = stiffness_layout%total_unknowns
        call pack_terpsichore_problem(stiffness_layout, stiffness, mass, &
            stiffness_blocks, mass_blocks, widths, permutation, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        call solve_terpsichore_lowest_negative(stiffness_blocks, mass_blocks, &
            result%eigenvalue, vector, result%residual, result%resolution, &
            result%certificate, result%negative_count, info, message)
    end subroutine solve_terpsichore_fixed_boundary_file

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

    subroutine assemble_fixed_problem(fixture, stiffness, mass, &
            stiffness_layout, mass_layout, info, message)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: stiffness(:, :), mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: stiffness_layout
        type(dynamic_family_layout_t), intent(out) :: mass_layout
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        info = terpsichore_fixed_spectrum_compute_error
        if (fixture%legacy_modelk /= 0) then
            message = "TERPSICHORE fixed-boundary solve requires MODELK=0"
            return
        end if
        if (fixture%parity /= 0.0_dp) then
            message = "TERPSICHORE fixed-boundary solve requires sine parity"
            return
        end if
        call assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
            fixture, stiffness, stiffness_layout, info)
        if (info /= terpsichore_noninteracting_ok) then
            message = "TERPSICHORE stiffness assembly failed"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        call assemble_terpsichore_fixture_reduced_mass(fixture, mass, &
            mass_layout, info)
        if (info /= terpsichore_reduced_adapter_ok) then
            message = "TERPSICHORE reduced mass assembly failed"
            info = terpsichore_fixed_spectrum_compute_error
            return
        end if
        if (.not. terpsichore_layouts_match(stiffness_layout, mass_layout)) then
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

    subroutine solve_terpsichore_lowest_negative(stiffness, mass, eigenvalue, &
            eigenvector, residual, resolution, certificate, negative_count, &
            info, message)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(out) :: eigenvalue, residual, resolution, certificate
        real(dp), allocatable, intent(out) :: eigenvector(:)
        integer, intent(out) :: negative_count, info
        character(len=*), intent(out) :: message
        real(dp) :: lower, upper, middle
        integer :: below, iteration

        call variable_generalized_inertia(stiffness, mass, 0.0_dp, &
            negative_count, info)
        if (info /= variable_generalized_ok) then
            call solve_failure("TERPSICHORE inertia failed", info, message)
            return
        end if
        if (negative_count <= 0) then
            call solve_failure("TERPSICHORE matrix has no negative eigenvalue", &
                info, message)
            return
        end if
        call bracket_lowest(stiffness, mass, lower, upper, info, message)
        if (info /= terpsichore_fixed_spectrum_ok) return
        do iteration = 1, 200
            middle = lower + 0.5_dp * (upper - lower)
            if (upper - lower <= 5.0e-5_dp * abs(middle) + 1.0e-16_dp) exit
            call variable_generalized_inertia(stiffness, mass, middle, below, &
                info)
            if (info /= variable_generalized_ok) then
                middle = nearest(middle, 1.0_dp)
                call variable_generalized_inertia(stiffness, mass, middle, &
                    below, info)
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
            eigenvalue, eigenvector, residual, resolution, info)
        if (info /= variable_generalized_ok) then
            call solve_failure("TERPSICHORE inverse iteration failed", info, &
                message)
            return
        end if
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine solve_terpsichore_lowest_negative

    subroutine bracket_lowest(stiffness, mass, lower, upper, info, message)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(out) :: lower, upper
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer :: below, iteration

        lower = -1.0e-8_dp
        do iteration = 1, 200
            call variable_generalized_inertia(stiffness, mass, lower, below, &
                info)
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
        upper = 0.0_dp
        info = terpsichore_fixed_spectrum_ok
        message = ""
    end subroutine bracket_lowest

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
