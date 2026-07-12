program gliss_terpsichore_eigen
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
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
        iterate_variable_generalized_eigenvalue, variable_generalized_ok, &
        variable_generalized_inertia
    implicit none

    type(terpsichore_matrix_fixture_t) :: fixture
    type(dynamic_family_layout_t) :: layout, mass_layout
    type(variable_block_tridiagonal_t) :: k_blocks, m_blocks
    real(dp), allocatable :: stiffness(:, :), mass(:, :), vector(:)
    integer, allocatable :: permutation(:), widths(:)
    character(len=1024) :: filename, vacuum_argument
    real(dp) :: eigenvalue, residual, resolution, lower, upper, middle
    real(dp) :: certificate
    integer :: info, io_status, unit, vacuum_intervals, count, iteration
    integer :: negatives

    call get_command_argument(1, filename, status=io_status)
    if (io_status /= 0 .or. len_trim(filename) == 0) then
        write (error_unit, "(a)") &
            "usage: gliss_terpsichore_eigen FORT.23 IVAC"
        error stop 2
    end if
    call get_command_argument(2, vacuum_argument, status=io_status)
    if (io_status /= 0) error stop 2
    read (vacuum_argument, *, iostat=io_status) vacuum_intervals
    if (io_status /= 0 .or. vacuum_intervals /= 0) then
        write (error_unit, "(a)") "fixed-boundary input requires IVAC=0"
        error stop 2
    end if
    open (newunit=unit, file=trim(filename), status="old", action="read", &
        access="sequential", form="unformatted", iostat=io_status)
    if (io_status /= 0) then
        write (error_unit, "(a)") "cannot open TERPSICHORE matrix fixture"
        error stop 1
    end if
    call read_terpsichore_fixed_boundary_potential_fixture(unit, &
        vacuum_intervals, fixture, info)
    close (unit)
    if (info /= terpsichore_matrix_fixture_ok) then
        write (error_unit, "(a,i0)") "TERPSICHORE fixture status: ", info
        error stop 1
    end if
    call assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
        fixture, stiffness, layout, info)
    if (info /= terpsichore_noninteracting_ok) then
        write (error_unit, "(a,i0)") "stiffness assembly status: ", info
        error stop 1
    end if
    call assemble_terpsichore_fixture_reduced_mass(fixture, mass, &
        mass_layout, info)
    if (info /= terpsichore_reduced_adapter_ok) then
        write (error_unit, "(a,i0)") "mass assembly status: ", info
        error stop 1
    end if
    if (mass_layout%total_unknowns /= layout%total_unknowns) then
        write (error_unit, "(a)") "stiffness and mass layouts differ"
        error stop 1
    end if
    call build_dynamic_block_permutation(layout, widths, permutation, info)
    if (info /= dynamic_layout_ok) then
        write (error_unit, "(a,i0)") "permutation status: ", info
        error stop 1
    end if
    call pack_permuted_variable_blocks(stiffness, permutation, widths, &
        k_blocks, info)
    if (info /= variable_block_ok) then
        write (error_unit, "(a,i0)") "stiffness packing status: ", info
        error stop 1
    end if
    deallocate (stiffness)
    call pack_permuted_variable_blocks(mass, permutation, widths, &
        m_blocks, info)
    if (info /= variable_block_ok) then
        write (error_unit, "(a,i0)") "mass packing status: ", info
        error stop 1
    end if
    deallocate (mass)

    call variable_generalized_inertia(k_blocks, m_blocks, 0.0_dp, &
        negatives, info)
    if (info /= variable_generalized_ok) then
        write (error_unit, "(a,i0)") "inertia status: ", info
        error stop 1
    end if
    write (error_unit, "(a,i0)") "negative count at zero shift: ", &
        negatives

    ! certified most-negative eigenvalue by inertia bisection in the
    ! TERPSICHORE normalization, resolved to half the manifest's 1e-4
    ! relative certificate, then inverse iteration at the bracket
    ! bottom.
    if (negatives > 0) then
        lower = -1.0e-8_dp
        do iteration = 1, 200
            call variable_generalized_inertia(k_blocks, m_blocks, lower, &
                count, info)
            if (info /= variable_generalized_ok) then
                lower = lower * (1.0_dp + 1.0e-8_dp)
                cycle
            end if
            write (error_unit, "(a,es12.4,a,i0)") "probe ", lower, &
                " below: ", count
            if (count == 0) exit
            lower = 2.0_dp * lower
        end do
        upper = 0.0_dp
        do iteration = 1, 200
            middle = 0.5_dp * (lower + upper)
            if (upper - lower <= 5.0e-5_dp * abs(middle) &
                + 1.0e-16_dp) exit
            call variable_generalized_inertia(k_blocks, m_blocks, middle, &
                count, info)
            if (info /= variable_generalized_ok) then
                middle = middle * (1.0_dp + 1.0e-8_dp)
                cycle
            end if
            write (error_unit, "(a,es12.4,a,i0)") "bisect ", middle, &
                " below: ", count
            if (count == 0) then
                lower = middle
            else
                upper = middle
            end if
        end do
        certificate = upper - lower
        call iterate_variable_generalized_eigenvalue(k_blocks, m_blocks, &
            lower, eigenvalue, vector, residual, resolution, info)
    else
        write (error_unit, "(a)") "no negative eigenvalue at IVAC=0"
        error stop 1
    end if
    if (info /= variable_generalized_ok) then
        write (error_unit, "(a,i0)") "inverse iteration status: ", info
        error stop 1
    end if

    write (*, "(a)") "unknowns,negative_count,eigenvalue,certificate," // &
        "residual,resolution"
    write (*, "(i0,a,i0,a,es24.16,a,es9.2,a,es9.2,a,es9.2)") &
        layout%total_unknowns, ",", negatives, ",", eigenvalue, ",", &
        certificate, ",", residual, ",", resolution
end program gliss_terpsichore_eigen
