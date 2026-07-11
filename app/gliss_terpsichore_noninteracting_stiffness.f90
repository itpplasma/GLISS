program gliss_terpsichore_noninteracting_stiffness
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: dynamic_family_layout_t
    use terpsichore_matrix_fixture, only: &
        read_terpsichore_fixed_boundary_potential_fixture, &
        terpsichore_matrix_fixture_ok, terpsichore_matrix_fixture_t
    use terpsichore_noninteracting_stiffness, only: &
        assemble_terpsichore_noninteracting_fixed_boundary_stiffness, &
        terpsichore_noninteracting_ok
    implicit none

    type(terpsichore_matrix_fixture_t) :: fixture
    type(dynamic_family_layout_t) :: layout
    real(dp), allocatable :: stiffness(:, :)
    character(len=1024) :: filename, vacuum_argument
    integer :: column, info, io_status, row, unit, vacuum_intervals

    call get_command_argument(1, filename, status=io_status)
    if (io_status /= 0 .or. len_trim(filename) == 0) then
        write (error_unit, "(a)") &
            "usage: gliss_terpsichore_noninteracting_stiffness FORT.23 IVAC"
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
        write (error_unit, "(a,i0)") "TERPSICHORE stiffness status: ", info
        error stop 1
    end if

    write (*, "(a)") "row,column,stiffness"
    do column = 1, layout%total_unknowns
        do row = 1, layout%total_unknowns
            write (*, "(i0,a,i0,a,es24.16e3)") row, ",", column, ",", &
                stiffness(row, column)
        end do
    end do
end program gliss_terpsichore_noninteracting_stiffness
