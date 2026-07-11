program gliss_terpsichore_reduced_mass
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: dynamic_family_layout_t
    use terpsichore_matrix_fixture, only: &
        read_terpsichore_fixed_boundary_fixture, &
        terpsichore_matrix_fixture_ok, terpsichore_matrix_fixture_t
    use terpsichore_reduced_mass_adapter, only: &
        assemble_terpsichore_fixture_reduced_mass, &
        terpsichore_reduced_adapter_ok
    implicit none

    type(terpsichore_matrix_fixture_t) :: fixture
    type(dynamic_family_layout_t) :: layout
    real(dp), allocatable :: mass(:, :)
    character(len=1024) :: filename, vacuum_argument
    integer :: column, info, io_status, row, unit, vacuum_intervals

    call get_command_argument(1, filename, status=io_status)
    if (io_status /= 0 .or. len_trim(filename) == 0) then
        write (error_unit, "(a)") &
            "usage: gliss_terpsichore_reduced_mass FORT.23 IVAC"
        error stop 2
    end if
    call get_command_argument(2, vacuum_argument, status=io_status)
    if (io_status /= 0) error stop 2
    read (vacuum_argument, *, iostat=io_status) vacuum_intervals
    if (io_status /= 0 .or. vacuum_intervals /= 0) then
        write (error_unit, "(a)") "reduced fixed-boundary input requires IVAC=0"
        error stop 2
    end if
    open (newunit=unit, file=trim(filename), status="old", action="read", &
        access="sequential", form="unformatted", iostat=io_status)
    if (io_status /= 0) then
        write (error_unit, "(a)") "cannot open TERPSICHORE matrix fixture"
        error stop 1
    end if
    call read_terpsichore_fixed_boundary_fixture(unit, vacuum_intervals, &
        fixture, info)
    close (unit)
    if (info /= terpsichore_matrix_fixture_ok) then
        write (error_unit, "(a,i0)") "TERPSICHORE fixture status: ", info
        error stop 1
    end if
    call assemble_terpsichore_fixture_reduced_mass(fixture, mass, layout, info)
    if (info /= terpsichore_reduced_adapter_ok) then
        write (error_unit, "(a,i0)") "TERPSICHORE adapter status: ", info
        error stop 1
    end if

    write (*, "(a)") "row,column,mass"
    do column = 1, layout%total_unknowns
        do row = 1, layout%total_unknowns
            write (*, "(i0,a,i0,a,es24.16e3)") row, ",", column, ",", &
                mass(row, column)
        end do
    end do
end program gliss_terpsichore_reduced_mass
