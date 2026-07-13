program gliss_terpsichore_pseudoplasma_stiffness
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use terpsichore_pseudoplasma_fixture, only: &
        pseudoplasma_fixture_ok, read_terpsichore_pseudoplasma_fixture, &
        terpsichore_pseudoplasma_fixture_t
    use terpsichore_pseudoplasma_stiffness, only: &
        assemble_terpsichore_pseudoplasma_stiffness, &
        pseudoplasma_stiffness_ok
    implicit none

    type(terpsichore_pseudoplasma_fixture_t) :: fixture
    real(dp), allocatable :: stiffness(:, :)
    character(len=1024) :: filename
    integer :: column, info, io_status, row, unit

    if (command_argument_count() /= 1) call fail_usage()
    call get_command_argument(1, filename, status=io_status)
    if (io_status /= 0 .or. len_trim(filename) == 0) call fail_usage()
    open (newunit=unit, file=trim(filename), status="old", action="read", &
        access="sequential", form="unformatted", iostat=io_status)
    if (io_status /= 0) call fail("cannot open pseudoplasma fixture")
    call read_terpsichore_pseudoplasma_fixture(unit, fixture, info)
    close (unit)
    if (info /= pseudoplasma_fixture_ok) call fail("invalid pseudoplasma fixture")
    call assemble_terpsichore_pseudoplasma_stiffness(fixture, stiffness, info)
    if (info /= pseudoplasma_stiffness_ok) &
        call fail("pseudoplasma stiffness assembly failed")
    write (*, "(a)") "row,column,stiffness"
    do column = 1, size(stiffness, 2)
        do row = 1, size(stiffness, 1)
            write (*, "(i0,',',i0,',',es24.16)") row, column, &
                stiffness(row, column)
        end do
    end do

contains

    subroutine fail_usage()
        write (error_unit, "(a)") &
            "usage: gliss_terpsichore_pseudoplasma_stiffness FORT.24"
        flush (error_unit)
        stop 2
    end subroutine fail_usage

    subroutine fail(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") &
            "gliss_terpsichore_pseudoplasma_stiffness: " // message
        flush (error_unit)
        stop 1
    end subroutine fail

end program gliss_terpsichore_pseudoplasma_stiffness
