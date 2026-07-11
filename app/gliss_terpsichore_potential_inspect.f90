program gliss_terpsichore_potential_inspect
    use, intrinsic :: iso_fortran_env, only: error_unit
    use terpsichore_matrix_fixture, only: &
        read_terpsichore_fixed_boundary_potential_fixture, &
        terpsichore_matrix_fixture_ok, terpsichore_matrix_fixture_t
    use terpsichore_model_policy, only: decode_terpsichore_model, &
        terpsichore_model_config_t, terpsichore_model_ok
    implicit none

    type(terpsichore_matrix_fixture_t) :: fixture
    type(terpsichore_model_config_t) :: model
    character(len=1024) :: filename, vacuum_argument
    integer :: info, io_status, unit, vacuum_intervals

    call get_command_argument(1, filename, status=io_status)
    call get_command_argument(2, vacuum_argument, status=info)
    if (io_status /= 0 .or. info /= 0 .or. len_trim(filename) == 0) then
        write (error_unit, "(a)") &
            "usage: gliss_terpsichore_potential_inspect FORT.23 IVAC"
        error stop 2
    end if
    read (vacuum_argument, *, iostat=io_status) vacuum_intervals
    if (io_status /= 0 .or. vacuum_intervals /= 0) then
        write (error_unit, "(a)") "potential fixture requires IVAC=0"
        error stop 2
    end if
    open (newunit=unit, file=trim(filename), status="old", action="read", &
        access="sequential", form="unformatted", iostat=io_status)
    if (io_status /= 0) then
        write (error_unit, "(a)") "cannot open TERPSICHORE potential fixture"
        error stop 1
    end if
    call read_terpsichore_fixed_boundary_potential_fixture(unit, &
        vacuum_intervals, fixture, info)
    close (unit)
    if (info /= terpsichore_matrix_fixture_ok) then
        write (error_unit, "(a,i0)") "TERPSICHORE fixture status: ", info
        error stop 1
    end if
    call decode_terpsichore_model(fixture%legacy_modelk, model, info)
    if (info /= terpsichore_model_ok) error stop 1

    write (*, "(a,i0)") "intervals=", fixture%intervals
    write (*, "(a,i0)") "angular_points=", &
        fixture%poloidal_points * fixture%toroidal_points
    write (*, "(a,i0)") "modes=", fixture%modes
    write (*, "(a,i0)") "legacy_modelk=", fixture%legacy_modelk
    write (*, "(a,i0)") "potential_model=", model%potential_model
    write (*, "(a,i0)") "kinetic_norm=", model%kinetic_norm
end program gliss_terpsichore_potential_inspect
