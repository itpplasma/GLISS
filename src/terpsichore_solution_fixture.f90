module terpsichore_solution_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: dynamic_family_layout_t, &
        eta_global_index, normal_global_index
    implicit none
    private

    integer, parameter, public :: terpsichore_solution_ok = 0
    integer, parameter, public :: terpsichore_solution_invalid = -1
    integer, parameter :: maximum_intervals = 996
    integer, parameter :: maximum_modes = 10000
    integer, parameter :: maximum_points = 1000000
    integer, parameter :: maximum_solution_values = 50000000
    integer, parameter :: maximum_mode_number = 1000000

    type, public :: terpsichore_solution_fixture_t
        integer :: plasma_intervals = 0
        integer :: vacuum_intervals = 0
        integer :: modes = 0
        integer :: poloidal_points = 0
        integer :: toroidal_points = 0
        real(dp) :: potential_energy = 0.0_dp
        real(dp) :: kinetic_energy = 0.0_dp
        integer, allocatable :: mode_m(:)
        integer, allocatable :: mode_n(:)
        real(dp), allocatable :: normal(:, :)
        real(dp), allocatable :: tangential(:, :)
    end type terpsichore_solution_fixture_t

    public :: build_terpsichore_plasma_solution
    public :: read_terpsichore_solution_fixture

contains

    subroutine read_terpsichore_solution_fixture(unit, vacuum_intervals, &
            fixture, info)
        integer, intent(in) :: unit, vacuum_intervals
        type(terpsichore_solution_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        real(dp), allocatable :: radial_power(:), compressional(:, :)
        integer :: allocation_status, field_periods, stability_periods
        integer :: total_intervals

        info = terpsichore_solution_invalid
        read (unit, iostat=allocation_status) fixture%plasma_intervals, &
            fixture%poloidal_points, fixture%toroidal_points, &
            stability_periods, field_periods, fixture%modes
        if (allocation_status /= 0) return
        fixture%vacuum_intervals = vacuum_intervals
        if (.not. valid_metadata(fixture, stability_periods, &
            field_periods)) return
        total_intervals = fixture%plasma_intervals + vacuum_intervals
        read (unit, iostat=allocation_status)
        if (allocation_status /= 0) return
        allocate (fixture%mode_m(fixture%modes), &
            fixture%mode_n(fixture%modes), radial_power(fixture%modes), &
            fixture%normal(fixture%modes, 0:total_intervals), &
            fixture%tangential(fixture%modes, total_intervals), &
            compressional(fixture%modes, total_intervals), &
            stat=allocation_status)
        if (allocation_status /= 0) return
        read (unit, iostat=allocation_status) fixture%mode_m, fixture%mode_n, &
            radial_power, fixture%normal, fixture%tangential, compressional
        if (allocation_status /= 0) return
        read (unit, iostat=allocation_status)
        if (allocation_status /= 0) return
        read (unit, iostat=allocation_status)
        if (allocation_status /= 0) return
        call read_energy_record(unit, fixture, allocation_status)
        if (allocation_status /= 0) return
        if (.not. fixture_is_valid(fixture)) return
        info = terpsichore_solution_ok
    end subroutine read_terpsichore_solution_fixture

    pure function valid_metadata(fixture, stability_periods, field_periods) &
            result(valid)
        type(terpsichore_solution_fixture_t), intent(in) :: fixture
        integer, intent(in) :: stability_periods, field_periods
        logical :: valid
        integer :: total_intervals

        valid = fixture%plasma_intervals >= 2 &
            .and. fixture%vacuum_intervals >= 0
        if (.not. valid) return
        total_intervals = fixture%plasma_intervals + fixture%vacuum_intervals
        valid = total_intervals <= maximum_intervals &
            .and. fixture%modes >= 1 .and. fixture%modes <= maximum_modes &
            .and. fixture%poloidal_points >= 1 &
            .and. fixture%toroidal_points >= 1 &
            .and. stability_periods >= 1 .and. field_periods >= 1
        if (.not. valid) return
        valid = fixture%poloidal_points <= maximum_points &
            / fixture%toroidal_points
        if (.not. valid) return
        valid = modulo(field_periods, stability_periods) == 0
        if (.not. valid) return
        valid = fixture%modes <= maximum_solution_values &
            / (2 * total_intervals + 1)
    end function valid_metadata

    subroutine read_energy_record(unit, fixture, io_status)
        integer, intent(in) :: unit
        type(terpsichore_solution_fixture_t), intent(inout) :: fixture
        integer, intent(out) :: io_status
        real(dp), allocatable :: cell_a(:), cell_b(:), mode_surface(:, :)
        real(dp), allocatable :: point_surface(:, :), surface_a(:)
        real(dp), allocatable :: surface_b(:)
        integer, allocatable :: equilibrium_m(:), equilibrium_n(:)
        integer :: equilibrium_modes, points

        read (unit, iostat=io_status) equilibrium_modes
        if (io_status /= 0) return
        if (equilibrium_modes < 1 .or. equilibrium_modes > maximum_modes) then
            io_status = 1
            return
        end if
        backspace (unit, iostat=io_status)
        if (io_status /= 0) return
        points = fixture%poloidal_points * fixture%toroidal_points
        allocate (equilibrium_m(equilibrium_modes), &
            equilibrium_n(equilibrium_modes), &
            mode_surface(equilibrium_modes, 0:fixture%plasma_intervals), &
            point_surface(points, 0:fixture%plasma_intervals), &
            surface_a(fixture%plasma_intervals + 1), &
            surface_b(fixture%plasma_intervals + 1), &
            cell_a(fixture%plasma_intervals), &
            cell_b(fixture%plasma_intervals), stat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) equilibrium_modes, equilibrium_m, &
            equilibrium_n, mode_surface, point_surface, surface_a, surface_b, &
            cell_a, cell_b, fixture%potential_energy, fixture%kinetic_energy
    end subroutine read_energy_record

    pure function fixture_is_valid(fixture) result(valid)
        type(terpsichore_solution_fixture_t), intent(in) :: fixture
        logical :: valid
        integer :: total_intervals

        valid = allocated(fixture%mode_m) .and. allocated(fixture%mode_n)
        if (.not. valid) return
        valid = allocated(fixture%normal) .and. &
            allocated(fixture%tangential)
        if (.not. valid) return
        total_intervals = fixture%plasma_intervals + fixture%vacuum_intervals
        valid = size(fixture%normal, 1) == fixture%modes &
            .and. size(fixture%normal, 2) == total_intervals + 1 &
            .and. size(fixture%tangential, 1) == fixture%modes &
            .and. size(fixture%tangential, 2) == total_intervals
        if (.not. valid) return
        valid = all(abs(real(fixture%mode_m, dp)) <= maximum_mode_number) &
            .and. all(abs(real(fixture%mode_n, dp)) <= maximum_mode_number)
        if (.not. valid) return
        valid = all(ieee_is_finite(fixture%normal)) &
            .and. all(ieee_is_finite(fixture%tangential)) &
            .and. ieee_is_finite(fixture%potential_energy) &
            .and. ieee_is_finite(fixture%kinetic_energy)
        if (.not. valid) return
        valid = fixture%kinetic_energy > 0.0_dp
    end function fixture_is_valid

    subroutine build_terpsichore_plasma_solution(fixture, layout, vector, info)
        type(terpsichore_solution_fixture_t), intent(in) :: fixture
        type(dynamic_family_layout_t), intent(in) :: layout
        real(dp), allocatable, intent(out) :: vector(:)
        integer, intent(out) :: info
        integer :: index, mode, radial

        info = terpsichore_solution_invalid
        if (.not. fixture_is_valid(fixture)) return
        if (layout%trials /= fixture%modes &
            .or. layout%intervals /= fixture%plasma_intervals) return
        allocate (vector(layout%total_unknowns), source=0.0_dp)
        do radial = 1, fixture%plasma_intervals
            do mode = 1, fixture%modes
                index = normal_global_index(layout, radial, mode)
                if (index > 0) vector(index) = fixture%normal(mode, radial)
                index = eta_global_index(layout, radial, mode)
                if (index > 0) vector(index) = fixture%tangential(mode, radial)
            end do
        end do
        if (.not. all(ieee_is_finite(vector))) return
        info = terpsichore_solution_ok
    end subroutine build_terpsichore_plasma_solution

end module terpsichore_solution_fixture
