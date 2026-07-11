module terpsichore_matrix_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: terpsichore_matrix_fixture_ok = 0
    integer, parameter, public :: terpsichore_matrix_fixture_invalid = -1
    integer, parameter :: maximum_intervals = 996
    integer, parameter :: maximum_angular_points = 999
    integer, parameter :: maximum_modes = 10000
    integer, parameter :: maximum_fixture_values = 50000000
    integer, parameter :: maximum_phase_values = 10000000
    integer, parameter :: maximum_dense_order = 4000

    type, public :: terpsichore_matrix_fixture_t
        integer :: intervals = 0
        integer :: poloidal_points = 0
        integer :: toroidal_points = 0
        integer :: stability_periods = 0
        integer :: field_periods = 0
        integer :: modes = 0
        real(dp) :: parity = 0.0_dp
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: flux_t_slope(:)
        real(dp), allocatable :: radial_power(:)
        real(dp), allocatable :: signed_bjac(:, :)
        integer, allocatable :: mode_m(:)
        integer, allocatable :: mode_n(:)
    end type terpsichore_matrix_fixture_t

    public :: read_terpsichore_fixed_boundary_fixture
    public :: terpsichore_fixed_fixture_is_valid

contains

    subroutine read_terpsichore_fixed_boundary_fixture(unit, vacuum_intervals, &
            fixture, info)
        integer, intent(in) :: unit, vacuum_intervals
        type(terpsichore_matrix_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        real(dp), allocatable :: bjac(:, :), flux_t_slope(:), radial_power(:)
        real(dp), allocatable :: radial_grid(:)
        integer :: io_status

        info = terpsichore_matrix_fixture_invalid
        if (vacuum_intervals /= 0) return
        read (unit, iostat=io_status) fixture%intervals, &
            fixture%poloidal_points, fixture%toroidal_points, &
            fixture%stability_periods, fixture%field_periods, fixture%modes
        if (io_status /= 0 .or. .not. valid_sizes(fixture)) return
        call allocate_record_arrays(fixture, radial_grid, flux_t_slope, &
            radial_power, bjac, io_status)
        if (io_status /= 0) return
        call read_radial_record(unit, fixture, radial_grid, flux_t_slope, &
            io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) fixture%mode_m, fixture%mode_n, &
            radial_power
        if (io_status /= 0) return
        read (unit, iostat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) bjac
        if (io_status /= 0) return
        call move_fixture_arrays(fixture, radial_grid, flux_t_slope, &
            radial_power, bjac)
        if (.not. terpsichore_fixed_fixture_is_valid(fixture)) return
        info = terpsichore_matrix_fixture_ok
    end subroutine read_terpsichore_fixed_boundary_fixture

    pure function valid_sizes(fixture) result(valid)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        logical :: valid

        valid = fixture%intervals >= 2 &
            .and. fixture%intervals <= maximum_intervals &
            .and. fixture%poloidal_points >= 1 &
            .and. fixture%poloidal_points <= maximum_angular_points &
            .and. fixture%toroidal_points >= 1 &
            .and. fixture%toroidal_points <= maximum_angular_points &
            .and. fixture%stability_periods >= 1 &
            .and. fixture%field_periods >= 1 .and. fixture%modes >= 1 &
            .and. fixture%modes <= maximum_modes
        if (.not. valid) return
        valid = modulo(fixture%field_periods, fixture%stability_periods) == 0
        if (.not. valid) return
        valid = fixture%poloidal_points <= maximum_fixture_values &
            / fixture%toroidal_points / (fixture%intervals + 1)
        if (.not. valid) return
        valid = fixture%modes <= maximum_phase_values &
            / fixture%poloidal_points / fixture%toroidal_points &
            / fixture%intervals / 2
        if (.not. valid) return
        valid = fixture%modes <= maximum_dense_order &
            / (2 * fixture%intervals - 1)
    end function valid_sizes

    subroutine allocate_record_arrays(fixture, radial_grid, flux_t_slope, &
            radial_power, bjac, allocation_status)
        type(terpsichore_matrix_fixture_t), intent(inout) :: fixture
        real(dp), allocatable, intent(out) :: radial_grid(:)
        real(dp), allocatable, intent(out) :: flux_t_slope(:)
        real(dp), allocatable, intent(out) :: radial_power(:)
        real(dp), allocatable, intent(out) :: bjac(:, :)
        integer, intent(out) :: allocation_status
        integer :: points

        allocation_status = 1
        points = fixture%poloidal_points * fixture%toroidal_points
        allocate (radial_grid(0:fixture%intervals), stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (flux_t_slope(fixture%intervals + 1), &
            stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (radial_power(fixture%modes), stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (bjac(points, 0:fixture%intervals), stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (fixture%mode_m(fixture%modes), stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (fixture%mode_n(fixture%modes), stat=allocation_status)
    end subroutine allocate_record_arrays

    subroutine read_radial_record(unit, fixture, radial_grid, &
            flux_t_slope, io_status)
        integer, intent(in) :: unit
        type(terpsichore_matrix_fixture_t), intent(inout) :: fixture
        real(dp), intent(out) :: radial_grid(0:)
        real(dp), intent(out) :: flux_t_slope(:)
        integer, intent(out) :: io_status
        real(dp) :: parity
        real(dp), allocatable :: interval_values(:), surface_values(:)

        allocate (interval_values(fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        allocate (surface_values(fixture%intervals + 1), stat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) radial_grid, interval_values, &
            surface_values, flux_t_slope, interval_values, interval_values, &
            interval_values, parity
        if (io_status == 0) fixture%parity = real(parity, dp)
    end subroutine read_radial_record

    subroutine move_fixture_arrays(fixture, radial_grid, flux_t_slope, &
            radial_power, bjac)
        type(terpsichore_matrix_fixture_t), intent(inout) :: fixture
        real(dp), allocatable, intent(inout) :: radial_grid(:)
        real(dp), allocatable, intent(inout) :: flux_t_slope(:)
        real(dp), allocatable, intent(inout) :: radial_power(:)
        real(dp), allocatable, intent(inout) :: bjac(:, :)

        call move_alloc(radial_grid, fixture%s)
        call move_alloc(flux_t_slope, fixture%flux_t_slope)
        call move_alloc(radial_power, fixture%radial_power)
        call move_alloc(bjac, fixture%signed_bjac)
    end subroutine move_fixture_arrays

    pure function terpsichore_fixed_fixture_is_valid(fixture) result(valid)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        logical :: valid
        real(dp), parameter :: endpoint_tolerance = 64.0_dp * epsilon(1.0_dp)
        integer :: points

        valid = .false.
        if (.not. valid_sizes(fixture)) return
        if (.not. allocated(fixture%s) &
            .or. .not. allocated(fixture%flux_t_slope)) return
        if (.not. allocated(fixture%radial_power) &
            .or. .not. allocated(fixture%signed_bjac)) return
        if (.not. allocated(fixture%mode_m) &
            .or. .not. allocated(fixture%mode_n)) return
        points = fixture%poloidal_points * fixture%toroidal_points
        if (lbound(fixture%s, 1) /= 0 &
            .or. ubound(fixture%s, 1) /= fixture%intervals) return
        if (size(fixture%flux_t_slope) /= fixture%intervals + 1) return
        if (size(fixture%radial_power) /= fixture%modes) return
        if (size(fixture%mode_m) /= fixture%modes &
            .or. size(fixture%mode_n) /= fixture%modes) return
        if (any(shape(fixture%signed_bjac) /= &
            [points, fixture%intervals + 1])) return
        if (.not. ieee_is_finite(fixture%parity)) return
        if (.not. all(ieee_is_finite(fixture%s))) return
        if (.not. all(ieee_is_finite(fixture%flux_t_slope))) return
        if (.not. all(ieee_is_finite(fixture%radial_power))) return
        if (.not. all(ieee_is_finite(fixture%signed_bjac))) return
        if (abs(fixture%s(0)) > endpoint_tolerance) return
        if (abs(fixture%s(fixture%intervals) - 1.0_dp) &
            > endpoint_tolerance) return
        if (any(fixture%s(1:) <= fixture%s(:fixture%intervals - 1))) return
        if (any(fixture%flux_t_slope(:fixture%intervals) == 0.0_dp)) return
        if (any(fixture%signed_bjac(:, 1:) >= 0.0_dp)) return
        valid = .true.
    end function terpsichore_fixed_fixture_is_valid

end module terpsichore_matrix_fixture
