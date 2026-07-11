module terpsichore_matrix_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use terpsichore_model_policy, only: decode_terpsichore_model, &
        terpsichore_model_config_t, terpsichore_model_ok
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
    integer, parameter :: maximum_potential_values = 20000000

    type, public :: terpsichore_matrix_fixture_t
        integer :: intervals = 0
        integer :: poloidal_points = 0
        integer :: toroidal_points = 0
        integer :: stability_periods = 0
        integer :: field_periods = 0
        integer :: modes = 0
        integer :: legacy_modelk = -1
        real(dp) :: parity = 0.0_dp
        real(dp) :: current_factor = 0.0_dp
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: flux_p_slope(:)
        real(dp), allocatable :: flux_t_slope(:)
        real(dp), allocatable :: flux_p_curve(:)
        real(dp), allocatable :: flux_t_curve(:)
        real(dp), allocatable :: current_i(:)
        real(dp), allocatable :: current_j(:)
        real(dp), allocatable :: pressure_slope(:)
        real(dp), allocatable :: radial_power(:)
        real(dp), allocatable :: signed_bjac(:, :)
        real(dp), allocatable :: signed_bjac_radial(:, :)
        real(dp), allocatable :: sigma_b_s(:, :)
        real(dp), allocatable :: metric_ss_over_jacobian(:, :)
        real(dp), allocatable :: metric_st_over_jacobian(:, :)
        real(dp), allocatable :: metric_tt_over_jacobian(:, :)
        real(dp), allocatable :: sigma_b(:, :)
        real(dp), allocatable :: parallel_current(:, :)
        integer, allocatable :: mode_m(:)
        integer, allocatable :: mode_n(:)
    end type terpsichore_matrix_fixture_t

    public :: read_terpsichore_fixed_boundary_fixture
    public :: read_terpsichore_fixed_boundary_potential_fixture
    public :: terpsichore_fixed_fixture_is_valid
    public :: terpsichore_potential_metadata_is_valid
    public :: terpsichore_potential_fixture_is_valid

contains

    subroutine read_terpsichore_fixed_boundary_fixture(unit, vacuum_intervals, &
            fixture, info)
        integer, intent(in) :: unit, vacuum_intervals
        type(terpsichore_matrix_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        integer :: io_status

        info = terpsichore_matrix_fixture_invalid
        if (vacuum_intervals /= 0) return
        call read_fixture_prefix(unit, .false., fixture, io_status)
        if (io_status /= 0) return
        if (.not. terpsichore_fixed_fixture_is_valid(fixture)) return
        info = terpsichore_matrix_fixture_ok
    end subroutine read_terpsichore_fixed_boundary_fixture

    subroutine read_terpsichore_fixed_boundary_potential_fixture(unit, &
            vacuum_intervals, fixture, info)
        integer, intent(in) :: unit, vacuum_intervals
        type(terpsichore_matrix_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        integer :: io_status

        info = terpsichore_matrix_fixture_invalid
        if (vacuum_intervals /= 0) return
        call read_fixture_prefix(unit, .true., fixture, io_status)
        if (io_status /= 0) return
        call read_potential_records(unit, fixture, io_status)
        if (io_status /= 0) return
        if (.not. terpsichore_potential_fixture_is_valid(fixture)) return
        info = terpsichore_matrix_fixture_ok
    end subroutine read_terpsichore_fixed_boundary_potential_fixture

    subroutine read_fixture_prefix(unit, retain_potential, fixture, info)
        integer, intent(in) :: unit
        logical, intent(in) :: retain_potential
        type(terpsichore_matrix_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        real(dp), allocatable :: bjac(:, :), flux_t_slope(:), radial_power(:)
        real(dp), allocatable :: radial_grid(:)

        info = 1
        read (unit, iostat=info) fixture%intervals, &
            fixture%poloidal_points, fixture%toroidal_points, &
            fixture%stability_periods, fixture%field_periods, fixture%modes
        if (info /= 0 .or. .not. valid_sizes(fixture)) return
        if (retain_potential .and. &
            .not. terpsichore_potential_metadata_is_valid(fixture)) return
        call allocate_record_arrays(fixture, radial_grid, flux_t_slope, &
            radial_power, bjac, info)
        if (info /= 0) return
        call read_radial_record(unit, fixture, radial_grid, flux_t_slope, &
            info)
        if (info /= 0) return
        read (unit, iostat=info) fixture%mode_m, fixture%mode_n, &
            radial_power
        if (info /= 0) return
        read (unit, iostat=info)
        if (info /= 0) return
        call read_geometry_record(unit, retain_potential, fixture, bjac, &
            info)
        if (info /= 0) return
        call move_fixture_arrays(fixture, radial_grid, flux_t_slope, &
            radial_power, bjac)
        info = 0
    end subroutine read_fixture_prefix

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
        real(dp), allocatable :: interval_values(:)

        allocate (interval_values(fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%flux_p_slope(fixture%intervals + 1), &
            stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%current_i(fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%current_j(fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%pressure_slope(fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) radial_grid, interval_values, &
            fixture%flux_p_slope, flux_t_slope, fixture%current_j, &
            fixture%current_i, fixture%pressure_slope, parity
        if (io_status == 0) fixture%parity = real(parity, dp)
    end subroutine read_radial_record

    subroutine read_geometry_record(unit, retain_potential, fixture, bjac, &
            io_status)
        integer, intent(in) :: unit
        logical, intent(in) :: retain_potential
        type(terpsichore_matrix_fixture_t), intent(inout) :: fixture
        real(dp), intent(out) :: bjac(:, 0:)
        integer, intent(out) :: io_status
        integer :: points

        if (.not. retain_potential) then
            read (unit, iostat=io_status) bjac
            return
        end if
        points = fixture%poloidal_points * fixture%toroidal_points
        allocate (fixture%sigma_b_s(points, 0:fixture%intervals), &
            stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%metric_ss_over_jacobian(points, &
            0:fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%metric_st_over_jacobian(points, &
            0:fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%metric_tt_over_jacobian(points, &
            0:fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) bjac, fixture%sigma_b_s, &
            fixture%metric_ss_over_jacobian, &
            fixture%metric_st_over_jacobian, &
            fixture%metric_tt_over_jacobian
    end subroutine read_geometry_record

    subroutine read_potential_records(unit, fixture, io_status)
        integer, intent(in) :: unit
        type(terpsichore_matrix_fixture_t), intent(inout) :: fixture
        integer, intent(out) :: io_status
        real(dp), allocatable :: bjac_modes(:, :), point_cell(:, :)
        integer, allocatable :: equilibrium_m(:), equilibrium_n(:)
        integer :: equilibrium_modes, points

        read (unit, iostat=io_status) equilibrium_modes
        if (io_status /= 0 .or. equilibrium_modes < 1 &
            .or. equilibrium_modes > maximum_modes) return
        backspace (unit, iostat=io_status)
        if (io_status /= 0) return
        points = fixture%poloidal_points * fixture%toroidal_points
        allocate (equilibrium_m(equilibrium_modes), &
            equilibrium_n(equilibrium_modes), stat=io_status)
        if (io_status /= 0) return
        allocate (bjac_modes(equilibrium_modes, 0:fixture%intervals), &
            stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%signed_bjac_radial(points, 0:fixture%intervals), &
            stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%flux_t_curve(fixture%intervals + 1), &
            fixture%flux_p_curve(fixture%intervals + 1), stat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) equilibrium_modes, equilibrium_m, &
            equilibrium_n, bjac_modes, fixture%signed_bjac_radial, &
            fixture%flux_t_curve, fixture%flux_p_curve
        if (io_status /= 0) return
        allocate (point_cell(points, fixture%intervals), stat=io_status)
        if (io_status /= 0) return
        call read_potential_point_record(unit, fixture, point_cell, io_status)
    end subroutine read_potential_records

    subroutine read_potential_point_record(unit, fixture, scratch, io_status)
        integer, intent(in) :: unit
        type(terpsichore_matrix_fixture_t), intent(inout) :: fixture
        real(dp), intent(inout) :: scratch(:, :)
        integer, intent(out) :: io_status
        real(dp), allocatable :: gparp(:, :), gperp(:, :), taub(:, :)
        real(dp), allocatable :: parkur(:, :)
        integer :: points

        points = size(scratch, 1)
        allocate (gparp, mold=scratch, stat=io_status)
        if (io_status /= 0) return
        allocate (gperp, mold=scratch, stat=io_status)
        if (io_status /= 0) return
        allocate (taub, mold=scratch, stat=io_status)
        if (io_status /= 0) return
        allocate (parkur, mold=scratch, stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%sigma_b(points, fixture%intervals), &
            stat=io_status)
        if (io_status /= 0) return
        allocate (fixture%parallel_current(points, fixture%intervals), &
            stat=io_status)
        if (io_status /= 0) return
        read (unit, iostat=io_status) gparp, gperp, fixture%sigma_b, taub, &
            parkur, fixture%parallel_current, fixture%current_factor, &
            fixture%legacy_modelk
    end subroutine read_potential_point_record

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

    pure function terpsichore_potential_fixture_is_valid(fixture) &
            result(valid)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        logical :: valid
        integer :: cell_shape(2), points, surface_shape(2)
        integer :: model_status
        type(terpsichore_model_config_t) :: model

        valid = .false.
        if (.not. terpsichore_fixed_fixture_is_valid(fixture)) return
        if (.not. terpsichore_potential_metadata_is_valid(fixture)) return
        call decode_terpsichore_model(fixture%legacy_modelk, model, &
            model_status)
        if (model_status /= terpsichore_model_ok) return
        if (.not. allocated(fixture%flux_p_slope) &
            .or. .not. allocated(fixture%flux_t_curve)) return
        if (.not. allocated(fixture%flux_p_curve) &
            .or. .not. allocated(fixture%current_i)) return
        if (.not. allocated(fixture%current_j) &
            .or. .not. allocated(fixture%pressure_slope)) return
        if (.not. allocated(fixture%signed_bjac_radial) &
            .or. .not. allocated(fixture%sigma_b_s)) return
        if (.not. allocated(fixture%metric_ss_over_jacobian) &
            .or. .not. allocated(fixture%metric_st_over_jacobian)) return
        if (.not. allocated(fixture%metric_tt_over_jacobian) &
            .or. .not. allocated(fixture%sigma_b)) return
        if (.not. allocated(fixture%parallel_current)) return
        points = fixture%poloidal_points * fixture%toroidal_points
        surface_shape = [points, fixture%intervals + 1]
        cell_shape = [points, fixture%intervals]
        if (size(fixture%flux_p_slope) /= fixture%intervals + 1) return
        if (size(fixture%flux_t_curve) /= fixture%intervals + 1 &
            .or. size(fixture%flux_p_curve) /= fixture%intervals + 1) return
        if (size(fixture%current_i) /= fixture%intervals &
            .or. size(fixture%current_j) /= fixture%intervals) return
        if (size(fixture%pressure_slope) /= fixture%intervals) return
        if (any(shape(fixture%signed_bjac_radial) /= surface_shape)) return
        if (any(shape(fixture%sigma_b_s) /= surface_shape)) return
        if (any(shape(fixture%metric_ss_over_jacobian) /= surface_shape)) &
            return
        if (any(shape(fixture%metric_st_over_jacobian) /= surface_shape)) &
            return
        if (any(shape(fixture%metric_tt_over_jacobian) /= surface_shape)) &
            return
        if (any(shape(fixture%sigma_b) /= cell_shape)) return
        if (any(shape(fixture%parallel_current) /= cell_shape)) return
        if (.not. all(ieee_is_finite(fixture%flux_p_slope))) return
        if (.not. all(ieee_is_finite(fixture%flux_t_curve))) return
        if (.not. all(ieee_is_finite(fixture%flux_p_curve))) return
        if (.not. all(ieee_is_finite(fixture%current_i))) return
        if (.not. all(ieee_is_finite(fixture%current_j))) return
        if (.not. all(ieee_is_finite(fixture%pressure_slope))) return
        if (.not. all(ieee_is_finite(fixture%signed_bjac_radial))) return
        if (.not. all(ieee_is_finite(fixture%sigma_b_s))) return
        if (.not. all(ieee_is_finite(fixture%metric_ss_over_jacobian))) return
        if (.not. all(ieee_is_finite(fixture%metric_st_over_jacobian))) return
        if (.not. all(ieee_is_finite(fixture%metric_tt_over_jacobian))) return
        if (.not. all(ieee_is_finite(fixture%sigma_b))) return
        if (.not. all(ieee_is_finite(fixture%parallel_current))) return
        if (.not. ieee_is_finite(fixture%current_factor)) return
        if (any(fixture%sigma_b == 0.0_dp)) return
        if (any(fixture%current_j * fixture%flux_p_slope(:fixture%intervals) &
            - fixture%current_i * fixture%flux_t_slope(:fixture%intervals) &
            == 0.0_dp)) return
        valid = .true.
    end function terpsichore_potential_fixture_is_valid

    pure function terpsichore_potential_metadata_is_valid(fixture) &
            result(valid)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        logical :: valid
        integer :: points

        valid = valid_sizes(fixture)
        if (.not. valid) return
        points = fixture%poloidal_points * fixture%toroidal_points
        valid = points <= maximum_potential_values &
            / (13 * (fixture%intervals + 1))
    end function terpsichore_potential_metadata_is_valid

end module terpsichore_matrix_fixture
