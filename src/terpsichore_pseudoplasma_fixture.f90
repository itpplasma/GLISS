module terpsichore_pseudoplasma_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: pseudoplasma_fixture_ok = 0
    integer, parameter, public :: pseudoplasma_fixture_invalid = -1
    integer, parameter :: fixture_magic = int(z'47565031')
    integer, parameter :: fixture_schema = 1
    integer, parameter :: maximum_intervals = 996
    integer, parameter :: maximum_modes = 10000
    integer, parameter :: maximum_coefficients = 50000000
    integer, parameter :: maximum_matrix_entries = 16000000
    integer, parameter :: maximum_mode_number = 1000000

    type, public :: terpsichore_pseudoplasma_fixture_t
        integer :: plasma_intervals = 0
        integer :: vacuum_intervals = 0
        integer :: modes = 0
        real(dp) :: flux_t_slope = 0.0_dp
        real(dp) :: flux_p_slope = 0.0_dp
        real(dp), allocatable :: s(:)
        integer, allocatable :: mode_m(:)
        integer, allocatable :: mode_n(:)
        real(dp), allocatable :: coefficient(:, :, :, :)
    end type terpsichore_pseudoplasma_fixture_t

    public :: read_terpsichore_pseudoplasma_fixture
    public :: terpsichore_pseudoplasma_fixture_is_valid

contains

    subroutine read_terpsichore_pseudoplasma_fixture(unit, fixture, info)
        integer, intent(in) :: unit
        type(terpsichore_pseudoplasma_fixture_t), intent(out) :: fixture
        integer, intent(out) :: info
        integer :: allocation_status, magic, schema

        info = pseudoplasma_fixture_invalid
        read (unit, iostat=allocation_status) magic, schema, &
            fixture%plasma_intervals, fixture%vacuum_intervals, fixture%modes
        if (allocation_status /= 0) return
        if (magic /= fixture_magic .or. schema /= fixture_schema) return
        if (.not. metadata_is_valid(fixture)) return
        allocate (fixture%s(0:fixture%vacuum_intervals), &
            fixture%mode_m(fixture%modes), fixture%mode_n(fixture%modes), &
            fixture%coefficient(6, fixture%modes, fixture%modes, &
            fixture%vacuum_intervals), stat=allocation_status)
        if (allocation_status /= 0) return
        read (unit, iostat=allocation_status) fixture%s, &
            fixture%flux_t_slope, fixture%flux_p_slope
        if (allocation_status /= 0) return
        read (unit, iostat=allocation_status) fixture%mode_m, fixture%mode_n
        if (allocation_status /= 0) return
        read (unit, iostat=allocation_status) fixture%coefficient(1, :, :, :), &
            fixture%coefficient(2, :, :, :), &
            fixture%coefficient(3, :, :, :), &
            fixture%coefficient(4, :, :, :), &
            fixture%coefficient(5, :, :, :), &
            fixture%coefficient(6, :, :, :)
        if (allocation_status /= 0) return
        if (.not. terpsichore_pseudoplasma_fixture_is_valid(fixture)) return
        info = pseudoplasma_fixture_ok
    end subroutine read_terpsichore_pseudoplasma_fixture

    pure function metadata_is_valid(fixture) result(valid)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        logical :: valid
        integer :: element_order, matrix_order

        valid = fixture%plasma_intervals >= 2 &
            .and. fixture%plasma_intervals <= maximum_intervals &
            .and. fixture%vacuum_intervals >= 1 &
            .and. fixture%vacuum_intervals <= maximum_intervals &
            .and. fixture%modes >= 1 .and. fixture%modes <= maximum_modes
        if (.not. valid) return
        valid = fixture%modes <= maximum_coefficients / fixture%modes &
            / fixture%vacuum_intervals / 6
        if (.not. valid) return
        matrix_order = 2 * fixture%modes * fixture%vacuum_intervals
        valid = matrix_order <= maximum_matrix_entries / matrix_order
        if (.not. valid) return
        element_order = 3 * fixture%modes
        valid = element_order <= maximum_matrix_entries / element_order
    end function metadata_is_valid

    pure function terpsichore_pseudoplasma_fixture_is_valid(fixture) &
            result(valid)
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: fixture
        logical :: valid

        valid = metadata_is_valid(fixture)
        if (.not. valid) return
        if (.not. allocated(fixture%s)) then
            valid = .false.
            return
        end if
        if (.not. allocated(fixture%mode_m) .or. &
            .not. allocated(fixture%mode_n)) then
            valid = .false.
            return
        end if
        if (.not. allocated(fixture%coefficient)) then
            valid = .false.
            return
        end if
        valid = lbound(fixture%s, 1) == 0 &
            .and. ubound(fixture%s, 1) == fixture%vacuum_intervals
        valid = valid .and. size(fixture%mode_m) == fixture%modes &
            .and. size(fixture%mode_n) == fixture%modes
        valid = valid .and. all(shape(fixture%coefficient) == &
            [6, fixture%modes, fixture%modes, fixture%vacuum_intervals])
        if (.not. valid) return
        valid = all(abs(real(fixture%mode_m, dp)) <= maximum_mode_number) &
            .and. all(abs(real(fixture%mode_n, dp)) <= maximum_mode_number)
        if (.not. valid) return
        valid = all(ieee_is_finite(fixture%s)) &
            .and. ieee_is_finite(fixture%flux_t_slope) &
            .and. ieee_is_finite(fixture%flux_p_slope) &
            .and. all(ieee_is_finite(fixture%coefficient))
        if (.not. valid) return
        valid = all(fixture%s(1:) > fixture%s(:fixture%vacuum_intervals - 1))
    end function terpsichore_pseudoplasma_fixture_is_valid

end module terpsichore_pseudoplasma_fixture
