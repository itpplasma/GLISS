module terpsichore_topology
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: terpsichore_topology_ok = 0
    integer, parameter, public :: terpsichore_topology_invalid = -1
    integer, parameter, public :: parity_xi_cosine = 1
    integer, parameter, public :: parity_xi_sine = 2

    type, public :: terpsichore_topology_config_t
        integer :: equilibrium_periods = 0
        integer :: field_periods_per_stability_period = 0
        integer :: poloidal_shift = 0
        real(dp) :: parfac = 0.0_dp
        real(dp) :: qn = 0.0_dp
    end type terpsichore_topology_config_t

    type, public :: terpsichore_mode_mask_t
        integer :: poloidal_min = 0
        integer :: toroidal_min = 0
        logical, allocatable :: selected(:, :)
    end type terpsichore_mode_mask_t

    type, public :: terpsichore_mode_selection_t
        integer :: field_periods = 0
        integer :: parity_class = 0
        integer, allocatable :: poloidal(:)
        integer, allocatable :: toroidal(:)
        real(dp), allocatable :: stored_variable_power(:)
    end type terpsichore_mode_selection_t

    public :: convert_terpsichore_mask

contains

    pure subroutine convert_terpsichore_mask(config, mask, selection, info)
        type(terpsichore_topology_config_t), intent(in) :: config
        type(terpsichore_mode_mask_t), intent(in) :: mask
        type(terpsichore_mode_selection_t), intent(out) :: selection
        integer, intent(out) :: info
        integer :: selected_count, i, j, mode, ratio, shifted_m, deck_n

        info = terpsichore_topology_invalid
        if (.not. valid_config(config)) return
        if (.not. allocated(mask%selected)) return
        selected_count = count(mask%selected)
        if (selected_count == 0) return
        ratio = config%equilibrium_periods &
            / config%field_periods_per_stability_period
        selection%field_periods = config%equilibrium_periods
        selection%parity_class = parity_from(config%parfac)
        allocate (selection%poloidal(selected_count), &
            selection%toroidal(selected_count), &
            selection%stored_variable_power(selected_count))
        mode = 0
        do j = 1, size(mask%selected, 2)
            deck_n = mask%toroidal_min + j - 1
            do i = 1, size(mask%selected, 1)
                if (.not. mask%selected(i, j)) cycle
                shifted_m = mask%poloidal_min + i - 1 &
                    + config%poloidal_shift
                if (shifted_m < 0) return
                mode = mode + 1
                selection%poloidal(mode) = shifted_m
                selection%toroidal(mode) = deck_n * ratio
                selection%stored_variable_power(mode) = 0.0_dp
                if (shifted_m == 1) then
                    selection%stored_variable_power(mode) = config%qn
                end if
            end do
        end do
        info = terpsichore_topology_ok
    end subroutine convert_terpsichore_mask

    pure function valid_config(config) result(valid)
        type(terpsichore_topology_config_t), intent(in) :: config
        logical :: valid
        real(dp), parameter :: tolerance = 64.0_dp * epsilon(1.0_dp)

        valid = config%equilibrium_periods > 0
        valid = valid .and. config%field_periods_per_stability_period > 0
        if (.not. valid) return
        valid = mod(config%equilibrium_periods, &
            config%field_periods_per_stability_period) == 0
        valid = valid .and. ieee_is_finite(config%parfac)
        valid = valid .and. ieee_is_finite(config%qn)
        valid = valid .and. (abs(config%parfac) <= tolerance .or. &
            abs(config%parfac - 0.5_dp) <= tolerance)
    end function valid_config

    pure function parity_from(parfac) result(parity)
        real(dp), intent(in) :: parfac
        integer :: parity

        parity = parity_xi_sine
        if (parfac > 0.25_dp) parity = parity_xi_cosine
    end function parity_from

end module terpsichore_topology
