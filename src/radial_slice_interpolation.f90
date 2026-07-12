module radial_slice_interpolation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    ! Two-point Gauss nodes on the unit cell.
    real(dp), parameter, public :: gauss_two_lower = &
        0.5_dp - 0.5_dp / sqrt(3.0_dp)
    real(dp), parameter, public :: gauss_two_upper = &
        0.5_dp + 0.5_dp / sqrt(3.0_dp)

    public :: interpolate_slice_pair

contains

    ! Linear interpolation of per-cell surface slices to an intra-cell
    ! evaluation coordinate: slices live at cell midpoints, so the
    ! target between midpoints blends the two neighbours; the first
    ! and last half cells clamp to their own slice.
    pure subroutine interpolate_slice_pair(interval, intervals, &
            local_coordinate, left_interval, right_interval, &
            right_fraction)
        integer, intent(in) :: interval, intervals
        real(dp), intent(in) :: local_coordinate
        integer, intent(out) :: left_interval, right_interval
        real(dp), intent(out) :: right_fraction

        if (local_coordinate >= 0.5_dp) then
            left_interval = interval
            right_interval = min(interval + 1, intervals)
            right_fraction = local_coordinate - 0.5_dp
        else
            left_interval = max(interval - 1, 1)
            right_interval = interval
            right_fraction = local_coordinate + 0.5_dp
        end if
        if (left_interval == right_interval) right_fraction = 0.0_dp
    end subroutine interpolate_slice_pair

end module radial_slice_interpolation
