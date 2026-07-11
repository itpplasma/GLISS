module helical_cylinder_limit
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter, public :: helical_limit_invalid_elongation = 1
    integer, parameter, public :: helical_limit_ok = 0

    public :: benchmark_helical_vertical_margin
    public :: elliptical_wall_radius_ratio
    public :: helical_vertical_margin
    public :: helical_vertical_threshold

contains

    pure function helical_vertical_margin(elongation, external_fraction) &
            result(margin)
        real(dp), intent(in) :: elongation, external_fraction
        real(dp) :: margin

        margin = 1.0_dp + elongation - (elongation**2 + 1.0_dp) &
            * (1.0_dp - external_fraction)
    end function helical_vertical_margin

    pure function helical_vertical_threshold(elongation) result(fraction)
        real(dp), intent(in) :: elongation
        real(dp) :: fraction

        fraction = (elongation**2 - elongation) &
            / (elongation**2 + 1.0_dp)
    end function helical_vertical_threshold

    pure subroutine elliptical_wall_radius_ratio(elongation, ratio, info)
        real(dp), intent(in) :: elongation
        real(dp), intent(out) :: ratio
        integer, intent(out) :: info

        if (.not. ieee_is_finite(elongation)) then
            ratio = 0.0_dp
            info = helical_limit_invalid_elongation
            return
        end if
        if (elongation <= 1.0_dp) then
            ratio = 0.0_dp
            info = helical_limit_invalid_elongation
            return
        end if
        ratio = sqrt((elongation + 1.0_dp) / (elongation - 1.0_dp))
        info = helical_limit_ok
    end subroutine elliptical_wall_radius_ratio

    pure function benchmark_helical_vertical_margin(external_fraction) &
            result(margin) bind(c, name="gliss_benchmark_helical_vertical_margin")
        real(c_double), intent(in), value :: external_fraction
        real(c_double) :: margin

        margin = helical_vertical_margin(2.0_dp, external_fraction)
    end function benchmark_helical_vertical_margin

end module helical_cylinder_limit
