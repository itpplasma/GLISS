module cas3d_coefficient_mass
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: cas3d2mn_envelope_mass_scale
    public :: cas3d2_physical_mass_scale

contains

    pure function cas3d2_physical_mass_scale(radial_intervals, &
            poloidal_points, toroidal_points, reference_length) result(scale)
        integer, intent(in) :: radial_intervals, poloidal_points
        integer, intent(in) :: toroidal_points
        real(dp), intent(in) :: reference_length
        real(dp) :: normalized_length, scale

        scale = 0.0_dp
        if (radial_intervals < 1) return
        if (poloidal_points < 1 .or. toroidal_points < 1) return
        if (.not. ieee_is_finite(reference_length)) return
        if (reference_length <= 0.0_dp) return
        normalized_length = reference_length / real(radial_intervals, dp)
        scale = normalized_length &
            * reference_length * reference_length &
            / real(poloidal_points, dp) / real(toroidal_points, dp)
        if (.not. ieee_is_finite(scale)) scale = 0.0_dp
    end function cas3d2_physical_mass_scale

    pure function cas3d2mn_envelope_mass_scale(radial_intervals, &
            poloidal_points, toroidal_points, reference_length) result(scale)
        integer, intent(in) :: radial_intervals, poloidal_points
        integer, intent(in) :: toroidal_points
        real(dp), intent(in) :: reference_length
        real(dp) :: scale

        scale = 0.5_dp * cas3d2_physical_mass_scale(radial_intervals, &
            poloidal_points, toroidal_points, reference_length)
    end function cas3d2mn_envelope_mass_scale

end module cas3d_coefficient_mass
