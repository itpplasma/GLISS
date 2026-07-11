module terpsichore_normalization
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: terpsichore_normalization_ok = 0
    integer, parameter, public :: terpsichore_normalization_invalid = -1

    public :: map_gliss_export_flip_pol_cell_to_terpsichore
    public :: map_gliss_internal_flip_pol_cell_to_terpsichore

contains

    pure subroutine map_gliss_export_flip_pol_cell_to_terpsichore( &
            field_periods, reference_intervals, s_left, s_right, &
            signed_jacobian, exported_phi_slope, exported_chi_slope, &
            signed_bjac, terpsichore_ftp, terpsichore_fpp, &
            normalized_radial_weight, info)
        integer, intent(in) :: field_periods, reference_intervals
        real(dp), intent(in) :: s_left, s_right, signed_jacobian(:)
        real(dp), intent(in) :: exported_phi_slope, exported_chi_slope
        real(dp), intent(out) :: signed_bjac(:), terpsichore_ftp
        real(dp), intent(out) :: terpsichore_fpp
        real(dp), intent(out) :: normalized_radial_weight
        integer, intent(out) :: info

        call map_angular_cell(field_periods, reference_intervals, s_left, &
            s_right, signed_jacobian, exported_phi_slope, exported_chi_slope, &
            signed_bjac, terpsichore_ftp, terpsichore_fpp, &
            normalized_radial_weight, info)
    end subroutine map_gliss_export_flip_pol_cell_to_terpsichore

    pure subroutine map_gliss_internal_flip_pol_cell_to_terpsichore( &
            field_periods, reference_intervals, s_left, s_right, &
            signed_jacobian, flux_t_slope, flux_p_slope, signed_bjac, &
            terpsichore_ftp, terpsichore_fpp, normalized_radial_weight, info)
        integer, intent(in) :: field_periods, reference_intervals
        real(dp), intent(in) :: s_left, s_right, signed_jacobian(:)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(out) :: signed_bjac(:), terpsichore_ftp
        real(dp), intent(out) :: terpsichore_fpp
        real(dp), intent(out) :: normalized_radial_weight
        integer, intent(out) :: info

        call map_angular_cell(field_periods, reference_intervals, s_left, &
            s_right, signed_jacobian, -flux_t_slope, &
            -real(field_periods, dp) * flux_p_slope, signed_bjac, &
            terpsichore_ftp, terpsichore_fpp, normalized_radial_weight, info)
    end subroutine map_gliss_internal_flip_pol_cell_to_terpsichore

    pure subroutine map_angular_cell(field_periods, reference_intervals, &
            s_left, s_right, signed_jacobian, ftp_input, fpp_input, &
            signed_bjac, terpsichore_ftp, terpsichore_fpp, &
            normalized_radial_weight, info)
        integer, intent(in) :: field_periods, reference_intervals
        real(dp), intent(in) :: s_left, s_right, signed_jacobian(:)
        real(dp), intent(in) :: ftp_input, fpp_input
        real(dp), intent(out) :: signed_bjac(:), terpsichore_ftp
        real(dp), intent(out) :: terpsichore_fpp
        real(dp), intent(out) :: normalized_radial_weight
        integer, intent(out) :: info
        real(dp) :: angular_scale

        info = terpsichore_normalization_invalid
        if (field_periods < 1 .or. reference_intervals < 1) return
        if (size(signed_jacobian) < 1) return
        if (size(signed_bjac) /= size(signed_jacobian)) return
        if (.not. ieee_is_finite(s_left) .or. .not. ieee_is_finite(s_right)) &
            return
        if (s_left < 0.0_dp .or. s_right > 1.0_dp &
            .or. s_right <= s_left) return
        if (.not. all(ieee_is_finite(signed_jacobian))) return
        if (any(signed_jacobian >= 0.0_dp)) return
        if (.not. ieee_is_finite(ftp_input) .or. ftp_input == 0.0_dp) return
        if (.not. ieee_is_finite(fpp_input)) return
        angular_scale = real(field_periods, dp) / (4.0_dp * acos(-1.0_dp)**2)
        signed_bjac = angular_scale * signed_jacobian
        terpsichore_ftp = ftp_input
        terpsichore_fpp = fpp_input
        normalized_radial_weight = real(reference_intervals, dp) &
            * (s_right - s_left)
        info = terpsichore_normalization_ok
    end subroutine map_angular_cell

end module terpsichore_normalization
