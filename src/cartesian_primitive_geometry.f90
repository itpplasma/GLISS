module cartesian_primitive_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: primitive_geometry_ok = 0
    integer, parameter, public :: primitive_geometry_invalid = -1

    type, public :: primitive_geometry_point_t
        real(dp) :: metric(3, 3) = 0.0_dp
        real(dp) :: signed_jacobian = 0.0_dp
        real(dp) :: jacobian_s = 0.0_dp
        real(dp) :: jacobian_theta = 0.0_dp
        real(dp) :: jacobian_zeta = 0.0_dp
        real(dp) :: b_contravariant(2) = 0.0_dp
        real(dp) :: b_covariant(2) = 0.0_dp
        real(dp) :: mod_b = 0.0_dp
        real(dp) :: second_form(2, 2) = 0.0_dp
    end type primitive_geometry_point_t

    public :: build_primitive_geometry_point

contains

    ! Share the generated value leaves in this compilation unit so the native
    ! compiler can inline the pointwise build without link-time optimization.
    include 'generated/gliss_geometry_primal.inc'

    subroutine build_primitive_geometry_point(position_s, position_theta, &
            position_zeta, position_ss, position_s_theta, position_s_zeta, &
            position_theta_theta, position_theta_zeta, position_zeta_zeta, &
            field_periods, toroidal_flux_slope, poloidal_flux_slope, &
            geometry, info)
        real(dp), intent(in) :: position_s(3), position_theta(3)
        real(dp), intent(in) :: position_zeta(3), position_ss(3)
        real(dp), intent(in) :: position_s_theta(3), position_s_zeta(3)
        real(dp), intent(in) :: position_theta_theta(3)
        real(dp), intent(in) :: position_theta_zeta(3)
        real(dp), intent(in) :: position_zeta_zeta(3)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: toroidal_flux_slope, poloidal_flux_slope
        type(primitive_geometry_point_t), intent(out) :: geometry
        integer, intent(out) :: info
        type(primitive_geometry_point_t) :: candidate
        real(dp) :: polynomial_input(27), polynomial_value(13)
        real(dp) :: field_input(7), field_value(5)
        real(dp) :: second_input(12), second_value(3)
        real(dp) :: angular_cross(3), area, frame_scale, field_squared

        geometry = primitive_geometry_point_t()
        info = primitive_geometry_invalid
        if (field_periods < 1) return
        if (.not. inputs_are_finite(position_s, position_theta, &
            position_zeta, position_ss, position_s_theta, position_s_zeta, &
            position_theta_theta, position_theta_zeta, position_zeta_zeta, &
            toroidal_flux_slope, poloidal_flux_slope)) return
        polynomial_input(1:3) = position_s
        polynomial_input(4:6) = position_theta
        polynomial_input(7:9) = position_zeta
        polynomial_input(10:12) = position_ss
        polynomial_input(13:15) = position_s_theta
        polynomial_input(16:18) = position_s_zeta
        polynomial_input(19:21) = position_theta_theta
        polynomial_input(22:24) = position_theta_zeta
        polynomial_input(25:27) = position_zeta_zeta
        call generated_geometry_polynomial(polynomial_input, polynomial_value)
        angular_cross = polynomial_value(11:13)
        area = sqrt(dot_product(angular_cross, angular_cross))
        frame_scale = vector_norm(position_s) * vector_norm(position_theta) &
            * vector_norm(position_zeta)
        if (.not. ieee_is_finite(area) .or. &
            .not. ieee_is_finite(frame_scale)) return
        if (area <= tiny(1.0_dp) .or. frame_scale <= tiny(1.0_dp)) return
        candidate%signed_jacobian = polynomial_value(7)
        if (abs(candidate%signed_jacobian) <= 128.0_dp * epsilon(1.0_dp) &
            * frame_scale) return
        candidate%metric(1, 1) = polynomial_value(1)
        candidate%metric(1, 2) = polynomial_value(2)
        candidate%metric(2, 1) = polynomial_value(2)
        candidate%metric(1, 3) = polynomial_value(3)
        candidate%metric(3, 1) = polynomial_value(3)
        candidate%metric(2, 2) = polynomial_value(4)
        candidate%metric(2, 3) = polynomial_value(5)
        candidate%metric(3, 2) = polynomial_value(5)
        candidate%metric(3, 3) = polynomial_value(6)
        candidate%jacobian_s = polynomial_value(8)
        candidate%jacobian_theta = polynomial_value(9)
        candidate%jacobian_zeta = polynomial_value(10)
        field_input(1:3) = polynomial_value(4:6)
        field_input(4) = candidate%signed_jacobian
        field_input(5) = toroidal_flux_slope
        field_input(6) = poloidal_flux_slope
        field_input(7) = 1.0_dp / real(field_periods, dp)
        call generated_geometry_field(field_input, field_value)
        candidate%b_contravariant = field_value(1:2)
        candidate%b_covariant = field_value(3:4)
        field_squared = dot_product(candidate%b_contravariant, &
            candidate%b_covariant)
        if (.not. ieee_is_finite(field_squared) .or. &
            field_squared <= 0.0_dp) return
        candidate%mod_b = field_value(5)
        second_input(1:3) = angular_cross
        second_input(4:6) = position_theta_theta
        second_input(7:9) = position_theta_zeta
        second_input(10:12) = position_zeta_zeta
        call generated_geometry_second_form(second_input, &
            sign(1.0_dp, candidate%signed_jacobian), second_value)
        candidate%second_form(1, 1) = second_value(1)
        candidate%second_form(1, 2) = second_value(2)
        candidate%second_form(2, 1) = second_value(2)
        candidate%second_form(2, 2) = second_value(3)
        if (.not. candidate_is_finite(candidate)) return
        if (candidate%mod_b <= 0.0_dp) return
        geometry = candidate
        info = primitive_geometry_ok
    end subroutine build_primitive_geometry_point

    pure function vector_norm(vector) result(norm)
        real(dp), intent(in) :: vector(3)
        real(dp) :: norm

        norm = sqrt(dot_product(vector, vector))
    end function vector_norm

    pure function inputs_are_finite(position_s, position_theta, &
            position_zeta, position_ss, position_s_theta, position_s_zeta, &
            position_theta_theta, position_theta_zeta, position_zeta_zeta, &
            toroidal_flux_slope, poloidal_flux_slope) result(finite)
        real(dp), intent(in) :: position_s(3), position_theta(3)
        real(dp), intent(in) :: position_zeta(3), position_ss(3)
        real(dp), intent(in) :: position_s_theta(3), position_s_zeta(3)
        real(dp), intent(in) :: position_theta_theta(3)
        real(dp), intent(in) :: position_theta_zeta(3)
        real(dp), intent(in) :: position_zeta_zeta(3)
        real(dp), intent(in) :: toroidal_flux_slope, poloidal_flux_slope
        logical :: finite

        finite = all(ieee_is_finite(position_s)) &
            .and. all(ieee_is_finite(position_theta)) &
            .and. all(ieee_is_finite(position_zeta)) &
            .and. all(ieee_is_finite(position_ss)) &
            .and. all(ieee_is_finite(position_s_theta)) &
            .and. all(ieee_is_finite(position_s_zeta)) &
            .and. all(ieee_is_finite(position_theta_theta)) &
            .and. all(ieee_is_finite(position_theta_zeta)) &
            .and. all(ieee_is_finite(position_zeta_zeta)) &
            .and. ieee_is_finite(toroidal_flux_slope) &
            .and. ieee_is_finite(poloidal_flux_slope)
    end function inputs_are_finite

    pure function candidate_is_finite(candidate) result(finite)
        type(primitive_geometry_point_t), intent(in) :: candidate
        logical :: finite

        finite = all(ieee_is_finite(candidate%metric)) &
            .and. ieee_is_finite(candidate%signed_jacobian) &
            .and. ieee_is_finite(candidate%jacobian_s) &
            .and. ieee_is_finite(candidate%jacobian_theta) &
            .and. ieee_is_finite(candidate%jacobian_zeta) &
            .and. all(ieee_is_finite(candidate%b_contravariant)) &
            .and. all(ieee_is_finite(candidate%b_covariant)) &
            .and. ieee_is_finite(candidate%mod_b) &
            .and. all(ieee_is_finite(candidate%second_form))
    end function candidate_is_finite

end module cartesian_primitive_geometry
