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
        real(dp) :: angular_cross(3), area, frame_scale, normal(3)
        real(dp) :: field_squared

        geometry = primitive_geometry_point_t()
        info = primitive_geometry_invalid
        if (field_periods < 1) return
        if (.not. inputs_are_finite(position_s, position_theta, &
            position_zeta, position_ss, position_s_theta, position_s_zeta, &
            position_theta_theta, position_theta_zeta, position_zeta_zeta, &
            toroidal_flux_slope, poloidal_flux_slope)) return
        call cross_product3(position_theta, position_zeta, angular_cross)
        area = sqrt(dot_product(angular_cross, angular_cross))
        frame_scale = vector_norm(position_s) * vector_norm(position_theta) &
            * vector_norm(position_zeta)
        if (.not. ieee_is_finite(area) .or. &
            .not. ieee_is_finite(frame_scale)) return
        if (area <= tiny(1.0_dp) .or. frame_scale <= tiny(1.0_dp)) return
        candidate%signed_jacobian = dot_product(position_s, angular_cross)
        if (abs(candidate%signed_jacobian) <= 128.0_dp * epsilon(1.0_dp) &
            * frame_scale) return
        call build_metric(position_s, position_theta, position_zeta, &
            candidate%metric)
        call build_jacobian_derivatives(position_s, position_theta, &
            position_zeta, position_ss, position_s_theta, position_s_zeta, &
            position_theta_theta, position_theta_zeta, position_zeta_zeta, &
            candidate%jacobian_s, candidate%jacobian_theta, &
            candidate%jacobian_zeta)
        candidate%b_contravariant(1) = -poloidal_flux_slope &
            / (real(field_periods, dp) * candidate%signed_jacobian)
        candidate%b_contravariant(2) = -toroidal_flux_slope &
            / candidate%signed_jacobian
        candidate%b_covariant(1) = dot_product(candidate%metric(2, 2:3), &
            candidate%b_contravariant)
        candidate%b_covariant(2) = dot_product(candidate%metric(3, 2:3), &
            candidate%b_contravariant)
        field_squared = dot_product(candidate%b_contravariant, &
            candidate%b_covariant)
        if (.not. ieee_is_finite(field_squared) .or. &
            field_squared <= 0.0_dp) return
        candidate%mod_b = sqrt(field_squared)
        normal = sign(1.0_dp, candidate%signed_jacobian) &
            * angular_cross / area
        candidate%second_form(1, 1) = dot_product(normal, &
            position_theta_theta)
        candidate%second_form(1, 2) = dot_product(normal, &
            position_theta_zeta)
        candidate%second_form(2, 1) = candidate%second_form(1, 2)
        candidate%second_form(2, 2) = dot_product(normal, &
            position_zeta_zeta)
        if (.not. candidate_is_finite(candidate)) return
        geometry = candidate
        info = primitive_geometry_ok
    end subroutine build_primitive_geometry_point

    pure subroutine build_metric(position_s, position_theta, position_zeta, &
            metric)
        real(dp), intent(in) :: position_s(3), position_theta(3)
        real(dp), intent(in) :: position_zeta(3)
        real(dp), intent(out) :: metric(3, 3)
        real(dp) :: basis(3, 3)
        integer :: first, second

        basis(:, 1) = position_s
        basis(:, 2) = position_theta
        basis(:, 3) = position_zeta
        do second = 1, 3
            do first = 1, second
                metric(first, second) = dot_product(basis(:, first), &
                    basis(:, second))
                metric(second, first) = metric(first, second)
            end do
        end do
    end subroutine build_metric

    pure subroutine build_jacobian_derivatives(position_s, position_theta, &
            position_zeta, position_ss, position_s_theta, position_s_zeta, &
            position_theta_theta, position_theta_zeta, position_zeta_zeta, &
            jacobian_s, jacobian_theta, jacobian_zeta)
        real(dp), intent(in) :: position_s(3), position_theta(3)
        real(dp), intent(in) :: position_zeta(3), position_ss(3)
        real(dp), intent(in) :: position_s_theta(3), position_s_zeta(3)
        real(dp), intent(in) :: position_theta_theta(3)
        real(dp), intent(in) :: position_theta_zeta(3)
        real(dp), intent(in) :: position_zeta_zeta(3)
        real(dp), intent(out) :: jacobian_s, jacobian_theta, jacobian_zeta
        real(dp) :: product(3)

        call cross_product3(position_theta, position_zeta, product)
        jacobian_s = dot_product(position_ss, product)
        jacobian_theta = dot_product(position_s_theta, product)
        jacobian_zeta = dot_product(position_s_zeta, product)
        call cross_product3(position_s_theta, position_zeta, product)
        jacobian_s = jacobian_s + dot_product(position_s, product)
        call cross_product3(position_theta, position_s_zeta, product)
        jacobian_s = jacobian_s + dot_product(position_s, product)
        call cross_product3(position_theta_theta, position_zeta, product)
        jacobian_theta = jacobian_theta + dot_product(position_s, product)
        call cross_product3(position_theta, position_theta_zeta, product)
        jacobian_theta = jacobian_theta + dot_product(position_s, product)
        call cross_product3(position_theta_zeta, position_zeta, product)
        jacobian_zeta = jacobian_zeta + dot_product(position_s, product)
        call cross_product3(position_theta, position_zeta_zeta, product)
        jacobian_zeta = jacobian_zeta + dot_product(position_s, product)
    end subroutine build_jacobian_derivatives

    pure subroutine cross_product3(left, right, product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp), intent(out) :: product(3)

        product(1) = left(2) * right(3) - left(3) * right(2)
        product(2) = left(3) * right(1) - left(1) * right(3)
        product(3) = left(1) * right(2) - left(2) * right(1)
    end subroutine cross_product3

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
