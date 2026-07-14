program test_cartesian_primitive_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cartesian_primitive_geometry, only: build_primitive_geometry_point, &
        primitive_geometry_invalid, primitive_geometry_ok, &
        primitive_geometry_point_t
    implicit none

    real(dp), parameter :: rs(3) = [1.0_dp, 2.0_dp, -1.0_dp]
    real(dp), parameter :: rt(3) = [0.0_dp, 2.0_dp, 1.0_dp]
    real(dp), parameter :: rz(3) = [3.0_dp, -1.0_dp, 2.0_dp]
    real(dp), parameter :: rss(3) = [2.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter :: rst(3) = [1.0_dp, -1.0_dp, 0.0_dp]
    real(dp), parameter :: rsz(3) = [0.0_dp, 2.0_dp, -1.0_dp]
    real(dp), parameter :: rtt(3) = [1.0_dp, 0.0_dp, 2.0_dp]
    real(dp), parameter :: rtz(3) = [-1.0_dp, 2.0_dp, 1.0_dp]
    real(dp), parameter :: rzz(3) = [2.0_dp, -2.0_dp, 0.0_dp]
    real(dp), parameter :: phi_slope = 7.0_dp, chi_slope = -3.0_dp
    integer, parameter :: periods = 5
    type(primitive_geometry_point_t) :: geometry
    integer :: info

    call build_fixture(rs, rt, rz, rss, rst, rsz, rtt, rtz, rzz, &
        phi_slope, chi_slope, geometry, info)
    call require(info == primitive_geometry_ok, "point construction failed")
    call check_mathematica_fixture(geometry)
    call check_invariants(geometry)
    call check_orientation_change(geometry)
    call check_second_form_tangent_invariance(geometry)
    call check_invalid_inputs()
    write (*, "(a)") "PASS"

contains

    subroutine check_mathematica_fixture(actual)
        type(primitive_geometry_point_t), intent(in) :: actual
        real(dp), parameter :: expected_metric(3, 3) = reshape([ &
            6.0_dp, 3.0_dp, -1.0_dp, 3.0_dp, 5.0_dp, 0.0_dp, &
            -1.0_dp, 0.0_dp, 14.0_dp], [3, 3])
        real(dp), parameter :: expected_second(2, 2) = reshape([ &
            -7.0_dp, -5.0_dp, -5.0_dp, 4.0_dp], [2, 2]) / sqrt(70.0_dp)

        call require(close_matrix(actual%metric, expected_metric), &
            "Mathematica metric fixture differs")
        call require(close(actual%signed_jacobian, 17.0_dp), &
            "Mathematica Jacobian fixture differs")
        call require(close(actual%jacobian_s, 7.0_dp) &
            .and. close(actual%jacobian_theta, 9.0_dp) &
            .and. close(actual%jacobian_zeta, 42.0_dp), &
            "Mathematica Jacobian derivative fixture differs")
        call require(close_vector(actual%b_contravariant, &
            [3.0_dp / 85.0_dp, -7.0_dp / 17.0_dp]), &
            "Mathematica contravariant-field fixture differs")
        call require(close_vector(actual%b_covariant, &
            [3.0_dp / 17.0_dp, -98.0_dp / 17.0_dp]), &
            "Mathematica covariant-field fixture differs")
        call require(close(actual%mod_b, sqrt(3439.0_dp / 1445.0_dp)), &
            "Mathematica magnetic-magnitude fixture differs")
        call require(close_matrix(actual%second_form, expected_second), &
            "Mathematica second-form fixture differs")
    end subroutine check_mathematica_fixture

    subroutine check_invariants(actual)
        type(primitive_geometry_point_t), intent(in) :: actual
        real(dp) :: determinant, field_squared

        determinant = actual%metric(1, 1) * (actual%metric(2, 2) &
            * actual%metric(3, 3) - actual%metric(2, 3)**2) &
            - actual%metric(1, 2) * (actual%metric(1, 2) &
            * actual%metric(3, 3) - actual%metric(1, 3) &
            * actual%metric(2, 3)) + actual%metric(1, 3) &
            * (actual%metric(1, 2) * actual%metric(2, 3) &
            - actual%metric(1, 3) * actual%metric(2, 2))
        field_squared = dot_product(actual%b_contravariant, &
            actual%b_covariant)
        call require(all(actual%metric == transpose(actual%metric)), &
            "metric is not exactly symmetric")
        call require(close(determinant, actual%signed_jacobian**2), &
            "metric determinant is not Jacobian squared")
        call require(close(-actual%signed_jacobian &
            * actual%b_contravariant(2), phi_slope), &
            "toroidal flux convention differs")
        call require(close(-real(periods, dp) * actual%signed_jacobian &
            * actual%b_contravariant(1), chi_slope), &
            "poloidal flux convention differs")
        call require(close(field_squared, actual%mod_b**2) &
            .and. actual%mod_b > 0.0_dp, &
            "magnetic magnitude identity differs")
        call require(actual%second_form(1, 2) &
            == actual%second_form(2, 1), "second form is not symmetric")
    end subroutine check_invariants

    subroutine check_orientation_change(reference)
        type(primitive_geometry_point_t), intent(in) :: reference
        type(primitive_geometry_point_t) :: swapped
        integer :: status

        call build_fixture(rs, rz, rt, rss, rsz, rst, rzz, rtz, rtt, &
            -chi_slope / real(periods, dp), &
            -real(periods, dp) * phi_slope, swapped, status)
        call require(status == primitive_geometry_ok, &
            "orientation-swapped point failed")
        call require(swapped%signed_jacobian == -reference%signed_jacobian, &
            "orientation swap did not reverse the Jacobian")
        call require(close(swapped%jacobian_s, -reference%jacobian_s) &
            .and. close(swapped%jacobian_theta, -reference%jacobian_zeta) &
            .and. close(swapped%jacobian_zeta, -reference%jacobian_theta), &
            "orientation-swapped Jacobian derivatives differ")
        call require(close(swapped%mod_b, reference%mod_b), &
            "orientation swap changed magnetic magnitude")
        call require(close(swapped%second_form(1, 1), &
            reference%second_form(2, 2)) &
            .and. close(swapped%second_form(1, 2), &
            reference%second_form(1, 2)) &
            .and. close(swapped%second_form(2, 2), &
            reference%second_form(1, 1)), &
            "orientation-swapped second form differs")
    end subroutine check_orientation_change

    subroutine check_second_form_tangent_invariance(reference)
        type(primitive_geometry_point_t), intent(in) :: reference
        type(primitive_geometry_point_t) :: shifted
        integer :: status

        call build_fixture(rs, rt, rz, rss, rst, rsz, &
            rtt + 2.0_dp * rt - 3.0_dp * rz, &
            rtz - rt + 4.0_dp * rz, rzz + 5.0_dp * rt + rz, &
            phi_slope, chi_slope, shifted, status)
        call require(status == primitive_geometry_ok, &
            "tangent-shifted point failed")
        call require(close_matrix(shifted%second_form, &
            reference%second_form), &
            "tangential second derivatives changed the second form")
    end subroutine check_second_form_tangent_invariance

    subroutine check_invalid_inputs()
        type(primitive_geometry_point_t) :: invalid
        real(dp) :: nonfinite(3)
        integer :: status

        call build_fixture(rs, rt, rz, rss, rst, rsz, rtt, rtz, rzz, &
            phi_slope, chi_slope, invalid, status, field_periods=0)
        call require(status == primitive_geometry_invalid, &
            "zero field periods were accepted")
        nonfinite = rs
        nonfinite(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        call build_fixture(nonfinite, rt, rz, rss, rst, rsz, rtt, rtz, rzz, &
            phi_slope, chi_slope, invalid, status)
        call require(status == primitive_geometry_invalid, &
            "nonfinite tangent was accepted")
        call build_fixture(rs, rt, 2.0_dp * rt, rss, rst, rsz, rtt, rtz, &
            rzz, phi_slope, chi_slope, invalid, status)
        call require(status == primitive_geometry_invalid, &
            "singular surface frame was accepted")
        call build_fixture(rt + rz, rt, rz, rss, rst, rsz, rtt, rtz, rzz, &
            phi_slope, chi_slope, invalid, status)
        call require(status == primitive_geometry_invalid, &
            "zero signed Jacobian was accepted")
        call build_fixture(rs, rt, rz, rss, rst, rsz, rtt, rtz, rzz, &
            0.0_dp, 0.0_dp, invalid, status)
        call require(status == primitive_geometry_invalid, &
            "zero magnetic field was accepted")
        call require(all(invalid%metric == 0.0_dp) &
            .and. invalid%signed_jacobian == 0.0_dp &
            .and. invalid%mod_b == 0.0_dp, &
            "failed construction returned partial geometry")
    end subroutine check_invalid_inputs

    subroutine build_fixture(local_rs, local_rt, local_rz, local_rss, &
            local_rst, local_rsz, local_rtt, local_rtz, local_rzz, &
            local_phi, local_chi, result, status, field_periods)
        real(dp), intent(in) :: local_rs(3), local_rt(3), local_rz(3)
        real(dp), intent(in) :: local_rss(3), local_rst(3), local_rsz(3)
        real(dp), intent(in) :: local_rtt(3), local_rtz(3), local_rzz(3)
        real(dp), intent(in) :: local_phi, local_chi
        type(primitive_geometry_point_t), intent(out) :: result
        integer, intent(out) :: status
        integer, intent(in), optional :: field_periods
        integer :: local_periods

        local_periods = periods
        if (present(field_periods)) local_periods = field_periods
        call build_primitive_geometry_point(local_rs, local_rt, local_rz, &
            local_rss, local_rst, local_rsz, local_rtt, local_rtz, &
            local_rzz, local_periods, local_phi, local_chi, result, status)
    end subroutine build_fixture

    elemental function close(actual, expected) result(matches)
        real(dp), intent(in) :: actual, expected
        logical :: matches

        matches = abs(actual - expected) &
            <= 5.0e-14_dp * max(1.0_dp, abs(expected))
    end function close

    function close_vector(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:), expected(:)
        logical :: matches

        matches = size(actual) == size(expected) &
            .and. all(close(actual, expected))
    end function close_vector

    function close_matrix(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:, :), expected(:, :)
        logical :: matches

        matches = all(shape(actual) == shape(expected)) &
            .and. all(close(actual, expected))
    end function close_matrix

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_cartesian_primitive_geometry
