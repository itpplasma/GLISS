program test_mercier_gradient
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use mercier_diagnostic, only: mercier_d_terms, mercier_d_terms_gradient, &
        mercier_d_terms_t, mercier_gradient_t
    implicit none

    real(dp), parameter :: h = 1.0e-7_dp
    real(dp), parameter :: tol = 1.0e-8_dp

    call check_iota_slope_gradient()
    call check_pressure_slope_gradient()
    call check_zero_structure()
    write (*, "(a)") "PASS"

contains

    pure function d_mercier_at(iota_slope, pressure_slope, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared) result(value)
        real(dp), intent(in) :: iota_slope, pressure_slope, psi_slope
        real(dp), intent(in) :: covariant_zeta, integral_xi, d2v_dpsi2
        real(dp), intent(in) :: integral_inverse, integral_bsq
        real(dp), intent(in) :: integral_jb, integral_jb_squared
        real(dp) :: value
        type(mercier_d_terms_t) :: d_terms

        d_terms = mercier_d_terms(iota_slope, pressure_slope, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared)
        value = d_terms%mercier
    end function d_mercier_at

    subroutine check_iota_slope_gradient()
        real(dp), parameter :: pressure_slope = -0.9_dp, psi_slope = 1.3_dp
        real(dp), parameter :: covariant_zeta = -2.1_dp, integral_xi = 0.6_dp
        real(dp), parameter :: d2v_dpsi2 = 1.7_dp, integral_inverse = 0.4_dp
        real(dp), parameter :: integral_bsq = 2.5_dp, integral_jb = 1.1_dp
        real(dp), parameter :: integral_jb_squared = 0.8_dp
        real(dp), parameter :: iota_slope = 0.7_dp
        real(dp) :: fd
        type(mercier_gradient_t) :: grad

        fd = (d_mercier_at(iota_slope + h, pressure_slope, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared) &
            - d_mercier_at(iota_slope - h, pressure_slope, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared)) / (2.0_dp * h)
        grad = mercier_d_terms_gradient(iota_slope, pressure_slope, &
            psi_slope, covariant_zeta, integral_xi, d2v_dpsi2, &
            integral_inverse, integral_bsq)
        call require(abs(fd - grad%d_iota_slope) <= &
            tol * abs(grad%d_iota_slope), &
            "iota_slope gradient disagrees with finite difference")
    end subroutine check_iota_slope_gradient

    subroutine check_pressure_slope_gradient()
        real(dp), parameter :: iota_slope = 0.0_dp, psi_slope = 1.3_dp
        real(dp), parameter :: covariant_zeta = -2.1_dp, integral_xi = 0.6_dp
        real(dp), parameter :: d2v_dpsi2 = 1.7_dp, integral_inverse = 0.4_dp
        real(dp), parameter :: integral_bsq = 2.5_dp, integral_jb = 0.0_dp
        real(dp), parameter :: integral_jb_squared = 0.0_dp
        real(dp), parameter :: pressure_slope = -0.9_dp
        real(dp) :: fd
        type(mercier_gradient_t) :: grad

        fd = (d_mercier_at(iota_slope, pressure_slope + h, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared) &
            - d_mercier_at(iota_slope, pressure_slope - h, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared)) / (2.0_dp * h)
        grad = mercier_d_terms_gradient(iota_slope, pressure_slope, &
            psi_slope, covariant_zeta, integral_xi, d2v_dpsi2, &
            integral_inverse, integral_bsq)
        call require(abs(fd - grad%d_pressure_slope) <= &
            tol * abs(grad%d_pressure_slope), &
            "pressure_slope gradient disagrees with finite difference")
    end subroutine check_pressure_slope_gradient

    subroutine check_zero_structure()
        real(dp), parameter :: psi_slope = 1.3_dp, covariant_zeta = -2.1_dp
        real(dp), parameter :: integral_xi = 0.6_dp, d2v_dpsi2 = 1.7_dp
        real(dp), parameter :: integral_inverse = 0.4_dp
        type(mercier_gradient_t) :: no_well, low_iota, high_iota

        no_well = mercier_d_terms_gradient(0.7_dp, -0.9_dp, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            0.0_dp)
        call require(no_well%d_pressure_slope == 0.0_dp, &
            "pressure gradient must vanish when integral_bsq is zero")

        low_iota = mercier_d_terms_gradient(0.7_dp, -0.9_dp, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            2.5_dp)
        high_iota = mercier_d_terms_gradient(4.2_dp, -0.9_dp, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            2.5_dp)
        call require(low_iota%d_pressure_slope == high_iota%d_pressure_slope, &
            "pressure gradient must not depend on iota_slope")
    end subroutine check_zero_structure

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_mercier_gradient
