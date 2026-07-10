program test_mercier_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: compute_mercier, mercier_ok, &
        mercier_result_t
    use cylinder_fixture, only: b_poloidal, create_cylinder_fixture, &
        radius_of
    implicit none

    integer, parameter :: ns = 33
    character(len=*), parameter :: fixture = "mercier_cylinder.nc"

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: mu0 = 4.0e-7_dp * pi
    real(dp), parameter :: minor_radius = 0.5_dp
    real(dp), parameter :: period_length = 6.0_dp * pi
    real(dp), parameter :: b_axial = 1.0_dp
    real(dp), parameter :: b_linear = 0.3_dp
    real(dp), parameter :: b_cubic = 0.4_dp

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(mercier_result_t) :: result
    integer :: info, i
    real(dp) :: s_value, ratio, expected_ratio, expected_shear

    call create_cylinder_fixture(fixture)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "cylinder fixture was rejected")
    call compute_mercier(equilibrium, 32, 16, result, info)
    call require(info == mercier_ok, "mercier computation failed")

    do i = 1, ns
        s_value = equilibrium%s(i)
        expected_shear = shear_term_expected()
        call require(abs(result%d_shear(i) - expected_shear) < &
            1.0e-8_dp * expected_shear, "shear term is wrong")
        ratio = result%d_mercier(i) / result%d_shear(i)
        expected_ratio = suydam_ratio_expected(radius_of(s_value))
        call require(abs(ratio - expected_ratio) < 1.0e-3_dp * &
            max(1.0_dp, abs(expected_ratio)), &
            "Mercier does not match the Suydam ratio")
        call require(abs(result%d_geodesic(i)) < 1.0e-6_dp * &
            (abs(result%d_well(i)) + expected_shear), &
            "geodesic term does not vanish")
        call require(result%iota_deviation(i) < 1.0e-10_dp, &
            "iota convention check failed")
        call require(result%boozer_deviation(i) < 1.0e-10_dp, &
            "covariant components are not flux functions")
        call require(result%force_balance_residual(i) < 1.0e-3_dp, &
            "force balance residual is too large")
        call require(result%jacobian_identity_deviation(i) < 1.0e-10_dp, &
            "Boozer Jacobian identity is violated")
        call require(result%beta_chart_deviation(i) < 1.0e-6_dp, &
            "position and spectral beta disagree on the cylinder")
    end do
    call require(any(result%d_mercier < 0.0_dp), &
        "expected Mercier-unstable surfaces were not found")

    open (unit=13, file=fixture, status="old")
    close (13, status="delete")
    write (*, "(a)") "PASS"

contains






    pure function shear_term_expected() result(value)
        real(dp) :: value
        real(dp) :: iota_slope, psi_slope

        iota_slope = period_length * b_cubic * minor_radius**2 &
            / (2.0_dp * pi * b_axial)
        psi_slope = -minor_radius**2 * b_axial / 2.0_dp
        value = (iota_slope / psi_slope)**2 / (16.0_dp * pi**2)
    end function shear_term_expected

    pure function suydam_ratio_expected(radius) result(value)
        real(dp), intent(in) :: radius
        real(dp) :: value
        real(dp) :: pressure_slope, shear_log

        pressure_slope = -b_poloidal(radius) &
            * (2.0_dp * b_linear * radius + 4.0_dp * b_cubic * radius**3) &
            / (mu0 * radius)
        shear_log = 2.0_dp * b_cubic * radius &
            / (b_linear + b_cubic * radius**2)
        value = 1.0_dp + 8.0_dp * mu0 * pressure_slope &
            / (radius * b_axial**2 * shear_log**2)
    end function suydam_ratio_expected



    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require


end program test_mercier_diagnostic
