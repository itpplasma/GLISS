program test_field_periodic_cartesian
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cartesian_harmonic_spline, only: cartesian_jet_grid_t
    use field_periodic_cartesian, only: convert_field_periodic_jet, &
        field_periodic_cartesian_invalid, field_periodic_cartesian_ok
    implicit none

    type(cartesian_jet_grid_t) :: jet
    real(dp), parameter :: zeta(1) = [1.0_dp]
    integer :: info

    call build_polynomial_jet(jet)
    call convert_field_periodic_jet(zeta, 4, 1, jet, info)
    call require(info == field_periodic_cartesian_ok, &
        "valid rotating-frame jet was rejected")
    call check_exact_jet(jet)
    call check_zero_winding()
    call check_invalid_inputs()
    write (*, "(a)") "PASS"

contains

    subroutine build_polynomial_jet(actual)
        type(cartesian_jet_grid_t), intent(out) :: actual

        call allocate_jet(actual, 1, 1)
        actual%value(1, 1, :) = [11.0_dp, 29.0_dp, 41.0_dp]
        actual%radial(1, 1, :) = [1.0_dp, 12.0_dp, 6.0_dp]
        actual%poloidal(1, 1, :) = [2.0_dp, 8.0_dp, 7.0_dp]
        actual%toroidal(1, 1, :) = [6.0_dp, 5.0_dp, 37.0_dp]
        actual%radial_radial(1, 1, :) = 0.0_dp
        actual%radial_poloidal(1, 1, :) = [0.0_dp, 4.0_dp, 0.0_dp]
        actual%radial_toroidal(1, 1, :) = 0.0_dp
        actual%poloidal_poloidal(1, 1, :) = 0.0_dp
        actual%poloidal_toroidal(1, 1, :) = [0.0_dp, 0.0_dp, 7.0_dp]
        actual%toroidal_toroidal(1, 1, :) = [6.0_dp, 0.0_dp, 16.0_dp]
    end subroutine build_polynomial_jet

    subroutine check_exact_jet(actual)
        type(cartesian_jet_grid_t), intent(in) :: actual
        real(dp) :: expected(3), q

        q = acos(-1.0_dp) / 2.0_dp
        call require(close(actual%value(1, 1, :), &
            [-29.0_dp, 11.0_dp, 41.0_dp]), "position rotation differs")
        call require(close(actual%radial(1, 1, :), &
            [-12.0_dp, 1.0_dp, 6.0_dp]), "radial rotation differs")
        call require(close(actual%poloidal(1, 1, :), &
            [-8.0_dp, 2.0_dp, 7.0_dp]), "poloidal rotation differs")
        expected(1) = -5.0_dp - 11.0_dp * q
        expected(2) = 6.0_dp - 29.0_dp * q
        expected(3) = 37.0_dp
        call require(close(actual%toroidal(1, 1, :), expected), &
            "toroidal derivative rotation differs")
        call require(close(actual%radial_radial(1, 1, :), &
            [0.0_dp, 0.0_dp, 0.0_dp]), "second radial rotation differs")
        call require(close(actual%radial_poloidal(1, 1, :), &
            [-4.0_dp, 0.0_dp, 0.0_dp]), &
            "radial-poloidal rotation differs")
        expected(1) = -q
        expected(2) = -12.0_dp * q
        expected(3) = 0.0_dp
        call require(close(actual%radial_toroidal(1, 1, :), expected), &
            "radial-toroidal rotation differs")
        call require(close(actual%poloidal_poloidal(1, 1, :), &
            [0.0_dp, 0.0_dp, 0.0_dp]), "second poloidal rotation differs")
        expected(1) = -2.0_dp * q
        expected(2) = -8.0_dp * q
        expected(3) = 7.0_dp
        call require(close(actual%poloidal_toroidal(1, 1, :), expected), &
            "poloidal-toroidal rotation differs")
        expected(1) = -12.0_dp * q + 29.0_dp * q**2
        expected(2) = 6.0_dp - 10.0_dp * q - 11.0_dp * q**2
        expected(3) = 16.0_dp
        call require(close(actual%toroidal_toroidal(1, 1, :), expected), &
            "second toroidal rotation differs")
    end subroutine check_exact_jet

    subroutine check_zero_winding()
        type(cartesian_jet_grid_t) :: actual, expected
        integer :: status

        call build_polynomial_jet(actual)
        expected = actual
        call convert_field_periodic_jet(zeta, 5, 0, actual, status)
        call require(status == field_periodic_cartesian_ok, &
            "zero winding was rejected")
        call require(jets_close(actual, expected), &
            "zero winding changed the Cartesian jet")
    end subroutine check_zero_winding

    subroutine check_invalid_inputs()
        type(cartesian_jet_grid_t) :: actual, expected
        real(dp), allocatable :: empty(:)
        real(dp) :: invalid_zeta(1)
        integer :: status

        call build_polynomial_jet(actual)
        expected = actual
        call convert_field_periodic_jet(zeta, 0, 1, actual, status)
        call require(status == field_periodic_cartesian_invalid &
            .and. jets_close(actual, expected), &
            "zero field periods did not fail without mutation")
        allocate (empty(0))
        call convert_field_periodic_jet(empty, 5, 1, actual, status)
        call require(status == field_periodic_cartesian_invalid &
            .and. jets_close(actual, expected), &
            "empty toroidal grid did not fail without mutation")
        invalid_zeta = ieee_value(0.0_dp, ieee_quiet_nan)
        call convert_field_periodic_jet(invalid_zeta, 5, 1, actual, status)
        call require(status == field_periodic_cartesian_invalid &
            .and. jets_close(actual, expected), &
            "nonfinite toroidal grid did not fail without mutation")
        actual%value(1, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call convert_field_periodic_jet(zeta, 5, 1, actual, status)
        call require(status == field_periodic_cartesian_invalid, &
            "nonfinite jet was accepted")
        call build_polynomial_jet(actual)
        deallocate (actual%radial_toroidal)
        call convert_field_periodic_jet(zeta, 5, 1, actual, status)
        call require(status == field_periodic_cartesian_invalid, &
            "incomplete jet was accepted")
    end subroutine check_invalid_inputs

    subroutine allocate_jet(actual, n_theta, n_zeta)
        type(cartesian_jet_grid_t), intent(out) :: actual
        integer, intent(in) :: n_theta, n_zeta

        allocate (actual%value(n_theta, n_zeta, 3), &
            actual%radial(n_theta, n_zeta, 3), &
            actual%poloidal(n_theta, n_zeta, 3), &
            actual%toroidal(n_theta, n_zeta, 3), &
            actual%radial_radial(n_theta, n_zeta, 3), &
            actual%radial_poloidal(n_theta, n_zeta, 3), &
            actual%radial_toroidal(n_theta, n_zeta, 3), &
            actual%poloidal_poloidal(n_theta, n_zeta, 3), &
            actual%poloidal_toroidal(n_theta, n_zeta, 3), &
            actual%toroidal_toroidal(n_theta, n_zeta, 3))
    end subroutine allocate_jet

    function jets_close(actual, expected) result(matches)
        type(cartesian_jet_grid_t), intent(in) :: actual, expected
        logical :: matches

        matches = close(actual%value, expected%value) &
            .and. close(actual%radial, expected%radial) &
            .and. close(actual%poloidal, expected%poloidal) &
            .and. close(actual%toroidal, expected%toroidal) &
            .and. close(actual%radial_radial, expected%radial_radial) &
            .and. close(actual%radial_poloidal, expected%radial_poloidal) &
            .and. close(actual%radial_toroidal, expected%radial_toroidal) &
            .and. close(actual%poloidal_poloidal, &
            expected%poloidal_poloidal) &
            .and. close(actual%poloidal_toroidal, &
            expected%poloidal_toroidal) &
            .and. close(actual%toroidal_toroidal, &
            expected%toroidal_toroidal)
    end function jets_close

    function close(actual, expected) result(matches)
        real(dp), intent(in) :: actual(..), expected(..)
        logical :: matches

        select rank (actual)
            rank (1)
            select rank (expected)
                rank (1)
                matches = all(shape(actual) == shape(expected)) &
                    .and. all(abs(actual - expected) <= 3.0e-12_dp &
                    * max(1.0_dp, abs(expected)))
                rank default
                matches = .false.
            end select
            rank (3)
            select rank (expected)
                rank (3)
                matches = all(shape(actual) == shape(expected)) &
                    .and. all(abs(actual - expected) <= 3.0e-12_dp &
                    * max(1.0_dp, abs(expected)))
                rank default
                matches = .false.
            end select
            rank default
            matches = .false.
        end select
    end function close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_field_periodic_cartesian
