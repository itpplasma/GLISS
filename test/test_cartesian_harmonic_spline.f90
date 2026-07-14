program test_cartesian_harmonic_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cartesian_harmonic_spline, only: cartesian_harmonic_invalid, &
        cartesian_harmonic_ok, cartesian_harmonic_spline_t, &
        cartesian_jet_grid_t, evaluate_cartesian_harmonic_spline, &
        fit_cartesian_harmonic_spline
    use gvec_cas3d_types, only: harmonic_pair_t
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        radial_cubic_spline_grid_t, radial_cubic_spline_ok
    implicit none

    real(dp), parameter :: nodes(6) = [0.04_dp, 0.17_dp, 0.36_dp, &
        0.58_dp, 0.79_dp, 0.96_dp]
    integer, parameter :: modes_m(2) = [0, 1]
    integer, parameter :: modes_n(3) = [-1, 0, 1]
    real(dp), parameter :: major_radius = 3.0_dp, minor_radius = 0.7_dp
    real(dp), parameter :: theta(2) = [0.17_dp, 0.31_dp]
    real(dp), parameter :: zeta(2) = [0.23_dp, 0.44_dp]
    real(dp), parameter :: query_s = 0.36_dp
    type(radial_cubic_spline_grid_t) :: grid
    type(harmonic_pair_t) :: x, y, z
    type(cartesian_harmonic_spline_t) :: spline
    type(cartesian_jet_grid_t) :: jet
    integer :: info

    call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
    call require(info == radial_cubic_spline_ok, "grid construction failed")
    call build_torus_coefficients(x, y, z)
    call fit_cartesian_harmonic_spline(grid, modes_m, modes_n, x, y, z, &
        spline, info)
    call require(info == cartesian_harmonic_ok, "Cartesian fit failed")
    call evaluate_cartesian_harmonic_spline(grid, spline, query_s, theta, &
        zeta, jet, info)
    call require(info == cartesian_harmonic_ok, &
        "Cartesian evaluation failed")
    call check_manufactured_torus(jet)
    call check_invalid_inputs(grid, x, y, z, spline)
    write (*, "(a)") "PASS"

contains

    subroutine build_torus_coefficients(x_pair, y_pair, z_pair)
        type(harmonic_pair_t), intent(out) :: x_pair, y_pair, z_pair
        real(dp) :: radius
        integer :: surface

        allocate (x_pair%cosine(size(nodes), 2, 3), &
            x_pair%sine(size(nodes), 2, 3), &
            y_pair%cosine(size(nodes), 2, 3), &
            y_pair%sine(size(nodes), 2, 3), &
            z_pair%cosine(size(nodes), 2, 3), &
            z_pair%sine(size(nodes), 2, 3))
        x_pair%cosine = 0.0_dp
        x_pair%sine = 0.0_dp
        y_pair%cosine = 0.0_dp
        y_pair%sine = 0.0_dp
        z_pair%cosine = 0.0_dp
        z_pair%sine = 0.0_dp
        do surface = 1, size(nodes)
            radius = minor_radius * sqrt(nodes(surface))
            x_pair%cosine(surface, 1, 3) = major_radius
            x_pair%cosine(surface, 2, 1) = 0.5_dp * radius
            x_pair%cosine(surface, 2, 3) = 0.5_dp * radius
            y_pair%sine(surface, 1, 3) = -major_radius
            y_pair%sine(surface, 2, 1) = 0.5_dp * radius
            y_pair%sine(surface, 2, 3) = -0.5_dp * radius
            z_pair%sine(surface, 2, 2) = radius
        end do
    end subroutine build_torus_coefficients

    subroutine check_manufactured_torus(actual)
        type(cartesian_jet_grid_t), intent(in) :: actual
        real(dp) :: expected(3, 10)
        integer :: j, k

        do k = 1, size(zeta)
            do j = 1, size(theta)
                call torus_jet(query_s, theta(j), zeta(k), expected)
                call require(close_vector(actual%value(j, k, :), &
                    expected(:, 1)), "position differs")
                call require(close_vector(actual%radial(j, k, :), &
                    expected(:, 2)), "radial derivative differs")
                call require(close_vector(actual%poloidal(j, k, :), &
                    expected(:, 3)), "poloidal derivative differs")
                call require(close_vector(actual%toroidal(j, k, :), &
                    expected(:, 4)), "toroidal derivative differs")
                call require(close_vector(actual%radial_radial(j, k, :), &
                    expected(:, 5)), "second radial derivative differs")
                call require(close_vector(actual%radial_poloidal(j, k, :), &
                    expected(:, 6)), "radial-poloidal derivative differs")
                call require(close_vector(actual%radial_toroidal(j, k, :), &
                    expected(:, 7)), "radial-toroidal derivative differs")
                call require(close_vector(actual%poloidal_poloidal(j, k, :), &
                    expected(:, 8)), "second poloidal derivative differs")
                call require(close_vector(actual%poloidal_toroidal(j, k, :), &
                    expected(:, 9)), "mixed angular derivative differs")
                call require(close_vector(actual%toroidal_toroidal(j, k, :), &
                    expected(:, 10)), "second toroidal derivative differs")
            end do
        end do
    end subroutine check_manufactured_torus

    subroutine torus_jet(s, theta_value, zeta_value, expected)
        real(dp), intent(in) :: s, theta_value, zeta_value
        real(dp), intent(out) :: expected(3, 10)
        real(dp) :: angle_t, angle_z, radius, radius_s, radius_ss
        real(dp) :: surface_radius, surface_theta
        real(dp) :: ct, st, cz, sz, two_pi

        two_pi = 2.0_dp * acos(-1.0_dp)
        angle_t = two_pi * theta_value
        angle_z = two_pi * zeta_value
        ct = cos(angle_t)
        st = sin(angle_t)
        cz = cos(angle_z)
        sz = sin(angle_z)
        radius = minor_radius * sqrt(s)
        radius_s = minor_radius / (2.0_dp * sqrt(s))
        radius_ss = -minor_radius / (4.0_dp * s**1.5_dp)
        surface_radius = major_radius + radius * ct
        surface_theta = -two_pi * radius * st
        expected(:, 1) = [surface_radius * cz, surface_radius * sz, &
            radius * st]
        expected(:, 2) = radius_s * [ct * cz, ct * sz, st]
        expected(:, 3) = two_pi * radius * [-st * cz, -st * sz, ct]
        expected(:, 4) = two_pi * surface_radius * [-sz, cz, 0.0_dp]
        expected(:, 5) = radius_ss * [ct * cz, ct * sz, st]
        expected(:, 6) = two_pi * radius_s * [-st * cz, -st * sz, ct]
        expected(:, 7) = two_pi * radius_s * [-ct * sz, ct * cz, 0.0_dp]
        expected(:, 8) = -two_pi**2 * radius * [ct * cz, ct * sz, st]
        expected(:, 9) = two_pi * surface_theta * [-sz, cz, 0.0_dp]
        expected(:, 10) = -two_pi**2 * surface_radius * [cz, sz, 0.0_dp]
    end subroutine torus_jet

    subroutine check_invalid_inputs(valid_grid, valid_x, valid_y, valid_z, &
            valid_spline)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        type(harmonic_pair_t), intent(in) :: valid_x, valid_y, valid_z
        type(cartesian_harmonic_spline_t), intent(in) :: valid_spline
        type(harmonic_pair_t) :: invalid_pair
        type(cartesian_harmonic_spline_t) :: invalid_spline
        type(cartesian_jet_grid_t) :: invalid_jet
        real(dp) :: invalid_theta(1)
        integer :: status

        allocate (invalid_pair%cosine(5, 2, 3), &
            invalid_pair%sine(5, 2, 3))
        invalid_pair%cosine = 0.0_dp
        invalid_pair%sine = 0.0_dp
        call fit_cartesian_harmonic_spline(valid_grid, modes_m, modes_n, &
            valid_x, invalid_pair, valid_z, invalid_spline, status)
        call require(status == cartesian_harmonic_invalid, &
            "mismatched position shape was accepted")
        invalid_pair = valid_y
        invalid_pair%cosine(2, 1, 1) = &
            ieee_value(0.0_dp, ieee_quiet_nan)
        call fit_cartesian_harmonic_spline(valid_grid, modes_m, modes_n, &
            valid_x, invalid_pair, valid_z, invalid_spline, status)
        call require(status == cartesian_harmonic_invalid, &
            "nonfinite position coefficient was accepted")
        call evaluate_cartesian_harmonic_spline(valid_grid, valid_spline, &
            0.0_dp, theta, zeta, invalid_jet, status)
        call require(status == cartesian_harmonic_invalid, &
            "singular axis jet was accepted")
        call require(.not. allocated(invalid_jet%value), &
            "failed axis evaluation retained output")
        invalid_theta = ieee_value(0.0_dp, ieee_quiet_nan)
        call evaluate_cartesian_harmonic_spline(valid_grid, valid_spline, &
            query_s, invalid_theta, zeta, invalid_jet, status)
        call require(status == cartesian_harmonic_invalid, &
            "nonfinite angle was accepted")
        call evaluate_cartesian_harmonic_spline(valid_grid, invalid_spline, &
            query_s, theta, zeta, invalid_jet, status)
        call require(status == cartesian_harmonic_invalid, &
            "invalid spline was accepted")
    end subroutine check_invalid_inputs

    function close_vector(actual, expected) result(matches)
        real(dp), intent(in) :: actual(:), expected(:)
        logical :: matches

        matches = size(actual) == size(expected) &
            .and. all(abs(actual - expected) <= 3.0e-11_dp &
            * max(1.0_dp, abs(expected)))
    end function close_vector

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_cartesian_harmonic_spline
