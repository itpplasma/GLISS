program test_radial_bspline_basis
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use radial_bspline_basis, only: bspline_invalid_degree, &
        bspline_invalid_knots, bspline_ok, bspline_outside_domain, &
        evaluate_bspline_basis, &
        evaluate_fixed_boundary_basis
    implicit none

    real(dp), parameter :: knots(7) = &
        [0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), allocatable :: basis(:), derivative(:)
    real(dp), allocatable :: constrained(:), constrained_derivative(:)
    real(dp) :: coordinate, step
    real(dp), allocatable :: basis_minus(:), basis_plus(:), ignored(:)
    integer :: info

    call evaluate_bspline_basis(knots, 2, 0.25_dp, basis, derivative, info)
    call require(info == bspline_ok, "quadratic basis evaluation failed")
    call require(maxval(abs(basis - [0.25_dp, 0.625_dp, 0.125_dp, &
        0.0_dp])) < 1.0e-14_dp, "quadratic basis values are wrong")
    call require(maxval(abs(derivative - [-2.0_dp, 1.0_dp, 1.0_dp, &
        0.0_dp])) < 1.0e-14_dp, "quadratic derivatives are wrong")

    call verify_basis_identities([0.0_dp, 1.0_dp], 0)
    call verify_basis_identities([0.0_dp, 0.0_dp, 0.3_dp, 1.0_dp, 1.0_dp], 1)
    call verify_basis_identities(knots, 2)
    call verify_basis_identities([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, &
        0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 3)

    call evaluate_bspline_basis(knots, 2, 0.0_dp, basis, derivative, info)
    call require(all(basis == [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]), &
        "axis endpoint interpolation failed")
    call evaluate_bspline_basis(knots, 2, 1.0_dp, basis, derivative, info)
    call require(all(basis == [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]), &
        "boundary endpoint interpolation failed")
    call require(all(derivative == [0.0_dp, 0.0_dp, -4.0_dp, 4.0_dp]), &
        "boundary endpoint derivative failed")

    call evaluate_fixed_boundary_basis(knots, 2, 1.0_dp, constrained, &
        constrained_derivative, info)
    call require(info == bspline_ok, "fixed-boundary basis failed")
    call require(all(constrained == 0.0_dp), &
        "fixed-boundary basis does not vanish at the surface")

    coordinate = 0.37_dp
    step = 1.0e-6_dp
    call evaluate_bspline_basis(knots, 2, coordinate, basis, derivative, info)
    call evaluate_bspline_basis(knots, 2, coordinate - step, basis_minus, &
        ignored, info)
    call evaluate_bspline_basis(knots, 2, coordinate + step, basis_plus, &
        ignored, info)
    call require(maxval(abs(derivative - (basis_plus - basis_minus) / &
        (2.0_dp * step))) < 1.0e-9_dp, "derivative finite difference failed")

    call evaluate_bspline_basis([0.0_dp, 0.0_dp, 0.5_dp, 0.4_dp, 1.0_dp, &
        1.0_dp], 1, 0.5_dp, basis, derivative, info)
    call require(info == bspline_invalid_knots, &
        "decreasing knot vector was accepted")
    call evaluate_bspline_basis(knots, 2, 1.1_dp, basis, derivative, info)
    call require(info == bspline_outside_domain, &
        "out-of-domain coordinate was accepted")
    call evaluate_bspline_basis(knots, -1, 0.5_dp, basis, derivative, info)
    call require(info == bspline_invalid_degree, &
        "negative degree was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine verify_basis_identities(knot_vector, degree)
        real(dp), intent(in) :: knot_vector(:)
        integer, intent(in) :: degree
        real(dp), allocatable :: values(:), derivatives(:)
        real(dp) :: sample_coordinate
        integer :: sample, status

        do sample = 0, 16
            sample_coordinate = knot_vector(1) + &
                (knot_vector(size(knot_vector)) - knot_vector(1)) * &
                real(sample, dp) / 16.0_dp
            call evaluate_bspline_basis(knot_vector, degree, &
                sample_coordinate, values, derivatives, status)
            call require(status == bspline_ok, "basis sample was rejected")
            call require(abs(sum(values) - 1.0_dp) < 1.0e-14_dp, &
                "partition of unity failed")
            call require(abs(sum(derivatives)) < 1.0e-13_dp, &
                "derivative partition failed")
            call require(minval(values) >= 0.0_dp, "basis became negative")
            call require(count(values > 0.0_dp) <= degree + 1, &
                "local support was lost")
        end do
    end subroutine verify_basis_identities

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_radial_bspline_basis
