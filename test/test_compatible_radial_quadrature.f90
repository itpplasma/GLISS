program test_compatible_radial_quadrature
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_radial_quadrature, only: accurate_nodes, &
        accurate_weights, build_constraint_quadrature, &
        compatible_quadrature_invalid, compatible_quadrature_ok
    implicit none

    integer :: degree

    call verify_rule(accurate_nodes, accurate_weights, 9)
    do degree = 1, 4
        call verify_constraint_rule(degree)
    end do
    call verify_rejections()
    write (*, "(a)") "PASS"

contains

    subroutine verify_constraint_rule(degree)
        integer, intent(in) :: degree
        real(dp), allocatable :: nodes(:), weights(:)
        integer :: info

        call build_constraint_quadrature(degree, nodes, weights, info)
        call require(info == compatible_quadrature_ok, &
            "valid constraint quadrature was rejected")
        call require(size(nodes) == degree .and. size(weights) == degree, &
            "constraint quadrature size differs")
        call verify_rule(nodes, weights, 2 * degree - 1)
    end subroutine verify_constraint_rule

    subroutine verify_rule(nodes, weights, exact_degree)
        real(dp), intent(in) :: nodes(:), weights(:)
        integer, intent(in) :: exact_degree
        real(dp) :: actual, expected
        integer :: power

        call require(all(weights > 0.0_dp), "quadrature weight is not positive")
        call require(maxval(abs(nodes + nodes(size(nodes):1:-1))) &
            < 2.0e-15_dp, "quadrature nodes are not symmetric")
        do power = 0, exact_degree
            actual = sum(weights * nodes**power)
            if (modulo(power, 2) == 0) then
                expected = 2.0_dp / real(power + 1, dp)
            else
                expected = 0.0_dp
            end if
            call require(abs(actual - expected) < 3.0e-15_dp, &
                "quadrature polynomial exactness failed")
        end do
    end subroutine verify_rule

    subroutine verify_rejections()
        real(dp), allocatable :: nodes(:), weights(:)
        integer :: info

        call build_constraint_quadrature(0, nodes, weights, info)
        call require(info == compatible_quadrature_invalid &
            .and. .not. allocated(nodes) .and. .not. allocated(weights), &
            "degree zero quadrature was accepted")
        call build_constraint_quadrature(5, nodes, weights, info)
        call require(info == compatible_quadrature_invalid &
            .and. .not. allocated(nodes) .and. .not. allocated(weights), &
            "degree five quadrature was accepted")
    end subroutine verify_rejections

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_compatible_radial_quadrature
