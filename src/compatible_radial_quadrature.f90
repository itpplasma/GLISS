module compatible_radial_quadrature
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: compatible_quadrature_ok = 0
    integer, parameter, public :: compatible_quadrature_invalid = -1
    real(dp), parameter, public :: accurate_nodes(5) = &
        [-0.9061798459386640_dp, -0.5384693101056831_dp, 0.0_dp, &
        0.5384693101056831_dp, 0.9061798459386640_dp]
    real(dp), parameter, public :: accurate_weights(5) = &
        [0.2369268850561891_dp, 0.4786286704993665_dp, &
        0.5688888888888889_dp, 0.4786286704993665_dp, &
        0.2369268850561891_dp]

    public :: build_constraint_quadrature

contains

    subroutine build_constraint_quadrature(degree, nodes, weights, info)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: nodes(:), weights(:)
        integer, intent(out) :: info

        info = compatible_quadrature_invalid
        if (degree < 1 .or. degree > 4) return
        allocate (nodes(degree), weights(degree))
        select case (degree)
        case (1)
            nodes = [0.0_dp]
            weights = [2.0_dp]
        case (2)
            nodes = [-0.5773502691896258_dp, 0.5773502691896258_dp]
            weights = [1.0_dp, 1.0_dp]
        case (3)
            nodes = [-0.7745966692414834_dp, 0.0_dp, &
                0.7745966692414834_dp]
            weights = [0.5555555555555556_dp, 0.8888888888888889_dp, &
                0.5555555555555556_dp]
        case (4)
            nodes = [-0.8611363115940526_dp, -0.3399810435848563_dp, &
                0.3399810435848563_dp, 0.8611363115940526_dp]
            weights = [0.3478548451374539_dp, 0.6521451548625461_dp, &
                0.6521451548625461_dp, 0.3478548451374539_dp]
        end select
        info = compatible_quadrature_ok
    end subroutine build_constraint_quadrature

end module compatible_radial_quadrature
