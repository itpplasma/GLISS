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

    public :: axis_quadrature_points
    public :: build_axis_quadrature
    public :: build_constraint_quadrature
    public :: build_gauss_legendre

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

    pure function axis_quadrature_points(degree, maximum_poloidal_mode) &
            result(points)
        integer, intent(in) :: degree, maximum_poloidal_mode
        integer :: points

        ! On the axis element a regular trial function of poloidal number m
        ! and degree p is s^(|m|/2) times a polynomial of degree p - 1 in s,
        ! a polynomial of degree |m| + 2p - 2 in u = sqrt(s/h). A product of
        ! two such functions with ds = 2 h u du has degree 2|m| + 4p - 3;
        ! the Gauss rule is exact from |m| + 2p - 1 points. Four further
        ! points absorb the smooth, non-polynomial geometry factors.
        points = max(size(accurate_nodes), &
            max(0, maximum_poloidal_mode) + 2 * degree + 3)
    end function axis_quadrature_points

    subroutine build_axis_quadrature(width, points, nodes, weights, info)
        real(dp), intent(in) :: width
        integer, intent(in) :: points
        real(dp), allocatable, intent(out) :: nodes(:), weights(:)
        integer, intent(out) :: info
        real(dp), allocatable :: x(:), w(:)
        real(dp) :: u
        integer :: point

        ! Rule on the axis element [0, width] in u = sqrt(s / width):
        ! s = width u^2 and ds = 2 width u du. Integrands that are smooth in
        ! sqrt(s), as all regular-axis energy densities are, become smooth
        ! in u; a Gauss rule in s converges only fractionally for them.
        info = compatible_quadrature_invalid
        if (.not. (width > 0.0_dp)) return
        call build_gauss_legendre(points, x, w, info)
        if (info /= compatible_quadrature_ok) return
        allocate (nodes(points), weights(points))
        do point = 1, points
            u = 0.5_dp * (x(point) + 1.0_dp)
            nodes(point) = width * u**2
            weights(point) = width * u * w(point)
        end do
        info = compatible_quadrature_ok
    end subroutine build_axis_quadrature

    subroutine build_gauss_legendre(points, nodes, weights, info)
        integer, intent(in) :: points
        real(dp), allocatable, intent(out) :: nodes(:), weights(:)
        integer, intent(out) :: info
        real(dp) :: derivative, previous, root, step, value
        integer :: iteration, node

        ! Nodes on [-1, 1] by Newton iteration on the Legendre recurrence
        ! from the Chebyshev-type initial guess, ascending order.
        info = compatible_quadrature_invalid
        if (points < 1 .or. points > 256) return
        allocate (nodes(points), weights(points))
        do node = 1, points
            root = cos(acos(-1.0_dp) * (real(node, dp) - 0.25_dp) &
                / (real(points, dp) + 0.5_dp))
            do iteration = 1, 100
                call legendre(points, root, value, previous)
                derivative = real(points, dp) * (root * value - previous) &
                    / (root**2 - 1.0_dp)
                step = value / derivative
                root = root - step
                if (abs(step) <= 4.0_dp * epsilon(1.0_dp)) exit
            end do
            call legendre(points, root, value, previous)
            derivative = real(points, dp) * (root * value - previous) &
                / (root**2 - 1.0_dp)
            nodes(points + 1 - node) = root
            weights(points + 1 - node) = 2.0_dp &
                / ((1.0_dp - root**2) * derivative**2)
        end do
        info = compatible_quadrature_ok
    contains
        pure subroutine legendre(order, x, current, below)
            integer, intent(in) :: order
            real(dp), intent(in) :: x
            real(dp), intent(out) :: current, below
            real(dp) :: older
            integer :: degree

            current = 1.0_dp
            below = 0.0_dp
            do degree = 1, order
                older = below
                below = current
                current = (real(2 * degree - 1, dp) * x * below &
                    - real(degree - 1, dp) * older) / real(degree, dp)
            end do
        end subroutine legendre
    end subroutine build_gauss_legendre

end module compatible_radial_quadrature
