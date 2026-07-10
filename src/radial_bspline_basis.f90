module radial_bspline_basis
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: bspline_ok = 0
    integer, parameter, public :: bspline_invalid_degree = 1
    integer, parameter, public :: bspline_invalid_knots = 2
    integer, parameter, public :: bspline_outside_domain = 3

    public :: evaluate_bspline_basis
    public :: evaluate_fixed_boundary_basis

contains

    subroutine evaluate_bspline_basis(knots, degree, coordinate, basis, &
            derivative, info)
        real(dp), intent(in) :: knots(:), coordinate
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: basis(:), derivative(:)
        integer, intent(out) :: info
        real(dp), allocatable :: level(:), lower_level(:)
        integer :: basis_count, level_degree

        call validate_inputs(knots, degree, coordinate, basis_count, info)
        if (info /= bspline_ok) return
        allocate (basis(basis_count), derivative(basis_count))
        if (coordinate == knots(size(knots))) then
            call evaluate_right_endpoint(knots, degree, basis, derivative)
            return
        end if
        call initialize_degree_zero(knots, coordinate, level)
        if (degree == 0) then
            basis = level(1:basis_count)
            derivative = 0.0_dp
            return
        end if
        do level_degree = 1, degree
            if (level_degree == degree) lower_level = level
            call elevate_degree(knots, coordinate, level_degree, level)
        end do
        basis = level(1:basis_count)
        call evaluate_derivative(knots, degree, lower_level, derivative)
    end subroutine evaluate_bspline_basis

    subroutine evaluate_fixed_boundary_basis(knots, degree, coordinate, &
            basis, derivative, info)
        real(dp), intent(in) :: knots(:), coordinate
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: basis(:), derivative(:)
        integer, intent(out) :: info
        real(dp), allocatable :: unconstrained(:), unconstrained_derivative(:)
        integer :: constrained_count

        call evaluate_bspline_basis(knots, degree, coordinate, unconstrained, &
            unconstrained_derivative, info)
        if (info /= bspline_ok) return
        constrained_count = size(unconstrained) - 1
        allocate (basis(constrained_count), derivative(constrained_count))
        basis = unconstrained(1:constrained_count)
        derivative = unconstrained_derivative(1:constrained_count)
    end subroutine evaluate_fixed_boundary_basis

    subroutine validate_inputs(knots, degree, coordinate, basis_count, info)
        real(dp), intent(in) :: knots(:), coordinate
        integer, intent(in) :: degree
        integer, intent(out) :: basis_count, info
        integer :: index

        basis_count = 0
        info = bspline_invalid_degree
        if (degree < 0) return
        if (size(knots) < 2 * degree + 2) return
        info = bspline_invalid_knots
        if (.not. all(ieee_is_finite(knots))) return
        if (.not. ieee_is_finite(coordinate)) return
        if (any(knots(2:) < knots(:size(knots) - 1))) return
        if (any(knots(1:degree + 1) /= knots(1))) return
        if (any(knots(size(knots) - degree:) /= knots(size(knots)))) return
        if (knots(degree + 2) <= knots(1)) return
        if (knots(size(knots) - degree - 1) >= knots(size(knots))) return
        do index = degree + 2, size(knots) - degree - 1
            if (knots(index) == knots(index - degree - 1)) return
        end do
        basis_count = size(knots) - degree - 1
        info = bspline_outside_domain
        if (coordinate < knots(1)) return
        if (coordinate > knots(size(knots))) return
        info = bspline_ok
    end subroutine validate_inputs

    subroutine initialize_degree_zero(knots, coordinate, level)
        real(dp), intent(in) :: knots(:), coordinate
        real(dp), allocatable, intent(out) :: level(:)
        integer :: index

        allocate (level(size(knots) - 1))
        level = 0.0_dp
        do index = 1, size(level)
            if (coordinate >= knots(index) .and. &
                coordinate < knots(index + 1)) level(index) = 1.0_dp
        end do
    end subroutine initialize_degree_zero

    subroutine elevate_degree(knots, coordinate, degree, level)
        real(dp), intent(in) :: knots(:), coordinate
        integer, intent(in) :: degree
        real(dp), allocatable, intent(inout) :: level(:)
        real(dp), allocatable :: elevated(:)
        real(dp) :: denominator
        integer :: index, level_size

        level_size = size(knots) - degree - 1
        allocate (elevated(level_size))
        elevated = 0.0_dp
        do index = 1, level_size
            denominator = knots(index + degree) - knots(index)
            if (denominator > 0.0_dp) elevated(index) = &
                (coordinate - knots(index)) * level(index) / denominator
            denominator = knots(index + degree + 1) - knots(index + 1)
            if (denominator > 0.0_dp) elevated(index) = elevated(index) + &
                (knots(index + degree + 1) - coordinate) * &
                level(index + 1) / denominator
        end do
        call move_alloc(elevated, level)
    end subroutine elevate_degree

    subroutine evaluate_derivative(knots, degree, lower_level, derivative)
        real(dp), intent(in) :: knots(:), lower_level(:)
        integer, intent(in) :: degree
        real(dp), intent(out) :: derivative(:)
        real(dp) :: denominator
        integer :: index

        derivative = 0.0_dp
        do index = 1, size(derivative)
            denominator = knots(index + degree) - knots(index)
            if (denominator > 0.0_dp) derivative(index) = &
                real(degree, dp) * lower_level(index) / denominator
            denominator = knots(index + degree + 1) - knots(index + 1)
            if (denominator > 0.0_dp) derivative(index) = derivative(index) - &
                real(degree, dp) * lower_level(index + 1) / denominator
        end do
    end subroutine evaluate_derivative

    subroutine evaluate_right_endpoint(knots, degree, basis, derivative)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), intent(out) :: basis(:), derivative(:)
        real(dp) :: endpoint_span

        basis = 0.0_dp
        derivative = 0.0_dp
        basis(size(basis)) = 1.0_dp
        if (degree == 0) return
        endpoint_span = knots(size(knots)) - knots(size(basis))
        derivative(size(derivative) - 1) = -real(degree, dp) / endpoint_span
        derivative(size(derivative)) = real(degree, dp) / endpoint_span
    end subroutine evaluate_right_endpoint

end module radial_bspline_basis
