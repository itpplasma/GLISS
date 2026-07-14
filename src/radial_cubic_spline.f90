module radial_cubic_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: radial_cubic_spline_ok = 0
    integer, parameter, public :: radial_cubic_spline_invalid = -1
    integer, parameter, public :: radial_cubic_spline_allocation_error = -2

    type, public :: radial_cubic_spline_grid_t
        real(dp) :: domain_min = 0.0_dp
        real(dp) :: domain_max = 0.0_dp
        real(dp), allocatable :: nodes(:)
        real(dp), allocatable :: intervals(:)
        real(dp), allocatable :: lower_factor(:)
        real(dp), allocatable :: diagonal_factor(:)
        real(dp), allocatable :: upper(:)
    end type radial_cubic_spline_grid_t

    type, public :: radial_cubic_spline_t
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: second_derivatives(:)
    end type radial_cubic_spline_t

    public :: build_radial_cubic_spline_grid
    public :: evaluate_radial_cubic_spline
    public :: fit_radial_cubic_spline

contains

    subroutine build_radial_cubic_spline_grid(nodes, domain_min, domain_max, &
            grid, info)
        real(dp), intent(in) :: nodes(:), domain_min, domain_max
        type(radial_cubic_spline_grid_t), intent(out) :: grid
        integer, intent(out) :: info
        real(dp), allocatable :: diagonal(:), lower(:), upper(:)
        integer :: allocation_status, interior

        grid = radial_cubic_spline_grid_t()
        info = radial_cubic_spline_invalid
        if (size(nodes) < 4) return
        if (.not. all(ieee_is_finite(nodes))) return
        if (.not. ieee_is_finite(domain_min)) return
        if (.not. ieee_is_finite(domain_max)) return
        if (domain_min >= domain_max) return
        if (any(nodes(2:) <= nodes(:size(nodes) - 1))) return
        if (nodes(1) < domain_min .or. nodes(size(nodes)) > domain_max) return
        interior = size(nodes) - 2
        allocate (grid%nodes(size(nodes)), grid%intervals(size(nodes) - 1), &
            diagonal(interior), lower(interior - 1), upper(interior - 1), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = radial_cubic_spline_allocation_error
            return
        end if
        grid%domain_min = domain_min
        grid%domain_max = domain_max
        grid%nodes = nodes
        grid%intervals = nodes(2:) - nodes(:size(nodes) - 1)
        call build_reduced_system(grid%intervals, lower, diagonal, upper)
        call factor_tridiagonal(lower, diagonal, upper, info)
        if (info /= radial_cubic_spline_ok) return
        call move_alloc(lower, grid%lower_factor)
        call move_alloc(diagonal, grid%diagonal_factor)
        call move_alloc(upper, grid%upper)
        if (.not. grid_is_valid(grid)) then
            info = radial_cubic_spline_invalid
            return
        end if
        info = radial_cubic_spline_ok
    end subroutine build_radial_cubic_spline_grid

    subroutine fit_radial_cubic_spline(grid, values, spline, info)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        real(dp), intent(in) :: values(:)
        type(radial_cubic_spline_t), intent(out) :: spline
        integer, intent(out) :: info
        real(dp), allocatable :: right_hand_side(:)
        real(dp) :: first, last
        integer :: allocation_status, i, interior

        spline = radial_cubic_spline_t()
        info = radial_cubic_spline_invalid
        if (.not. grid_is_valid(grid)) return
        if (size(values) /= size(grid%nodes)) return
        if (.not. all(ieee_is_finite(values))) return
        interior = size(values) - 2
        allocate (spline%values(size(values)), &
            spline%second_derivatives(size(values)), &
            right_hand_side(interior), stat=allocation_status)
        if (allocation_status /= 0) then
            info = radial_cubic_spline_allocation_error
            return
        end if
        spline%values = values
        do i = 1, interior
            right_hand_side(i) = 6.0_dp * ((values(i + 2) - values(i + 1)) &
                / grid%intervals(i + 1) - (values(i + 1) - values(i)) &
                / grid%intervals(i))
        end do
        call solve_factored_system(grid, right_hand_side)
        spline%second_derivatives(2:size(values) - 1) = right_hand_side
        first = ((grid%intervals(1) + grid%intervals(2)) &
            * right_hand_side(1) - grid%intervals(1) &
            * right_hand_side(2)) / grid%intervals(2)
        last = ((grid%intervals(size(values) - 2) &
            + grid%intervals(size(values) - 1)) * right_hand_side(interior) &
            - grid%intervals(size(values) - 1) &
            * right_hand_side(interior - 1)) &
            / grid%intervals(size(values) - 2)
        spline%second_derivatives(1) = first
        spline%second_derivatives(size(values)) = last
        if (.not. spline_is_valid(grid, spline)) return
        info = radial_cubic_spline_ok
    end subroutine fit_radial_cubic_spline

    subroutine evaluate_radial_cubic_spline(grid, spline, coordinate, value, &
            derivative, info, second_derivative)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        type(radial_cubic_spline_t), intent(in) :: spline
        real(dp), intent(in) :: coordinate
        real(dp), intent(out) :: value, derivative
        integer, intent(out) :: info
        real(dp), intent(out), optional :: second_derivative
        real(dp) :: distance_left, distance_right, interval
        integer :: left

        info = radial_cubic_spline_invalid
        if (.not. grid_is_valid(grid)) return
        if (.not. spline_is_valid(grid, spline)) return
        if (.not. ieee_is_finite(coordinate)) return
        if (coordinate < grid%domain_min .or. coordinate > grid%domain_max) &
            return
        left = find_interval(grid%nodes, coordinate)
        interval = grid%intervals(left)
        distance_left = coordinate - grid%nodes(left)
        distance_right = grid%nodes(left + 1) - coordinate
        value = spline%second_derivatives(left) * distance_right**3 &
            / (6.0_dp * interval) &
            + spline%second_derivatives(left + 1) * distance_left**3 &
            / (6.0_dp * interval) &
            + (spline%values(left) - spline%second_derivatives(left) &
            * interval**2 / 6.0_dp) * distance_right / interval &
            + (spline%values(left + 1) - spline%second_derivatives(left + 1) &
            * interval**2 / 6.0_dp) * distance_left / interval
        derivative = -spline%second_derivatives(left) * distance_right**2 &
            / (2.0_dp * interval) &
            + spline%second_derivatives(left + 1) * distance_left**2 &
            / (2.0_dp * interval) &
            + (spline%values(left + 1) - spline%values(left)) / interval &
            - interval * (spline%second_derivatives(left + 1) &
            - spline%second_derivatives(left)) / 6.0_dp
        if (present(second_derivative)) then
            second_derivative = (spline%second_derivatives(left) &
                * distance_right + spline%second_derivatives(left + 1) &
                * distance_left) / interval
        end if
        if (.not. ieee_is_finite(value)) return
        if (.not. ieee_is_finite(derivative)) return
        if (present(second_derivative)) then
            if (.not. ieee_is_finite(second_derivative)) return
        end if
        info = radial_cubic_spline_ok
    end subroutine evaluate_radial_cubic_spline

    subroutine build_reduced_system(intervals, lower, diagonal, upper)
        real(dp), intent(in) :: intervals(:)
        real(dp), intent(out) :: lower(:), diagonal(:), upper(:)
        integer :: i, interior

        interior = size(diagonal)
        diagonal(1) = (intervals(1) + intervals(2)) &
            * (2.0_dp + intervals(1) / intervals(2))
        upper(1) = intervals(2) - intervals(1)**2 / intervals(2)
        do i = 2, interior - 1
            lower(i - 1) = intervals(i)
            diagonal(i) = 2.0_dp * (intervals(i) + intervals(i + 1))
            upper(i) = intervals(i + 1)
        end do
        lower(interior - 1) = intervals(size(intervals) - 1) &
            - intervals(size(intervals))**2 &
            / intervals(size(intervals) - 1)
        diagonal(interior) = (intervals(size(intervals) - 1) &
            + intervals(size(intervals))) * (2.0_dp &
            + intervals(size(intervals)) / intervals(size(intervals) - 1))
    end subroutine build_reduced_system

    subroutine factor_tridiagonal(lower, diagonal, upper, info)
        real(dp), intent(inout) :: lower(:), diagonal(:)
        real(dp), intent(in) :: upper(:)
        integer, intent(out) :: info
        real(dp) :: scale
        integer :: i

        info = radial_cubic_spline_invalid
        scale = max(1.0_dp, maxval(abs(diagonal)), maxval(abs(lower)), &
            maxval(abs(upper)))
        if (abs(diagonal(1)) <= epsilon(1.0_dp) * scale) return
        do i = 2, size(diagonal)
            lower(i - 1) = lower(i - 1) / diagonal(i - 1)
            diagonal(i) = diagonal(i) - lower(i - 1) * upper(i - 1)
            if (abs(diagonal(i)) <= epsilon(1.0_dp) * scale) return
        end do
        if (.not. all(ieee_is_finite(lower))) return
        if (.not. all(ieee_is_finite(diagonal))) return
        info = radial_cubic_spline_ok
    end subroutine factor_tridiagonal

    subroutine solve_factored_system(grid, values)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        real(dp), intent(inout) :: values(:)
        integer :: i

        do i = 2, size(values)
            values(i) = values(i) - grid%lower_factor(i - 1) * values(i - 1)
        end do
        values(size(values)) = values(size(values)) &
            / grid%diagonal_factor(size(values))
        do i = size(values) - 1, 1, -1
            values(i) = (values(i) - grid%upper(i) * values(i + 1)) &
                / grid%diagonal_factor(i)
        end do
    end subroutine solve_factored_system

    function find_interval(nodes, coordinate) result(left)
        real(dp), intent(in) :: nodes(:), coordinate
        integer :: left
        integer :: high, low, middle

        if (coordinate <= nodes(1)) then
            left = 1
            return
        end if
        if (coordinate >= nodes(size(nodes))) then
            left = size(nodes) - 1
            return
        end if
        low = 1
        high = size(nodes)
        do while (high - low > 1)
            middle = low + (high - low) / 2
            if (coordinate < nodes(middle)) then
                high = middle
            else
                low = middle
            end if
        end do
        left = low
    end function find_interval

    function grid_is_valid(grid) result(valid)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        logical :: valid
        integer :: interior, nodes

        valid = .false.
        if (.not. allocated(grid%nodes)) return
        if (.not. allocated(grid%intervals)) return
        if (.not. allocated(grid%lower_factor)) return
        if (.not. allocated(grid%diagonal_factor)) return
        if (.not. allocated(grid%upper)) return
        nodes = size(grid%nodes)
        interior = nodes - 2
        if (nodes < 4) return
        if (size(grid%intervals) /= nodes - 1) return
        if (size(grid%lower_factor) /= interior - 1) return
        if (size(grid%diagonal_factor) /= interior) return
        if (size(grid%upper) /= interior - 1) return
        if (.not. all(ieee_is_finite(grid%nodes))) return
        if (.not. all(ieee_is_finite(grid%intervals))) return
        if (.not. all(ieee_is_finite(grid%lower_factor))) return
        if (.not. all(ieee_is_finite(grid%diagonal_factor))) return
        if (.not. all(ieee_is_finite(grid%upper))) return
        if (any(grid%intervals <= 0.0_dp)) return
        if (any(grid%nodes(2:) <= grid%nodes(:nodes - 1))) return
        if (any(grid%intervals /= &
            grid%nodes(2:) - grid%nodes(:nodes - 1))) return
        if (.not. ieee_is_finite(grid%domain_min)) return
        if (.not. ieee_is_finite(grid%domain_max)) return
        if (any(grid%diagonal_factor == 0.0_dp)) return
        if (grid%domain_min > grid%nodes(1)) return
        if (grid%domain_max < grid%nodes(nodes)) return
        valid = grid%domain_min < grid%domain_max
    end function grid_is_valid

    function spline_is_valid(grid, spline) result(valid)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        type(radial_cubic_spline_t), intent(in) :: spline
        logical :: valid

        valid = .false.
        if (.not. allocated(spline%values)) return
        if (.not. allocated(spline%second_derivatives)) return
        if (size(spline%values) /= size(grid%nodes)) return
        if (size(spline%second_derivatives) /= size(grid%nodes)) return
        if (.not. all(ieee_is_finite(spline%values))) return
        if (.not. all(ieee_is_finite(spline%second_derivatives))) return
        valid = .true.
    end function spline_is_valid

end module radial_cubic_spline
