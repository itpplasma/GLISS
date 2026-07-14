module radial_feec_complex
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use radial_bspline_basis, only: bspline_ok, evaluate_bspline_basis
    implicit none
    private

    integer, parameter, public :: radial_feec_ok = 0
    integer, parameter, public :: radial_feec_invalid = -1
    integer, parameter, public :: radial_feec_allocation_error = -2

    type, public :: radial_feec_complex_t
        integer :: h1_degree = 0
        integer :: l2_degree = -1
        integer :: h1_dofs = 0
        integer :: l2_dofs = 0
        logical :: left_trace = .false.
        logical :: right_trace = .false.
        real(dp), allocatable :: h1_knots(:)
        real(dp), allocatable :: l2_knots(:)
        real(dp), allocatable :: derivative(:, :)
        integer, allocatable :: h1_basis_index(:)
    end type radial_feec_complex_t

    public :: build_radial_feec_complex
    public :: evaluate_radial_feec_complex

contains

    subroutine build_radial_feec_complex(breaks, degree, left_trace, &
            right_trace, complex, info)
        real(dp), intent(in) :: breaks(:)
        integer, intent(in) :: degree
        logical, intent(in) :: left_trace, right_trace
        type(radial_feec_complex_t), intent(out) :: complex
        integer, intent(out) :: info
        integer :: active, allocation_status, first, full_h1, last

        complex = radial_feec_complex_t()
        info = radial_feec_invalid
        if (degree < 1 .or. size(breaks) < 2) return
        if (degree > (huge(degree) - size(breaks)) / 2) return
        if (.not. all(ieee_is_finite(breaks))) return
        if (any(breaks(2:) <= breaks(:size(breaks) - 1))) return
        call build_knot_vectors(breaks, degree, complex%h1_knots, &
            complex%l2_knots, info)
        if (info /= radial_feec_ok) return
        full_h1 = size(complex%h1_knots) - degree - 1
        first = 1 + merge(1, 0, left_trace)
        last = full_h1 - merge(1, 0, right_trace)
        complex%h1_degree = degree
        complex%l2_degree = degree - 1
        complex%h1_dofs = max(0, last - first + 1)
        complex%l2_dofs = full_h1 - 1
        complex%left_trace = left_trace
        complex%right_trace = right_trace
        allocate (complex%h1_basis_index(complex%h1_dofs), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = radial_feec_allocation_error
            return
        end if
        complex%h1_basis_index = [(active, active=first, last)]
        call build_derivative_map(complex, info)
    end subroutine build_radial_feec_complex

    subroutine evaluate_radial_feec_complex(complex, coordinate, h1_basis, &
            h1_derivative, l2_basis, info)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: coordinate
        real(dp), allocatable, intent(out) :: h1_basis(:)
        real(dp), allocatable, intent(out) :: h1_derivative(:)
        real(dp), allocatable, intent(out) :: l2_basis(:)
        integer, intent(out) :: info
        real(dp), allocatable :: full_basis(:), full_derivative(:)
        real(dp), allocatable :: ignored(:)
        integer :: basis_info

        info = radial_feec_invalid
        if (.not. complex_is_valid(complex)) return
        call evaluate_bspline_basis(complex%h1_knots, complex%h1_degree, &
            coordinate, full_basis, full_derivative, basis_info)
        if (basis_info /= bspline_ok) return
        call evaluate_bspline_basis(complex%l2_knots, complex%l2_degree, &
            coordinate, l2_basis, ignored, basis_info)
        if (basis_info /= bspline_ok) return
        allocate (h1_basis(complex%h1_dofs), &
            h1_derivative(complex%h1_dofs))
        h1_basis = full_basis(complex%h1_basis_index)
        h1_derivative = full_derivative(complex%h1_basis_index)
        info = radial_feec_ok
    end subroutine evaluate_radial_feec_complex

    subroutine build_knot_vectors(breaks, degree, h1_knots, l2_knots, info)
        real(dp), intent(in) :: breaks(:)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: h1_knots(:), l2_knots(:)
        integer, intent(out) :: info
        integer :: allocation_status, internal_end, knot_count

        info = radial_feec_allocation_error
        knot_count = size(breaks) + 2 * degree
        allocate (h1_knots(knot_count), stat=allocation_status)
        if (allocation_status /= 0) return
        h1_knots(1:degree + 1) = breaks(1)
        internal_end = degree + size(breaks) - 1
        if (size(breaks) > 2) &
            h1_knots(degree + 2:internal_end) = breaks(2:size(breaks) - 1)
        h1_knots(internal_end + 1:) = breaks(size(breaks))
        allocate (l2_knots(knot_count - 2), stat=allocation_status)
        if (allocation_status /= 0) return
        l2_knots = h1_knots(2:knot_count - 1)
        info = radial_feec_ok
    end subroutine build_knot_vectors

    subroutine build_derivative_map(complex, info)
        type(radial_feec_complex_t), intent(inout) :: complex
        integer, intent(out) :: info
        real(dp) :: denominator, weight
        integer :: active, allocation_status, full, full_h1

        info = radial_feec_invalid
        full_h1 = complex%l2_dofs + 1
        allocate (complex%derivative(complex%l2_dofs, complex%h1_dofs), &
            source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) then
            info = radial_feec_allocation_error
            return
        end if
        do active = 1, complex%h1_dofs
            full = complex%h1_basis_index(active)
            if (full > 1) then
                denominator = complex%h1_knots(full + complex%h1_degree) &
                    - complex%h1_knots(full)
                if (denominator <= 0.0_dp) return
                weight = real(complex%h1_degree, dp) / denominator
                complex%derivative(full - 1, active) = weight
            end if
            if (full < full_h1) then
                denominator = complex%h1_knots(full + complex%h1_degree + 1) &
                    - complex%h1_knots(full + 1)
                if (denominator <= 0.0_dp) return
                weight = real(complex%h1_degree, dp) / denominator
                complex%derivative(full, active) = -weight
            end if
        end do
        if (.not. complex_is_valid(complex)) return
        info = radial_feec_ok
    end subroutine build_derivative_map

    function complex_is_valid(complex) result(valid)
        type(radial_feec_complex_t), intent(in) :: complex
        logical :: valid
        integer :: full_h1

        valid = .false.
        if (complex%h1_degree < 1) return
        if (complex%l2_degree /= complex%h1_degree - 1) return
        if (complex%h1_dofs < 0 .or. complex%l2_dofs < 1) return
        if (.not. allocated(complex%h1_knots)) return
        if (.not. allocated(complex%l2_knots)) return
        if (.not. allocated(complex%h1_basis_index)) return
        if (.not. allocated(complex%derivative)) return
        full_h1 = size(complex%h1_knots) - complex%h1_degree - 1
        if (complex%l2_dofs /= full_h1 - 1) return
        if (size(complex%l2_knots) /= size(complex%h1_knots) - 2) return
        if (size(complex%h1_basis_index) /= complex%h1_dofs) return
        if (any(shape(complex%derivative) /= &
            [complex%l2_dofs, complex%h1_dofs])) return
        if (any(complex%h1_basis_index < 1)) return
        if (any(complex%h1_basis_index > full_h1)) return
        if (.not. all(ieee_is_finite(complex%derivative))) return
        valid = .true.
    end function complex_is_valid

end module radial_feec_complex
