program test_radial_feec_complex
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use radial_feec_complex, only: build_radial_feec_complex, &
        evaluate_radial_feec_complex, radial_feec_complex_t, &
        radial_feec_invalid, radial_feec_ok
    implicit none

    integer :: degree

    do degree = 1, 4
        call verify_complex(degree, .false., .false.)
        call verify_complex(degree, .true., .false.)
        call verify_complex(degree, .false., .true.)
        call verify_complex(degree, .true., .true.)
    end do
    call verify_rejections()
    write (*, "(a)") "PASS"

contains

    subroutine verify_complex(degree, left_trace, right_trace)
        integer, intent(in) :: degree
        logical, intent(in) :: left_trace, right_trace
        real(dp), parameter :: breaks(5) = &
            [0.0_dp, 0.17_dp, 0.52_dp, 0.81_dp, 1.0_dp]
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: h1(:), h1_derivative(:), l2(:)
        real(dp), allocatable :: coefficients(:), mapped(:)
        real(dp) :: coordinate, direct
        integer :: expected_h1, expected_l2, info, sample

        call build_radial_feec_complex(breaks, degree, left_trace, &
            right_trace, complex, info)
        call require(info == radial_feec_ok, "valid complex was rejected")
        expected_h1 = (size(breaks) - 1) * degree + 1 &
            - merge(1, 0, left_trace) - merge(1, 0, right_trace)
        expected_l2 = (size(breaks) - 1) * degree
        call require(complex%h1_dofs == expected_h1, &
            "H1 dimension is wrong")
        call require(complex%l2_dofs == expected_l2, &
            "L2 dimension is wrong")
        call require(all(shape(complex%derivative) == &
            [expected_l2, expected_h1]), "derivative map shape is wrong")
        call require(numerical_rank(complex%derivative) == &
            min(expected_h1, expected_l2), "derivative rank is wrong")
        if (.not. left_trace .and. .not. right_trace) then
            call require(maxval(abs(matmul(complex%derivative, &
                [(1.0_dp, sample=1, expected_h1)]))) < 1.0e-13_dp, &
                "constant is not the derivative kernel")
        end if

        allocate (coefficients(expected_h1), mapped(expected_l2))
        coefficients = [(real(sample * sample - sample, dp) + 0.25_dp, &
            sample=1, expected_h1)]
        mapped = matmul(complex%derivative, coefficients)
        do sample = 0, 20
            coordinate = breaks(1) + (breaks(size(breaks)) - breaks(1)) &
                * real(sample, dp) / 20.0_dp
            call evaluate_radial_feec_complex(complex, coordinate, h1, &
                h1_derivative, l2, info)
            call require(info == radial_feec_ok, &
                "complex evaluation failed")
            call require(maxval(abs(h1_derivative &
                - matmul(transpose(complex%derivative), l2))) &
                < 2.0e-13_dp, "commuting derivative identity failed")
            call require(abs(sum(l2) - 1.0_dp) < 2.0e-14_dp, &
                "L2 partition of unity failed")
            direct = dot_product(h1_derivative, coefficients)
            call require(abs(direct - dot_product(l2, mapped)) &
                < 2.0e-12_dp, "coefficient derivative does not commute")
            if (.not. left_trace .and. .not. right_trace) then
                call require(abs(sum(h1) - 1.0_dp) < 2.0e-14_dp, &
                    "H1 partition of unity failed")
            end if
        end do
        call verify_traces(complex, left_trace, right_trace)
        call verify_knot_multiplicity(complex, breaks)
        call verify_fundamental_theorem(complex, breaks)
        call verify_local_l2_space(complex, breaks)
    end subroutine verify_complex

    subroutine verify_knot_multiplicity(complex, breaks)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:)
        integer :: boundary, expected

        call require(size(complex%h1_knots) &
            == complex%h1_degree * size(breaks) + 2, &
            "H1 knot count differs")
        do boundary = 1, size(breaks)
            expected = complex%h1_degree
            if (boundary == 1 .or. boundary == size(breaks)) &
                expected = expected + 1
            call require(count(complex%h1_knots == breaks(boundary)) &
                == expected, "H1 knot multiplicity differs")
            call require(count(complex%l2_knots == breaks(boundary)) &
                == complex%h1_degree, "L2 knot multiplicity differs")
        end do
    end subroutine verify_knot_multiplicity

    subroutine verify_local_l2_space(complex, breaks)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:)
        real(dp), allocatable :: collocation(:, :), h1(:), h1_derivative(:)
        real(dp), allocatable :: l2(:), nodes(:)
        integer, allocatable :: active(:), indices(:)
        real(dp) :: coordinate, half_width, midpoint
        integer :: cell, info, point

        nodes = local_gauss_nodes(complex%h1_degree)
        allocate (collocation(complex%h1_degree, complex%h1_degree))
        indices = [(point, point=1, complex%l2_dofs)]
        do cell = 1, size(breaks) - 1
            midpoint = 0.5_dp * (breaks(cell) + breaks(cell + 1))
            half_width = 0.5_dp * (breaks(cell + 1) - breaks(cell))
            call evaluate_radial_feec_complex(complex, midpoint, h1, &
                h1_derivative, l2, info)
            call require(info == radial_feec_ok, "midpoint evaluation failed")
            active = pack(indices, abs(l2) > 1.0e-14_dp)
            call require(size(active) == complex%h1_degree, &
                "L2 cell has the wrong number of active basis functions")
            call require(all(active == [(complex%h1_degree * (cell - 1) &
                + point, point=1, complex%h1_degree)]), &
                "L2 basis support crosses a cell boundary")
            do point = 1, complex%h1_degree
                coordinate = midpoint + half_width * nodes(point)
                call evaluate_radial_feec_complex(complex, coordinate, h1, &
                    h1_derivative, l2, info)
                call require(info == radial_feec_ok, &
                    "Gauss-node evaluation failed")
                collocation(point, :) = l2(active)
            end do
            call require(numerical_rank(collocation) == complex%h1_degree, &
                "L2 basis is singular at matched Gauss nodes")
        end do
    end subroutine verify_local_l2_space

    function local_gauss_nodes(degree) result(nodes)
        integer, intent(in) :: degree
        real(dp), allocatable :: nodes(:)

        allocate (nodes(degree))
        select case (degree)
        case (1)
            nodes = [0.0_dp]
        case (2)
            nodes = [-0.5773502691896258_dp, 0.5773502691896258_dp]
        case (3)
            nodes = [-0.7745966692414834_dp, 0.0_dp, &
                0.7745966692414834_dp]
        case (4)
            nodes = [-0.8611363115940526_dp, -0.3399810435848563_dp, &
                0.3399810435848563_dp, 0.8611363115940526_dp]
        end select
    end function local_gauss_nodes

    subroutine verify_fundamental_theorem(complex, breaks)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:)
        real(dp), parameter :: nodes(4) = [-0.8611363115940526_dp, &
            -0.3399810435848563_dp, 0.3399810435848563_dp, &
            0.8611363115940526_dp]
        real(dp), parameter :: weights(4) = [0.3478548451374539_dp, &
            0.6521451548625461_dp, 0.6521451548625461_dp, &
            0.3478548451374539_dp]
        real(dp), allocatable :: h1(:), h1_derivative(:), l2(:)
        real(dp), allocatable :: integral(:), left(:), right(:)
        real(dp) :: coordinate, half_width, midpoint
        integer :: cell, info, point

        allocate (integral(complex%h1_dofs), source=0.0_dp)
        do cell = 1, size(breaks) - 1
            midpoint = 0.5_dp * (breaks(cell) + breaks(cell + 1))
            half_width = 0.5_dp * (breaks(cell + 1) - breaks(cell))
            do point = 1, size(nodes)
                coordinate = midpoint + half_width * nodes(point)
                call evaluate_radial_feec_complex(complex, coordinate, h1, &
                    h1_derivative, l2, info)
                call require(info == radial_feec_ok, &
                    "quadrature evaluation failed")
                integral = integral + half_width * weights(point) &
                    * matmul(transpose(complex%derivative), l2)
            end do
        end do
        call evaluate_radial_feec_complex(complex, breaks(1), left, &
            h1_derivative, l2, info)
        call evaluate_radial_feec_complex(complex, breaks(size(breaks)), &
            right, h1_derivative, l2, info)
        call require(maxval(abs(integral - (right - left))) < 2.0e-13_dp, &
            "discrete fundamental theorem failed")
    end subroutine verify_fundamental_theorem

    subroutine verify_traces(complex, left_trace, right_trace)
        type(radial_feec_complex_t), intent(in) :: complex
        logical, intent(in) :: left_trace, right_trace
        real(dp), allocatable :: h1(:), h1_derivative(:), l2(:)
        integer :: info

        call evaluate_radial_feec_complex(complex, &
            complex%h1_knots(1), h1, h1_derivative, l2, info)
        call require(info == radial_feec_ok, "left trace evaluation failed")
        if (left_trace) then
            call require(all(h1 == 0.0_dp), "left trace was not eliminated")
        else
            call require(maxval(h1) == 1.0_dp, "left trace was lost")
        end if
        call evaluate_radial_feec_complex(complex, &
            complex%h1_knots(size(complex%h1_knots)), h1, h1_derivative, &
            l2, info)
        call require(info == radial_feec_ok, "right trace evaluation failed")
        if (right_trace) then
            call require(all(h1 == 0.0_dp), "right trace was not eliminated")
        else
            call require(maxval(h1) == 1.0_dp, "right trace was lost")
        end if
    end subroutine verify_traces

    subroutine verify_rejections()
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: h1(:), h1_derivative(:), l2(:)
        real(dp) :: nan
        integer :: info

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call build_radial_feec_complex([0.0_dp, 1.0_dp], 0, .false., &
            .false., complex, info)
        call require(info == radial_feec_invalid, "degree zero was accepted")
        call build_radial_feec_complex([0.0_dp, 1.0_dp], huge(0), .false., &
            .false., complex, info)
        call require(info == radial_feec_invalid, &
            "overflowing degree was accepted")
        call build_radial_feec_complex([0.0_dp], 2, .false., .false., &
            complex, info)
        call require(info == radial_feec_invalid, "one break was accepted")
        call build_radial_feec_complex([0.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], &
            2, .false., .false., complex, info)
        call require(info == radial_feec_invalid, &
            "repeated break was accepted")
        call build_radial_feec_complex([0.0_dp, nan, 1.0_dp], 2, .false., &
            .false., complex, info)
        call require(info == radial_feec_invalid, &
            "nonfinite break was accepted")
        call build_radial_feec_complex([0.0_dp, 1.0_dp], 2, .false., &
            .false., complex, info)
        call require(info == radial_feec_ok, "rejection fixture failed")
        call evaluate_radial_feec_complex(complex, -0.1_dp, h1, &
            h1_derivative, l2, info)
        call require(info == radial_feec_invalid, &
            "outside coordinate was accepted")
    end subroutine verify_rejections

    function numerical_rank(matrix) result(rank)
        real(dp), intent(in) :: matrix(:, :)
        integer :: rank
        real(dp), allocatable :: work(:, :), row_buffer(:)
        real(dp) :: pivot_value, tolerance
        integer :: column, pivot, row

        work = matrix
        allocate (row_buffer(size(matrix, 2)))
        tolerance = 1.0e-12_dp * max(1.0_dp, maxval(abs(matrix)))
        rank = 0
        do column = 1, size(work, 2)
            if (rank == size(work, 1)) exit
            pivot = rank + maxloc(abs(work(rank + 1:, column)), dim=1)
            if (abs(work(pivot, column)) <= tolerance) cycle
            rank = rank + 1
            if (pivot /= rank) then
                row_buffer = work(rank, :)
                work(rank, :) = work(pivot, :)
                work(pivot, :) = row_buffer
            end if
            pivot_value = work(rank, column)
            work(rank, :) = work(rank, :) / pivot_value
            do row = rank + 1, size(work, 1)
                work(row, :) = work(row, :) &
                    - work(row, column) * work(rank, :)
            end do
        end do
    end function numerical_rank

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_radial_feec_complex
