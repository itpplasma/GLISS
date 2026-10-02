! Manufactured convergence and exactness of the radial FEEC complex.
!
! For degree p the H1 space holds C0 piecewise polynomials of degree p and
! the L2 space their discontinuous degree-(p-1) derivatives. L2 projections of a
! smooth function converge in L2 like h^(p+1) and h^p, and d/ds of the H1
! projection converges like h^p. Monomials up to the space degree are
! reproduced to roundoff (the element oracle), and the derivative map is
! exact: D applied to the H1 coefficients of s^k gives the L2 coefficients
! of k s^(k-1). A corrupted derivative entry must break that identity.
program test_radial_feec_convergence
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use radial_feec_complex, only: build_radial_feec_complex, &
        evaluate_radial_feec_complex, radial_feec_complex_t, radial_feec_ok
    implicit none
    integer, parameter :: meshes(3) = [8, 16, 32]
    real(dp), parameter :: gauss_x(5) = [ &
        -0.9061798459386640_dp, -0.5384693101056831_dp, 0.0_dp, &
        0.5384693101056831_dp, 0.9061798459386640_dp]
    real(dp), parameter :: gauss_w(5) = [ &
        0.2369268850561891_dp, 0.4786286704993665_dp, &
        0.5688888888888889_dp, 0.4786286704993665_dp, &
        0.2369268850561891_dp]
    real(dp) :: errors(3, 3), rate
    integer :: degree, mesh, kind

    do degree = 1, 4
        do mesh = 1, size(meshes)
            call projection_errors(degree, meshes(mesh), errors(:, mesh))
        end do
        do kind = 1, 3
            rate = log(errors(kind, 2) / errors(kind, 3)) / log(2.0_dp)
            write (*, "(a,i0,a,i0,a,3es10.2,a,f6.3)") "degree ", degree, &
                " kind ", kind, " errors ", errors(kind, :), " rate ", rate
            call require(abs(rate - expected_rate(degree, kind)) <= 0.15_dp, &
                "projection does not converge at the optimal rate")
        end do
        call check_monomials(degree)
    end do
    call check_corrupted_derivative()
    write (*, "(a)") "PASS"

contains

    pure function expected_rate(degree, kind) result(rate)
        integer, intent(in) :: degree, kind
        real(dp) :: rate

        ! kind 1: H1 projection in L2; 2: L2 projection in L2;
        ! 3: derivative of the H1 projection in L2.
        rate = real(degree, dp)
        if (kind == 1) rate = real(degree + 1, dp)
    end function expected_rate

    pure function exact(s) result(value)
        real(dp), intent(in) :: s
        real(dp) :: value

        value = sin(2.3_dp * s) + exp(-s) * s**2
    end function exact

    pure function exact_derivative(s) result(value)
        real(dp), intent(in) :: s
        real(dp) :: value

        value = 2.3_dp * cos(2.3_dp * s) + exp(-s) * (2.0_dp * s - s**2)
    end function exact_derivative

    subroutine uniform_breaks(elements, breaks)
        integer, intent(in) :: elements
        real(dp), allocatable, intent(out) :: breaks(:)
        integer :: i

        allocate (breaks(elements + 1))
        do i = 0, elements
            ! A smooth graded mesh keeps the elements quasi-uniform.
            breaks(i + 1) = real(i, dp) / elements
            breaks(i + 1) = breaks(i + 1) + 0.1_dp * breaks(i + 1) &
                * (1.0_dp - breaks(i + 1))
        end do
    end subroutine uniform_breaks

    subroutine projection_errors(degree, elements, errors)
        integer, intent(in) :: degree, elements
        real(dp), intent(out) :: errors(3)
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: breaks(:), h1_coefficients(:)
        real(dp), allocatable :: l2_coefficients(:), derivative(:)
        integer :: info

        call uniform_breaks(elements, breaks)
        call build_radial_feec_complex(breaks, degree, .false., .false., &
            complex, info)
        call require(info == radial_feec_ok, "complex construction failed")
        call project(complex, breaks, .true., h1_coefficients)
        call project(complex, breaks, .false., l2_coefficients)
        derivative = matmul(complex%derivative, h1_coefficients)
        errors(1) = l2_error(complex, breaks, h1_coefficients, 1)
        errors(2) = l2_error(complex, breaks, l2_coefficients, 2)
        errors(3) = l2_error(complex, breaks, derivative, 3)
    end subroutine projection_errors

    subroutine project(complex, breaks, h1_space, coefficients, power)
        ! L2 projection of the manufactured function (or of s^power).
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:)
        logical, intent(in) :: h1_space
        real(dp), allocatable, intent(out) :: coefficients(:)
        integer, intent(in), optional :: power
        real(dp), allocatable :: mass(:, :), load(:), basis(:)
        real(dp) :: s, weight, value
        integer :: n, element, point, info, column

        n = merge(complex%h1_dofs, complex%l2_dofs, h1_space)
        allocate (mass(n, n), load(n))
        mass = 0.0_dp
        load = 0.0_dp
        do element = 1, size(breaks) - 1
            do point = 1, size(gauss_x)
                call quadrature_point(breaks, element, point, 1, 1, s, weight)
                call space_basis(complex, s, h1_space, basis)
                value = exact(s)
                if (present(power)) value = s**power
                do column = 1, n
                    mass(:, column) = mass(:, column) &
                        + weight * basis(column) * basis
                end do
                load = load + weight * value * basis
            end do
        end do
        call solve_spd(n, mass, load, info)
        call require(info == 0, "projection mass is not positive definite")
        coefficients = load
    end subroutine project

    function l2_error(complex, breaks, coefficients, kind) result(error)
        ! Integrated on four subcells per element so that the error norm is
        ! not evaluated at the projection's own quadrature points.
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: breaks(:), coefficients(:)
        integer, intent(in) :: kind
        real(dp) :: error
        real(dp), allocatable :: basis(:)
        real(dp) :: s, weight, reference
        integer :: element, point, cell

        error = 0.0_dp
        do element = 1, size(breaks) - 1
            do cell = 1, 4
                do point = 1, size(gauss_x)
                    call quadrature_point(breaks, element, point, cell, 4, s, &
                        weight)
                    call space_basis(complex, s, kind == 1, basis)
                    reference = exact(s)
                    if (kind == 3) reference = exact_derivative(s)
                    error = error + weight &
                        * (dot_product(basis, coefficients) - reference)**2
                end do
            end do
        end do
        error = sqrt(error)
    end function l2_error

    subroutine check_monomials(degree)
        ! Element oracle: s^k for k <= p is in the H1 space and k s^(k-1)
        ! in the L2 space, so projection and derivative map are exact.
        integer, intent(in) :: degree
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: breaks(:), h1(:), l2(:), mapped(:)
        real(dp), allocatable :: basis(:)
        real(dp) :: s, worst
        integer :: power, info, sample

        call uniform_breaks(5, breaks)
        call build_radial_feec_complex(breaks, degree, .false., .false., &
            complex, info)
        call require(info == radial_feec_ok, "complex construction failed")
        do power = 0, degree
            call project(complex, breaks, .true., h1, power)
            mapped = matmul(complex%derivative, h1)
            worst = 0.0_dp
            do sample = 0, 40
                s = real(sample, dp) / 40.0_dp
                call space_basis(complex, s, .true., basis)
                worst = max(worst, abs(dot_product(basis, h1) - s**power))
                call space_basis(complex, s, .false., basis)
                if (power > 0) then
                    worst = max(worst, abs(dot_product(basis, mapped) &
                        - power * s**(power - 1)))
                else
                    worst = max(worst, abs(dot_product(basis, mapped)))
                end if
            end do
            call require(worst <= 1.0e-12_dp, &
                "monomial is not reproduced by the FEEC complex")
            if (power < degree) then
                call project(complex, breaks, .false., l2, power)
                do sample = 0, 40
                    s = real(sample, dp) / 40.0_dp
                    call space_basis(complex, s, .false., basis)
                    call require(abs(dot_product(basis, l2) - s**power) &
                        <= 1.0e-12_dp, "L2 monomial is not reproduced")
                end do
            end if
        end do
    end subroutine check_monomials

    subroutine check_corrupted_derivative()
        ! Control: a single corrupted derivative-map entry must break the
        ! exact sequence identity that check_monomials relies on.
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: breaks(:), h1(:), mapped(:), basis(:)
        real(dp) :: s, worst
        integer :: info, sample

        call uniform_breaks(5, breaks)
        call build_radial_feec_complex(breaks, 2, .false., .false., complex, &
            info)
        call require(info == radial_feec_ok, "complex construction failed")
        complex%derivative(3, 4) = complex%derivative(3, 4) + 1.0e-3_dp
        call project(complex, breaks, .true., h1, 2)
        mapped = matmul(complex%derivative, h1)
        worst = 0.0_dp
        do sample = 0, 40
            s = real(sample, dp) / 40.0_dp
            call space_basis(complex, s, .false., basis)
            worst = max(worst, abs(dot_product(basis, mapped) - 2.0_dp * s))
        end do
        call require(worst > 1.0e-6_dp, &
            "corrupted derivative map was not detected")
    end subroutine check_corrupted_derivative

    subroutine quadrature_point(breaks, element, point, cell, cells, s, &
            weight)
        real(dp), intent(in) :: breaks(:)
        integer, intent(in) :: element, point, cell, cells
        real(dp), intent(out) :: s, weight
        real(dp) :: left, width

        width = (breaks(element + 1) - breaks(element)) / cells
        left = breaks(element) + (cell - 1) * width
        s = left + 0.5_dp * width * (gauss_x(point) + 1.0_dp)
        weight = 0.5_dp * width * gauss_w(point)
    end subroutine quadrature_point

    subroutine space_basis(complex, s, h1_space, basis)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: s
        logical, intent(in) :: h1_space
        real(dp), allocatable, intent(out) :: basis(:)
        real(dp), allocatable :: h1(:), h1_derivative(:), l2(:)
        integer :: info

        call evaluate_radial_feec_complex(complex, s, h1, h1_derivative, l2, &
            info)
        call require(info == radial_feec_ok, "basis evaluation failed")
        if (h1_space) then
            call move_alloc(h1, basis)
        else
            call move_alloc(l2, basis)
        end if
    end subroutine space_basis

    subroutine solve_spd(n, matrix, rhs, info)
        integer, intent(in) :: n
        real(dp), intent(inout) :: matrix(n, n), rhs(n)
        integer, intent(out) :: info
        interface
            subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
                import dp
                character(len=1), intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: info
            end subroutine dposv
        end interface

        call dposv("U", n, 1, matrix, n, rhs, n, info)
    end subroutine solve_spd

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_radial_feec_convergence
