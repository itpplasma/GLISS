program test_nonuniform_derivative
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use nonuniform_derivative, only: first_derivative_nonuniform
    implicit none

    integer, parameter :: n = 17
    real(dp) :: x(n), values(n), derivative(n), exact(n)
    integer :: i

    do i = 1, n
        x(i) = (real(i - 1, dp) / real(n - 1, dp))**2
    end do
    values = 3.0_dp * x**2 - 2.0_dp * x + 5.0_dp
    exact = 6.0_dp * x - 2.0_dp
    call first_derivative_nonuniform(x, values, derivative)
    if (maxval(abs(derivative - exact)) > 1.0e-12_dp) then
        write (error_unit, "(a)") "quadratic is not differentiated exactly"
        error stop 1
    end if

    values = sin(x)
    call first_derivative_nonuniform(x, values, derivative)
    if (maxval(abs(derivative - cos(x))) > 1.0e-2_dp) then
        write (error_unit, "(a)") "smooth derivative error is too large"
        error stop 1
    end if
    write (*, "(a)") "PASS"
end program test_nonuniform_derivative
