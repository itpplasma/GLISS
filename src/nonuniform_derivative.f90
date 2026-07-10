module nonuniform_derivative
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: first_derivative_nonuniform

contains

    pure subroutine first_derivative_nonuniform(x, values, derivative)
        real(dp), intent(in) :: x(:), values(:)
        real(dp), intent(out) :: derivative(:)
        integer :: i, n
        real(dp) :: left, right

        n = size(x)
        do i = 2, n - 1
            left = x(i) - x(i - 1)
            right = x(i + 1) - x(i)
            derivative(i) = (values(i + 1) * left**2 &
                - values(i - 1) * right**2 &
                + values(i) * (right**2 - left**2)) &
                / (left * right * (left + right))
        end do
        left = x(2) - x(1)
        right = x(3) - x(2)
        derivative(1) = -values(1) * (2.0_dp * left + right) &
            / (left * (left + right)) &
            + values(2) * (left + right) / (left * right) &
            - values(3) * left / (right * (left + right))
        left = x(n - 1) - x(n - 2)
        right = x(n) - x(n - 1)
        derivative(n) = values(n - 2) * right / (left * (left + right)) &
            - values(n - 1) * (left + right) / (left * right) &
            + values(n) * (2.0_dp * right + left) &
            / (right * (left + right))
    end subroutine first_derivative_nonuniform

end module nonuniform_derivative
