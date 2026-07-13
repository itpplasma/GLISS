module stable_reduction
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: stable_dot_product

contains

    pure function stable_dot_product(first, second) result(product)
        real(dp), intent(in) :: first(:), second(:)
        real(dp) :: product, correction, term, updated
        integer :: i

        product = 0.0_dp
        correction = 0.0_dp
        do i = 1, size(first)
            term = first(i) * second(i)
            updated = product + term
            if (abs(product) >= abs(term)) then
                correction = correction + (product - updated) + term
            else
                correction = correction + (term - updated) + product
            end if
            product = updated
        end do
        product = product + correction
    end function stable_dot_product

end module stable_reduction
