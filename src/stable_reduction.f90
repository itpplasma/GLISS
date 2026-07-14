module stable_reduction
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: stable_dot_product
    public :: stable_norm2

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

    pure function stable_norm2(values) result(norm)
        real(dp), intent(in) :: values(:)
        real(dp) :: norm, scale, scaled(size(values)), squared

        if (size(values) == 0) then
            norm = 0.0_dp
            return
        end if
        scale = maxval(abs(values))
        if (scale == 0.0_dp) then
            norm = 0.0_dp
            return
        end if
        scaled = values / scale
        squared = stable_dot_product(scaled, scaled)
        norm = scale * sqrt(max(0.0_dp, squared))
    end function stable_norm2

end module stable_reduction
