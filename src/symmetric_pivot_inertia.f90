module symmetric_pivot_inertia
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: pivot_negative_count

contains

    pure function pivot_negative_count(factored, pivots) result(count)
        real(dp), intent(in) :: factored(:, :)
        integer, intent(in) :: pivots(:)
        integer :: count
        real(dp) :: determinant, trace, scale, a, b, c
        integer :: j

        count = 0
        j = 1
        do while (j <= size(pivots))
            if (pivots(j) > 0) then
                if (factored(j, j) < 0.0_dp) count = count + 1
                j = j + 1
            else
                ! Positive scaling preserves inertia and keeps the determinant
                ! products finite even for very small or large physical units.
                scale = max(abs(factored(j, j)), abs(factored(j, j + 1)), &
                    abs(factored(j + 1, j + 1)))
                if (scale == 0.0_dp) then
                    j = j + 2
                    cycle
                end if
                a = factored(j, j) / scale
                b = factored(j, j + 1) / scale
                c = factored(j + 1, j + 1) / scale
                determinant = a * c - b * b
                trace = a + c
                if (determinant < 0.0_dp) then
                    count = count + 1
                else if (trace < 0.0_dp) then
                    count = count + 2
                end if
                j = j + 2
            end if
        end do
    end function pivot_negative_count

end module symmetric_pivot_inertia
