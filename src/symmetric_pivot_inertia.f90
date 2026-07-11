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
        real(dp) :: determinant, trace
        integer :: j

        count = 0
        j = 1
        do while (j <= size(pivots))
            if (pivots(j) > 0) then
                if (factored(j, j) < 0.0_dp) count = count + 1
                j = j + 1
            else
                determinant = factored(j, j) * factored(j + 1, j + 1) &
                    - factored(j, j + 1)**2
                trace = factored(j, j) + factored(j + 1, j + 1)
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
