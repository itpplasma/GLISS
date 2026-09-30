module period_averaged_assembly
    ! Angular quadrature of bilinear forms averaged over field periods, as
    ! matrix products.
    !
    ! Column a responds at angular point p as R_a = C_a cos(phi_a) +
    ! S_a sin(phi_a), with phi_a the phase of its trial on one field period.
    ! Writing R_a = Re(U_a), U_a = (C_a - i S_a) exp(i phi_a), the average
    ! of R_a R_b over the N_FP periods keeps (1/2) Re(U_a conj(U_b)) when
    ! n_a - n_b = 0 (mod N_FP) and (1/2) Re(U_a U_b) when n_a + n_b = 0
    ! (mod N_FP). With X = Re(U) and Y = Im(U) the weighted sum over points is
    !     1/2 [(d- + d+) X^T W X + (d- - d+) Y^T W Y],
    ! the same values as the pairwise phase-product sum, in two GEMMs.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: accumulate_period_averaged, period_masks

    interface
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta
            real(dp), intent(in) :: a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
    end interface

contains

    ! Column weights of the two period averages: plus(a, b) multiplies
    ! X^T W X and minus(a, b) multiplies Y^T W Y, from same = (n_a - n_b = 0)
    ! and opposite = (n_a + n_b = 0) modulo the field periods of the trials
    ! of columns a and b (column a has trial modulo(a - 1, trials) + 1);
    ! mixed is false when minus vanishes, and Y then drops out.
    pure subroutine period_masks(trial_n, field_periods, columns, plus, &
            minus, mixed)
        integer, intent(in) :: trial_n(:), field_periods, columns
        real(dp), intent(out) :: plus(:, :), minus(:, :)
        logical, intent(out) :: mixed
        real(dp) :: trial_plus(size(trial_n), size(trial_n))
        real(dp) :: trial_minus(size(trial_n), size(trial_n))
        real(dp) :: same, opposite
        integer :: a, b, trials

        trials = size(trial_n)
        do b = 1, trials
            do a = 1, trials
                same = merge(1.0_dp, 0.0_dp, &
                    modulo(trial_n(a) - trial_n(b), field_periods) == 0)
                opposite = merge(1.0_dp, 0.0_dp, &
                    modulo(trial_n(a) + trial_n(b), field_periods) == 0)
                trial_plus(a, b) = 0.5_dp * (same + opposite)
                trial_minus(a, b) = 0.5_dp * (same - opposite)
            end do
        end do
        mixed = any(trial_minus /= 0.0_dp)
        do b = 1, columns
            do a = 1, columns
                plus(a, b) = trial_plus(modulo(a - 1, trials) + 1, &
                    modulo(b - 1, trials) + 1)
                minus(a, b) = trial_minus(modulo(a - 1, trials) + 1, &
                    modulo(b - 1, trials) + 1)
            end do
        end do
    end subroutine period_masks

    ! cosine_part(p, a), sine_part(p, a): C_a and S_a at point p;
    ! cosine_phase(p, a), sine_phase(p, a): cos and sin of phi_a;
    ! weight(p): the signed quadrature weight of this channel.
    ! target(a, b) += sum_p weight(p) <R_a R_b>_periods.
    subroutine accumulate_period_averaged(cosine_part, sine_part, &
            cosine_phase, sine_phase, weight, plus, minus, mixed, target)
        real(dp), contiguous, intent(in) :: cosine_part(:, :), sine_part(:, :)
        real(dp), contiguous, intent(in) :: cosine_phase(:, :)
        real(dp), contiguous, intent(in) :: sine_phase(:, :), weight(:)
        real(dp), intent(in) :: plus(:, :), minus(:, :)
        logical, intent(in) :: mixed
        real(dp), intent(inout) :: target(:, :)
        real(dp), allocatable :: real_part(:, :), weighted(:, :)
        real(dp), allocatable :: product(:, :)
        integer :: a, b, columns, p, points

        points = size(cosine_part, 1)
        columns = size(cosine_part, 2)
        allocate (real_part(points, columns), weighted(points, columns), &
            product(columns, columns))
        do a = 1, columns
            do p = 1, points
                real_part(p, a) = cosine_part(p, a) * cosine_phase(p, a) &
                    + sine_part(p, a) * sine_phase(p, a)
                weighted(p, a) = weight(p) * real_part(p, a)
            end do
        end do
        call dgemm("T", "N", columns, columns, points, 1.0_dp, real_part, &
            points, weighted, points, 0.0_dp, product, columns)
        do b = 1, columns
            do a = 1, columns
                target(a, b) = target(a, b) + plus(a, b) * product(a, b)
            end do
        end do
        if (.not. mixed) return
        do a = 1, columns
            do p = 1, points
                real_part(p, a) = cosine_part(p, a) * sine_phase(p, a) &
                    - sine_part(p, a) * cosine_phase(p, a)
                weighted(p, a) = weight(p) * real_part(p, a)
            end do
        end do
        call dgemm("T", "N", columns, columns, points, 1.0_dp, real_part, &
            points, weighted, points, 0.0_dp, product, columns)
        do b = 1, columns
            do a = 1, columns
                target(a, b) = target(a, b) + minus(a, b) * product(a, b)
            end do
        end do
    end subroutine accumulate_period_averaged

end module period_averaged_assembly
