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
    ! the same values as the pairwise phase-product sum. Each product is
    ! symmetric: with rows scaled by sqrt|w| and grouped by the sign of w it
    ! is the difference of two symmetric rank-k updates, half the flops of a
    ! general matrix product.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: accumulate_period_averaged, accumulate_period_averaged_tangent
    public :: period_masks

    interface
        subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: uplo, trans
            integer, intent(in) :: n, k, lda, ldc
            real(dp), intent(in) :: alpha, beta
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dsyrk
        subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: uplo, trans
            integer, intent(in) :: n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta
            real(dp), intent(in) :: a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dsyr2k
    end interface

contains

    ! Trial weights of the two period averages: plus(a, b) multiplies
    ! X^T W X and minus(a, b) multiplies Y^T W Y for columns of trials a and
    ! b, from same = (n_a - n_b = 0) and opposite = (n_a + n_b = 0) modulo the
    ! field periods; mixed is false when minus vanishes, and Y then drops out.
    ! Column a has trial modulo(a - 1, trials) + 1.
    pure subroutine period_masks(trial_n, field_periods, plus, minus, mixed)
        integer, intent(in) :: trial_n(:), field_periods
        real(dp), intent(out) :: plus(:, :), minus(:, :)
        logical, intent(out) :: mixed
        real(dp) :: same, opposite
        integer :: a, b

        do b = 1, size(trial_n)
            do a = 1, size(trial_n)
                same = merge(1.0_dp, 0.0_dp, &
                    modulo(trial_n(a) - trial_n(b), field_periods) == 0)
                opposite = merge(1.0_dp, 0.0_dp, &
                    modulo(trial_n(a) + trial_n(b), field_periods) == 0)
                plus(a, b) = 0.5_dp * (same + opposite)
                minus(a, b) = 0.5_dp * (same - opposite)
            end do
        end do
        mixed = any(minus /= 0.0_dp)
    end subroutine period_masks

    ! cosine_part(p, a), sine_part(p, a): C_a and S_a at point p;
    ! cosine_phase(p, a), sine_phase(p, a): cos and sin of phi_a;
    ! weight(p): the signed quadrature weight of this channel.
    ! target(a, b) += sum_p weight(p) <R_a R_b>_periods.
    !
    ! rotated_target, with rotation_sign: the same sum for the columns whose
    ! trial basis is turned a quarter period, U' = -i s U with s = +-1 per
    ! column. Then X' = s Y and Y' = -s X, so the rotated form is
    ! s_a s_b [plus Y^T W Y + minus X^T W X](a, b) from the same products.
    subroutine accumulate_period_averaged(cosine_part, sine_part, &
            cosine_phase, sine_phase, weight, plus, minus, mixed, target, &
            rotated_target, rotation_sign)
        real(dp), contiguous, intent(in) :: cosine_part(:, :), sine_part(:, :)
        real(dp), contiguous, intent(in) :: cosine_phase(:, :)
        real(dp), contiguous, intent(in) :: sine_phase(:, :), weight(:)
        real(dp), intent(in) :: plus(:, :), minus(:, :)
        logical, intent(in) :: mixed
        real(dp), intent(inout) :: target(:, :)
        real(dp), intent(inout), optional :: rotated_target(:, :)
        real(dp), intent(in), optional :: rotation_sign(:)
        real(dp), allocatable :: scaled(:, :), product(:, :)
        integer :: column_trial(size(cosine_part, 2))
        integer :: order(size(weight))
        integer :: a, columns, negative, p, positive, trials

        columns = size(cosine_part, 2)
        trials = size(plus, 1)
        do a = 1, columns
            column_trial(a) = modulo(a - 1, trials) + 1
        end do
        ! Points with positive weights first, then negative ones; zero
        ! weights drop out. Rows scaled by sqrt|w| turn each weighted
        ! product into the difference of two symmetric rank-k updates.
        positive = 0
        do p = 1, size(weight)
            if (weight(p) > 0.0_dp) then
                positive = positive + 1
                order(positive) = p
            end if
        end do
        negative = 0
        do p = 1, size(weight)
            if (weight(p) < 0.0_dp) then
                negative = negative + 1
                order(positive + negative) = p
            end if
        end do
        if (positive + negative == 0) return
        allocate (scaled(positive + negative, columns), &
            product(columns, columns))
        call real_rows(.true.)
        call weighted_gram(positive, negative, columns, scaled, product)
        call add_masked(plus, product, column_trial, target)
        if (present(rotated_target)) call add_masked(minus, product, &
            column_trial, rotated_target, rotation_sign)
        if (.not. mixed .and. .not. present(rotated_target)) return
        call real_rows(.false.)
        call weighted_gram(positive, negative, columns, scaled, product)
        call add_masked(minus, product, column_trial, target)
        if (present(rotated_target)) call add_masked(plus, product, &
            column_trial, rotated_target, rotation_sign)
    contains
        ! X = Re(U) = C cos + S sin, or Y = Im(U) = C sin - S cos, at the
        ! ordered points, scaled by sqrt|w|.
        subroutine real_rows(real_part)
            logical, intent(in) :: real_part
            real(dp) :: root
            integer :: column, point, row

            do column = 1, columns
                do row = 1, positive + negative
                    point = order(row)
                    root = sqrt(abs(weight(point)))
                    if (real_part) then
                        scaled(row, column) = root &
                            * (cosine_part(point, column) &
                            * cosine_phase(point, column) &
                            + sine_part(point, column) &
                            * sine_phase(point, column))
                    else
                        scaled(row, column) = root &
                            * (cosine_part(point, column) &
                            * sine_phase(point, column) &
                            - sine_part(point, column) &
                            * cosine_phase(point, column))
                    end if
                end do
            end do
        end subroutine real_rows
    end subroutine accumulate_period_averaged

    ! Exact product rule at fixed weights, phases and field-period masks:
    ! target += sum_p w_p <dR_a R_b + R_a dR_b>_periods.
    subroutine accumulate_period_averaged_tangent(cosine_part, sine_part, &
            cosine_tangent, sine_tangent, cosine_phase, sine_phase, weight, &
            plus, minus, mixed, target)
        real(dp), contiguous, intent(in) :: cosine_part(:, :), sine_part(:, :)
        real(dp), contiguous, intent(in) :: cosine_tangent(:, :), sine_tangent(:, :)
        real(dp), contiguous, intent(in) :: cosine_phase(:, :), sine_phase(:, :)
        real(dp), contiguous, intent(in) :: weight(:)
        real(dp), intent(in) :: plus(:, :), minus(:, :)
        logical, intent(in) :: mixed
        real(dp), intent(inout) :: target(:, :)
        real(dp), allocatable :: scaled(:, :), tangent(:, :), product(:, :)
        integer :: column_trial(size(cosine_part, 2)), order(size(weight))
        integer :: positive, negative, rows, columns, p, a

        columns = size(cosine_part, 2)
        do a = 1, columns
            column_trial(a) = modulo(a - 1, size(plus, 1)) + 1
        end do
        positive = 0
        do p = 1, size(weight)
            if (weight(p) <= 0.0_dp) cycle
            positive = positive + 1
            order(positive) = p
        end do
        negative = 0
        do p = 1, size(weight)
            if (weight(p) >= 0.0_dp) cycle
            negative = negative + 1
            order(positive + negative) = p
        end do
        rows = positive + negative
        if (rows == 0) return
        allocate (scaled(rows, columns), tangent(rows, columns), &
            product(columns, columns))
        call product_rule(.true.)
        call add_masked(plus, product, column_trial, target)
        if (.not. mixed) return
        call product_rule(.false.)
        call add_masked(minus, product, column_trial, target)

    contains

        subroutine product_rule(real_part)
            logical, intent(in) :: real_part
            real(dp) :: root, phase_c, phase_s
            integer :: row, point, column

            do column = 1, columns
                do row = 1, rows
                    point = order(row)
                    root = sqrt(abs(weight(point)))
                    if (real_part) then
                        phase_c = cosine_phase(point, column)
                        phase_s = sine_phase(point, column)
                    else
                        phase_c = sine_phase(point, column)
                        phase_s = -cosine_phase(point, column)
                    end if
                    scaled(row, column) = root * (cosine_part(point, column) &
                        * phase_c + sine_part(point, column) * phase_s)
                    tangent(row, column) = root * (cosine_tangent(point, column) &
                        * phase_c + sine_tangent(point, column) * phase_s)
                end do
            end do
            product = 0.0_dp
            if (positive > 0) call dsyr2k('U', 'T', columns, positive, 1.0_dp, &
                scaled, rows, tangent, rows, 0.0_dp, product, columns)
            if (negative > 0) call dsyr2k('U', 'T', columns, negative, -1.0_dp, &
                scaled(positive + 1, 1), rows, tangent(positive + 1, 1), rows, &
                1.0_dp, product, columns)
        end subroutine product_rule

    end subroutine accumulate_period_averaged_tangent

    ! Upper triangle of A_+^T A_+ - A_-^T A_- for the first positive and
    ! the next negative rows of scaled.
    subroutine weighted_gram(positive, negative, columns, scaled, product)
        integer, intent(in) :: positive, negative, columns
        real(dp), intent(in) :: scaled(positive + negative, columns)
        real(dp), intent(out) :: product(columns, columns)
        integer :: rows

        rows = positive + negative
        product = 0.0_dp
        if (positive > 0) call dsyrk("U", "T", columns, positive, 1.0_dp, &
            scaled, rows, 0.0_dp, product, columns)
        if (negative > 0) call dsyrk("U", "T", columns, negative, -1.0_dp, &
            scaled(positive + 1, 1), rows, 1.0_dp, product, columns)
    end subroutine weighted_gram

    ! target(a, b) += sign(a) sign(b) mask(trial a, trial b) product(a, b)
    ! from the upper triangle of the symmetric product.
    pure subroutine add_masked(mask, product, column_trial, target, sign)
        real(dp), intent(in) :: mask(:, :), product(:, :)
        integer, intent(in) :: column_trial(:)
        real(dp), intent(inout) :: target(:, :)
        real(dp), intent(in), optional :: sign(:)
        real(dp) :: value
        integer :: a, b

        do b = 1, size(product, 2)
            do a = 1, b - 1
                value = mask(column_trial(a), column_trial(b)) * product(a, b)
                if (present(sign)) value = sign(a) * sign(b) * value
                target(a, b) = target(a, b) + value
                target(b, a) = target(b, a) + value
            end do
            target(b, b) = target(b, b) + mask(column_trial(b), &
                column_trial(b)) * product(b, b)
        end do
    end subroutine add_masked

end module period_averaged_assembly
