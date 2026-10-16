program test_parallel_blas
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_dynamic, omp_get_max_threads, &
        omp_get_num_threads, omp_set_dynamic, omp_set_num_threads
    implicit none
    integer, parameter :: rows = 256, columns = 72, jobs = 8
    real(dp) :: input(rows, columns, jobs)
    real(dp) :: reference(columns, columns, jobs), output(columns, columns, jobs)
    complex(dp) :: complex_input(rows, columns, jobs)
    complex(dp) :: complex_reference(columns, columns, jobs)
    complex(dp) :: complex_output(columns, columns, jobs)
    real(dp) :: real_scale, complex_scale, real_error, complex_error
    integer :: previous_threads, actual_threads, job, row, column, repeat
    integer :: seed, total, team_size
    logical :: previous_dynamic, valid_team, finite_results

    interface
        subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: uplo, trans
            integer, intent(in) :: n, k, lda, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dsyrk

        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            complex(dp), intent(inout) :: c(ldc, *)
        end subroutine zgemm
    end interface

    ! Each job owns disjoint contiguous buffers. This independently checks
    ! concurrent BLAS calls from the same OpenMP runtime as the application.
    previous_threads = omp_get_max_threads()
    previous_dynamic = omp_get_dynamic()
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(1)
    reference = 0.0_dp
    output = 0.0_dp
    do job = 1, jobs
        do column = 1, columns
            do row = 1, rows
                input(row, column, job) = sin(0.15_dp * row + 0.19_dp * column + job)
                complex_input(row, column, job) = cmplx(input(row, column, job), &
                    cos(0.17_dp * row + column + job), dp)
            end do
        end do
        call dsyrk('U', 'T', columns, rows, 1.0_dp, input(:, :, job), rows, &
            0.0_dp, reference(:, :, job), columns)
        call zgemm('C', 'N', columns, columns, rows, cmplx(1.0_dp, 0.0_dp, dp), &
            complex_input(:, :, job), rows, complex_input(:, :, job), rows, &
            cmplx(0.0_dp, 0.0_dp, dp), complex_reference(:, :, job), columns)
    end do
    real_scale = max(1.0_dp, maxval(abs(reference)))
    complex_scale = max(1.0_dp, maxval(abs(complex_reference)))
    finite_results = all(ieee_is_finite(reference)) &
        .and. all(ieee_is_finite(abs(complex_reference)))
    call omp_set_num_threads(2)
    seed = 7
    total = 0
    actual_threads = 0
    valid_team = .true.
    !$omp parallel do collapse(2) schedule(dynamic,1) default(none) &
    !$omp firstprivate(seed) private(team_size) &
    !$omp reduction(+:total) reduction(max:actual_threads) &
    !$omp reduction(.and.:valid_team)
    do job = 1, jobs
        do row = 1, 16
            total = total + seed + job + row
            team_size = omp_get_num_threads()
            actual_threads = max(actual_threads, team_size)
            valid_team = valid_team .and. team_size == 2
        end do
    end do
    !$omp end parallel do
    real_error = 0.0_dp
    complex_error = 0.0_dp
    do repeat = 1, 50
        !$omp parallel do default(none) &
        !$omp shared(input, output, complex_input, complex_output)
        do job = 1, jobs
            call dsyrk('U', 'T', columns, rows, 1.0_dp, input(:, :, job), rows, &
                0.0_dp, output(:, :, job), columns)
            call zgemm('C', 'N', columns, columns, rows, cmplx(1.0_dp, 0.0_dp, dp), &
                complex_input(:, :, job), rows, complex_input(:, :, job), rows, &
                cmplx(0.0_dp, 0.0_dp, dp), complex_output(:, :, job), columns)
        end do
        !$omp end parallel do
        real_error = maxval(abs(output - reference)) / real_scale
        complex_error = maxval(abs(complex_output - complex_reference)) / complex_scale
        finite_results = finite_results .and. all(ieee_is_finite(output)) &
            .and. all(ieee_is_finite(abs(complex_output))) &
            .and. ieee_is_finite(real_error) .and. ieee_is_finite(complex_error)
        if (.not. finite_results) exit
        if (max(real_error, complex_error) > 1.0e-12_dp) exit
    end do
    call omp_set_num_threads(previous_threads)
    call omp_set_dynamic(previous_dynamic)
    if (actual_threads /= 2) error stop 'parallel BLAS oracle requires two workers'
    if (total /= 2560 .or. .not. valid_team) &
        error stop 'OpenMP collapsed loop or reduction admission failed'
    if (.not. finite_results) error stop 'parallel BLAS produced nonfinite values'
    if (max(real_error, complex_error) > 1.0e-12_dp) then
        print '(a,2es14.5)', 'concurrent DSYRK/ZGEMM relative errors: ', &
            real_error, complex_error
        error stop 'concurrent BLAS differs from serial reference; check OpenMP runtimes'
    end if
    print *, 'PASS'
end program test_parallel_blas
