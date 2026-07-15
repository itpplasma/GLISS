module dense_generalized_inverse_iteration
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t, &
        valid_fixed_boundary_solver_controls
    use symmetric_eigensolver, only: solve_symmetric_generalized, &
        symmetric_eigensolver_ok
    implicit none
    private

    integer, parameter, public :: dense_inverse_ok = 0
    integer, parameter, public :: dense_inverse_invalid = -1
    integer, parameter, public :: dense_inverse_mass_not_spd = -2
    integer, parameter, public :: dense_inverse_allocation = -3
    integer, parameter, public :: dense_inverse_factorization = -4
    integer, parameter, public :: dense_inverse_no_convergence = -5

    public :: solve_dense_generalized_near_shift
    public :: solve_dense_generalized_subspace_near_shift

    interface
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
        subroutine dpocon(uplo, n, a, lda, anorm, rcond, work, iwork, info)
            import dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(in) :: a(lda, *), anorm
            real(dp), intent(out) :: rcond, work(*)
            integer, intent(out) :: iwork(*), info
        end subroutine dpocon
        subroutine dpotrf(uplo, n, a, lda, info)
            import dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, &
                x, ldx, ferr, berr, work, iwork, info)
            import dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldaf, ldb, ldx
            real(dp), intent(in) :: a(lda, *), af(ldaf, *), b(ldb, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: x(ldx, *)
            real(dp), intent(out) :: ferr(*), berr(*), work(*)
            integer, intent(out) :: iwork(*), info
        end subroutine dsyrfs
        subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
            import dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: work(*)
        end subroutine dsytrf
        subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dsytrs
        subroutine dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
            import dp
            character(len=1), intent(in) :: uplo, trans, diag
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dtrtrs
        subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, &
                b, ldb)
            import dp
            character(len=1), intent(in) :: side, uplo, transa, diag
            integer, intent(in) :: m, n, lda, ldb
            real(dp), intent(in) :: alpha, a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
        end subroutine dtrsm
    end interface

contains

    subroutine solve_dense_generalized_near_shift(stiffness, mass, shift, &
            eigenvalue, vector, action_residual, equilibrated_action_residual, &
            backward_error, standard_residual_norm, &
            mass_reciprocal_condition, iterations, info, controls)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :), shift
        real(dp), intent(out) :: eigenvalue, action_residual
        real(dp), intent(out) :: equilibrated_action_residual, backward_error
        real(dp), intent(out) :: standard_residual_norm
        real(dp), intent(out) :: mass_reciprocal_condition
        real(dp), allocatable, intent(out) :: vector(:)
        integer, intent(out) :: iterations, info
        type(fixed_boundary_solver_controls_t), intent(in), optional :: controls
        type(fixed_boundary_solver_controls_t) :: stopping
        integer, allocatable :: integer_work(:), pivots(:)
        real(dp), allocatable :: factor(:, :), mass_factor(:, :)
        real(dp), allocatable :: mass_scaled(:, :), right_hand_side(:, :)
        real(dp), allocatable :: scales(:), shifted_scaled(:, :)
        real(dp), allocatable :: source_right_hand_side(:, :), work(:)
        real(dp), allocatable :: best_vector(:), iterate(:)
        real(dp) :: backward_solve_error(1), forward_solve_error(1)
        real(dp) :: best_action_residual, best_backward_error
        real(dp) :: best_eigenvalue, best_equilibrated_action_residual
        real(dp) :: best_standard_residual_norm
        real(dp) :: mass_norm, previous, quotient
        real(dp) :: previous_best_standard_residual
        integer :: allocation_status, column, lapack_info, n, row
        integer :: stagnation, work_size
        integer(int64) :: requested_work_size

        info = dense_inverse_invalid
        eigenvalue = 0.0_dp
        action_residual = huge(1.0_dp)
        equilibrated_action_residual = huge(1.0_dp)
        backward_error = huge(1.0_dp)
        standard_residual_norm = huge(1.0_dp)
        mass_reciprocal_condition = 0.0_dp
        iterations = 0
        stopping = fixed_boundary_solver_controls_t()
        if (present(controls)) stopping = controls
        if (.not. valid_fixed_boundary_solver_controls(stopping)) return
        if (.not. valid_problem(stiffness, mass, shift)) return
        n = size(stiffness, 1)
        requested_work_size = 64_int64 * int(n, int64)
        if (requested_work_size > int(huge(work_size), int64)) return
        work_size = int(requested_work_size)
        allocate (factor(n, n), mass_factor(n, n), mass_scaled(n, n), &
            shifted_scaled(n, n), right_hand_side(n, 1), scales(n), &
            best_vector(n), iterate(n), vector(n), &
            source_right_hand_side(n, 1), pivots(n), integer_work(n), &
            work(work_size), stat=allocation_status)
        if (allocation_status /= 0) then
            info = dense_inverse_allocation
            return
        end if
        do row = 1, n
            if (mass(row, row) <= 0.0_dp) then
                info = dense_inverse_mass_not_spd
                return
            end if
            scales(row) = 1.0_dp / sqrt(mass(row, row))
        end do
        do column = 1, n
            do row = 1, n
                mass_scaled(row, column) = scales(row) &
                    * mass(row, column) * scales(column)
                mass_factor(row, column) = mass_scaled(row, column)
                shifted_scaled(row, column) = scales(row) &
                    * (stiffness(row, column) - shift * mass(row, column)) &
                    * scales(column)
                factor(row, column) = shifted_scaled(row, column)
            end do
        end do
        mass_norm = matrix_one_norm(mass_factor)
        call dpotrf("U", n, mass_factor, n, lapack_info)
        if (lapack_info /= 0) then
            info = dense_inverse_mass_not_spd
            return
        end if
        call dpocon("U", n, mass_factor, n, mass_norm, &
            mass_reciprocal_condition, work, integer_work, lapack_info)
        if (lapack_info /= 0 &
            .or. .not. ieee_is_finite(mass_reciprocal_condition)) then
            info = dense_inverse_factorization
            return
        end if
        call dsytrf("U", n, factor, n, pivots, work, work_size, lapack_info)
        if (lapack_info /= 0) then
            info = dense_inverse_factorization
            return
        end if
        do row = 1, n
            iterate(row) = 1.0_dp + real(row, dp) / real(n, dp)
        end do
        call normalize_scaled_mass(iterate, mass_scaled, mass_norm)
        if (mass_norm <= 0.0_dp) then
            info = dense_inverse_mass_not_spd
            return
        end if
        previous = huge(1.0_dp)
        best_standard_residual_norm = huge(1.0_dp)
        best_backward_error = huge(1.0_dp)
        stagnation = 0
        do iterations = 1, stopping%inverse_iteration_limit
            call matrix_vector(mass_scaled, iterate, &
                source_right_hand_side(:, 1))
            do row = 1, n
                right_hand_side(row, 1) = source_right_hand_side(row, 1)
            end do
            call dsytrs("U", n, 1, factor, n, pivots, right_hand_side, n, &
                lapack_info)
            if (lapack_info /= 0) then
                info = dense_inverse_factorization
                return
            end if
            call dsyrfs("U", n, 1, shifted_scaled, n, factor, n, pivots, &
                source_right_hand_side, n, right_hand_side, n, &
                forward_solve_error, backward_solve_error, work, integer_work, &
                lapack_info)
            if (lapack_info /= 0) then
                info = dense_inverse_factorization
                return
            end if
            call normalize_scaled_mass(right_hand_side(:, 1), mass_scaled, &
                mass_norm)
            if (mass_norm <= 0.0_dp) then
                info = dense_inverse_mass_not_spd
                return
            end if
            do row = 1, n
                iterate(row) = right_hand_side(row, 1)
                vector(row) = scales(row) * iterate(row)
            end do
            call generalized_diagnostics(stiffness, mass, shifted_scaled, &
                mass_scaled, mass_factor, shift, vector, iterate, &
                quotient, action_residual, equilibrated_action_residual, &
                backward_error, standard_residual_norm, lapack_info)
            if (lapack_info /= 0) then
                info = dense_inverse_factorization
                return
            end if
            eigenvalue = quotient
            previous_best_standard_residual = best_standard_residual_norm
            if (standard_residual_norm < best_standard_residual_norm) then
                best_eigenvalue = eigenvalue
                best_action_residual = action_residual
                best_equilibrated_action_residual = &
                    equilibrated_action_residual
                best_backward_error = backward_error
                best_standard_residual_norm = standard_residual_norm
                do row = 1, n
                    best_vector(row) = vector(row)
                end do
            end if
            if (standard_residual_norm &
                < (1.0_dp - 1.0e-3_dp) &
                * previous_best_standard_residual) then
                stagnation = 0
            else
                stagnation = stagnation + 1
            end if
            if (iterations > 1) then
                if (abs(eigenvalue - previous) <= &
                    stopping%eigenvalue_relative &
                    * max(1.0_dp, abs(eigenvalue)) &
                    .and. best_backward_error &
                    <= stopping%residual_relative) then
                    call restore_best
                    info = dense_inverse_ok
                    return
                end if
                if (iterations >= 8 .and. stagnation >= 8 &
                    .and. best_backward_error &
                    <= stopping%residual_relative) then
                    call restore_best
                    info = dense_inverse_ok
                    return
                end if
            end if
            previous = eigenvalue
        end do
        iterations = stopping%inverse_iteration_limit
        if (best_backward_error <= stopping%residual_relative) then
            call restore_best
            info = dense_inverse_ok
        else
            info = dense_inverse_no_convergence
        end if
    contains
        subroutine restore_best
            eigenvalue = best_eigenvalue
            action_residual = best_action_residual
            equilibrated_action_residual = &
                best_equilibrated_action_residual
            backward_error = best_backward_error
            standard_residual_norm = best_standard_residual_norm
            do row = 1, n
                vector(row) = best_vector(row)
            end do
        end subroutine restore_best
    end subroutine solve_dense_generalized_near_shift

    subroutine solve_dense_generalized_subspace_near_shift(stiffness, mass, &
            shift, initial, iteration_limit, eigenvalues, vectors, residuals, &
            iterations, info)
        real(dp), contiguous, intent(in) :: stiffness(:, :), mass(:, :)
        real(dp), contiguous, intent(in) :: initial(:, :)
        real(dp), intent(in) :: shift
        integer, intent(in) :: iteration_limit
        real(dp), allocatable, intent(out) :: eigenvalues(:), vectors(:, :)
        real(dp), allocatable, intent(out) :: residuals(:)
        integer, intent(out) :: iterations, info
        real(dp), allocatable :: candidate(:, :), coefficients(:, :)
        real(dp), allocatable :: current(:, :), factor(:, :), image(:, :)
        real(dp), allocatable :: mass_reduced(:, :), scales(:)
        real(dp), allocatable :: stiffness_reduced(:, :), work(:)
        integer, allocatable :: pivots(:)
        integer :: allocation_status, column, k, lapack_info, n, row
        integer :: work_size
        integer(int64) :: requested_work_size

        info = dense_inverse_invalid
        iterations = 0
        if (.not. valid_problem(stiffness, mass, shift)) return
        n = size(stiffness, 1)
        k = size(initial, 2)
        if (size(initial, 1) /= n .or. k < 1 .or. k >= n) return
        if (iteration_limit < 1 .or. iteration_limit > 1000) return
        if (.not. all(ieee_is_finite(initial))) return
        requested_work_size = 64_int64 * int(n, int64)
        if (requested_work_size > int(huge(work_size), int64)) return
        work_size = int(requested_work_size)
        allocate (candidate(n, k), current(n, k), factor(n, n), image(n, k), &
            mass_reduced(k, k), residuals(k), scales(n), &
            stiffness_reduced(k, k), vectors(n, k), pivots(n), &
            work(work_size), stat=allocation_status)
        if (allocation_status /= 0) then
            info = dense_inverse_allocation
            return
        end if
        do row = 1, n
            if (mass(row, row) <= 0.0_dp) then
                info = dense_inverse_mass_not_spd
                return
            end if
            scales(row) = 1.0_dp / sqrt(mass(row, row))
        end do
        do column = 1, n
            do row = 1, n
                factor(row, column) = scales(row) &
                    * (stiffness(row, column) - shift * mass(row, column)) &
                    * scales(column)
            end do
        end do
        call dsytrf("U", n, factor, n, pivots, work, work_size, lapack_info)
        if (lapack_info /= 0) then
            info = dense_inverse_factorization
            return
        end if
        current = initial
        call orthonormalize_subspace(current, mass, image, mass_reduced, &
            lapack_info)
        if (lapack_info /= 0) then
            info = dense_inverse_mass_not_spd
            return
        end if
        do iterations = 1, iteration_limit
            call contract_dense_subspace(mass, current, image, mass_reduced)
            do column = 1, k
                do row = 1, n
                    candidate(row, column) = scales(row) * image(row, column)
                end do
            end do
            call dsytrs("U", n, k, factor, n, pivots, candidate, n, &
                lapack_info)
            if (lapack_info /= 0) then
                info = dense_inverse_factorization
                return
            end if
            do column = 1, k
                do row = 1, n
                    candidate(row, column) = scales(row) &
                        * candidate(row, column)
                end do
            end do
            call orthonormalize_subspace(candidate, mass, image, mass_reduced, &
                lapack_info)
            if (lapack_info /= 0) then
                info = dense_inverse_mass_not_spd
                return
            end if
            call contract_dense_subspace(stiffness, candidate, image, &
                stiffness_reduced)
            call contract_dense_subspace(mass, candidate, image, mass_reduced)
            call solve_symmetric_generalized(stiffness_reduced, mass_reduced, &
                eigenvalues, coefficients, lapack_info)
            if (lapack_info /= symmetric_eigensolver_ok) then
                info = dense_inverse_factorization
                return
            end if
            call dgemm("N", "N", n, k, k, 1.0_dp, candidate, n, &
                coefficients, k, 0.0_dp, vectors, n)
            current = vectors
        end do
        iterations = iteration_limit
        call contract_dense_subspace(stiffness, current, image, &
            stiffness_reduced)
        call contract_dense_subspace(mass, current, candidate, mass_reduced)
        call calculate_subspace_residuals(image, candidate, eigenvalues, &
            residuals)
        vectors = current
        info = dense_inverse_ok
    end subroutine solve_dense_generalized_subspace_near_shift

    subroutine orthonormalize_subspace(vectors, mass, image, gram, info)
        real(dp), contiguous, intent(inout) :: vectors(:, :)
        real(dp), contiguous, intent(in) :: mass(:, :)
        real(dp), contiguous, intent(out) :: image(:, :), gram(:, :)
        integer, intent(out) :: info

        call contract_dense_subspace(mass, vectors, image, gram)
        call dpotrf("U", size(gram, 1), gram, size(gram, 1), info)
        if (info /= 0) return
        call dtrsm("R", "U", "N", "N", size(vectors, 1), &
            size(vectors, 2), 1.0_dp, gram, size(gram, 1), vectors, &
            size(vectors, 1))
    end subroutine orthonormalize_subspace

    subroutine contract_dense_subspace(matrix, vectors, image, reduced)
        real(dp), contiguous, intent(in) :: matrix(:, :), vectors(:, :)
        real(dp), contiguous, intent(out) :: image(:, :), reduced(:, :)
        integer :: first, row, second

        call dgemm("N", "N", size(matrix, 1), size(vectors, 2), &
            size(matrix, 2), 1.0_dp, matrix, size(matrix, 1), vectors, &
            size(vectors, 1), 0.0_dp, image, size(image, 1))
        do second = 1, size(vectors, 2)
            do first = 1, second
                reduced(first, second) = 0.0_dp
                do row = 1, size(vectors, 1)
                    reduced(first, second) = reduced(first, second) &
                        + vectors(row, first) * image(row, second)
                end do
                reduced(second, first) = reduced(first, second)
            end do
        end do
    end subroutine contract_dense_subspace

    subroutine calculate_subspace_residuals(stiffness_images, mass_images, &
            eigenvalues, residuals)
        real(dp), intent(in) :: stiffness_images(:, :), mass_images(:, :)
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(out) :: residuals(:)
        real(dp) :: difference, difference_norm, mass_norm, stiffness_norm
        integer :: column, row

        do column = 1, size(eigenvalues)
            difference_norm = 0.0_dp
            mass_norm = 0.0_dp
            stiffness_norm = 0.0_dp
            do row = 1, size(stiffness_images, 1)
                difference = stiffness_images(row, column) &
                    - eigenvalues(column) * mass_images(row, column)
                difference_norm = difference_norm + difference**2
                stiffness_norm = stiffness_norm &
                    + stiffness_images(row, column)**2
                mass_norm = mass_norm + mass_images(row, column)**2
            end do
            residuals(column) = sqrt(difference_norm) / &
                max(tiny(1.0_dp), sqrt(stiffness_norm) &
                + abs(eigenvalues(column)) * sqrt(mass_norm))
        end do
    end subroutine calculate_subspace_residuals

    subroutine normalize_scaled_mass(vector, mass, norm)
        real(dp), intent(inout) :: vector(:)
        real(dp), intent(in) :: mass(:, :)
        real(dp), intent(out) :: norm
        real(dp) :: image(size(vector)), squared
        integer :: row

        call matrix_vector(mass, vector, image)
        squared = vector_dot(vector, image)
        if (.not. ieee_is_finite(squared) .or. squared <= 0.0_dp) then
            norm = -1.0_dp
            return
        end if
        norm = sqrt(squared)
        do row = 1, size(vector)
            vector(row) = vector(row) / norm
        end do
    end subroutine normalize_scaled_mass

    subroutine generalized_diagnostics(stiffness, mass, shifted_scaled, &
            mass_scaled, mass_factor, shift, vector, scaled_vector, &
            quotient, action_residual, equilibrated_action_residual, &
            backward_error, standard_residual_norm, info)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :)
        real(dp), intent(in) :: shifted_scaled(:, :), mass_scaled(:, :)
        real(dp), intent(in) :: mass_factor(:, :), shift
        real(dp), intent(in) :: vector(:), scaled_vector(:)
        real(dp), intent(out) :: quotient, action_residual
        real(dp), intent(out) :: equilibrated_action_residual, backward_error
        real(dp), intent(out) :: standard_residual_norm
        integer, intent(out) :: info
        real(dp) :: action(size(vector)), difference, difference_squared
        real(dp) :: equilibrated_difference_squared
        real(dp) :: equilibrated_mass_norm_squared
        real(dp) :: equilibrated_operator_norm_squared
        real(dp) :: mass_action(size(vector)), mass_scaled_action(size(vector))
        real(dp) :: operator_norm_squared, mass_norm_squared
        real(dp) :: shifted_action(size(vector)), shifted_norm_squared
        real(dp) :: shifted_potential, scaled_norm_squared
        real(dp) :: mass_matrix_norm_squared, standard_residual(size(vector), 1)
        integer :: column, row

        call matrix_vector(stiffness, vector, action)
        call matrix_vector(mass, vector, mass_action)
        call matrix_vector(shifted_scaled, scaled_vector, shifted_action)
        call matrix_vector(mass_scaled, scaled_vector, mass_scaled_action)
        shifted_potential = vector_dot(scaled_vector, shifted_action)
        quotient = shift + shifted_potential &
            / vector_dot(scaled_vector, mass_scaled_action)
        difference_squared = 0.0_dp
        equilibrated_difference_squared = 0.0_dp
        equilibrated_operator_norm_squared = 0.0_dp
        equilibrated_mass_norm_squared = 0.0_dp
        operator_norm_squared = 0.0_dp
        mass_norm_squared = 0.0_dp
        shifted_norm_squared = 0.0_dp
        mass_matrix_norm_squared = 0.0_dp
        scaled_norm_squared = vector_dot(scaled_vector, scaled_vector)
        do column = 1, size(vector)
            do row = 1, size(vector)
                shifted_norm_squared = shifted_norm_squared &
                    + shifted_scaled(row, column)**2
                mass_matrix_norm_squared = mass_matrix_norm_squared &
                    + mass_scaled(row, column)**2
            end do
        end do
        do row = 1, size(vector)
            difference = action(row) - quotient * mass_action(row)
            difference_squared = difference_squared + difference**2
            operator_norm_squared = operator_norm_squared + action(row)**2
            mass_norm_squared = mass_norm_squared + mass_action(row)**2
            difference = shifted_action(row) &
                - (quotient - shift) * mass_scaled_action(row)
            standard_residual(row, 1) = difference
            equilibrated_difference_squared = equilibrated_difference_squared &
                + difference**2
            equilibrated_operator_norm_squared = &
                equilibrated_operator_norm_squared &
                + (shifted_action(row) &
                + shift * mass_scaled_action(row))**2
            equilibrated_mass_norm_squared = equilibrated_mass_norm_squared &
                + mass_scaled_action(row)**2
        end do
        action_residual = sqrt(difference_squared) / &
            max(tiny(1.0_dp), sqrt(operator_norm_squared) &
            + abs(quotient) * sqrt(mass_norm_squared))
        equilibrated_action_residual = sqrt( &
            equilibrated_difference_squared) / max(tiny(1.0_dp), &
            sqrt(equilibrated_operator_norm_squared) &
            + abs(quotient) * sqrt(equilibrated_mass_norm_squared))
        backward_error = sqrt(equilibrated_difference_squared) / &
            max(tiny(1.0_dp), (sqrt(shifted_norm_squared) &
            + abs(quotient - shift) * sqrt(mass_matrix_norm_squared)) &
            * sqrt(scaled_norm_squared))
        call dtrtrs("U", "T", "N", size(vector), 1, mass_factor, &
            size(vector), standard_residual, size(vector), info)
        if (info /= 0) return
        standard_residual_norm = sqrt( &
            vector_dot(standard_residual(:, 1), standard_residual(:, 1)))
        if (.not. ieee_is_finite(quotient) &
            .or. .not. ieee_is_finite(action_residual) &
            .or. .not. ieee_is_finite(equilibrated_action_residual) &
            .or. .not. ieee_is_finite(backward_error) &
            .or. .not. ieee_is_finite(standard_residual_norm)) info = -1
    end subroutine generalized_diagnostics

    subroutine matrix_vector(matrix, vector, result)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp), intent(out) :: result(:)
        integer :: column, row

        do row = 1, size(vector)
            result(row) = 0.0_dp
            do column = 1, size(vector)
                result(row) = result(row) + matrix(row, column) * vector(column)
            end do
        end do
    end subroutine matrix_vector

    function vector_dot(left, right) result(value)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: value
        integer :: row

        value = 0.0_dp
        do row = 1, size(left)
            value = value + left(row) * right(row)
        end do
    end function vector_dot

    function matrix_one_norm(matrix) result(value)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: column_sum, value
        integer :: column, row

        value = 0.0_dp
        do column = 1, size(matrix, 2)
            column_sum = 0.0_dp
            do row = 1, size(matrix, 1)
                column_sum = column_sum + abs(matrix(row, column))
            end do
            value = max(value, column_sum)
        end do
    end function matrix_one_norm

    function valid_problem(stiffness, mass, shift) result(valid)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :), shift
        logical :: valid
        integer :: column, n, row

        valid = .false.
        n = size(stiffness, 1)
        if (n < 1 .or. size(stiffness, 2) /= n) return
        if (size(mass, 1) /= n .or. size(mass, 2) /= n) return
        if (.not. ieee_is_finite(shift)) return
        do column = 1, n
            do row = 1, n
                if (.not. ieee_is_finite(stiffness(row, column))) return
                if (.not. ieee_is_finite(mass(row, column))) return
            end do
        end do
        if (.not. symmetric_matrix(stiffness)) return
        if (.not. symmetric_matrix(mass)) return
        valid = .true.
    end function valid_problem

    function symmetric_matrix(matrix) result(symmetric)
        real(dp), intent(in) :: matrix(:, :)
        logical :: symmetric
        real(dp) :: scale, tolerance
        integer :: column, row

        scale = 1.0_dp
        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                scale = max(scale, abs(matrix(row, column)))
            end do
        end do
        tolerance = 64.0_dp * epsilon(1.0_dp) &
            * real(size(matrix, 1), dp) * scale
        symmetric = .true.
        do column = 1, size(matrix, 2)
            do row = column + 1, size(matrix, 1)
                if (abs(matrix(row, column) - matrix(column, row)) &
                    > tolerance) then
                    symmetric = .false.
                    return
                end if
            end do
        end do
    end function symmetric_matrix

end module dense_generalized_inverse_iteration
