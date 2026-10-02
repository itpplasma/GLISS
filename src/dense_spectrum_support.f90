module dense_spectrum_support
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_boundary_eigen_bracket, only: bounded_inertia_probe, &
        fixed_boundary_bracket_ok
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use variable_block_tridiagonal, only: apply_variable_block_tridiagonal, &
        factorize_variable_shifted, variable_block_factor_t, &
        variable_block_ok, variable_block_to_dense, &
        variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, pencil_roundoff, &
        variable_eigenvalue_bound, variable_generalized_diagnostics, &
        variable_generalized_inertia, variable_generalized_ok, &
        validate_variable_pencil
    implicit none
    private

    integer, parameter, public :: dense_spectrum_ok = 0
    integer, parameter, public :: dense_spectrum_invalid = -1
    integer, parameter, public :: dense_spectrum_allocation = -2

    public :: certify_dense_spectrum_inertia
    public :: certify_dense_spectrum_orthogonality
    public :: diagnose_dense_spectrum
    public :: dense_spectrum_is_certified
    public :: refine_dense_eigenpair
    public :: refine_dense_spectrum
    public :: unpermute_dense_vectors

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dsygst(itype, uplo, n, a, lda, b, ldb, info)
            import :: dp
            integer, intent(in) :: itype, n, lda, ldb
            character(len=1), intent(in) :: uplo
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(in) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dsygst
        subroutine dsytrd(uplo, n, a, lda, d, e, tau, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: d(*), e(*), tau(*), work(*)
            integer, intent(out) :: info
        end subroutine dsytrd
    end interface

contains

    subroutine certify_dense_spectrum_orthogonality(mass, eigenvectors, info)
        type(variable_block_tridiagonal_t), intent(in) :: mass
        real(dp), contiguous, intent(in) :: eigenvectors(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: gram(:, :), mass_images(:, :)
        integer :: allocation_status, index, n

        info = dense_spectrum_invalid
        n = size(eigenvectors, 1)
        if (n < 1 .or. size(eigenvectors, 2) /= n) return
        if (.not. all(ieee_is_finite(eigenvectors))) return
        info = dense_spectrum_allocation
        allocate (mass_images(n, n), gram(n, n), stat=allocation_status)
        if (allocation_status /= 0) return
        do index = 1, n
            call apply_variable_block_tridiagonal(mass, &
                eigenvectors(:, index), mass_images(:, index), info)
            if (info /= variable_block_ok) then
                info = dense_spectrum_invalid
                return
            end if
        end do
        gram = matmul(transpose(eigenvectors), mass_images)
        do index = 1, n
            gram(index, index) = gram(index, index) - 1.0_dp
        end do
        if (maxval(abs(gram)) > 64.0_dp * sqrt(epsilon(1.0_dp))) then
            info = dense_spectrum_invalid
            return
        end if
        info = dense_spectrum_ok
    end subroutine certify_dense_spectrum_orthogonality

    ! Inertia at every resolved gap midpoint.  The Cholesky congruence
    ! K - sigma M = U^T (C - sigma I) U, C = U^(-T) K U^(-1), preserves the
    ! count, and one orthogonal tridiagonal reduction of C gives Sturm counts
    ! in O(n) per shift instead of one dense factorization per gap.
    subroutine certify_dense_spectrum_inertia(stiffness, mass, eigenvalues, &
            residuals, resolutions, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: eigenvalues(:), residuals(:), resolutions(:)
        integer, intent(out) :: info
        real(dp), allocatable :: diagonal(:), offdiagonal(:)
        real(dp) :: gap, shift, uncertainty, roundoff
        integer :: index

        info = dense_spectrum_invalid
        if (size(eigenvalues) < 1) return
        if (size(residuals) /= size(eigenvalues)) return
        if (size(resolutions) /= size(eigenvalues)) return
        if (.not. all(ieee_is_finite(eigenvalues))) return
        if (.not. all(ieee_is_finite(residuals))) return
        if (.not. all(ieee_is_finite(resolutions))) return
        if (any(residuals < 0.0_dp) .or. any(resolutions < 0.0_dp)) return
        if (any(eigenvalues(2:) < eigenvalues(:size(eigenvalues) - 1))) return
        call congruent_tridiagonal(stiffness, mass, diagonal, offdiagonal, &
            info)
        if (info /= dense_spectrum_ok) return
        info = dense_spectrum_invalid
        if (.not. allocated(diagonal) .or. .not. allocated(offdiagonal)) return
        ! The Sturm count of the congruent tridiagonal matrix T is exact for
        ! a matrix within about n eps ||T|| of it; closer eigenvalue pairs,
        ! such as the parity-degenerate pairs of a coupled axisymmetric
        ! operator, form one cluster and are not separated by a probe.
        roundoff = 8.0_dp * real(size(diagonal), dp) * epsilon(1.0_dp) &
            * (maxval(abs(diagonal)) + 2.0_dp * maxval(abs(offdiagonal)))
        do index = 1, size(eigenvalues) - 1
            if (eigenvalues(index) == eigenvalues(index + 1)) cycle
            gap = eigenvalues(index + 1) - eigenvalues(index)
            uncertainty = residuals(index) + resolutions(index) &
                + residuals(index + 1) + resolutions(index + 1) + roundoff
            if (gap <= uncertainty) cycle
            shift = eigenvalues(index) &
                + 0.5_dp * gap
            if (shift <= eigenvalues(index)) cycle
            if (shift >= eigenvalues(index + 1)) cycle
            if (sturm_count(diagonal, offdiagonal, shift) /= index) return
        end do
        info = dense_spectrum_ok
    end subroutine certify_dense_spectrum_inertia

    subroutine congruent_tridiagonal(stiffness, mass, diagonal, offdiagonal, &
            info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), allocatable, intent(out) :: diagonal(:), offdiagonal(:)
        integer, intent(out) :: info
        real(dp), allocatable :: dense_k(:, :), dense_m(:, :), tau(:), work(:)
        integer :: allocation_status, n, status

        info = dense_spectrum_invalid
        call variable_block_to_dense(stiffness, dense_k, status)
        if (status /= variable_block_ok) return
        call variable_block_to_dense(mass, dense_m, status)
        if (status /= variable_block_ok) return
        n = size(dense_k, 1)
        info = dense_spectrum_allocation
        allocate (diagonal(n), offdiagonal(max(1, n - 1)), tau(max(1, n - 1)), &
            work(64 * n), stat=allocation_status)
        if (allocation_status /= 0) return
        ! DSYTRD leaves E untouched for n = 1; keep its dummy entry defined.
        offdiagonal = 0.0_dp
        info = dense_spectrum_invalid
        call dpotrf("U", n, dense_m, n, status)
        if (status /= 0) return
        call dsygst(1, "U", n, dense_k, n, dense_m, n, status)
        if (status /= 0) return
        call dsytrd("U", n, dense_k, n, diagonal, offdiagonal, tau, work, &
            size(work), status)
        if (status /= 0) return
        if (.not. all(ieee_is_finite(diagonal))) return
        if (.not. all(ieee_is_finite(offdiagonal))) return
        info = dense_spectrum_ok
    end subroutine congruent_tridiagonal

    ! Number of eigenvalues of the symmetric tridiagonal matrix below shift,
    ! from the signs of the LDL^T pivots with the usual tiny-pivot guard.
    pure function sturm_count(diagonal, offdiagonal, shift) result(count)
        real(dp), intent(in) :: diagonal(:), offdiagonal(:), shift
        integer :: count
        real(dp) :: pivot, guard
        integer :: i

        guard = tiny(1.0_dp) / epsilon(1.0_dp)
        count = 0
        pivot = diagonal(1) - shift
        if (abs(pivot) < guard) pivot = -guard
        if (pivot < 0.0_dp) count = 1
        do i = 2, size(diagonal)
            pivot = diagonal(i) - shift - offdiagonal(i - 1)**2 / pivot
            if (abs(pivot) < guard) pivot = -guard
            if (pivot < 0.0_dp) count = count + 1
        end do
    end function sturm_count

    subroutine refine_dense_spectrum(stiffness, mass, controls, eigenvalues, &
            eigenvectors, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        real(dp), intent(inout) :: eigenvalues(:)
        real(dp), contiguous, intent(inout) :: eigenvectors(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: seeds(:), vector(:), images(:, :)
        real(dp) :: eigenvalue, quotient, residual, resolution, roundoff
        real(dp) :: tolerance
        integer :: allocation_status, index, other
        logical :: keep

        info = dense_spectrum_invalid
        if (size(eigenvalues) < 1) return
        if (size(eigenvectors, 1) /= size(eigenvalues)) return
        if (size(eigenvectors, 2) /= size(eigenvalues)) return
        if (.not. all(ieee_is_finite(eigenvalues))) return
        if (.not. all(ieee_is_finite(eigenvectors))) return
        if (any(eigenvalues(2:) < eigenvalues(:size(eigenvalues) - 1))) return
        info = dense_spectrum_allocation
        allocate (seeds, source=eigenvalues, stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (images(size(eigenvalues), size(eigenvalues)), &
            stat=allocation_status)
        if (allocation_status /= 0) return
        tolerance = 64.0_dp * sqrt(epsilon(1.0_dp))
        roundoff = pencil_roundoff(stiffness, mass)
        do index = 1, size(eigenvalues)
            ! A dense pair is kept when it meets the inverse-iteration
            ! residual criterion and is mass orthonormal to every pair kept
            ! before it; the others are bracketed by index and refined.
            call variable_generalized_diagnostics(stiffness, mass, &
                eigenvectors(:, index), eigenvalues(index), eigenvalue, &
                residual, resolution, info, validated=index > 1)
            if (info /= variable_generalized_ok) then
                info = dense_spectrum_invalid
                return
            end if
            call apply_variable_block_tridiagonal(mass, &
                eigenvectors(:, index), images(:, index), info)
            if (info /= variable_block_ok) then
                info = dense_spectrum_invalid
                return
            end if
            keep = residual <= max(controls%residual_relative &
                * abs(eigenvalues(index)), roundoff, resolution)
            if (keep) keep = abs(dot_product(eigenvectors(:, index), &
                images(:, index)) - 1.0_dp) <= tolerance
            do other = 1, index - 1
                if (.not. keep) exit
                keep = abs(dot_product(eigenvectors(:, other), &
                    images(:, index))) <= tolerance
            end do
            if (keep) then
                ! The Rayleigh quotient of a kept vector is more accurate than
                ! the dense eigenvalue (quadratic in the vector error).
                eigenvalues(index) = eigenvalue
                cycle
            end if
            call refine_dense_eigenpair(stiffness, mass, controls, seeds, &
                index, eigenvectors(:, index), eigenvalue, vector, residual, &
                resolution, info)
            if (info /= dense_spectrum_ok) return
            ! Inverse iteration inside a degenerate or roundoff-split cluster
            ! can drift toward a partner already accepted. Eigenvectors are
            ! mass orthogonal, so the earlier ones are projected out (a
            ! roundoff-level change for separated eigenvalues) and the pair
            ! is diagnosed again.
            call orthonormalize_against(mass, eigenvectors, images, index, &
                vector, info)
            if (info /= dense_spectrum_ok) return
            call variable_generalized_diagnostics(stiffness, mass, vector, &
                eigenvalue, quotient, residual, resolution, info, &
                validated=.true.)
            if (info /= variable_generalized_ok) then
                info = dense_spectrum_invalid
                return
            end if
            eigenvalues(index) = quotient
            eigenvectors(:, index) = vector
            call apply_variable_block_tridiagonal(mass, &
                eigenvectors(:, index), images(:, index), info)
            if (info /= variable_block_ok) then
                info = dense_spectrum_invalid
                return
            end if
        end do
        ! Rayleigh quotients of a degenerate cluster (the two parities of a
        ! coupled axisymmetric operator) may leave the dense order by
        ! roundoff; restore ascending order for the inertia certificate.
        call sort_eigenpairs(eigenvalues, eigenvectors, vector)
        info = dense_spectrum_ok
    end subroutine refine_dense_spectrum

    ! Mass-orthogonalize vector against the first index - 1 columns of
    ! vectors (with their mass images) and normalize it in the mass norm;
    ! a vector lying in their span is rejected.
    subroutine orthonormalize_against(mass, vectors, images, index, vector, &
            info)
        type(variable_block_tridiagonal_t), intent(in) :: mass
        real(dp), intent(in) :: vectors(:, :), images(:, :)
        integer, intent(in) :: index
        real(dp), contiguous, intent(inout) :: vector(:)
        integer, intent(out) :: info
        real(dp) :: image(size(vector)), norm_before, squared
        integer :: other, pass

        info = dense_spectrum_invalid
        call apply_variable_block_tridiagonal(mass, vector, image, info)
        if (info /= variable_block_ok) then
            info = dense_spectrum_invalid
            return
        end if
        norm_before = sqrt(dot_product(vector, image))
        if (.not. (norm_before > 0.0_dp)) then
            info = dense_spectrum_invalid
            return
        end if
        ! Two passes of classical Gram-Schmidt reach orthogonality to working
        ! precision.
        do pass = 1, 2
            do other = 1, index - 1
                vector = vector - dot_product(images(:, other), vector) &
                    * vectors(:, other)
            end do
        end do
        call apply_variable_block_tridiagonal(mass, vector, image, info)
        if (info /= variable_block_ok) then
            info = dense_spectrum_invalid
            return
        end if
        squared = dot_product(vector, image)
        info = dense_spectrum_invalid
        if (.not. (squared > (sqrt(epsilon(1.0_dp)) * norm_before)**2)) return
        vector = vector / sqrt(squared)
        info = dense_spectrum_ok
    end subroutine orthonormalize_against

    pure subroutine sort_eigenpairs(eigenvalues, eigenvectors, column)
        real(dp), intent(inout) :: eigenvalues(:)
        real(dp), intent(inout) :: eigenvectors(:, :)
        real(dp), allocatable, intent(inout) :: column(:)
        real(dp) :: key
        integer :: index, position

        if (.not. allocated(column)) allocate (column(size(eigenvectors, 1)))
        if (size(column) /= size(eigenvectors, 1)) then
            deallocate (column)
            allocate (column(size(eigenvectors, 1)))
        end if
        do index = 2, size(eigenvalues)
            if (eigenvalues(index) >= eigenvalues(index - 1)) cycle
            key = eigenvalues(index)
            column = eigenvectors(:, index)
            position = index - 1
            do while (position >= 1)
                if (eigenvalues(position) <= key) exit
                eigenvalues(position + 1) = eigenvalues(position)
                eigenvectors(:, position + 1) = eigenvectors(:, position)
                position = position - 1
            end do
            eigenvalues(position + 1) = key
            eigenvectors(:, position + 1) = column
        end do
    end subroutine sort_eigenpairs

    subroutine refine_dense_eigenpair(stiffness, mass, controls, seeds, &
            target, seed_vector, eigenvalue, vector, residual, resolution, &
            info, interval)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        real(dp), intent(in) :: seeds(:), seed_vector(:)
        integer, intent(in) :: target
        real(dp), intent(out) :: eigenvalue, residual, resolution
        real(dp), allocatable, intent(out) :: vector(:)
        integer, intent(out) :: info
        real(dp), intent(out), optional :: interval
        real(dp), allocatable :: initial(:)
        real(dp) :: initial_scale, shift, width
        integer :: allocation_status, entry

        info = dense_spectrum_invalid
        if (target < 1 .or. target > size(seeds)) return
        if (size(seed_vector) /= size(seeds)) return
        if (.not. all(ieee_is_finite(seeds))) return
        if (.not. all(ieee_is_finite(seed_vector))) return
        call validate_variable_pencil(stiffness, mass, info)
        if (info /= variable_generalized_ok) then
            info = dense_spectrum_invalid
            return
        end if
        info = dense_spectrum_allocation
        allocate (initial, source=seed_vector, stat=allocation_status)
        if (allocation_status /= 0) return
        ! Every shift below factors this validated pencil.
        call bracket_indexed_eigenvalue(stiffness, mass, seeds, target, &
            controls, shift, info, width)
        if (info /= dense_spectrum_ok) return
        if (present(interval)) interval = width
        initial_scale = sqrt(epsilon(1.0_dp)) * norm2(seed_vector) &
            / sqrt(real(size(initial), dp))
        do entry = 1, size(initial)
            initial(entry) = initial(entry) + initial_scale &
                * (1.0_dp + 0.1_dp * real(entry, dp))
        end do
        call iterate_variable_generalized_eigenvalue(stiffness, mass, shift, &
            eigenvalue, vector, residual, resolution, info, controls, &
            initial=initial, validated=.true.)
        if (info /= variable_generalized_ok) then
            info = dense_spectrum_invalid
            return
        end if
        info = dense_spectrum_ok
    end subroutine refine_dense_eigenpair

    subroutine bracket_indexed_eigenvalue(stiffness, mass, seeds, target, &
            controls, shift, info, width)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: seeds(:)
        integer, intent(in) :: target
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        real(dp), intent(out) :: shift
        integer, intent(out) :: info
        real(dp), intent(out) :: width
        real(dp) :: center, lower, roundoff, step, tolerance, upper
        integer :: count, iteration, lower_count, upper_count

        info = dense_spectrum_invalid
        center = seeds(target)
        roundoff = pencil_roundoff(stiffness, mass)
        step = max(sqrt(epsilon(1.0_dp)) * abs(center), roundoff)
        if (target > 1) step = max(step, center - seeds(target - 1))
        if (target < size(seeds)) &
            step = max(step, seeds(target + 1) - center)
        do iteration = 1, controls%bracket_iteration_limit
            lower = center - step
            upper = center + step
            if (.not. ieee_is_finite(lower)) then
                info = dense_spectrum_invalid
                return
            end if
            if (.not. ieee_is_finite(upper)) then
                info = dense_spectrum_invalid
                return
            end if
            call directed_indexed_probe(stiffness, mass, lower, -1.0_dp, &
                lower_count, info)
            if (info /= dense_spectrum_ok) return
            call directed_indexed_probe(stiffness, mass, upper, 1.0_dp, &
                upper_count, info)
            if (info /= dense_spectrum_ok) return
            if (lower_count < target .and. upper_count >= target) exit
            step = 2.0_dp * step
        end do
        if (iteration > controls%bracket_iteration_limit) then
            info = dense_spectrum_invalid
            return
        end if
        do iteration = 1, controls%bracket_iteration_limit
            shift = lower + 0.5_dp * (upper - lower)
            tolerance = max(controls%eigenvalue_relative * abs(shift), &
                roundoff)
            if (upper - lower <= tolerance) exit
            call bounded_inertia_probe(stiffness, mass, lower, upper, &
                shift, count, info, validated=.true.)
            if (info /= fixed_boundary_bracket_ok) then
                info = dense_spectrum_invalid
                return
            end if
            if (count < target) then
                lower = shift
            else
                upper = shift
            end if
        end do
        if (iteration > controls%bracket_iteration_limit) then
            info = dense_spectrum_invalid
            return
        end if
        ! The lower endpoint has a successful factorization and fewer than
        ! target eigenvalues below it. An unprobed midpoint can be singular.
        shift = lower
        width = upper - lower
        info = dense_spectrum_ok
    end subroutine bracket_indexed_eigenvalue

    subroutine directed_indexed_probe(stiffness, mass, shift, direction, &
            count, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(inout) :: shift
        real(dp), intent(in) :: direction
        integer, intent(out) :: count, info
        real(dp) :: candidate, delta, origin
        integer :: attempt

        origin = shift
        delta = max(16.0_dp * epsilon(1.0_dp) * abs(origin), &
            pencil_roundoff(stiffness, mass))
        do attempt = 0, 15
            candidate = origin
            if (attempt > 0) candidate = origin + direction * delta
            if (.not. ieee_is_finite(candidate)) exit
            call variable_generalized_inertia(stiffness, mass, candidate, &
                count, info, validated=.true.)
            if (info == variable_generalized_ok) then
                shift = candidate
                info = dense_spectrum_ok
                return
            end if
            if (attempt > 0) delta = 16.0_dp * delta
        end do
        info = dense_spectrum_invalid
    end subroutine directed_indexed_probe

    subroutine diagnose_dense_spectrum(stiffness, mass, eigenvalues, &
            eigenvectors, rayleigh_quotients, residuals, resolutions, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), contiguous, intent(in) :: eigenvectors(:, :)
        real(dp), allocatable, intent(out) :: rayleigh_quotients(:)
        real(dp), allocatable, intent(out) :: residuals(:), resolutions(:)
        integer, intent(out) :: info
        type(variable_block_factor_t) :: mass_factor
        integer :: allocation_status, index, count, status

        info = dense_spectrum_invalid
        count = size(eigenvalues)
        if (count < 1) return
        if (size(eigenvectors, 1) /= count &
            .or. size(eigenvectors, 2) /= count) return
        info = dense_spectrum_allocation
        allocate (rayleigh_quotients(count), residuals(count), &
            resolutions(count), stat=allocation_status)
        if (allocation_status /= 0) return
        info = dense_spectrum_invalid
        call factorize_variable_shifted(mass, 0.0_dp, mass_factor, status)
        if (status /= variable_block_ok) return
        if (mass_factor%negative_count /= 0) return
        do index = 1, count
            call variable_generalized_diagnostics(stiffness, mass, &
                eigenvectors(:, index), eigenvalues(index), &
                rayleigh_quotients(index), residuals(index), &
                resolutions(index), info, validated=index > 1)
            if (info /= variable_generalized_ok) then
                info = dense_spectrum_invalid
                return
            end if
            call variable_eigenvalue_bound(stiffness, mass, &
                eigenvectors(:, index), eigenvalues(index), residuals(index), &
                info, mass_factor)
            if (info /= variable_generalized_ok) then
                info = dense_spectrum_invalid
                return
            end if
        end do
        info = dense_spectrum_ok
    end subroutine diagnose_dense_spectrum

    function dense_spectrum_is_certified(eigenvalues, rayleigh_quotients, &
            zero_floor, negative_count, floor_count, has_active, &
            lowest_active, certificate) result(valid)
        real(dp), intent(in) :: eigenvalues(:), rayleigh_quotients(:)
        real(dp), intent(in) :: zero_floor, lowest_active, certificate
        integer, intent(in) :: negative_count, floor_count
        logical, intent(in) :: has_active
        logical :: valid
        real(dp) :: scale, tolerance
        integer :: active

        valid = .false.
        if (size(eigenvalues) < 1) return
        if (size(rayleigh_quotients) /= size(eigenvalues)) return
        if (any(eigenvalues(2:) < eigenvalues(:size(eigenvalues) - 1))) return
        if (count(eigenvalues < -zero_floor) /= negative_count) return
        if (count(abs(eigenvalues) <= zero_floor) /= floor_count) return
        ! The dense spectrum spans the pencil, so its largest modulus sets
        ! the roundoff floor of every eigenvalue.
        scale = maxval(abs(eigenvalues))
        if (any(abs(eigenvalues - rayleigh_quotients) &
            > max(sqrt(epsilon(1.0_dp)) * max(abs(eigenvalues), &
            abs(rayleigh_quotients)), 1024.0_dp * epsilon(1.0_dp) * scale))) &
            return
        if (.not. has_active) then
            valid = .true.
            return
        end if
        active = 1
        if (negative_count == 0) active = floor_count + 1
        if (active > size(eigenvalues)) return
        tolerance = certificate + 16.0_dp * epsilon(1.0_dp) &
            * max(scale, abs(lowest_active))
        valid = abs(eigenvalues(active) - lowest_active) <= tolerance
    end function dense_spectrum_is_certified

    subroutine unpermute_dense_vectors(vectors, permutation, info)
        real(dp), intent(inout) :: vectors(:, :)
        integer, intent(in) :: permutation(:)
        integer, intent(out) :: info
        logical, allocatable :: visited(:)
        real(dp), allocatable :: temporary(:)
        real(dp) :: value
        integer :: allocation_status, column, current, destination, start

        info = dense_spectrum_invalid
        if (size(vectors, 1) /= size(permutation)) return
        info = dense_spectrum_allocation
        allocate (visited(size(permutation)), source=.false., &
            stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (temporary(size(vectors, 2)), stat=allocation_status)
        if (allocation_status /= 0) return
        do start = 1, size(permutation)
            destination = permutation(start)
            if (destination < 1 .or. destination > size(permutation)) then
                info = dense_spectrum_invalid
                return
            end if
            if (visited(destination)) then
                info = dense_spectrum_invalid
                return
            end if
            visited(destination) = .true.
        end do
        visited = .false.
        do start = 1, size(permutation)
            if (visited(start)) cycle
            temporary = vectors(start, :)
            current = start
            do
                destination = permutation(current)
                do column = 1, size(vectors, 2)
                    value = vectors(destination, column)
                    vectors(destination, column) = temporary(column)
                    temporary(column) = value
                end do
                visited(current) = .true.
                current = destination
                if (current == start) exit
            end do
        end do
        info = dense_spectrum_ok
    end subroutine unpermute_dense_vectors

end module dense_spectrum_support
