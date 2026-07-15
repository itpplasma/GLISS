module dense_spectrum_support
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use variable_block_tridiagonal, only: apply_variable_block_tridiagonal, &
        variable_block_ok, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, &
        variable_generalized_diagnostics, variable_generalized_inertia, &
        variable_generalized_ok
    implicit none
    private

    integer, parameter, public :: dense_spectrum_ok = 0
    integer, parameter, public :: dense_spectrum_invalid = -1
    integer, parameter, public :: dense_spectrum_allocation = -2

    public :: certify_dense_spectrum_inertia
    public :: certify_dense_spectrum_orthogonality
    public :: diagnose_dense_spectrum
    public :: dense_spectrum_is_certified
    public :: refine_dense_spectrum
    public :: unpermute_dense_vectors

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

    subroutine certify_dense_spectrum_inertia(stiffness, mass, eigenvalues, &
            residuals, resolutions, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: eigenvalues(:), residuals(:), resolutions(:)
        integer, intent(out) :: info
        real(dp) :: gap, shift, uncertainty
        integer :: count, index

        info = dense_spectrum_invalid
        if (size(eigenvalues) < 1) return
        if (size(residuals) /= size(eigenvalues)) return
        if (size(resolutions) /= size(eigenvalues)) return
        if (.not. all(ieee_is_finite(eigenvalues))) return
        if (.not. all(ieee_is_finite(residuals))) return
        if (.not. all(ieee_is_finite(resolutions))) return
        if (any(residuals < 0.0_dp) .or. any(resolutions < 0.0_dp)) return
        do index = 1, size(eigenvalues) - 1
            if (eigenvalues(index) > eigenvalues(index + 1)) return
            if (eigenvalues(index) == eigenvalues(index + 1)) cycle
            gap = eigenvalues(index + 1) - eigenvalues(index)
            uncertainty = residuals(index) + resolutions(index) &
                + residuals(index + 1) + resolutions(index + 1)
            if (gap <= uncertainty) cycle
            shift = eigenvalues(index) &
                + 0.5_dp * gap
            if (shift <= eigenvalues(index)) cycle
            if (shift >= eigenvalues(index + 1)) cycle
            call variable_generalized_inertia(stiffness, mass, shift, count, &
                info)
            if (info /= variable_generalized_ok) then
                info = dense_spectrum_invalid
                return
            end if
            if (count /= index) then
                info = dense_spectrum_invalid
                return
            end if
        end do
        info = dense_spectrum_ok
    end subroutine certify_dense_spectrum_inertia

    subroutine refine_dense_spectrum(stiffness, mass, controls, eigenvalues, &
            eigenvectors, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        real(dp), intent(inout) :: eigenvalues(:), eigenvectors(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: initial(:), seeds(:), vector(:)
        real(dp) :: eigenvalue, initial_scale, residual, resolution, shift
        integer :: allocation_status, entry, index

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
        allocate (initial(size(eigenvalues)), stat=allocation_status)
        if (allocation_status /= 0) return
        do index = 1, size(eigenvalues)
            call bracket_indexed_eigenvalue(stiffness, mass, seeds, index, &
                controls, shift, info)
            if (info /= dense_spectrum_ok) return
            initial = eigenvectors(:, index)
            initial_scale = sqrt(epsilon(1.0_dp)) &
                * norm2(eigenvectors(:, index)) &
                / sqrt(real(size(initial), dp))
            do entry = 1, size(initial)
                initial(entry) = initial(entry) + initial_scale &
                    * (1.0_dp + 0.1_dp * real(entry, dp))
            end do
            call iterate_variable_generalized_eigenvalue(stiffness, mass, &
                shift, eigenvalue, vector, residual, resolution, info, &
                controls, initial=initial)
            if (info /= variable_generalized_ok) then
                info = dense_spectrum_invalid
                return
            end if
            eigenvalues(index) = eigenvalue
            eigenvectors(:, index) = vector
        end do
        info = dense_spectrum_ok
    end subroutine refine_dense_spectrum

    subroutine bracket_indexed_eigenvalue(stiffness, mass, seeds, target, &
            controls, shift, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: seeds(:)
        integer, intent(in) :: target
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        real(dp), intent(out) :: shift
        integer, intent(out) :: info
        real(dp) :: center, lower, step, tolerance, upper
        integer :: count, iteration, lower_count, upper_count

        info = dense_spectrum_invalid
        center = seeds(target)
        step = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(center))
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
            call variable_generalized_inertia(stiffness, mass, lower, &
                lower_count, info)
            if (info /= variable_generalized_ok) return
            call variable_generalized_inertia(stiffness, mass, upper, &
                upper_count, info)
            if (info /= variable_generalized_ok) return
            if (lower_count < target .and. upper_count >= target) exit
            step = 2.0_dp * step
        end do
        if (iteration > controls%bracket_iteration_limit) then
            info = dense_spectrum_invalid
            return
        end if
        do iteration = 1, controls%bracket_iteration_limit
            shift = lower + 0.5_dp * (upper - lower)
            tolerance = controls%eigenvalue_relative &
                * max(1.0_dp, abs(shift))
            if (upper - lower <= tolerance) exit
            call variable_generalized_inertia(stiffness, mass, shift, count, &
                info)
            if (info /= variable_generalized_ok) return
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
        info = dense_spectrum_ok
    end subroutine bracket_indexed_eigenvalue

    subroutine diagnose_dense_spectrum(stiffness, mass, eigenvalues, &
            eigenvectors, rayleigh_quotients, residuals, resolutions, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), contiguous, intent(in) :: eigenvectors(:, :)
        real(dp), allocatable, intent(out) :: rayleigh_quotients(:)
        real(dp), allocatable, intent(out) :: residuals(:), resolutions(:)
        integer, intent(out) :: info
        integer :: allocation_status, index, count

        info = dense_spectrum_invalid
        count = size(eigenvalues)
        if (count < 1) return
        if (size(eigenvectors, 1) /= count &
            .or. size(eigenvectors, 2) /= count) return
        info = dense_spectrum_allocation
        allocate (rayleigh_quotients(count), residuals(count), &
            resolutions(count), stat=allocation_status)
        if (allocation_status /= 0) return
        do index = 1, count
            call variable_generalized_diagnostics(stiffness, mass, &
                eigenvectors(:, index), eigenvalues(index), &
                rayleigh_quotients(index), residuals(index), &
                resolutions(index), info)
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
        if (any(abs(eigenvalues - rayleigh_quotients) &
            > sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(eigenvalues), &
            abs(rayleigh_quotients)))) return
        if (.not. has_active) then
            valid = .true.
            return
        end if
        active = 1
        if (negative_count == 0) active = floor_count + 1
        if (active > size(eigenvalues)) return
        scale = max(1.0_dp, abs(lowest_active))
        tolerance = certificate + 16.0_dp * epsilon(1.0_dp) * scale
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
