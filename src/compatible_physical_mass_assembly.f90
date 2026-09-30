module compatible_physical_mass_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use period_averaged_assembly, only: accumulate_period_averaged, &
        period_masks
    use phase_factor_topology, only: phase_cosine, phase_sine
    use perpendicular_kinetic_kernel, only: perpendicular_kinetic_matrix
    use physical_mass_kernel, only: physical_mass_matrix
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: assemble_compatible_perpendicular_mass_surface
    public :: assemble_compatible_physical_mass_surface

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

contains

    subroutine assemble_compatible_physical_mass_surface(fields, &
            density_kg_m3, trial_m, trial_n, trial_parity, field_periods, &
            h1_values, eta_values, l2_values, radial_weight, phase_assembly, &
            mass, info)
        real(dp), intent(in) :: fields(:, :, :), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods, phase_assembly
        real(dp), intent(in) :: h1_values(:, :), eta_values(:, :)
        real(dp), intent(in) :: l2_values(:, :), radial_weight
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: angular_weight
        integer :: j, k, period

        call validate_inputs(fields, density_kg_m3, trial_m, trial_n, &
            trial_parity, field_periods, h1_values, l2_values, &
            radial_weight, phase_assembly, 2, mass, info)
        if (info /= 0) return
        info = -1
        if (any(shape(eta_values) /= shape(l2_values))) return
        if (.not. all(ieee_is_finite(eta_values))) return
        if (phase_assembly == phase_assembly_direct) then
            angular_weight = radial_weight / real(size(fields, 1) &
                * size(fields, 2) * field_periods, dp)
            do period = 0, field_periods - 1
                do k = 1, size(fields, 2)
                    do j = 1, size(fields, 1)
                        call accumulate_direct(fields, j, k, density_kg_m3, &
                            trial_m, trial_n, trial_parity, &
                            field_periods, h1_values, eta_values, l2_values, &
                            real(j - 1, dp) / real(size(fields, 1), dp), &
                            real(k - 1, dp) / real(size(fields, 2), dp) &
                            + real(period, dp), angular_weight, mass)
                    end do
                end do
            end do
        else
            call assemble_transformed(fields, density_kg_m3, trial_m, &
                trial_n, trial_parity, field_periods, h1_values, eta_values, &
                l2_values, radial_weight, 3, mass, info)
            if (info /= 0) return
        end if
        info = 0
    end subroutine assemble_compatible_physical_mass_surface

    subroutine assemble_compatible_perpendicular_mass_surface(fields, &
            density_kg_m3, trial_m, trial_n, trial_parity, field_periods, &
            h1_values, l2_values, radial_weight, phase_assembly, mass, info)
        real(dp), intent(in) :: fields(:, :, :), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        integer, intent(in) :: field_periods, phase_assembly
        real(dp), intent(in) :: h1_values(:, :), l2_values(:, :)
        real(dp), intent(in) :: radial_weight
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: angular_weight
        integer :: j, k, period

        call validate_inputs(fields, density_kg_m3, trial_m, trial_n, &
            trial_parity, field_periods, h1_values, l2_values, &
            radial_weight, phase_assembly, 1, mass, info)
        if (info /= 0) return
        if (phase_assembly == phase_assembly_direct) then
            angular_weight = radial_weight / real(size(fields, 1) &
                * size(fields, 2) * field_periods, dp)
            do period = 0, field_periods - 1
                do k = 1, size(fields, 2)
                    do j = 1, size(fields, 1)
                        call accumulate_perpendicular_direct(fields, j, k, &
                            density_kg_m3, trial_m, trial_n, trial_parity, &
                            field_periods, h1_values, l2_values, &
                            real(j - 1, dp) / real(size(fields, 1), dp), &
                            real(k - 1, dp) / real(size(fields, 2), dp) &
                            + real(period, dp), angular_weight, mass)
                    end do
                end do
            end do
        else
            call assemble_transformed(fields, density_kg_m3, trial_m, &
                trial_n, trial_parity, field_periods, h1_values, l2_values, &
                l2_values, radial_weight, 2, mass, info)
            if (info /= 0) return
        end if
        info = 0
    end subroutine assemble_compatible_perpendicular_mass_surface

    ! Period-averaged angular quadrature as matrix products: a factor
    ! point mass = V diag(lambda) V^T gives one weighted channel per column
    ! of V, whose response is V(:, c)^T of the basis coefficients.
    ! components = 3 is the physical (xi^s, eta, nu) form, 2 the
    ! perpendicular (xi^s, eta) form, which ignores eta.
    subroutine assemble_transformed(fields, density, trial_m, trial_n, &
            parity, field_periods, h1, eta, l2, radial_weight, components, &
            mass, info)
        real(dp), intent(in) :: fields(:, :, :), density
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:)
        integer, intent(in) :: field_periods, components
        real(dp), intent(in) :: h1(:, :), eta(:, :), l2(:, :), radial_weight
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        integer, parameter :: chunk_limit = 256
        real(dp), allocatable :: cosine_part(:, :, :), sine_part(:, :, :)
        real(dp), allocatable :: cosine_phase(:, :), sine_phase(:, :)
        real(dp), allocatable :: weight(:, :)
        real(dp) :: coefficients(components, 2, size(mass, 1))
        real(dp) :: point_mass(components, components), eigenvalues(components)
        real(dp) :: work(8 * components), angular_weight, phase, theta, zeta
        real(dp) :: trial_cosine(size(trial_m)), trial_sine(size(trial_m))
        real(dp), allocatable :: plus(:, :), minus(:, :)
        logical :: mixed
        integer :: channel, chunk, column, columns, count, first, j, k
        integer :: lapack_info, point, points, trial, trials

        info = -1
        columns = size(mass, 1)
        trials = size(trial_m)
        points = size(fields, 1) * size(fields, 2)
        angular_weight = radial_weight / real(points, dp)
        chunk = min(chunk_limit, points)
        allocate (cosine_part(chunk, columns, components), &
            sine_part(chunk, columns, components), &
            cosine_phase(chunk, columns), sine_phase(chunk, columns), &
            weight(chunk, components), plus(trials, trials), &
            minus(trials, trials))
        if (components == 3) then
            call build_basis_coefficients(parity, h1, eta, l2, coefficients)
        else
            call build_perpendicular_basis(parity, h1, l2, coefficients)
        end if
        call period_masks(trial_n, field_periods, plus, minus, mixed)
        do first = 1, points, chunk
            count = min(chunk, points - first + 1)
            ! Rows past the last point carry zero weight and zero response.
            cosine_part(count + 1:, :, :) = 0.0_dp
            sine_part(count + 1:, :, :) = 0.0_dp
            cosine_phase(count + 1:, :) = 0.0_dp
            sine_phase(count + 1:, :) = 0.0_dp
            weight(count + 1:, :) = 0.0_dp
            do point = 1, count
                j = modulo(first + point - 2, size(fields, 1)) + 1
                k = (first + point - 2) / size(fields, 1) + 1
                if (components == 3) then
                    call point_mass_matrix(fields, j, k, density, point_mass)
                else
                    call perpendicular_point_mass(fields, j, k, density, &
                        point_mass)
                end if
                call factor_point_mass(components, point_mass, eigenvalues, &
                    work, size(work), lapack_info)
                if (lapack_info /= 0) return
                weight(point, :) = angular_weight * eigenvalues
                theta = real(j - 1, dp) / real(size(fields, 1), dp)
                zeta = real(k - 1, dp) / real(size(fields, 2), dp)
                do trial = 1, trials
                    phase = two_pi * (real(trial_m(trial), dp) * theta &
                        - real(trial_n(trial), dp) * zeta &
                        / real(field_periods, dp))
                    trial_cosine(trial) = cos(phase)
                    trial_sine(trial) = sin(phase)
                end do
                do column = 1, columns
                    trial = modulo(column - 1, trials) + 1
                    cosine_phase(point, column) = trial_cosine(trial)
                    sine_phase(point, column) = trial_sine(trial)
                    do channel = 1, components
                        cosine_part(point, column, channel) = dot_product( &
                            point_mass(:, channel), &
                            coefficients(:, phase_cosine, column))
                        sine_part(point, column, channel) = dot_product( &
                            point_mass(:, channel), &
                            coefficients(:, phase_sine, column))
                    end do
                end do
            end do
            do channel = 1, components
                call accumulate_period_averaged(cosine_part(:, :, channel), &
                    sine_part(:, :, channel), cosine_phase, sine_phase, &
                    weight(:, channel), plus, minus, mixed, mass)
            end do
        end do
        info = 0
    end subroutine assemble_transformed

    subroutine accumulate_perpendicular_direct(fields, j, k, density, &
            trial_m, trial_n, parity, field_periods, h1, l2, theta, zeta, &
            weight, mass)
        real(dp), intent(in) :: fields(:, :, :), density, h1(:, :), l2(:, :)
        integer, intent(in) :: j, k
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: coefficients(2, 2, size(mass, 1))
        real(dp) :: basis(2, size(mass, 1)), point_mass(2, 2)
        real(dp) :: phase, cosine, sine
        integer :: column, trial, trials

        call build_perpendicular_basis(parity, h1, l2, coefficients)
        call perpendicular_point_mass(fields, j, k, density, point_mass)
        trials = size(trial_m)
        do column = 1, size(mass, 1)
            trial = modulo(column - 1, trials) + 1
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            cosine = cos(phase)
            sine = sin(phase)
            basis(:, column) = coefficients(:, phase_cosine, column) * cosine &
                + coefficients(:, phase_sine, column) * sine
        end do
        call rank_update(basis, point_mass, weight, mass)
    end subroutine accumulate_perpendicular_direct

    pure subroutine build_perpendicular_basis(parity, h1, l2, coefficients)
        integer, intent(in) :: parity(:)
        real(dp), intent(in) :: h1(:, :), l2(:, :)
        real(dp), intent(out) :: coefficients(:, :, :)
        integer :: basis, column, kind, trial, trials

        coefficients = 0.0_dp
        trials = size(parity)
        do basis = 1, size(h1, 1)
            do trial = 1, trials
                column = (basis - 1) * trials + trial
                coefficients(1, parity(trial), column) = h1(basis, trial)
            end do
        end do
        do basis = 1, size(l2, 1)
            do trial = 1, trials
                kind = phase_cosine
                if (parity(trial) == phase_cosine) kind = phase_sine
                column = size(h1, 1) * trials &
                    + (basis - 1) * trials + trial
                coefficients(2, kind, column) = l2(basis, trial)
            end do
        end do
    end subroutine build_perpendicular_basis

    pure subroutine perpendicular_point_mass(fields, j, k, density, mass)
        real(dp), intent(in) :: fields(:, :, :), density
        integer, intent(in) :: j, k
        real(dp), intent(out) :: mass(2, 2)

        call perpendicular_kinetic_matrix(fields(j, k, 7), &
            fields(j, k, 8), fields(j, k, 9), fields(j, k, 12), density, mass)
    end subroutine perpendicular_point_mass

    subroutine accumulate_direct(fields, j, k, density, trial_m, trial_n, &
            parity, field_periods, h1, eta, l2, theta, zeta, weight, mass)
        real(dp), intent(in) :: fields(:, :, :), density, h1(:, :), l2(:, :)
        real(dp), intent(in) :: eta(:, :)
        integer, intent(in) :: j, k
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: coefficients(3, 2, size(mass, 1))
        real(dp) :: basis(3, size(mass, 1)), point_mass(3, 3)
        real(dp) :: phase, cosine, sine
        integer :: column, trial, trials

        call build_basis_coefficients(parity, h1, eta, l2, coefficients)
        call point_mass_matrix(fields, j, k, density, point_mass)
        trials = size(trial_m)
        do column = 1, size(mass, 1)
            trial = modulo(column - 1, trials) + 1
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            cosine = cos(phase)
            sine = sin(phase)
            basis(:, column) = coefficients(:, phase_cosine, column) * cosine &
                + coefficients(:, phase_sine, column) * sine
        end do
        call rank_update(basis, point_mass, weight, mass)
    end subroutine accumulate_direct

    pure subroutine build_basis_coefficients(parity, h1, eta, l2, &
            coefficients)
        integer, intent(in) :: parity(:)
        real(dp), intent(in) :: h1(:, :), eta(:, :), l2(:, :)
        real(dp), intent(out) :: coefficients(:, :, :)
        integer :: basis, column, kind, trial, trials

        coefficients = 0.0_dp
        trials = size(parity)
        do basis = 1, size(h1, 1)
            do trial = 1, trials
                column = (basis - 1) * trials + trial
                coefficients(1, parity(trial), column) = h1(basis, trial)
            end do
        end do
        do basis = 1, size(l2, 1)
            do trial = 1, trials
                if (parity(trial) == phase_cosine) then
                    kind = phase_sine
                else
                    kind = phase_cosine
                end if
                column = size(h1, 1) * trials &
                    + (basis - 1) * trials + trial
                coefficients(2, kind, column) = eta(basis, trial)
                column = (size(h1, 1) + size(l2, 1)) * trials &
                    + (basis - 1) * trials + trial
                coefficients(3, kind, column) = l2(basis, trial)
            end do
        end do
    end subroutine build_basis_coefficients

    pure subroutine point_mass_matrix(fields, j, k, density, mass)
        real(dp), intent(in) :: fields(:, :, :), density
        integer, intent(in) :: j, k
        real(dp), intent(out) :: mass(3, 3)

        call physical_mass_matrix(fields(j, k, 1), fields(j, k, 2), &
            fields(j, k, 5), fields(j, k, 6), fields(j, k, 7), &
            fields(j, k, 8), fields(j, k, 9), fields(j, k, 12), &
            fields(j, k, 13), density, mass)
        call apply_parallel_unknown(fields(j, k, 2) / fields(j, k, 1) &
            * fields(j, k, 7), mass)
    end subroutine point_mass_matrix

    ! The third unknown is nu=mu-c eta with c=(FP'/FT') sqrtg, which is
    ! proportional to sqrtg xi^zeta and regular at the axis.  The kinetic
    ! form in (xi^s, eta, nu) is T^T M T with mu=nu+c eta.
    pure subroutine apply_parallel_unknown(c, mass)
        real(dp), intent(in) :: c
        real(dp), intent(inout) :: mass(3, 3)

        mass(:, 2) = mass(:, 2) + c * mass(:, 3)
        mass(2, :) = mass(2, :) + c * mass(3, :)
    end subroutine apply_parallel_unknown

    subroutine rank_update(basis, point_mass, weight, mass)
        real(dp), intent(in) :: basis(:, :), point_mass(:, :), weight
        real(dp), intent(inout) :: mass(:, :)
        integer :: a, b

        do b = 1, size(mass, 2)
            do a = 1, size(mass, 1)
                mass(a, b) = mass(a, b) + weight &
                    * bilinear(basis(:, a), point_mass, basis(:, b))
            end do
        end do
    end subroutine rank_update

    pure function bilinear(first, matrix, second) result(value)
        real(dp), intent(in) :: first(:), matrix(:, :), second(:)
        real(dp) :: value
        integer :: i, j

        value = 0.0_dp
        do j = 1, size(second)
            do i = 1, size(first)
                value = value + first(i) * matrix(i, j) * second(j)
            end do
        end do
    end function bilinear

    ! point_mass = V diag(lambda) V^T with V overwriting point_mass: the
    ! Cholesky factor (lambda = 1) of the positive definite kinetic form,
    ! or its eigenpairs when it is only semidefinite.
    subroutine factor_point_mass(n, point_mass, lambda, work, lwork, info)
        integer, intent(in) :: n, lwork
        real(dp), intent(inout) :: point_mass(n, n)
        real(dp), intent(out) :: lambda(n), work(lwork)
        integer, intent(out) :: info
        real(dp) :: factor(n, n), pivot
        integer :: i, j

        factor = 0.0_dp
        do j = 1, n
            pivot = point_mass(j, j) - sum(factor(j, 1:j - 1)**2)
            if (pivot <= 1.0e-12_dp * maxval(abs(point_mass))) then
                call dsyev("V", "U", n, point_mass, n, lambda, work, lwork, &
                    info)
                return
            end if
            factor(j, j) = sqrt(pivot)
            do i = j + 1, n
                factor(i, j) = (point_mass(i, j) &
                    - sum(factor(i, 1:j - 1) * factor(j, 1:j - 1))) / factor(j, j)
            end do
        end do
        point_mass = factor
        lambda = 1.0_dp
        info = 0
    end subroutine factor_point_mass

    subroutine validate_inputs(fields, density, trial_m, trial_n, parity, &
            field_periods, h1, l2, radial_weight, phase_assembly, &
            tangential_components, mass, info)
        real(dp), intent(in) :: fields(:, :, :), density, h1(:, :), l2(:, :)
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:)
        integer, intent(in) :: field_periods, phase_assembly
        integer, intent(in) :: tangential_components
        real(dp), intent(in) :: radial_weight, mass(:, :)
        integer, intent(out) :: info
        integer :: expected, trials

        info = -1
        trials = size(trial_m)
        if (trials < 1 .or. field_periods < 1) return
        if (size(trial_n) /= trials .or. size(parity) /= trials) return
        if (any(trial_m < 0) .or. any(parity < 1) .or. any(parity > 2)) return
        if (size(h1, 1) < 1 .or. size(l2, 1) < 1) return
        if (size(h1, 2) /= trials .or. size(l2, 2) /= trials) return
        if (tangential_components < 1 .or. tangential_components > 2) return
        expected = trials * (size(h1, 1) &
            + tangential_components * size(l2, 1))
        if (any(shape(mass) /= expected)) return
        if (size(fields, 1) < 1 .or. size(fields, 2) < 1 &
            .or. size(fields, 3) < 13) return
        if (.not. all(ieee_is_finite(fields(:, :, 1:13))) &
            .or. .not. all(ieee_is_finite(h1)) &
            .or. .not. all(ieee_is_finite(l2))) return
        if (.not. ieee_is_finite(density) .or. density <= 0.0_dp) return
        if (.not. ieee_is_finite(radial_weight) .or. radial_weight <= 0.0_dp) &
            return
        if (phase_assembly /= phase_assembly_direct .and. &
            phase_assembly /= phase_assembly_transformed) return
        if (any(fields(:, :, 7) == 0.0_dp) &
            .or. any(fields(:, :, 8) <= 0.0_dp) &
            .or. any(fields(:, :, 9) <= 0.0_dp)) return
        if (any(fields(:, :, 1)**2 + fields(:, :, 2)**2 <= 0.0_dp)) return
        info = 0
    end subroutine validate_inputs

end module compatible_physical_mass_assembly
