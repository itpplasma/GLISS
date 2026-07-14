module compatible_physical_mass_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use phase_factor_topology, only: phase_cosine, &
        phase_product_coefficients, phase_sine
    use perpendicular_kinetic_kernel, only: perpendicular_kinetic_matrix
    use physical_mass_kernel, only: physical_mass_matrix
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: assemble_compatible_perpendicular_mass_surface
    public :: assemble_compatible_physical_mass_surface

contains

    subroutine assemble_compatible_physical_mass_surface(fields, &
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
            radial_weight, phase_assembly, 2, mass, info)
        if (info /= 0) return
        if (phase_assembly == phase_assembly_direct) then
            angular_weight = radial_weight / real(size(fields, 1) &
                * size(fields, 2) * field_periods, dp)
            do period = 0, field_periods - 1
                do k = 1, size(fields, 2)
                    do j = 1, size(fields, 1)
                        call accumulate_direct(fields, j, k, density_kg_m3, &
                            trial_m, trial_n, trial_parity, &
                            field_periods, h1_values, l2_values, &
                            real(j - 1, dp) / real(size(fields, 1), dp), &
                            real(k - 1, dp) / real(size(fields, 2), dp) &
                            + real(period, dp), angular_weight, mass)
                    end do
                end do
            end do
        else
            angular_weight = radial_weight / real(size(fields, 1) &
                * size(fields, 2), dp)
            do k = 1, size(fields, 2)
                do j = 1, size(fields, 1)
                    call accumulate_transformed(fields, j, k, density_kg_m3, &
                        trial_m, trial_n, trial_parity, &
                        field_periods, h1_values, l2_values, &
                        real(j - 1, dp) / real(size(fields, 1), dp), &
                        real(k - 1, dp) / real(size(fields, 2), dp), &
                        angular_weight, mass)
                end do
            end do
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
            angular_weight = radial_weight / real(size(fields, 1) &
                * size(fields, 2), dp)
            do k = 1, size(fields, 2)
                do j = 1, size(fields, 1)
                    call accumulate_perpendicular_transformed(fields, j, k, &
                        density_kg_m3, trial_m, trial_n, trial_parity, &
                        field_periods, h1_values, l2_values, &
                        real(j - 1, dp) / real(size(fields, 1), dp), &
                        real(k - 1, dp) / real(size(fields, 2), dp), &
                        angular_weight, mass)
                end do
            end do
        end if
        info = 0
    end subroutine assemble_compatible_perpendicular_mass_surface

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

    subroutine accumulate_perpendicular_transformed(fields, j, k, density, &
            trial_m, trial_n, parity, field_periods, h1, l2, theta, zeta, &
            weight, mass)
        real(dp), intent(in) :: fields(:, :, :), density, h1(:, :), l2(:, :)
        integer, intent(in) :: j, k
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: coefficients(2, 2, size(mass, 1)), point_mass(2, 2)
        real(dp) :: cosine(size(trial_m)), sine(size(trial_m)), phase
        integer :: trial

        call build_perpendicular_basis(parity, h1, l2, coefficients)
        call perpendicular_point_mass(fields, j, k, density, point_mass)
        do trial = 1, size(trial_m)
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            cosine(trial) = cos(phase)
            sine(trial) = sin(phase)
        end do
        call transformed_rank_update(coefficients, point_mass, cosine, sine, &
            trial_n, field_periods, weight, mass)
    end subroutine accumulate_perpendicular_transformed

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
            parity, field_periods, h1, l2, theta, zeta, weight, mass)
        real(dp), intent(in) :: fields(:, :, :), density, h1(:, :), l2(:, :)
        integer, intent(in) :: j, k
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: coefficients(3, 2, size(mass, 1))
        real(dp) :: basis(3, size(mass, 1)), point_mass(3, 3)
        real(dp) :: phase, cosine, sine
        integer :: column, trial, trials

        call build_basis_coefficients(parity, h1, l2, coefficients)
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

    subroutine accumulate_transformed(fields, j, k, density, trial_m, &
            trial_n, parity, field_periods, h1, l2, theta, zeta, weight, mass)
        real(dp), intent(in) :: fields(:, :, :), density, h1(:, :), l2(:, :)
        integer, intent(in) :: j, k
        integer, intent(in) :: trial_m(:), trial_n(:), parity(:), field_periods
        real(dp), intent(in) :: theta, zeta, weight
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: coefficients(3, 2, size(mass, 1)), point_mass(3, 3)
        real(dp) :: cosine(size(trial_m)), sine(size(trial_m)), phase
        integer :: trial

        call build_basis_coefficients(parity, h1, l2, coefficients)
        call point_mass_matrix(fields, j, k, density, point_mass)
        do trial = 1, size(trial_m)
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            cosine(trial) = cos(phase)
            sine(trial) = sin(phase)
        end do
        call transformed_rank_update(coefficients, point_mass, cosine, sine, &
            trial_n, field_periods, weight, mass)
    end subroutine accumulate_transformed

    pure subroutine build_basis_coefficients(parity, h1, l2, coefficients)
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
                if (parity(trial) == phase_cosine) then
                    kind = phase_sine
                else
                    kind = phase_cosine
                end if
                column = size(h1, 1) * trials &
                    + (basis - 1) * trials + trial
                coefficients(2, kind, column) = l2(basis, trial)
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
    end subroutine point_mass_matrix

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

    subroutine transformed_rank_update(coefficients, point_mass, cosine, &
            sine, trial_n, field_periods, weight, mass)
        real(dp), intent(in) :: coefficients(:, :, :), point_mass(:, :)
        real(dp), intent(in) :: cosine(:), sine(:), weight
        integer, intent(in) :: trial_n(:), field_periods
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: products(2, 2), contribution
        integer :: a, b, kind_a, kind_b, trial_a, trial_b, trials

        trials = size(trial_n)
        do b = 1, size(mass, 2)
            trial_b = modulo(b - 1, trials) + 1
            do a = 1, size(mass, 1)
                trial_a = modulo(a - 1, trials) + 1
                call phase_product_coefficients(cosine(trial_a), &
                    sine(trial_a), cosine(trial_b), sine(trial_b), &
                    trial_n(trial_a), trial_n(trial_b), field_periods, products)
                contribution = 0.0_dp
                do kind_b = phase_cosine, phase_sine
                    do kind_a = phase_cosine, phase_sine
                        contribution = contribution &
                            + products(kind_a, kind_b) * bilinear( &
                            coefficients(:, kind_a, a), &
                            point_mass, coefficients(:, kind_b, b))
                    end do
                end do
                mass(a, b) = mass(a, b) + weight * contribution
            end do
        end do
    end subroutine transformed_rank_update

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
