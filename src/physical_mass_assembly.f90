module physical_mass_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mass_density_policy, only: evaluate_mass_density, mass_density_ok, &
        mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use phase_factor_topology, only: phase_cosine, &
        phase_product_coefficients, phase_sine
    use physical_mass_kernel, only: physical_mass_matrix
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: assemble_physical_mass_surface
    public :: assemble_physical_mass_surface_resolved

contains

    subroutine assemble_physical_mass_surface(fields, density_profile, &
            trial_m, trial_n, trial_parity, normal_stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            phase_assembly, mass, info)
        real(dp), intent(in) :: fields(:, :, :)
        type(mass_density_profile_t), intent(in) :: density_profile
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: density_kg_m3

        call evaluate_mass_density(density_profile, radial_coordinate, &
            density_kg_m3, info)
        if (info /= mass_density_ok) then
            info = -1
            return
        end if
        allocate (mass(4 * size(trial_m), 4 * size(trial_m)))
        call assemble_physical_mass_surface_resolved(fields, density_kg_m3, &
            trial_m, trial_n, trial_parity, normal_stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            phase_assembly, mass, info)
    end subroutine assemble_physical_mass_surface

    subroutine assemble_physical_mass_surface_resolved(fields, &
            density_kg_m3, trial_m, trial_n, trial_parity, &
            normal_stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, phase_assembly, mass, info)
        real(dp), intent(in) :: fields(:, :, :), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: info

        call validate_inputs(fields, density_kg_m3, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, &
            radial_step, phase_assembly, info)
        if (info /= 0) return
        if (size(mass, 1) /= 4 * size(trial_m)) then
            info = -1
            return
        end if
        if (size(mass, 2) /= 4 * size(trial_m)) then
            info = -1
            return
        end if
        mass = 0.0_dp
        if (phase_assembly == phase_assembly_direct) then
            call assemble_direct(fields, density_kg_m3, trial_m, trial_n, &
                trial_parity, normal_stored_power, field_periods, &
                radial_space, radial_coordinate, radial_step, mass, info)
        else
            call assemble_transformed(fields, density_kg_m3, trial_m, &
                trial_n, trial_parity, normal_stored_power, field_periods, &
                radial_space, radial_coordinate, radial_step, mass, info)
        end if
    end subroutine assemble_physical_mass_surface_resolved

    subroutine assemble_direct(fields, density_kg_m3, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, mass, info)
        real(dp), intent(in) :: fields(:, :, :), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: theta, zeta, weight
        integer :: j, k, period

        weight = radial_step / real(size(fields, 1) * size(fields, 2) &
            * field_periods, dp)
        do period = 0, field_periods - 1
            do k = 1, size(fields, 2)
                zeta = (real(k - 1, dp) / real(size(fields, 2), dp)) &
                    + real(period, dp)
                do j = 1, size(fields, 1)
                    theta = real(j - 1, dp) / real(size(fields, 1), dp)
                    call accumulate_direct_point(fields(j, k, :), &
                        density_kg_m3, trial_m, trial_n, trial_parity, &
                        stored_power, field_periods, radial_space, &
                        radial_coordinate, radial_step, theta, zeta, weight, &
                        mass, info)
                    if (info /= 0) return
                end do
            end do
        end do
        info = 0
    end subroutine assemble_direct

    subroutine accumulate_direct_point(fields, density_kg_m3, trial_m, &
            trial_n, trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, theta, zeta, weight, mass, info)
        real(dp), intent(in) :: fields(:), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step, theta, zeta
        real(dp), intent(in) :: weight
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: basis(3, 4 * size(trial_m)), point_mass(3, 3)
        real(dp) :: normal_values(2), normal_derivatives(2), phase
        real(dp) :: xi_phase, tangential_phase
        integer :: trial, trials

        call point_mass_matrix(fields, density_kg_m3, point_mass)
        trials = size(trial_m)
        basis = 0.0_dp
        do trial = 1, trials
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, 0.5_dp, normal_values, &
                normal_derivatives, info, stored_power(trial))
            if (info /= radial_space_ok) return
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta &
                / real(field_periods, dp))
            call phase_values(phase, trial_parity(trial), xi_phase, &
                tangential_phase)
            basis(1, trial) = normal_values(1) * xi_phase
            basis(1, trials + trial) = normal_values(2) * xi_phase
            basis(2, 2 * trials + trial) = tangential_phase
            basis(3, 3 * trials + trial) = tangential_phase
        end do
        call accumulate_basis_mass(basis, point_mass, weight, mass)
        info = 0
    end subroutine accumulate_direct_point

    subroutine assemble_transformed(fields, density_kg_m3, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, mass, info)
        real(dp), intent(in) :: fields(:, :, :), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: theta, zeta, weight
        integer :: j, k

        weight = radial_step / real(size(fields, 1) * size(fields, 2), dp)
        do k = 1, size(fields, 2)
            zeta = real(k - 1, dp) / real(size(fields, 2), dp)
            do j = 1, size(fields, 1)
                theta = real(j - 1, dp) / real(size(fields, 1), dp)
                call accumulate_transformed_point(fields(j, k, :), &
                    density_kg_m3, trial_m, trial_n, trial_parity, &
                    stored_power, field_periods, radial_space, &
                    radial_coordinate, radial_step, theta, zeta, weight, &
                    mass, info)
                if (info /= 0) return
            end do
        end do
        info = 0
    end subroutine assemble_transformed

    subroutine accumulate_transformed_point(fields, density_kg_m3, trial_m, &
            trial_n, trial_parity, stored_power, field_periods, radial_space, &
            radial_coordinate, radial_step, theta, zeta, weight, mass, info)
        real(dp), intent(in) :: fields(:), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step, theta, zeta
        real(dp), intent(in) :: weight
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: point_mass(3, 3), cosine(size(trial_m))
        real(dp) :: sine(size(trial_m)), factors(4, size(trial_m))

        call point_mass_matrix(fields, density_kg_m3, point_mass)
        call prepare_transformed_trials(trial_m, trial_n, stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            theta, zeta, factors, cosine, sine, info)
        if (info /= 0) return
        call accumulate_transformed_pairs(point_mass, factors, cosine, sine, &
            trial_n, trial_parity, field_periods, weight, mass)
        info = 0
    end subroutine accumulate_transformed_point

    subroutine prepare_transformed_trials(trial_m, trial_n, stored_power, &
            field_periods, radial_space, radial_coordinate, radial_step, &
            theta, zeta, factors, cosine, sine, info)
        integer, intent(in) :: trial_m(:), trial_n(:), field_periods
        real(dp), intent(in) :: stored_power(:)
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step, theta, zeta
        real(dp), intent(out) :: factors(:, :), cosine(:), sine(:)
        integer, intent(out) :: info
        real(dp) :: normal_values(2), normal_derivatives(2), phase
        integer :: trial

        do trial = 1, size(trial_m)
            call evaluate_normal_basis(radial_space, trial_m(trial), &
                radial_coordinate, radial_step, 0.5_dp, normal_values, &
                normal_derivatives, info, stored_power(trial))
            if (info /= radial_space_ok) return
            factors(:, trial) = [normal_values(1), normal_values(2), &
                1.0_dp, 1.0_dp]
            phase = two_pi * (real(trial_m(trial), dp) * theta &
                - real(trial_n(trial), dp) * zeta / real(field_periods, dp))
            cosine(trial) = cos(phase)
            sine(trial) = sin(phase)
        end do
        info = 0
    end subroutine prepare_transformed_trials

    subroutine accumulate_transformed_pairs(point_mass, factors, cosine, &
            sine, trial_n, trial_parity, field_periods, weight, mass)
        real(dp), intent(in) :: point_mass(:, :), factors(:, :)
        real(dp), intent(in) :: cosine(:), sine(:), weight
        integer, intent(in) :: trial_n(:), trial_parity(:), field_periods
        real(dp), intent(inout) :: mass(:, :)
        real(dp) :: products(2, 2), product
        integer :: a, b, block_a, block_b, component_a, component_b
        integer :: kind_a, kind_b, trial_a, trial_b, trials

        trials = size(trial_n)
        do trial_b = 1, trials
            do trial_a = 1, trials
                call phase_product_coefficients(cosine(trial_a), &
                    sine(trial_a), cosine(trial_b), sine(trial_b), &
                    trial_n(trial_a), trial_n(trial_b), field_periods, &
                    products)
                do block_b = 1, 4
                    b = (block_b - 1) * trials + trial_b
                    component_b = physical_component(block_b)
                    kind_b = component_phase_kind(block_b, &
                        trial_parity(trial_b))
                    do block_a = 1, 4
                        a = (block_a - 1) * trials + trial_a
                        component_a = physical_component(block_a)
                        kind_a = component_phase_kind(block_a, &
                            trial_parity(trial_a))
                        product = products(kind_a, kind_b)
                        mass(a, b) = mass(a, b) + weight * product &
                            * factors(block_a, trial_a) &
                            * factors(block_b, trial_b) &
                            * point_mass(component_a, component_b)
                    end do
                end do
            end do
        end do
    end subroutine accumulate_transformed_pairs

    pure subroutine point_mass_matrix(fields, density_kg_m3, mass)
        real(dp), intent(in) :: fields(:), density_kg_m3
        real(dp), intent(out) :: mass(3, 3)

        call physical_mass_matrix(fields(1), fields(2), fields(5), fields(6), &
            fields(7), fields(8), fields(9), fields(12), fields(13), &
            density_kg_m3, mass)
    end subroutine point_mass_matrix

    pure subroutine phase_values(phase, parity, xi_phase, tangential_phase)
        real(dp), intent(in) :: phase
        integer, intent(in) :: parity
        real(dp), intent(out) :: xi_phase, tangential_phase

        if (parity == 1) then
            xi_phase = cos(phase)
            tangential_phase = sin(phase)
        else
            xi_phase = sin(phase)
            tangential_phase = cos(phase)
        end if
    end subroutine phase_values

    pure function physical_component(block) result(component)
        integer, intent(in) :: block
        integer :: component

        if (block <= 2) then
            component = 1
        else
            component = block - 1
        end if
    end function physical_component

    pure function component_phase_kind(block, parity) result(kind)
        integer, intent(in) :: block, parity
        integer :: kind

        if (block <= 2) then
            kind = parity
        else if (parity == phase_cosine) then
            kind = phase_sine
        else
            kind = phase_cosine
        end if
    end function component_phase_kind

    subroutine accumulate_basis_mass(basis, point_mass, weight, mass)
        real(dp), intent(in) :: basis(:, :), point_mass(:, :), weight
        real(dp), intent(inout) :: mass(:, :)
        integer :: a, b, first, second

        do b = 1, size(mass, 2)
            do a = 1, size(mass, 1)
                do second = 1, 3
                    do first = 1, 3
                        mass(a, b) = mass(a, b) + weight &
                            * basis(first, a) * point_mass(first, second) &
                            * basis(second, b)
                    end do
                end do
            end do
        end do
    end subroutine accumulate_basis_mass

    subroutine validate_inputs(fields, density_kg_m3, trial_m, trial_n, &
            trial_parity, stored_power, field_periods, radial_step, &
            phase_assembly, info)
        real(dp), intent(in) :: fields(:, :, :), density_kg_m3
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        real(dp), intent(in) :: radial_step
        integer, intent(out) :: info
        integer :: trials

        info = -1
        trials = size(trial_m)
        if (trials < 1 .or. size(trial_n) /= trials) return
        if (size(trial_parity) /= trials .or. size(stored_power) /= trials) &
            return
        if (any(trial_m < 0)) return
        if (any(trial_parity < 1) .or. any(trial_parity > 2)) return
        if (.not. all(ieee_is_finite(stored_power))) return
        if (.not. ieee_is_finite(density_kg_m3)) return
        if (density_kg_m3 <= 0.0_dp) return
        if (field_periods < 1) return
        if (phase_assembly /= phase_assembly_transformed .and. &
            phase_assembly /= phase_assembly_direct) return
        if (size(fields, 1) < 1 .or. size(fields, 2) < 1) return
        if (size(fields, 3) < 13) return
        if (.not. all(ieee_is_finite(fields(:, :, 1:13)))) return
        if (.not. ieee_is_finite(radial_step) .or. radial_step <= 0.0_dp) &
            return
        if (.not. all(fields(:, :, 7) /= 0.0_dp)) return
        if (.not. all(fields(:, :, 8) > 0.0_dp)) return
        if (.not. all(fields(:, :, 9) > 0.0_dp)) return
        if (.not. all(fields(:, :, 1)**2 + fields(:, :, 2)**2 > 0.0_dp)) &
            return
        info = 0
    end subroutine validate_inputs

end module physical_mass_assembly
