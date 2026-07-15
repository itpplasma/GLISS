program test_compressible_stiffness_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compressible_stiffness_assembly, only: &
        assemble_compressible_stiffness_surface
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use radial_space_policy, only: evaluate_normal_basis, &
        radial_space_config_t, radial_space_ok
    use three_component_kernel, only: three_component_density_value
    implicit none

    integer, parameter :: n_theta = 12, n_zeta = 10
    integer, parameter :: trial_m(3) = [1, 2, 1]
    integer, parameter :: trial_n(3) = [1, 4, -2]
    integer, parameter :: trial_parity(3) = [1, 2, 1]
    real(dp), parameter :: stored_power(3) = 0.0_dp
    real(dp), allocatable :: fields(:, :, :), drive(:, :)
    real(dp), allocatable :: jacobian_radial(:, :), jacobian_theta(:, :)
    real(dp), allocatable :: jacobian_zeta(:, :), gamma_pressure(:, :)
    real(dp), allocatable :: direct(:, :), transformed(:, :), zero(:, :)
    real(dp), allocatable :: doubled(:, :), invalid(:, :)
    real(dp), allocatable :: doubled_gamma(:, :), scaled_rows(:, :)
    real(dp), allocatable :: zero_gamma(:, :)
    type(radial_space_config_t) :: radial_space
    real(dp) :: probe(12), matrix_energy, oracle_energy
    integer :: info

    call build_fixture(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, gamma_pressure)
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, gamma_pressure, phase_assembly_direct, direct, info)
    call require(info == 0, "direct compressible stiffness assembly failed")
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, gamma_pressure, phase_assembly_transformed, &
        transformed, info)
    call require(info == 0, &
        "transformed compressible stiffness assembly failed")
    call require_close(direct, transformed, 1.0e-12_dp, &
        "one-period and all-period compressible stiffness matrices differ")
    call require_symmetric(transformed, 1.0e-13_dp, &
        "compressible stiffness element is not symmetric")
    call require(maxval(abs(transformed(10:12, :))) > 1.0e-8_dp, &
        "finite compression did not retain the mu block")

    allocate (zero_gamma, mold=gamma_pressure)
    zero_gamma = 0.0_dp
    doubled_gamma = 2.0_dp * gamma_pressure
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, zero_gamma, &
        phase_assembly_transformed, zero, info)
    call require(info == 0, "zero-compression stiffness assembly failed")
    call require(maxval(abs(zero(10:12, :))) < 1.0e-13_dp, &
        "gamma=0 stiffness retained a mu coupling")
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, doubled_gamma, &
        phase_assembly_transformed, doubled, info)
    call require(info == 0, "scaled-compression stiffness assembly failed")
    scaled_rows = 2.0_dp * transformed(10:12, :)
    call require_close(doubled(10:12, :), scaled_rows, &
        1.0e-13_dp, &
        "compressional stiffness is not linear in gamma pressure")

    probe = [0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.25_dp, 0.15_dp, &
        0.05_dp, -0.2_dp, 0.1_dp, -0.08_dp, 0.12_dp, 0.07_dp]
    matrix_energy = quadratic_form(direct, probe)
    call calculate_direct_energy(fields, drive, jacobian_radial, &
        jacobian_theta, jacobian_zeta, gamma_pressure, radial_space, probe, &
        oracle_energy, info)
    call require(info == 0, "independent direct-energy oracle failed")
    call require(abs(matrix_energy - oracle_energy) < 1.0e-12_dp &
        * max(1.0_dp, abs(matrix_energy), abs(oracle_energy)), &
        "assembled stiffness disagrees with point-energy quadrature")
    call check_analytic_product_derivatives(radial_space)
    call check_analytic_phase_parity(radial_space)

    fields(:, :, 7) = 0.0_dp
    call assemble(fields, drive, jacobian_radial, jacobian_theta, &
        jacobian_zeta, gamma_pressure, phase_assembly_transformed, invalid, &
        info)
    call require(info /= 0, "zero signed Jacobian was accepted")

    write (*, "(a)") "PASS"

contains

    pure function quadratic_form(matrix, vector) result(value)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp) :: value, image
        integer :: column, row

        value = 0.0_dp
        do row = 1, size(matrix, 1)
            image = 0.0_dp
            do column = 1, size(matrix, 2)
                image = image + matrix(row, column) * vector(column)
            end do
            value = value + vector(row) * image
        end do
    end function quadratic_form

    subroutine require_symmetric(matrix, tolerance, message)
        real(dp), intent(in) :: matrix(:, :), tolerance
        character(len=*), intent(in) :: message
        real(dp) :: error, scale
        integer :: column, row

        error = 0.0_dp
        do column = 1, size(matrix, 2)
            do row = 1, column - 1
                error = max(error, abs(matrix(row, column) &
                    - matrix(column, row)))
            end do
        end do
        scale = max(1.0_dp, maxval(abs(matrix)))
        call require(error <= tolerance * scale, message)
    end subroutine require_symmetric

    subroutine assemble(local_fields, local_drive, local_jacobian_radial, &
            local_jacobian_theta, local_jacobian_zeta, local_gamma_pressure, &
            phase_assembly, stiffness, info)
        real(dp), intent(in) :: local_fields(:, :, :), local_drive(:, :)
        real(dp), intent(in) :: local_jacobian_radial(:, :)
        real(dp), intent(in) :: local_jacobian_theta(:, :)
        real(dp), intent(in) :: local_jacobian_zeta(:, :)
        real(dp), intent(in) :: local_gamma_pressure(:, :)
        integer, intent(in) :: phase_assembly
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        integer, intent(out) :: info

        call assemble_compressible_stiffness_surface(local_fields, &
            local_drive, local_jacobian_radial, local_jacobian_theta, &
            local_jacobian_zeta, local_gamma_pressure, trial_m, trial_n, &
            trial_parity, stored_power, 3, radial_space, 0.375_dp, 0.25_dp, &
            phase_assembly, stiffness, info)
    end subroutine assemble

    subroutine build_fixture(local_fields, local_drive, local_jacobian_radial, &
            local_jacobian_theta, local_jacobian_zeta, local_gamma_pressure)
        real(dp), allocatable, intent(out) :: local_fields(:, :, :)
        real(dp), allocatable, intent(out) :: local_drive(:, :)
        real(dp), allocatable, intent(out) :: local_jacobian_radial(:, :)
        real(dp), allocatable, intent(out) :: local_jacobian_theta(:, :)
        real(dp), allocatable, intent(out) :: local_jacobian_zeta(:, :)
        real(dp), allocatable, intent(out) :: local_gamma_pressure(:, :)
        real(dp) :: theta, zeta, two_pi
        integer :: j, k

        two_pi = 2.0_dp * acos(-1.0_dp)
        allocate (local_fields(n_theta, n_zeta, 13), source=0.0_dp)
        allocate (local_drive(n_theta, n_zeta))
        allocate (local_jacobian_radial(n_theta, n_zeta))
        allocate (local_jacobian_theta(n_theta, n_zeta))
        allocate (local_jacobian_zeta(n_theta, n_zeta))
        allocate (local_gamma_pressure(n_theta, n_zeta))
        do k = 1, n_zeta
            zeta = two_pi * real(k - 1, dp) / real(n_zeta, dp)
            do j = 1, n_theta
                theta = two_pi * real(j - 1, dp) / real(n_theta, dp)
                local_fields(j, k, 1) = 1.2_dp + 0.05_dp * cos(theta)
                local_fields(j, k, 2) = 0.7_dp + 0.03_dp * sin(zeta)
                local_fields(j, k, 3) = 0.04_dp
                local_fields(j, k, 4) = -0.03_dp
                local_fields(j, k, 5) = 0.8_dp + 0.02_dp * cos(zeta)
                local_fields(j, k, 6) = 0.6_dp - 0.01_dp * sin(theta)
                local_fields(j, k, 7) = -1.1_dp &
                    - 0.02_dp * cos(theta) - 0.03_dp * sin(zeta)
                local_fields(j, k, 8) = 1.4_dp + 0.04_dp * sin(theta + zeta)
                local_fields(j, k, 9) = 1.3_dp + 0.03_dp * cos(zeta)
                local_fields(j, k, 10) = 0.2_dp - 0.01_dp * cos(theta)
                local_fields(j, k, 11) = -0.15_dp
                local_fields(j, k, 12) = 0.1_dp - 0.01_dp * sin(theta)
                local_fields(j, k, 13) = -0.12_dp + 0.02_dp * cos(zeta)
                local_drive(j, k) = 0.05_dp + 0.01_dp * cos(theta - zeta)
                local_jacobian_radial(j, k) = -0.08_dp &
                    + 0.01_dp * sin(theta + zeta)
                local_jacobian_theta(j, k) = 0.02_dp * sin(theta) * two_pi
                local_jacobian_zeta(j, k) = -0.03_dp * cos(zeta) * two_pi
                local_gamma_pressure(j, k) = 0.9_dp &
                    + 0.1_dp * cos(theta + zeta)
            end do
        end do
    end subroutine build_fixture

    subroutine check_analytic_product_derivatives(local_radial_space)
        type(radial_space_config_t), intent(in) :: local_radial_space
        integer, parameter :: modes(2) = [0, 0], parities(2) = [1, 2]
        real(dp), parameter :: powers(2) = 0.0_dp
        real(dp), parameter :: jacobian = -1.1_dp, jacobian_s = -0.08_dp
        real(dp), parameter :: jacobian_t = 0.03_dp, jacobian_z = -0.02_dp
        real(dp), parameter :: gamma_p = 2.0_dp, step = 0.25_dp
        real(dp) :: constant_fields(1, 1, 13), local_drive(1, 1)
        real(dp) :: local_js(1, 1), local_jt(1, 1), local_jz(1, 1)
        real(dp) :: local_gamma(1, 1), zero_local_gamma(1, 1)
        real(dp) :: expected(8, 8), difference(8, 8), slopes(2), eta_div
        real(dp), allocatable :: finite(:, :), zero_gamma(:, :)
        integer :: local_info

        constant_fields = 0.0_dp
        constant_fields(1, 1, 1:2) = [1.2_dp, 0.7_dp]
        constant_fields(1, 1, 7:9) = [jacobian, 1.4_dp, 1.3_dp]
        local_drive = 0.0_dp
        local_js = jacobian_s
        local_jt = jacobian_t
        local_jz = jacobian_z
        local_gamma = gamma_p
        zero_local_gamma = 0.0_dp
        call assemble_compressible_stiffness_surface(constant_fields, &
            local_drive, local_js, local_jt, local_jz, local_gamma, modes, &
            modes, parities, powers, 1, local_radial_space, 0.375_dp, step, &
            phase_assembly_transformed, finite, local_info)
        call require(local_info == 0, "analytic product gate failed to assemble")
        call assemble_compressible_stiffness_surface(constant_fields, &
            local_drive, local_js, local_jt, local_jz, zero_local_gamma, &
            modes, modes, parities, powers, 1, local_radial_space, 0.375_dp, &
            step, phase_assembly_transformed, zero_gamma, local_info)
        call require(local_info == 0, "analytic zero-gamma gate failed")
        slopes = jacobian_s / (2.0_dp * jacobian) + [-1.0_dp, 1.0_dp] / step
        eta_div = (1.2_dp * jacobian_t - 0.7_dp * jacobian_z) &
            / (jacobian * (1.2_dp**2 + 0.7_dp**2))
        expected = 0.0_dp
        expected(1, 1) = step * gamma_p * abs(jacobian) * slopes(1)**2
        expected(1, 3) = step * gamma_p * abs(jacobian) &
            * slopes(1) * slopes(2)
        expected(3, 1) = expected(1, 3)
        expected(3, 3) = step * gamma_p * abs(jacobian) * slopes(2)**2
        expected(1, 6) = step * gamma_p * abs(jacobian) &
            * slopes(1) * eta_div
        expected(6, 1) = expected(1, 6)
        expected(3, 6) = step * gamma_p * abs(jacobian) &
            * slopes(2) * eta_div
        expected(6, 3) = expected(3, 6)
        expected(6, 6) = step * gamma_p * abs(jacobian) * eta_div**2
        difference = finite - zero_gamma
        call require_close(difference, expected, 1.0e-13_dp, &
            "analytic signed-Jacobian product derivative is wrong")
    end subroutine check_analytic_product_derivatives

    subroutine check_analytic_phase_parity(local_radial_space)
        type(radial_space_config_t), intent(in) :: local_radial_space
        integer, parameter :: modes_m(2) = [1, 1], modes_n(2) = [1, 1]
        integer, parameter :: parities(2) = [1, 2]
        real(dp), parameter :: powers(2) = 0.0_dp
        real(dp), parameter :: ft = 1.2_dp, fp = 0.7_dp, jacobian = -1.1_dp
        real(dp), parameter :: gamma_p = 2.0_dp, step = 0.25_dp
        real(dp) :: constant_fields(4, 1, 13), zeros(4, 1), gamma(4, 1)
        real(dp) :: expected(2, 2), divergence_amplitude
        real(dp), allocatable :: stiffness(:, :)
        integer :: local_info

        constant_fields = 0.0_dp
        constant_fields(:, :, 1) = ft
        constant_fields(:, :, 2) = fp
        constant_fields(:, :, 7) = jacobian
        constant_fields(:, :, 8) = 1.4_dp
        constant_fields(:, :, 9) = 1.3_dp
        zeros = 0.0_dp
        gamma = gamma_p
        call assemble_compressible_stiffness_surface(constant_fields, zeros, &
            zeros, zeros, zeros, gamma, modes_m, modes_n, parities, powers, &
            3, local_radial_space, 0.375_dp, step, &
            phase_assembly_transformed, stiffness, local_info)
        call require(local_info == 0, "analytic phase gate failed to assemble")
        divergence_amplitude = 2.0_dp * acos(-1.0_dp) &
            * (fp - ft / 3.0_dp) / (jacobian * (ft**2 + fp**2))
        expected = 0.0_dp
        expected(1, 1) = 0.5_dp * step * gamma_p * abs(jacobian) &
            * divergence_amplitude**2
        expected(2, 2) = expected(1, 1)
        call require_close(stiffness(7:8, 7:8), expected, 1.0e-13_dp, &
            "analytic cosine/sine phase parity is wrong")
    end subroutine check_analytic_phase_parity

    subroutine calculate_direct_energy(local_fields, local_drive, &
            local_jacobian_radial, local_jacobian_theta, local_jacobian_zeta, &
            local_gamma_pressure, local_radial_space, coefficients, energy, &
            info)
        real(dp), intent(in) :: local_fields(:, :, :), local_drive(:, :)
        real(dp), intent(in) :: local_jacobian_radial(:, :)
        real(dp), intent(in) :: local_jacobian_theta(:, :)
        real(dp), intent(in) :: local_jacobian_zeta(:, :)
        real(dp), intent(in) :: local_gamma_pressure(:, :), coefficients(:)
        type(radial_space_config_t), intent(in) :: local_radial_space
        real(dp), intent(out) :: energy
        integer, intent(out) :: info
        real(dp) :: density, theta, zeta, values(9), weight
        integer :: j, k, period

        energy = 0.0_dp
        weight = 0.25_dp / real(n_theta * n_zeta * 3, dp)
        do period = 0, 2
            do k = 1, n_zeta
                zeta = real(k - 1, dp) / real(n_zeta, dp) + real(period, dp)
                do j = 1, n_theta
                    theta = real(j - 1, dp) / real(n_theta, dp)
                    call evaluate_displacement(theta, zeta, &
                        local_radial_space, coefficients, values, info)
                    if (info /= 0) return
                    density = point_density(local_fields(j, k, :), &
                        local_drive(j, k), local_jacobian_radial(j, k), &
                        local_jacobian_theta(j, k), &
                        local_jacobian_zeta(j, k), &
                        local_gamma_pressure(j, k), values)
                    energy = energy + weight * density
                end do
            end do
        end do
        info = 0
    end subroutine calculate_direct_energy

    subroutine evaluate_displacement(theta, zeta, local_radial_space, &
            coefficients, values, info)
        real(dp), intent(in) :: theta, zeta, coefficients(:)
        type(radial_space_config_t), intent(in) :: local_radial_space
        real(dp), intent(out) :: values(9)
        integer, intent(out) :: info
        real(dp) :: normal_values(2), normal_derivatives(2)
        real(dp) :: phase_values(6), normal, normal_radial
        integer :: trial

        values = 0.0_dp
        do trial = 1, size(trial_m)
            call evaluate_normal_basis(local_radial_space, trial_m(trial), &
                0.375_dp, 0.25_dp, 0.5_dp, normal_values, &
                normal_derivatives, info, stored_power(trial))
            if (info /= radial_space_ok) return
            call evaluate_phases(theta, zeta, trial_m(trial), trial_n(trial), &
                trial_parity(trial), phase_values)
            normal = coefficients(trial) * normal_values(1) &
                + coefficients(3 + trial) * normal_values(2)
            normal_radial = coefficients(trial) * normal_derivatives(1) &
                + coefficients(3 + trial) * normal_derivatives(2)
            values(1) = values(1) + normal * phase_values(1)
            values(2) = values(2) + normal_radial * phase_values(1)
            values(3) = values(3) + normal * phase_values(2)
            values(4) = values(4) + normal * phase_values(3)
            values(5) = values(5) + coefficients(6 + trial) * phase_values(4)
            values(6) = values(6) + coefficients(6 + trial) * phase_values(5)
            values(7) = values(7) + coefficients(6 + trial) * phase_values(6)
            values(8) = values(8) + coefficients(9 + trial) * phase_values(5)
            values(9) = values(9) + coefficients(9 + trial) * phase_values(6)
        end do
        info = 0
    end subroutine evaluate_displacement

    pure subroutine evaluate_phases(theta, zeta, mode_m, mode_n, parity, &
            values)
        real(dp), intent(in) :: theta, zeta
        integer, intent(in) :: mode_m, mode_n, parity
        real(dp), intent(out) :: values(6)
        real(dp) :: phase, xi, tangential, xi_derivative, tangential_derivative

        phase = 2.0_dp * acos(-1.0_dp) * (real(mode_m, dp) * theta &
            - real(mode_n, dp) * zeta / 3.0_dp)
        if (parity == 1) then
            xi = cos(phase)
            tangential = sin(phase)
            xi_derivative = -tangential
            tangential_derivative = xi
        else
            xi = sin(phase)
            tangential = cos(phase)
            xi_derivative = tangential
            tangential_derivative = -xi
        end if
        values(1) = xi
        values(2) = xi_derivative * 2.0_dp * acos(-1.0_dp) * mode_m
        values(3) = -xi_derivative * 2.0_dp * acos(-1.0_dp) * mode_n / 3.0_dp
        values(4) = tangential
        values(5) = tangential_derivative * 2.0_dp * acos(-1.0_dp) * mode_m
        values(6) = -tangential_derivative * 2.0_dp * acos(-1.0_dp) &
            * mode_n / 3.0_dp
    end subroutine evaluate_phases

    pure function point_density(local_fields, local_drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure, values) &
            result(density)
        real(dp), intent(in) :: local_fields(:), local_drive
        real(dp), intent(in) :: jacobian_radial, jacobian_theta
        real(dp), intent(in) :: jacobian_zeta, gamma_pressure, values(9)
        real(dp) :: density

        density = three_component_density_value(local_fields(1), &
            local_fields(2), local_fields(3), local_fields(4), &
            local_fields(5), local_fields(6), local_fields(7), &
            local_fields(8), local_fields(9), local_fields(10), &
            local_fields(11), local_fields(12), local_fields(13), &
            local_drive, gamma_pressure, values(1), values(2), values(3), &
            values(4), values(6), values(7), &
            jacobian_radial * values(1) + local_fields(7) * values(2), &
            jacobian_theta * values(5) + local_fields(7) * values(6), &
            jacobian_zeta * values(5) + local_fields(7) * values(7), &
            values(8), values(9))
    end function point_density

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

    subroutine require_close(first, second, tolerance, message)
        real(dp), intent(in) :: first(:, :), second(:, :), tolerance
        character(len=*), intent(in) :: message
        real(dp) :: difference, scale

        difference = maxval(abs(first - second))
        scale = max(1.0_dp, maxval(abs(first)), maxval(abs(second)))
        if (difference > tolerance * scale) then
            write (error_unit, "(a,2es24.16)") "FAIL: " // message // &
                "; difference and scale: ", difference, scale
            error stop 1
        end if
    end subroutine require_close

end program test_compressible_stiffness_assembly
