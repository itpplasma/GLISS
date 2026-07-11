module enzyme_compressible_stiffness_family_fixture
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compressible_stiffness_family_assembly, only: &
        assemble_compressible_family_stiffness_resolved
    use phase_assembly_policy, only: phase_assembly_transformed
    use radial_space_policy, only: radial_space_config_t
    implicit none
    private

    public :: global_compressible_stiffness_energy

contains

    function global_compressible_stiffness_energy(active) result(energy) &
            bind(c, name="global_compressible_stiffness_energy")
        real(c_double), value :: active
        real(c_double) :: energy
        integer, parameter :: trial_m(2) = [1, 2]
        integer, parameter :: trial_n(2) = [1, -2]
        integer, parameter :: parity(2) = [1, 2]
        type(radial_space_config_t) :: radial_space
        real(dp) :: fields(2, 2, 13, 3), drive(2, 2, 3)
        real(dp) :: jacobian_radial(2, 2, 3), jacobian_theta(2, 2, 3)
        real(dp) :: jacobian_zeta(2, 2, 3), gamma_pressure(2, 2, 3)
        real(dp) :: stiffness(16, 16), stored_power(2)
        integer :: info

        call build_active_fields(active, fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure)
        stored_power = [0.01_dp, -0.015_dp] * active
        call assemble_compressible_family_stiffness_resolved(fields, drive, &
            jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure, &
            trial_m, trial_n, parity, stored_power, 3, radial_space, &
            1.0_dp / 3.0_dp, phase_assembly_transformed, stiffness, info)
        energy = 0.0_c_double
        if (info /= 0) return
        energy = weighted_stiffness(stiffness)
    end function global_compressible_stiffness_energy

    pure subroutine build_active_fields(active, fields, drive, &
            jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure)
        real(c_double), intent(in) :: active
        real(dp), intent(out) :: fields(2, 2, 13, 3), drive(2, 2, 3)
        real(dp), intent(out) :: jacobian_radial(2, 2, 3)
        real(dp), intent(out) :: jacobian_theta(2, 2, 3)
        real(dp), intent(out) :: jacobian_zeta(2, 2, 3)
        real(dp), intent(out) :: gamma_pressure(2, 2, 3)
        real(dp) :: point_fields(13)
        integer :: i, j, interval

        fields = 0.0_dp
        do interval = 1, 3
            do j = 1, 2
                do i = 1, 2
                    call build_active_point(active, i, j, interval, &
                        point_fields, drive(i, j, interval), &
                        jacobian_radial(i, j, interval), &
                        jacobian_theta(i, j, interval), &
                        jacobian_zeta(i, j, interval), &
                        gamma_pressure(i, j, interval))
                    fields(i, j, :, interval) = point_fields
                end do
            end do
        end do
    end subroutine build_active_fields

    pure subroutine build_active_point(active, i, j, interval, fields, &
            drive, jacobian_radial, jacobian_theta, jacobian_zeta, &
            gamma_pressure)
        real(c_double), intent(in) :: active
        integer, intent(in) :: i, j, interval
        real(dp), intent(out) :: fields(13), drive, jacobian_radial
        real(dp), intent(out) :: jacobian_theta, jacobian_zeta, gamma_pressure

        fields = 0.0_dp
        fields(1) = 1.2_dp + 0.01_dp * interval + active * 0.002_dp * (i + j)
        fields(2) = 0.7_dp + 0.01_dp * j - active * 0.001_dp * (i + interval)
        fields(3) = 0.04_dp + active * 0.001_dp * i
        fields(4) = -0.03_dp + active * 0.001_dp * j
        fields(5) = 0.8_dp + active * 0.002_dp * interval
        fields(6) = 0.6_dp - active * 0.001_dp * j
        fields(7) = -1.1_dp + active * 0.001_dp * (i + j + interval)
        fields(8) = 1.4_dp + active * 0.002_dp * j
        fields(9) = 1.3_dp + active * 0.001_dp * i
        fields(10) = 0.2_dp - active * 0.001_dp * interval
        fields(11) = -0.15_dp + active * 0.001_dp * i
        fields(12) = 0.1_dp - active * 0.001_dp * j
        fields(13) = -0.12_dp + active * 0.002_dp * interval
        drive = 0.05_dp + active * 0.001_dp * (i - j + interval)
        jacobian_radial = -0.08_dp + active * 0.002_dp * i
        jacobian_theta = 0.03_dp - active * 0.001_dp * j
        jacobian_zeta = -0.02_dp + active * 0.001_dp * (i + interval)
        gamma_pressure = 0.9_dp + active * 0.003_dp * (i + j + interval)
    end subroutine build_active_point

    pure function weighted_stiffness(stiffness) result(energy)
        real(dp), intent(in) :: stiffness(:, :)
        real(dp) :: energy, weight
        integer :: row, column

        energy = 0.0_dp
        do column = 1, size(stiffness, 2)
            do row = 1, size(stiffness, 1)
                weight = real(row + 2 * column, dp) &
                    / real(3 * size(stiffness, 1), dp)
                energy = energy + weight * stiffness(row, column)
            end do
        end do
    end function weighted_stiffness

end module enzyme_compressible_stiffness_family_fixture
