module enzyme_compressible_stiffness_fixture
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compressible_stiffness_assembly, only: &
        assemble_compressible_stiffness_surface_resolved
    use phase_assembly_policy, only: phase_assembly_transformed
    use radial_space_policy, only: radial_space_config_t
    implicit none
    private

    public :: transformed_compressible_stiffness_energy
    public :: transformed_radial_coordinate_energy
    public :: transformed_radial_step_energy
    public :: transformed_stored_power_energy

contains

    function transformed_compressible_stiffness_energy(active) result(energy) &
            bind(c, name="transformed_compressible_stiffness_energy")
        real(c_double), value :: active
        real(c_double) :: energy

        energy = active_stiffness_energy(active, 1)
    end function transformed_compressible_stiffness_energy

    function transformed_radial_coordinate_energy(active) result(energy) &
            bind(c, name="transformed_radial_coordinate_energy")
        real(c_double), value :: active
        real(c_double) :: energy

        energy = active_stiffness_energy(active, 2)
    end function transformed_radial_coordinate_energy

    function transformed_radial_step_energy(active) result(energy) &
            bind(c, name="transformed_radial_step_energy")
        real(c_double), value :: active
        real(c_double) :: energy

        energy = active_stiffness_energy(active, 3)
    end function transformed_radial_step_energy

    function transformed_stored_power_energy(active) result(energy) &
            bind(c, name="transformed_stored_power_energy")
        real(c_double), value :: active
        real(c_double) :: energy

        energy = active_stiffness_energy(active, 4)
    end function transformed_stored_power_energy

    function active_stiffness_energy(active, direction) result(energy)
        real(c_double), intent(in) :: active
        integer, intent(in) :: direction
        real(c_double) :: energy
        integer, parameter :: trial_m(2) = [1, 2]
        integer, parameter :: trial_n(2) = [1, -2]
        integer, parameter :: parity(2) = [1, 2]
        type(radial_space_config_t) :: radial_space
        real(dp) :: fields(2, 2, 13), drive(2, 2), jacobian_radial(2, 2)
        real(dp) :: jacobian_theta(2, 2), jacobian_zeta(2, 2)
        real(dp) :: gamma_pressure(2, 2), stiffness(8, 8), stored_power(2)
        real(dp) :: field_active, radial_coordinate, radial_step
        integer :: info

        field_active = 0.0_dp
        if (direction == 1) field_active = active
        radial_coordinate = 0.375_dp
        if (direction == 2) radial_coordinate = radial_coordinate + 0.01_dp * active
        radial_step = 0.25_dp
        if (direction == 3) radial_step = radial_step + 0.01_dp * active
        stored_power = 0.0_dp
        if (direction == 2) stored_power = [0.02_dp, -0.015_dp]
        if (direction == 4) stored_power = [0.02_dp, -0.015_dp] * active
        call build_active_fields(field_active, fields, drive, jacobian_radial, &
            jacobian_theta, jacobian_zeta, gamma_pressure)
        call assemble_compressible_stiffness_surface_resolved(fields, drive, &
            jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure, &
            trial_m, trial_n, parity, stored_power, 3, radial_space, &
            radial_coordinate, radial_step, phase_assembly_transformed, &
            stiffness, info)
        energy = 0.0_c_double
        if (info /= 0) return
        energy = weighted_stiffness(stiffness)
    end function active_stiffness_energy

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

    pure subroutine build_active_fields(active, fields, drive, &
            jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure)
        real(c_double), intent(in) :: active
        real(dp), intent(out) :: fields(2, 2, 13), drive(2, 2)
        real(dp), intent(out) :: jacobian_radial(2, 2)
        real(dp), intent(out) :: jacobian_theta(2, 2)
        real(dp), intent(out) :: jacobian_zeta(2, 2)
        real(dp), intent(out) :: gamma_pressure(2, 2)
        real(dp) :: point_fields(13)
        integer :: i, j

        fields = 0.0_dp
        do j = 1, 2
            do i = 1, 2
                call build_active_point(active, i, j, point_fields, &
                    drive(i, j), jacobian_radial(i, j), &
                    jacobian_theta(i, j), jacobian_zeta(i, j), &
                    gamma_pressure(i, j))
                fields(i, j, :) = point_fields
            end do
        end do
    end subroutine build_active_fields

    pure subroutine build_active_point(active, i, j, fields, drive, &
            jacobian_radial, jacobian_theta, jacobian_zeta, gamma_pressure)
        real(c_double), intent(in) :: active
        integer, intent(in) :: i, j
        real(dp), intent(out) :: fields(13), drive, jacobian_radial
        real(dp), intent(out) :: jacobian_theta, jacobian_zeta, gamma_pressure

        fields = 0.0_dp
        fields(1) = 1.2_dp + 0.01_dp * i + active * 0.002_dp * (i + j)
        fields(2) = 0.7_dp + 0.01_dp * j - active * 0.001_dp * (2 * i + j)
        fields(3) = 0.04_dp + active * 0.001_dp * i
        fields(4) = -0.03_dp + active * 0.001_dp * j
        fields(5) = 0.8_dp + active * 0.002_dp * i
        fields(6) = 0.6_dp - active * 0.001_dp * j
        fields(7) = -1.1_dp + active * 0.001_dp * (i + j)
        fields(8) = 1.4_dp + active * 0.002_dp * j
        fields(9) = 1.3_dp + active * 0.001_dp * i
        fields(10) = 0.2_dp - active * 0.001_dp * j
        fields(11) = -0.15_dp + active * 0.001_dp * i
        fields(12) = 0.1_dp - active * 0.001_dp * j
        fields(13) = -0.12_dp + active * 0.002_dp * i
        drive = 0.05_dp + active * 0.001_dp * (i - j)
        jacobian_radial = -0.08_dp + active * 0.002_dp * i
        jacobian_theta = 0.03_dp - active * 0.001_dp * j
        jacobian_zeta = -0.02_dp + active * 0.001_dp * (i + j)
        gamma_pressure = 0.9_dp + active * 0.003_dp * (2 * i + j)
    end subroutine build_active_point

end module enzyme_compressible_stiffness_fixture
