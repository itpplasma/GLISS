module enzyme_phase_fixture
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use family_point_assembly, only: assemble_transformed_surface_resolved
    use radial_space_policy, only: radial_space_config_t
    implicit none
    private

    public :: transformed_phase_energy

contains

    function transformed_phase_energy(scale) result(energy) &
            bind(c, name="transformed_phase_energy")
        real(c_double), value :: scale
        real(c_double) :: energy
        integer, parameter :: trials = 4
        integer, parameter :: mode_m(trials) = [1, 1, 2, 2]
        integer, parameter :: mode_n(trials) = [1, 1, -1, -1]
        integer, parameter :: parity(trials) = [1, 2, 1, 2]
        real(dp), parameter :: stored_power(trials) = &
            [0.25_dp, 0.25_dp, 0.0_dp, 0.0_dp]
        type(radial_space_config_t) :: radial_space
        real(dp) :: fields(2, 2, 13), drive(2, 2)
        real(dp) :: full(3 * trials, 3 * trials), direction, weight
        integer :: i, j, k, row, column, info

        do k = 1, 13
            do j = 1, 2
                do i = 1, 2
                    fields(i, j, k) = 0.2_dp + 0.03_dp * real(k, dp) &
                        + 0.01_dp * real(i + 2 * j, dp)
                    direction = 0.002_dp * real(k * (i + j), dp)
                    fields(i, j, k) = fields(i, j, k) + scale * direction
                end do
            end do
        end do
        fields(:, :, 7) = fields(:, :, 7) + 1.5_dp
        fields(:, :, 8) = fields(:, :, 8) + 1.0_dp
        fields(:, :, 9) = fields(:, :, 9) + 0.8_dp
        do j = 1, 2
            do i = 1, 2
                drive(i, j) = 0.04_dp * real(i - j, dp) &
                    + scale * 0.003_dp * real(i + j, dp)
            end do
        end do
        full = 0.0_dp
        call assemble_transformed_surface_resolved(fields, drive, mode_m, &
            mode_n, parity, stored_power, 3, radial_space, 0.35_dp, 0.1_dp, &
            full, info)
        energy = 0.0_c_double
        if (info /= 0) return
        do column = 1, size(full, 2)
            do row = 1, size(full, 1)
                weight = real(row + 2 * column, dp) &
                    / real(3 * size(full, 1), dp)
                energy = energy + weight * full(row, column)
            end do
        end do
    end function transformed_phase_energy

end module enzyme_phase_fixture
