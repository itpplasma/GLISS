module enzyme_terpsichore_normalization_fixture
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use terpsichore_normalization, only: &
        map_gliss_export_flip_pol_cell_to_terpsichore
    use terpsichore_reduced_mass, only: terpsichore_reduced_element_energy
    implicit none
    private

    public :: normalized_terpsichore_reduced_energy

contains

    function normalized_terpsichore_reduced_energy(active) result(energy) &
            bind(c, name="normalized_terpsichore_reduced_energy")
        real(c_double), value :: active
        real(c_double) :: energy
        real(dp) :: signed_jacobian(4), signed_bjac(4), normal_phase(2, 4)
        real(dp) :: tangential_phase(2, 4), radial_factor(2), displacement(6)
        real(dp) :: ftp, fpp, radial_weight
        integer :: info

        signed_jacobian = [-28.0_dp - 0.1_dp * active, &
            -30.0_dp + 0.2_dp * active, -29.0_dp - 0.3_dp * active, &
            -31.0_dp + 0.4_dp * active]
        normal_phase = reshape([0.2_dp, -0.4_dp, 0.5_dp, 0.1_dp, &
            -0.3_dp, 0.6_dp, 0.4_dp, -0.2_dp], [2, 4]) + 0.01_dp * active
        tangential_phase = reshape([0.7_dp, 0.1_dp, -0.2_dp, 0.8_dp, &
            0.3_dp, -0.5_dp, 0.6_dp, 0.2_dp], [2, 4]) - 0.02_dp * active
        radial_factor = [0.9_dp + 0.01_dp * active, &
            1.1_dp - 0.02_dp * active]
        displacement = [0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, -0.2_dp, 0.5_dp] &
            + 0.03_dp * active
        call map_gliss_export_flip_pol_cell_to_terpsichore(3, 4, &
            0.25_dp - 0.01_dp * active, 0.5_dp + 0.01_dp * active, &
            signed_jacobian, -0.5_dp - 0.02_dp * active, &
            -0.335_dp + 0.01_dp * active, signed_bjac, ftp, fpp, &
            radial_weight, info)
        energy = 0.0_c_double
        if (info /= 0) return
        energy = terpsichore_reduced_element_energy(signed_bjac, ftp, &
            normal_phase, tangential_phase, radial_factor, radial_weight, &
            displacement) + 0.01_dp * fpp**2
    end function normalized_terpsichore_reduced_energy

end module enzyme_terpsichore_normalization_fixture
