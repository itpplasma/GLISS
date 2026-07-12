module terpsichore_reduced_mass
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: terpsichore_reduced_ok = 0
    integer, parameter, public :: terpsichore_reduced_invalid = -1

    public :: add_reduced_values
    public :: assemble_terpsichore_reduced_mass_element
    public :: assemble_terpsichore_reduced_mass_element_resolved
    public :: benchmark_terpsichore_reduced_mass_energy
    public :: terpsichore_reduced_element_energy

contains

    pure subroutine assemble_terpsichore_reduced_mass_element(signed_bjac, &
            flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, mass, info)
        ! TERPSICHORE inputs are BJAC, s**(-QL), and NI*Delta(s); output order
        ! is the GLISS element order [xi_left, xi_right, eta].
        real(dp), intent(in) :: signed_bjac(:), flux_t_slope
        real(dp), intent(in) :: normal_phase(:, :), tangential_phase(:, :)
        real(dp), intent(in) :: normal_radial_factor(:)
        real(dp), intent(in) :: normalized_radial_weight
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: info
        integer :: mode_count

        info = terpsichore_reduced_invalid
        if (.not. valid_reduced_inputs(signed_bjac, flux_t_slope, &
            normal_phase, tangential_phase, normal_radial_factor, &
            normalized_radial_weight)) return
        mode_count = size(normal_phase, 1)
        allocate (mass(3 * mode_count, 3 * mode_count))
        call assemble_terpsichore_reduced_mass_element_resolved(signed_bjac, &
            flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, mass, info)
    end subroutine assemble_terpsichore_reduced_mass_element

    pure subroutine assemble_terpsichore_reduced_mass_element_resolved( &
            signed_bjac, flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, mass, info)
        real(dp), intent(in) :: signed_bjac(:), flux_t_slope
        real(dp), intent(in) :: normal_phase(:, :), tangential_phase(:, :)
        real(dp), intent(in) :: normal_radial_factor(:)
        real(dp), intent(in) :: normalized_radial_weight
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: normal_value, tangential_value, weight
        integer :: first, mode_count, point, second

        info = terpsichore_reduced_invalid
        if (.not. valid_reduced_inputs(signed_bjac, flux_t_slope, &
            normal_phase, tangential_phase, normal_radial_factor, &
            normalized_radial_weight)) return
        mode_count = size(normal_phase, 1)
        if (any(shape(mass) /= 3 * mode_count)) return
        mass = 0.0_dp
        do point = 1, size(signed_bjac)
            weight = normalized_radial_weight * abs(signed_bjac(point)) &
                / real(size(signed_bjac), dp)
            do second = 1, mode_count
                do first = 1, mode_count
                    normal_value = 0.5_dp * weight &
                        * normal_phase(first, point) &
                        * normal_phase(second, point) &
                        * normal_radial_factor(first) &
                        * normal_radial_factor(second)
                    tangential_value = weight * tangential_phase(first, point) &
                        * tangential_phase(second, point) &
                        / (2.0_dp * flux_t_slope**2)
                    call add_reduced_values(mass, mode_count, first, second, &
                        normal_value, tangential_value)
                end do
            end do
        end do
        info = terpsichore_reduced_ok
    end subroutine assemble_terpsichore_reduced_mass_element_resolved

    pure function terpsichore_reduced_element_energy(signed_bjac, &
            flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, displacement) &
            result(energy)
        real(dp), intent(in) :: signed_bjac(:), flux_t_slope
        real(dp), intent(in) :: normal_phase(:, :), tangential_phase(:, :)
        real(dp), intent(in) :: normal_radial_factor(:)
        real(dp), intent(in) :: normalized_radial_weight, displacement(:)
        real(dp) :: energy, eta_value, normal_value, weight
        integer :: mode, mode_count, point

        energy = ieee_value(energy, ieee_quiet_nan)
        if (.not. valid_reduced_inputs(signed_bjac, flux_t_slope, &
            normal_phase, tangential_phase, normal_radial_factor, &
            normalized_radial_weight)) return
        mode_count = size(normal_phase, 1)
        if (size(displacement) /= 3 * mode_count) return
        if (.not. all(ieee_is_finite(displacement))) return
        energy = 0.0_dp
        do point = 1, size(signed_bjac)
            normal_value = 0.0_dp
            eta_value = 0.0_dp
            do mode = 1, mode_count
                normal_value = normal_value + 0.5_dp &
                    * (displacement(mode) + displacement(mode_count + mode)) &
                    * normal_phase(mode, point) * normal_radial_factor(mode)
                eta_value = eta_value + displacement(2 * mode_count + mode) &
                    * tangential_phase(mode, point)
            end do
            weight = normalized_radial_weight * abs(signed_bjac(point)) &
                / real(size(signed_bjac), dp)
            energy = energy + 2.0_dp * weight &
                * (normal_value**2 + eta_value**2 &
                / (4.0_dp * flux_t_slope**2))
        end do
        energy = 0.5_dp * energy
    end function terpsichore_reduced_element_energy

    pure subroutine add_reduced_values(mass, mode_count, first, second, &
            normal_value, tangential_value)
        real(dp), intent(inout) :: mass(:, :)
        integer, intent(in) :: mode_count, first, second
        real(dp), intent(in) :: normal_value, tangential_value

        mass(first, second) = mass(first, second) + normal_value
        mass(first, mode_count + second) = &
            mass(first, mode_count + second) + normal_value
        mass(mode_count + first, second) = &
            mass(mode_count + first, second) + normal_value
        mass(mode_count + first, mode_count + second) = &
            mass(mode_count + first, mode_count + second) + normal_value
        mass(2 * mode_count + first, 2 * mode_count + second) = &
            mass(2 * mode_count + first, 2 * mode_count + second) &
            + tangential_value
    end subroutine add_reduced_values

    pure function valid_reduced_inputs(signed_bjac, flux_t_slope, &
            normal_phase, tangential_phase, normal_radial_factor, &
            normalized_radial_weight) result(valid)
        real(dp), intent(in) :: signed_bjac(:), flux_t_slope
        real(dp), intent(in) :: normal_phase(:, :), tangential_phase(:, :)
        real(dp), intent(in) :: normal_radial_factor(:)
        real(dp), intent(in) :: normalized_radial_weight
        logical :: valid

        valid = size(signed_bjac) > 0
        if (.not. valid) return
        valid = size(normal_phase, 1) > 0
        if (.not. valid) return
        valid = size(normal_phase, 2) == size(signed_bjac)
        if (.not. valid) return
        valid = all(shape(tangential_phase) == shape(normal_phase))
        if (.not. valid) return
        valid = size(normal_radial_factor) == size(normal_phase, 1)
        if (.not. valid) return
        valid = ieee_is_finite(flux_t_slope) .and. flux_t_slope /= 0.0_dp
        if (.not. valid) return
        valid = ieee_is_finite(normalized_radial_weight) &
            .and. normalized_radial_weight > 0.0_dp
        if (.not. valid) return
        valid = all(ieee_is_finite(signed_bjac))
        if (.not. valid) return
        valid = all(ieee_is_finite(normal_radial_factor))
        if (.not. valid) return
        valid = all(ieee_is_finite(normal_phase)) &
            .and. all(ieee_is_finite(tangential_phase))
    end function valid_reduced_inputs

    pure function benchmark_terpsichore_reduced_mass_energy(active) &
            result(energy) bind(c, &
            name="gliss_benchmark_terpsichore_reduced_mass_energy")
        real(c_double), intent(in), value :: active
        real(c_double) :: energy
        real(dp) :: signed_bjac(4), normal_phase(2, 4)
        real(dp) :: normal_radial_factor(2)
        real(dp) :: tangential_phase(2, 4), displacement(6)

        signed_bjac = [-1.1_dp - 0.02_dp * active, &
            -1.2_dp + 0.01_dp * active, -0.9_dp - 0.03_dp * active, &
            -1.3_dp + 0.04_dp * active]
        normal_phase = reshape([0.2_dp, -0.4_dp, 0.5_dp, 0.1_dp, &
            -0.3_dp, 0.6_dp, 0.4_dp, -0.2_dp], [2, 4]) + 0.01_dp * active
        tangential_phase = reshape([0.7_dp, 0.1_dp, -0.2_dp, 0.8_dp, &
            0.3_dp, -0.5_dp, 0.6_dp, 0.2_dp], [2, 4]) - 0.02_dp * active
        displacement = [0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, -0.2_dp, 0.5_dp] &
            + 0.03_dp * active
        normal_radial_factor = [0.9_dp + 0.01_dp * active, &
            1.1_dp - 0.02_dp * active]
        energy = terpsichore_reduced_element_energy(signed_bjac, &
            1.3_dp + 0.02_dp * active, normal_phase, tangential_phase, &
            normal_radial_factor, 0.75_dp + 0.01_dp * active, displacement)
    end function benchmark_terpsichore_reduced_mass_energy

end module terpsichore_reduced_mass
