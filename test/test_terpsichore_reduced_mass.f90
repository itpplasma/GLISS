program test_terpsichore_reduced_mass
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fourier_phase_kind, only: phase_cosine, phase_sine
    use symmetric_eigensolver, only: solve_symmetric_generalized
    use terpsichore_reduced_mass, only: &
        assemble_terpsichore_reduced_mass_element, &
        terpsichore_reduced_element_energy, terpsichore_reduced_invalid, &
        terpsichore_reduced_ok
    implicit none

    integer, parameter :: modes = 3, points = 64
    real(dp), parameter :: flux_t_slope = 1.3_dp
    real(dp), parameter :: normalized_radial_weight = 0.75_dp
    real(dp), parameter :: normal_radial_factor(modes) = &
        [1.0_dp, 0.9_dp, 1.1_dp]
    real(dp) :: signed_bjac(points), reflected_bjac(points)
    real(dp) :: phase(modes, points)
    real(dp) :: normal_phase(modes, points), tangential_phase(modes, points)
    real(dp) :: displacement(3 * modes), energy, image(3 * modes)
    real(dp), allocatable :: mass(:, :), reflected(:, :)
    integer :: column, info, row

    call build_fixture(signed_bjac, phase)
    displacement = [0.2_dp, -0.3_dp, 0.4_dp, -0.1_dp, 0.5_dp, -0.2_dp, &
        0.3_dp, 0.1_dp, -0.4_dp]
    call check_parity(phase_sine, signed_bjac, phase)
    call check_parity(phase_cosine, signed_bjac, phase)

    normal_phase = sin(phase)
    tangential_phase = cos(phase)
    call assemble_terpsichore_reduced_mass_element(signed_bjac, &
        flux_t_slope, normal_phase, tangential_phase, normal_radial_factor, &
        normalized_radial_weight, mass, info)
    call require(info == terpsichore_reduced_ok, &
        "valid reduced mass fixture failed")
    call require(maxval(abs(mass - transpose(mass))) < 2.0e-14_dp, &
        "reduced mass is not symmetric")
    call require(maxval(abs(mass(1:2 * modes, 2 * modes + 1:))) == 0.0_dp, &
        "reduced mass has a C11 coupling")
    do row = 1, points
        reflected_bjac(row) = -signed_bjac(row)
    end do
    call assemble_terpsichore_reduced_mass_element(reflected_bjac, &
        flux_t_slope, normal_phase, tangential_phase, normal_radial_factor, &
        normalized_radial_weight, reflected, info)
    call require(info == terpsichore_reduced_ok, &
        "opposite-orientation reduced mass failed")
    call require(maxval(abs(reflected - mass)) == 0.0_dp, &
        "reduced mass depends on BJAC orientation")
    energy = terpsichore_reduced_element_energy(signed_bjac, flux_t_slope, &
        normal_phase, tangential_phase, normal_radial_factor, &
        normalized_radial_weight, displacement)
    do row = 1, size(image)
        image(row) = 0.0_dp
        do column = 1, size(displacement)
            image(row) = image(row) + mass(row, column) * displacement(column)
        end do
    end do
    call require(abs(energy - 0.5_dp * dot_product(displacement, image)) &
        < 2.0e-13_dp * max(1.0_dp, energy), &
        "reduced mass and energy disagree")
    call assemble_terpsichore_reduced_mass_element(signed_bjac, 0.0_dp, &
        normal_phase, tangential_phase, normal_radial_factor, &
        normalized_radial_weight, reflected, info)
    call require(info == terpsichore_reduced_invalid, &
        "zero toroidal flux slope was accepted")
    call assemble_terpsichore_reduced_mass_element(signed_bjac, &
        flux_t_slope, normal_phase, tangential_phase, &
        normal_radial_factor(:2), normalized_radial_weight, reflected, info)
    call require(info == terpsichore_reduced_invalid, &
        "wrong radial-factor shape was accepted")
    call require(ieee_is_nan(terpsichore_reduced_element_energy(signed_bjac, &
        flux_t_slope, normal_phase, tangential_phase, normal_radial_factor, &
        -normalized_radial_weight, displacement)), &
        "negative radial weight produced a finite energy")

    write (*, "(a)") "PASS"

contains

    subroutine build_fixture(local_bjac, local_phase)
        real(dp), intent(out) :: local_bjac(:), local_phase(:, :)
        real(dp) :: theta, zeta
        integer :: point

        do point = 1, size(local_bjac)
            theta = 2.0_dp * acos(-1.0_dp) * real(point - 1, dp) &
                / real(size(local_bjac), dp)
            zeta = 2.0_dp * acos(-1.0_dp) * real(modulo(5 * point, &
                size(local_bjac)), dp) / real(size(local_bjac), dp)
            local_bjac(point) = -(1.2_dp + 0.15_dp * cos(theta) &
                - 0.08_dp * cos(2.0_dp * theta - zeta))
            local_phase(1, point) = theta - zeta
            local_phase(2, point) = 2.0_dp * theta + zeta
            local_phase(3, point) = 3.0_dp * theta - 2.0_dp * zeta
        end do
    end subroutine build_fixture

    subroutine check_parity(normal_kind, local_bjac, local_phase)
        integer, intent(in) :: normal_kind
        real(dp), intent(in) :: local_bjac(:), local_phase(:, :)
        real(dp) :: expected(3 * modes, 3 * modes)
        real(dp), allocatable :: actual(:, :)
        integer :: info

        if (normal_kind == phase_sine) then
            normal_phase = sin(local_phase)
            tangential_phase = cos(local_phase)
        else
            normal_phase = cos(local_phase)
            tangential_phase = sin(local_phase)
        end if
        call assemble_terpsichore_reduced_mass_element(local_bjac, &
            flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, actual, info)
        call require(info == terpsichore_reduced_ok, &
            "reduced parity assembly failed")
        call source_formula_oracle(normal_kind, local_bjac, local_phase, &
            expected)
        call require(maxval(abs(actual - expected)) < 2.0e-13_dp &
            * max(1.0_dp, maxval(abs(expected))), &
            "reduced mass disagrees with the C8/C10 source formula")
        call check_psd_and_nullity(actual)
    end subroutine check_parity

    subroutine source_formula_oracle(normal_kind, local_bjac, local_phase, &
            expected)
        integer, intent(in) :: normal_kind
        real(dp), intent(in) :: local_bjac(:), local_phase(:, :)
        real(dp), intent(out) :: expected(:, :)
        real(dp) :: difference, summation, c8, c10, sign
        integer :: first, second

        expected = 0.0_dp
        sign = -1.0_dp
        if (normal_kind == phase_cosine) sign = 1.0_dp
        do second = 1, modes
            do first = 1, modes
                difference = weighted_cosine_pair(local_bjac, local_phase, &
                    first, second, -1.0_dp)
                summation = weighted_cosine_pair(local_bjac, local_phase, &
                    first, second, 1.0_dp)
                c8 = normalized_radial_weight &
                    * (difference + sign * summation)
                c10 = normalized_radial_weight &
                    * (difference - sign * summation) &
                    / (4.0_dp * flux_t_slope**2)
                call set_normal_blocks(expected, first, second, &
                    0.25_dp * c8 * normal_radial_factor(first) &
                    * normal_radial_factor(second))
                expected(2 * modes + first, 2 * modes + second) = c10
            end do
        end do
    end subroutine source_formula_oracle

    function weighted_cosine_pair(local_bjac, local_phase, first, second, &
            second_sign) result(average)
        real(dp), intent(in) :: local_bjac(:), local_phase(:, :)
        integer, intent(in) :: first, second
        real(dp), intent(in) :: second_sign
        real(dp) :: average
        integer :: point

        average = 0.0_dp
        do point = 1, size(local_bjac)
            average = average + abs(local_bjac(point)) &
                * cos(local_phase(first, point) &
                + second_sign * local_phase(second, point))
        end do
        average = average / real(size(local_bjac), dp)
    end function weighted_cosine_pair

    subroutine set_normal_blocks(matrix, first, second, value)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(in) :: first, second
        real(dp), intent(in) :: value

        matrix(first, second) = value
        matrix(first, modes + second) = value
        matrix(modes + first, second) = value
        matrix(modes + first, modes + second) = value
    end subroutine set_normal_blocks

    subroutine check_psd_and_nullity(matrix)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: identity(size(matrix, 1), size(matrix, 1)), scale, floor
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        integer :: i, info

        identity = 0.0_dp
        do i = 1, size(identity, 1)
            identity(i, i) = 1.0_dp
        end do
        call solve_symmetric_generalized(matrix, identity, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "reduced mass dense eigensolve failed")
        scale = max(1.0_dp, maxval(abs(eigenvalues)))
        floor = 512.0_dp * epsilon(1.0_dp) * scale
        call require(minval(eigenvalues) >= -floor, &
            "reduced element mass is not positive semidefinite")
        call require(count(abs(eigenvalues) <= floor) == modes, &
            "reduced element mass has the wrong endpoint nullity")
    end subroutine check_psd_and_nullity

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_reduced_mass
