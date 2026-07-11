program test_variable_spectrum_analysis
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_tridiagonal_t
    use variable_spectrum_analysis, only: analyze_variable_spectrum, &
        variable_spectrum_invalid, variable_spectrum_ok, &
        variable_spectrum_summary_t
    implicit none

    real(dp), parameter :: floor = 1.0e-8_dp
    real(dp) :: stiffness(7, 7), mass(7, 7)
    type(variable_block_tridiagonal_t) :: stiffness_blocks, mass_blocks
    type(variable_spectrum_summary_t) :: summary
    integer :: info

    call build_fixture(stiffness, mass)
    call pack_variable_blocks(stiffness, [2, 3, 2], stiffness_blocks, info)
    call require(info == 0, "spectrum stiffness packing failed")
    call pack_variable_blocks(mass, [2, 3, 2], mass_blocks, info)
    call require(info == 0, "spectrum mass packing failed")

    call analyze_variable_spectrum(stiffness_blocks, mass_blocks, floor, &
        summary, info)
    call require(info == variable_spectrum_ok, "spectrum analysis failed")
    call require(summary%negative_count == 1, "negative count is wrong")
    call require(summary%zero_count == 3, "zero-cluster count is wrong")
    call require(summary%has_positive, "positive spectrum was not found")
    call require(summary%first_positive_lower <= 3.0_dp .and. &
        summary%first_positive_upper > 3.0_dp, &
        "first positive eigenvalue is outside the bracket")
    call require(summary%first_positive_upper &
        - summary%first_positive_lower <= floor, &
        "first positive bracket exceeds the configured resolution")
    call require(summary%first_positive_cluster_count == 2, &
        "first positive resolution-cluster count is wrong")

    call analyze_variable_spectrum(stiffness_blocks, mass_blocks, 2.0_dp, &
        summary, info)
    call require(info == variable_spectrum_ok, &
        "singular midpoint spectrum analysis failed")
    call require(summary%negative_count == 0 .and. summary%zero_count == 4, &
        "floor-boundary spectrum counts are wrong")
    call require(summary%first_positive_lower < 3.0_dp .and. &
        summary%first_positive_upper > 3.0_dp, &
        "singular midpoint was not straddled")
    call require(summary%first_positive_cluster_count == 2, &
        "singular midpoint cluster count is wrong")

    call build_higher_singular_fixture(stiffness, mass)
    call pack_variable_blocks(stiffness, [2, 3, 2], stiffness_blocks, info)
    call pack_variable_blocks(mass, [2, 3, 2], mass_blocks, info)
    call analyze_variable_spectrum(stiffness_blocks, mass_blocks, 0.1_dp, &
        summary, info)
    call require(info == variable_spectrum_ok, &
        "higher singular probe spectrum analysis failed")
    call require(summary%first_positive_lower <= 0.82_dp .and. &
        summary%first_positive_upper > 0.82_dp, &
        "higher singular probe hid the first positive eigenvalue")

    call build_fixture(stiffness, mass)
    stiffness(5, 5) = 0.0_dp
    stiffness(6, 6) = 0.0_dp
    stiffness(7, 7) = 0.0_dp
    call pack_variable_blocks(stiffness, [2, 3, 2], stiffness_blocks, info)
    call pack_variable_blocks(mass, [2, 3, 2], mass_blocks, info)
    call analyze_variable_spectrum(stiffness_blocks, mass_blocks, floor, &
        summary, info)
    call require(info == variable_spectrum_ok, &
        "nonpositive spectrum analysis failed")
    call require(.not. summary%has_positive, &
        "nonpositive spectrum reports a positive gap")
    call require(summary%negative_count == 1 .and. summary%zero_count == 6, &
        "nonpositive spectrum counts are wrong")

    call analyze_variable_spectrum(stiffness_blocks, mass_blocks, 0.0_dp, &
        summary, info)
    call require(info == variable_spectrum_invalid, "zero floor was accepted")
    call analyze_variable_spectrum(stiffness_blocks, mass_blocks, &
        ieee_value(0.0_dp, ieee_quiet_nan), summary, info)
    call require(info == variable_spectrum_invalid, &
        "nonfinite floor was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine build_fixture(local_stiffness, local_mass)
        real(dp), intent(out) :: local_stiffness(:, :), local_mass(:, :)
        real(dp), parameter :: values(7) = &
            [-2.0_dp, -2.0e-10_dp, 0.0_dp, 4.0e-10_dp, &
            3.0_dp, 3.0_dp, 8.0_dp]
        integer :: i

        local_stiffness = 0.0_dp
        local_mass = 0.0_dp
        do i = 1, size(values)
            local_mass(i, i) = 1.0_dp + 0.25_dp * real(i, dp)
            local_stiffness(i, i) = values(i) * local_mass(i, i)
        end do
    end subroutine build_fixture

    subroutine build_higher_singular_fixture(local_stiffness, local_mass)
        real(dp), intent(out) :: local_stiffness(:, :), local_mass(:, :)
        real(dp), parameter :: values(7) = &
            [0.0_dp, 0.82_dp, 0.85_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        integer :: i

        local_stiffness = 0.0_dp
        local_mass = 0.0_dp
        do i = 1, size(values)
            local_mass(i, i) = 1.0_dp
            local_stiffness(i, i) = values(i)
        end do
    end subroutine build_higher_singular_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_variable_spectrum_analysis
