module terpsichore_reduced_mass_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_block_scatter, only: add_mapped_block_element, &
        allocate_dynamic_blocks, build_dynamic_block_map, &
        complete_dynamic_blocks, dynamic_block_map_t, dynamic_block_scatter_ok
    use dynamic_family_layout, only: add_mapped_dynamic_element, &
        dynamic_family_layout_t, dynamic_layout_ok
    use fourier_phase_kind, only: phase_sine
    use terpsichore_matrix_fixture, only: terpsichore_dense_order_is_valid, &
        terpsichore_fixed_fixture_is_valid, terpsichore_matrix_fixture_t
    use terpsichore_pair_average, only: terpsichore_pair_averages, &
        terpsichore_pair_ok
    use terpsichore_reduced_layout, only: &
        build_terpsichore_reduced_fixed_boundary_layout, &
        build_terpsichore_reduced_free_boundary_layout, &
        terpsichore_reduced_layout_ok
    use terpsichore_reduced_mass, only: add_reduced_values
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    implicit none
    private

    integer, parameter, public :: terpsichore_reduced_adapter_ok = 0
    integer, parameter, public :: terpsichore_reduced_adapter_invalid = -1

    public :: assemble_terpsichore_fixture_reduced_mass
    public :: assemble_terpsichore_fixture_reduced_mass_free_boundary
    public :: assemble_terpsichore_fixture_reduced_mass_blocks

contains

    ! Fixture-driven reduced mass through the pair-average transform:
    ! per interval one |BJAC| transform on the difference/sum mode
    ! table replaces the O(points modes^2) point-pair loop.  The
    ! phase-table family assembly remains the small-size oracle
    ! (test_terpsichore_reduced_mass_adapter compares both routes).
    subroutine assemble_terpsichore_fixture_reduced_mass(fixture, mass, &
            layout, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        call assemble_terpsichore_fixture_reduced_mass_with_boundary( &
            fixture, .false., layout, info, mass=mass)
    end subroutine assemble_terpsichore_fixture_reduced_mass

    ! The fixed-boundary reduced mass in the block-tridiagonal storage of
    ! map, assembled interval by interval without the dense matrix.
    subroutine assemble_terpsichore_fixture_reduced_mass_blocks(fixture, &
            blocks, layout, map, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        type(variable_block_tridiagonal_t), intent(out) :: blocks
        type(dynamic_family_layout_t), intent(out) :: layout
        type(dynamic_block_map_t), intent(out) :: map
        integer, intent(out) :: info

        call assemble_terpsichore_fixture_reduced_mass_with_boundary( &
            fixture, .false., layout, info, blocks=blocks, map=map)
    end subroutine assemble_terpsichore_fixture_reduced_mass_blocks

    subroutine assemble_terpsichore_fixture_reduced_mass_free_boundary( &
            fixture, mass, layout, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info

        call assemble_terpsichore_fixture_reduced_mass_with_boundary( &
            fixture, .true., layout, info, mass=mass)
    end subroutine assemble_terpsichore_fixture_reduced_mass_free_boundary

    ! Dense mass, or blocks and their map: exactly one of the two.
    subroutine assemble_terpsichore_fixture_reduced_mass_with_boundary( &
            fixture, retain_outer_normal, layout, info, mass, blocks, map)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        logical, intent(in) :: retain_outer_normal
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        real(dp), allocatable, optional, intent(out) :: mass(:, :)
        type(variable_block_tridiagonal_t), optional, intent(out) :: blocks
        type(dynamic_block_map_t), optional, intent(out) :: map
        real(dp), allocatable :: absolute_bjac(:)
        real(dp), allocatable :: radial_factor(:, :), radial_weight(:)
        real(dp), allocatable :: normal_normal(:, :), normal_tangent(:, :)
        real(dp), allocatable :: tangent_tangent(:, :), element(:, :)
        integer, allocatable :: element_to_global(:, :), parity(:)
        real(dp) :: normal_value, tangential_value
        integer :: allocation_status, interval, first, second, modes, point

        info = terpsichore_reduced_adapter_invalid
        if (.not. terpsichore_fixed_fixture_is_valid(fixture)) return
        if (fixture%parity /= 0.0_dp) return
        call build_radial_values(fixture, radial_factor, radial_weight, &
            allocation_status)
        if (allocation_status /= 0) return
        modes = fixture%modes
        allocate (parity(modes), source=phase_sine, stat=allocation_status)
        if (allocation_status /= 0) return
        if (retain_outer_normal) then
            call build_terpsichore_reduced_free_boundary_layout( &
                fixture%mode_m, fixture%mode_n, parity, fixture%intervals, &
                layout, element_to_global, info)
        else
            call build_terpsichore_reduced_fixed_boundary_layout( &
                fixture%mode_m, fixture%mode_n, parity, fixture%intervals, &
                layout, element_to_global, info)
        end if
        if (info /= terpsichore_reduced_layout_ok) then
            info = terpsichore_reduced_adapter_invalid
            return
        end if
        info = terpsichore_reduced_adapter_invalid
        if (present(mass)) then
            if (.not. terpsichore_dense_order_is_valid(layout%total_unknowns)) &
                return
            allocate (mass(layout%total_unknowns, layout%total_unknowns), &
                source=0.0_dp, stat=allocation_status)
            if (allocation_status /= 0) return
        else
            call build_dynamic_block_map(layout, map, info)
            if (info /= dynamic_block_scatter_ok) then
                info = terpsichore_reduced_adapter_invalid
                return
            end if
            call allocate_dynamic_blocks(map, blocks, info)
            if (info /= dynamic_block_scatter_ok) then
                info = terpsichore_reduced_adapter_invalid
                return
            end if
        end if
        allocate (absolute_bjac(size(fixture%signed_bjac, 1)), &
            normal_normal(modes, modes), normal_tangent(modes, modes), &
            tangent_tangent(modes, modes), element(3 * modes, 3 * modes), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = terpsichore_reduced_adapter_invalid
            return
        end if
        do interval = 1, fixture%intervals
            do point = 1, size(absolute_bjac)
                absolute_bjac(point) = abs(fixture%signed_bjac(point, interval))
            end do
            call terpsichore_pair_averages( &
                absolute_bjac, &
                fixture%poloidal_points, fixture%toroidal_points, &
                fixture%stability_periods, fixture%field_periods, &
                fixture%mode_m, fixture%mode_n, &
                normal_normal, normal_tangent, tangent_tangent, info)
            if (info /= terpsichore_pair_ok) then
                info = terpsichore_reduced_adapter_invalid
                return
            end if
            element = 0.0_dp
            do second = 1, modes
                do first = 1, modes
                    normal_value = 0.25_dp * radial_weight(interval) &
                        * radial_factor(first, interval) &
                        * radial_factor(second, interval) &
                        * normal_normal(first, second)
                    tangential_value = 0.25_dp * radial_weight(interval) &
                        * tangent_tangent(first, second) &
                        / fixture%flux_t_slope(interval)**2
                    call add_reduced_values(element, modes, first, second, &
                        normal_value, tangential_value)
                end do
            end do
            if (present(mass)) then
                call add_mapped_dynamic_element(element_to_global(:, &
                    interval), element, mass, info)
                if (info /= dynamic_layout_ok) then
                    info = terpsichore_reduced_adapter_invalid
                    return
                end if
            else
                call add_mapped_block_element(element_to_global(:, &
                    interval), element, map, blocks, info)
                if (info /= dynamic_block_scatter_ok) then
                    info = terpsichore_reduced_adapter_invalid
                    return
                end if
            end if
        end do
        if (.not. present(mass)) call complete_dynamic_blocks(blocks)
        info = terpsichore_reduced_adapter_ok
    end subroutine assemble_terpsichore_fixture_reduced_mass_with_boundary

    subroutine build_radial_values(fixture, factor, weight, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: factor(:, :), weight(:)
        integer, intent(out) :: info
        real(dp) :: midpoint
        integer :: interval

        allocate (factor(fixture%modes, fixture%intervals), stat=info)
        if (info /= 0) return
        allocate (weight(fixture%intervals), stat=info)
        if (info /= 0) return
        do interval = 1, fixture%intervals
            midpoint = 0.5_dp * (fixture%s(interval - 1) + fixture%s(interval))
            factor(:, interval) = midpoint**(-fixture%radial_power)
            weight(interval) = real(fixture%intervals, dp) &
                * (fixture%s(interval) - fixture%s(interval - 1))
        end do
    end subroutine build_radial_values

end module terpsichore_reduced_mass_adapter
