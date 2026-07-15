module terpsichore_reduced_mass_family_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: add_mapped_dynamic_element, &
        dynamic_family_layout_t, dynamic_layout_ok
    use terpsichore_reduced_layout, only: &
        build_terpsichore_reduced_fixed_boundary_layout, &
        terpsichore_reduced_layout_ok
    use terpsichore_reduced_mass, only: &
        assemble_terpsichore_reduced_mass_element_resolved, &
        terpsichore_reduced_ok
    implicit none
    private

    integer, parameter, public :: terpsichore_reduced_family_ok = 0
    integer, parameter, public :: terpsichore_reduced_family_invalid = -1

    public :: assemble_terpsichore_reduced_fixed_boundary_mass
    public :: assemble_terpsichore_reduced_family_mass_fixed_layout

contains

    subroutine assemble_terpsichore_reduced_fixed_boundary_mass(signed_bjac, &
            flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, trial_m, trial_n, &
            trial_parity, mass, layout, info)
        real(dp), intent(in) :: signed_bjac(:, :), flux_t_slope(:)
        real(dp), intent(in) :: normal_phase(:, :, :)
        real(dp), intent(in) :: tangential_phase(:, :, :)
        real(dp), intent(in) :: normal_radial_factor(:, :)
        real(dp), intent(in) :: normalized_radial_weight(:)
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), allocatable, intent(out) :: mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(out) :: info
        integer, allocatable :: element_to_global(:, :)
        integer :: allocation_status

        info = terpsichore_reduced_family_invalid
        call build_terpsichore_reduced_fixed_boundary_layout(trial_m, trial_n, &
            trial_parity, size(signed_bjac, 2), layout, element_to_global, info)
        if (info /= terpsichore_reduced_layout_ok) return
        allocate (mass(layout%total_unknowns, layout%total_unknowns), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = terpsichore_reduced_family_invalid
            return
        end if
        call assemble_terpsichore_reduced_family_mass_fixed_layout( &
            signed_bjac, flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, &
            element_to_global, mass, info)
    end subroutine assemble_terpsichore_reduced_fixed_boundary_mass

    subroutine assemble_terpsichore_reduced_family_mass_fixed_layout( &
            signed_bjac, flux_t_slope, normal_phase, tangential_phase, &
            normal_radial_factor, normalized_radial_weight, &
            element_to_global, mass, info)
        real(dp), intent(in) :: signed_bjac(:, :), flux_t_slope(:)
        real(dp), intent(in) :: normal_phase(:, :, :)
        real(dp), intent(in) :: tangential_phase(:, :, :)
        real(dp), intent(in) :: normal_radial_factor(:, :)
        real(dp), intent(in) :: normalized_radial_weight(:)
        integer, intent(in) :: element_to_global(:, :)
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: info
        real(dp) :: element(size(element_to_global, 1), &
            size(element_to_global, 1))
        integer :: interval

        info = terpsichore_reduced_family_invalid
        if (.not. valid_family_inputs(signed_bjac, flux_t_slope, &
            normal_phase, tangential_phase, normal_radial_factor, &
            normalized_radial_weight, element_to_global, mass)) return
        mass = 0.0_dp
        do interval = 1, size(signed_bjac, 2)
            call assemble_terpsichore_reduced_mass_element_resolved( &
                signed_bjac(:, interval), flux_t_slope(interval), &
                normal_phase(:, :, interval), &
                tangential_phase(:, :, interval), &
                normal_radial_factor(:, interval), &
                normalized_radial_weight(interval), element, info)
            if (info /= terpsichore_reduced_ok) then
                info = terpsichore_reduced_family_invalid
                return
            end if
            call add_mapped_dynamic_element(element_to_global(:, interval), &
                element, mass, info)
            if (info /= dynamic_layout_ok) then
                info = terpsichore_reduced_family_invalid
                return
            end if
        end do
        info = terpsichore_reduced_family_ok
    end subroutine assemble_terpsichore_reduced_family_mass_fixed_layout

    pure function valid_family_inputs(signed_bjac, flux_t_slope, &
            normal_phase, tangential_phase, normal_radial_factor, &
            normalized_radial_weight, element_to_global, mass) result(valid)
        real(dp), intent(in) :: signed_bjac(:, :), flux_t_slope(:)
        real(dp), intent(in) :: normal_phase(:, :, :)
        real(dp), intent(in) :: tangential_phase(:, :, :)
        real(dp), intent(in) :: normal_radial_factor(:, :)
        real(dp), intent(in) :: normalized_radial_weight(:)
        integer, intent(in) :: element_to_global(:, :)
        real(dp), intent(in) :: mass(:, :)
        logical :: valid
        integer :: global, intervals, trials

        valid = .false.
        trials = size(normal_phase, 1)
        intervals = size(signed_bjac, 2)
        if (trials < 1 .or. intervals < 2 .or. size(signed_bjac, 1) < 1) return
        if (any(shape(tangential_phase) /= shape(normal_phase))) return
        if (size(normal_phase, 2) /= size(signed_bjac, 1)) return
        if (size(normal_phase, 3) /= intervals) return
        if (size(normal_radial_factor, 1) /= trials &
            .or. size(normal_radial_factor, 2) /= intervals) return
        if (size(flux_t_slope) /= intervals) return
        if (size(normalized_radial_weight) /= intervals) return
        if (size(element_to_global, 1) /= 3 * trials &
            .or. size(element_to_global, 2) /= intervals) return
        if (size(mass, 1) /= size(mass, 2) .or. size(mass, 1) < 1) return
        if (any(element_to_global < 0) &
            .or. any(element_to_global > size(mass, 1))) return
        if (.not. all(ieee_is_finite(flux_t_slope)) &
            .or. any(flux_t_slope == 0.0_dp)) return
        if (.not. all(ieee_is_finite(normalized_radial_weight)) &
            .or. any(normalized_radial_weight <= 0.0_dp)) return
        if (.not. all(ieee_is_finite(signed_bjac))) return
        if (.not. all(ieee_is_finite(normal_phase))) return
        if (.not. all(ieee_is_finite(tangential_phase))) return
        if (.not. all(ieee_is_finite(normal_radial_factor))) return
        do global = 1, size(mass, 1)
            if (.not. any(element_to_global == global)) return
        end do
        valid = .true.
    end function valid_family_inputs

end module terpsichore_reduced_mass_family_assembly
