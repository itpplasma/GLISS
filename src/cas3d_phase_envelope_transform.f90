module cas3d_phase_envelope_transform
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fourier_phase_kind, only: phase_cosine, phase_sine
    implicit none
    private

    integer, parameter, public :: cas3d_phase_transform_ok = 0
    integer, parameter, public :: cas3d_phase_transform_invalid = -1
    integer, parameter, public :: cas3d_phase_transform_allocation = -2

    type, public :: cas3d_phase_envelope_map_t
        integer :: physical_mode_count = 0
        integer :: envelope_mode_count = 0
        integer, allocatable :: physical_row(:, :)
        real(dp), allocatable :: normal_weight(:, :)
        real(dp), allocatable :: eta_weight(:, :)
    end type cas3d_phase_envelope_map_t

    public :: apply_cas3d_phase_envelope_congruence
    public :: build_cas3d_phase_envelope_map

contains

    subroutine build_cas3d_phase_envelope_map(labeled_m, labeled_n, &
            orientation, physical_m, physical_n, parity, map, info)
        integer, intent(in) :: labeled_m(:), labeled_n(:), orientation(:)
        integer, intent(in) :: physical_m(:), physical_n(:), parity
        type(cas3d_phase_envelope_map_t), intent(out) :: map
        integer, intent(out) :: info
        integer :: allocation_status, column, first_row, second_row
        real(dp) :: first_eta_sign, first_normal_sign
        real(dp) :: second_eta_sign, second_normal_sign

        map = cas3d_phase_envelope_map_t()
        info = cas3d_phase_transform_invalid
        if (size(labeled_m) < 1 .or. modulo(size(labeled_m), 2) /= 1) return
        if (size(labeled_n) /= size(labeled_m)) return
        if (size(orientation) /= size(labeled_m)) return
        if (size(physical_m) < 1) return
        if (size(physical_n) /= size(physical_m)) return
        if (any(abs(orientation) /= 1)) return
        if (parity /= phase_cosine .and. parity /= phase_sine) return
        allocate (map%physical_row(2, size(labeled_m)), source=0, &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        allocate (map%normal_weight(2, size(labeled_m)), source=0.0_dp, &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        allocate (map%eta_weight(2, size(labeled_m)), source=0.0_dp, &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        first_row = find_physical_row(labeled_m(1), labeled_n(1), &
            physical_m, physical_n)
        if (first_row == 0) return
        map%physical_row(1, 1) = first_row
        map%normal_weight(1, 1) = phase_orientation_sign(parity, &
            orientation(1))
        map%eta_weight(1, 1) = phase_orientation_sign(3 - parity, &
            orientation(1))
        do column = 2, size(labeled_m), 2
            first_row = find_physical_row(labeled_m(column), &
                labeled_n(column), physical_m, physical_n)
            second_row = find_physical_row(labeled_m(column + 1), &
                labeled_n(column + 1), physical_m, physical_n)
            if (first_row == 0 .or. second_row == 0) return
            map%physical_row(1, column) = first_row
            map%physical_row(2, column) = second_row
            map%physical_row(1, column + 1) = first_row
            map%physical_row(2, column + 1) = second_row
            first_normal_sign = phase_orientation_sign(parity, &
                orientation(column))
            second_normal_sign = phase_orientation_sign(parity, &
                orientation(column + 1))
            first_eta_sign = phase_orientation_sign(3 - parity, &
                orientation(column))
            second_eta_sign = phase_orientation_sign(3 - parity, &
                orientation(column + 1))
            map%normal_weight(1, column) = 0.5_dp * first_normal_sign
            map%normal_weight(2, column) = 0.5_dp * second_normal_sign
            map%normal_weight(1, column + 1) = -0.5_dp * first_normal_sign
            map%normal_weight(2, column + 1) = 0.5_dp * second_normal_sign
            map%eta_weight(1, column) = 0.5_dp * first_eta_sign
            map%eta_weight(2, column) = 0.5_dp * second_eta_sign
            map%eta_weight(1, column + 1) = 0.5_dp * first_eta_sign
            map%eta_weight(2, column + 1) = -0.5_dp * second_eta_sign
        end do
        do column = 1, size(labeled_m)
            call canonicalize_support(map, column)
        end do
        map%physical_mode_count = size(physical_m)
        map%envelope_mode_count = size(labeled_m)
        if (.not. valid_map(map)) return
        info = cas3d_phase_transform_ok
    end subroutine build_cas3d_phase_envelope_map

    subroutine apply_cas3d_phase_envelope_congruence(map, h1_dofs, l2_dofs, &
            mass_scale, stiffness, stiffness_terms, mass, info)
        type(cas3d_phase_envelope_map_t), intent(in) :: map
        integer, intent(in) :: h1_dofs, l2_dofs
        real(dp), intent(in) :: mass_scale
        real(dp), allocatable, intent(inout) :: stiffness(:, :)
        real(dp), allocatable, intent(inout) :: stiffness_terms(:, :, :)
        real(dp), allocatable, intent(inout) :: mass(:, :)
        integer, intent(out) :: info
        integer, allocatable :: row(:, :)
        real(dp), allocatable :: transformed(:, :), transformed_terms(:, :, :)
        real(dp), allocatable :: transformed_mass(:, :), weight(:, :)
        integer :: allocation_status, basis, column, envelope_unknowns
        integer :: physical_unknowns

        info = cas3d_phase_transform_invalid
        if (.not. valid_map(map)) return
        if (h1_dofs < 0 .or. l2_dofs < 0) return
        if (.not. ieee_is_finite(mass_scale) .or. mass_scale <= 0.0_dp) return
        if (h1_dofs > huge(physical_unknowns) - l2_dofs) return
        if (h1_dofs + l2_dofs > huge(physical_unknowns) / &
            map%physical_mode_count) return
        if (h1_dofs + l2_dofs > huge(envelope_unknowns) / &
            map%envelope_mode_count) return
        physical_unknowns = (h1_dofs + l2_dofs) * &
            map%physical_mode_count
        envelope_unknowns = (h1_dofs + l2_dofs) * &
            map%envelope_mode_count
        if (physical_unknowns < 1 .or. envelope_unknowns < 1) return
        if (.not. valid_pencil_shape(stiffness, stiffness_terms, mass, &
            physical_unknowns)) return
        allocate (row(2, envelope_unknowns), source=0, &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        allocate (weight(2, envelope_unknowns), source=0.0_dp, &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        do basis = 1, h1_dofs
            do column = 1, map%envelope_mode_count
                row(:, (basis - 1) * map%envelope_mode_count + column) = &
                    (basis - 1) * map%physical_mode_count &
                    + map%physical_row(:, column)
                weight(:, (basis - 1) * map%envelope_mode_count + column) = &
                    map%normal_weight(:, column)
            end do
        end do
        do basis = 1, l2_dofs
            do column = 1, map%envelope_mode_count
                row(:, h1_dofs * map%envelope_mode_count &
                    + (basis - 1) * map%envelope_mode_count + column) = &
                    h1_dofs * map%physical_mode_count &
                    + (basis - 1) * map%physical_mode_count &
                    + map%physical_row(:, column)
                weight(:, h1_dofs * map%envelope_mode_count &
                    + (basis - 1) * map%envelope_mode_count + column) = &
                    map%eta_weight(:, column)
            end do
        end do
        allocate (transformed(envelope_unknowns, envelope_unknowns), &
            source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        allocate (transformed_terms(envelope_unknowns, envelope_unknowns, &
            size(stiffness_terms, 3)), source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        allocate (transformed_mass(envelope_unknowns, envelope_unknowns), &
            source=0.0_dp, stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        call sparse_congruence(row, weight, stiffness, stiffness_terms, &
            transformed, transformed_terms)
        do column = 1, envelope_unknowns
            transformed_mass(column, column) = mass_scale
        end do
        call move_alloc(transformed, stiffness)
        call move_alloc(transformed_terms, stiffness_terms)
        call move_alloc(transformed_mass, mass)
        info = cas3d_phase_transform_ok
    end subroutine apply_cas3d_phase_envelope_congruence

    subroutine sparse_congruence(row, weight, source, source_terms, &
            target, target_terms)
        integer, intent(in) :: row(:, :)
        real(dp), intent(in) :: weight(:, :), source(:, :)
        real(dp), intent(in) :: source_terms(:, :, :)
        real(dp), intent(out) :: target(:, :), target_terms(:, :, :)
        real(dp) :: factor
        integer :: first, first_support, second, second_support, term

        target = 0.0_dp
        target_terms = 0.0_dp
        do second = 1, size(target, 2)
            do first = 1, second
                do second_support = 1, 2
                    if (row(second_support, second) == 0) cycle
                    do first_support = 1, 2
                        if (row(first_support, first) == 0) cycle
                        factor = weight(first_support, first) &
                            * weight(second_support, second)
                        target(first, second) = target(first, second) &
                            + factor * source(row(first_support, first), &
                            row(second_support, second))
                        do term = 1, size(target_terms, 3)
                            target_terms(first, second, term) = &
                                target_terms(first, second, term) &
                                + factor * source_terms( &
                                row(first_support, first), &
                                row(second_support, second), term)
                        end do
                    end do
                end do
                target(second, first) = target(first, second)
                do term = 1, size(target_terms, 3)
                    target_terms(second, first, term) = &
                        target_terms(first, second, term)
                end do
            end do
        end do
    end subroutine sparse_congruence

    pure function find_physical_row(mode_m, mode_n, physical_m, physical_n) &
            result(row)
        integer, intent(in) :: mode_m, mode_n, physical_m(:), physical_n(:)
        integer :: row

        row = findloc(physical_m == mode_m .and. physical_n == mode_n, &
            .true., dim=1)
    end function find_physical_row

    pure function phase_orientation_sign(phase, orientation) result(sign)
        integer, intent(in) :: phase, orientation
        real(dp) :: sign

        sign = 1.0_dp
        if (phase == phase_sine .and. orientation == -1) sign = -1.0_dp
    end function phase_orientation_sign

    pure subroutine canonicalize_support(map, column)
        type(cas3d_phase_envelope_map_t), intent(inout) :: map
        integer, intent(in) :: column
        real(dp) :: temporary_weight
        integer :: temporary_row

        if (map%physical_row(2, column) == 0) return
        if (map%physical_row(1, column) > map%physical_row(2, column)) then
            temporary_row = map%physical_row(1, column)
            map%physical_row(1, column) = map%physical_row(2, column)
            map%physical_row(2, column) = temporary_row
            temporary_weight = map%normal_weight(1, column)
            map%normal_weight(1, column) = map%normal_weight(2, column)
            map%normal_weight(2, column) = temporary_weight
            temporary_weight = map%eta_weight(1, column)
            map%eta_weight(1, column) = map%eta_weight(2, column)
            map%eta_weight(2, column) = temporary_weight
        end if
        if (map%physical_row(1, column) &
            /= map%physical_row(2, column)) return
        map%normal_weight(1, column) = map%normal_weight(1, column) &
            + map%normal_weight(2, column)
        map%eta_weight(1, column) = map%eta_weight(1, column) &
            + map%eta_weight(2, column)
        map%physical_row(2, column) = 0
        map%normal_weight(2, column) = 0.0_dp
        map%eta_weight(2, column) = 0.0_dp
    end subroutine canonicalize_support

    pure function valid_map(map) result(valid)
        type(cas3d_phase_envelope_map_t), intent(in) :: map
        integer :: physical_row
        logical :: valid

        valid = map%physical_mode_count >= 1 &
            .and. map%envelope_mode_count >= 1
        if (.not. valid) return
        valid = allocated(map%physical_row) &
            .and. allocated(map%normal_weight) &
            .and. allocated(map%eta_weight)
        if (.not. valid) return
        valid = size(map%physical_row, 1) == 2 &
            .and. size(map%physical_row, 2) == map%envelope_mode_count
        if (.not. valid) return
        valid = all(map%physical_row(1, :) >= 1) &
            .and. all(map%physical_row >= 0) &
            .and. all(map%physical_row <= map%physical_mode_count)
        if (.not. valid) return
        valid = size(map%normal_weight, 1) == size(map%physical_row, 1) &
            .and. size(map%normal_weight, 2) == size(map%physical_row, 2) &
            .and. size(map%eta_weight, 1) == size(map%physical_row, 1) &
            .and. size(map%eta_weight, 2) == size(map%physical_row, 2)
        if (.not. valid) return
        valid = all(ieee_is_finite(map%normal_weight)) &
            .and. all(ieee_is_finite(map%eta_weight))
        if (.not. valid) return
        do physical_row = 1, map%physical_mode_count
            valid = any(map%physical_row == physical_row &
                .and. map%normal_weight /= 0.0_dp)
            if (.not. valid) return
            valid = any(map%physical_row == physical_row &
                .and. map%eta_weight /= 0.0_dp)
            if (.not. valid) return
        end do
    end function valid_map

    pure function valid_pencil_shape(stiffness, stiffness_terms, mass, &
            unknowns) result(valid)
        real(dp), intent(in) :: stiffness(:, :), stiffness_terms(:, :, :)
        real(dp), intent(in) :: mass(:, :)
        integer, intent(in) :: unknowns
        logical :: valid

        valid = size(stiffness, 1) == unknowns &
            .and. size(stiffness, 2) == unknowns &
            .and. size(mass, 1) == unknowns &
            .and. size(mass, 2) == unknowns &
            .and. size(stiffness_terms, 1) == unknowns &
            .and. size(stiffness_terms, 2) == unknowns &
            .and. size(stiffness_terms, 3) >= 1
    end function valid_pencil_shape

end module cas3d_phase_envelope_transform
