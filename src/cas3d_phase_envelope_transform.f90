module cas3d_phase_envelope_transform
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fourier_phase_kind, only: phase_cosine, phase_sine
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
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

    ! Dimensions of the labeled coefficient pencil. The envelope columns of
    ! one radial basis function map onto the physical modes of the same
    ! function (degree one), so the labeled pencil is the congruence of the
    ! physical one by a block-diagonal C of full row rank: its rank per
    ! component is the number of physical modes reached, the rest of the
    ! labels span an exact null space, and the physical inertia is the
    ! quotient inertia.
    type, public :: cas3d_quotient_t
        integer :: envelope_unknowns = 0
        integer :: physical_unknowns = 0
        integer :: quotient_rank = 0
        integer :: nullity = 0
        integer :: peak_block_width = 0
    end type cas3d_quotient_t

    public :: apply_cas3d_phase_envelope_block_congruence
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


    subroutine apply_cas3d_phase_envelope_block_congruence(map, h1_dofs, &
            l2_dofs, mass_scale, stiffness, mass, quotient, info)
        type(cas3d_phase_envelope_map_t), intent(in) :: map
        integer, intent(in) :: h1_dofs, l2_dofs
        real(dp), intent(in) :: mass_scale
        type(variable_block_tridiagonal_t), intent(inout) :: stiffness, mass
        type(cas3d_quotient_t), intent(out) :: quotient
        integer, intent(out) :: info
        type(variable_block_tridiagonal_t) :: labeled, labeled_mass
        integer, allocatable :: rows(:, :, :)
        real(dp), allocatable :: weights(:, :, :)
        integer :: allocation_status, block, blocks, column, envelope, width

        ! Degree-one block pencil of the physical problem: block b holds the
        ! normal unknowns of H1 function b (b <= h1_dofs), mode-minor, then
        ! the eta unknowns of L2 function b. The labeled pencil has the same
        ! blocks with envelope columns; diagonal blocks map as C_b^T A_b C_b
        ! and lower blocks as C_(b+1)^T L_b C_b. The coefficient mass is
        ! mass_scale times the identity, never a congruence.
        info = cas3d_phase_transform_invalid
        if (.not. valid_map(map)) return
        if (h1_dofs < 0 .or. l2_dofs /= h1_dofs + 1) return
        if (.not. ieee_is_finite(mass_scale) .or. mass_scale <= 0.0_dp) return
        blocks = l2_dofs
        if (.not. allocated(stiffness%widths)) return
        if (size(stiffness%widths) /= blocks) return
        if (.not. allocated(mass%widths)) return
        if (any(mass%widths /= stiffness%widths)) return
        do block = 1, blocks
            if (stiffness%widths(block) /= physical_width(block)) return
        end do
        envelope = map%envelope_mode_count
        allocate (rows(2, 2 * envelope, blocks), &
            weights(2, 2 * envelope, blocks), stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        rows = 0
        weights = 0.0_dp
        allocate (labeled%widths(blocks), labeled%diagonal(blocks), &
            labeled%lower(max(0, blocks - 1)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = cas3d_phase_transform_allocation
            return
        end if
        do block = 1, blocks
            width = 0
            if (block <= h1_dofs) then
                do column = 1, envelope
                    rows(:, width + column, block) = map%physical_row(:, column)
                    weights(:, width + column, block) = &
                        map%normal_weight(:, column)
                end do
                width = envelope
            end if
            do column = 1, envelope
                where (map%physical_row(:, column) > 0)
                    rows(:, width + column, block) = &
                        physical_normal_width(block) &
                        + map%physical_row(:, column)
                end where
                weights(:, width + column, block) = map%eta_weight(:, column)
            end do
            labeled%widths(block) = width + envelope
        end do
        do block = 1, blocks
            allocate (labeled%diagonal(block)%values(labeled%widths(block), &
                labeled%widths(block)), stat=allocation_status)
            if (allocation_status /= 0) then
                info = cas3d_phase_transform_allocation
                return
            end if
            call block_congruence(rows(:, :labeled%widths(block), block), &
                weights(:, :labeled%widths(block), block), &
                rows(:, :labeled%widths(block), block), &
                weights(:, :labeled%widths(block), block), &
                stiffness%diagonal(block)%values, &
                labeled%diagonal(block)%values)
            if (block == blocks) cycle
            allocate (labeled%lower(block)%values(labeled%widths(block + 1), &
                labeled%widths(block)), stat=allocation_status)
            if (allocation_status /= 0) then
                info = cas3d_phase_transform_allocation
                return
            end if
            call block_congruence( &
                rows(:, :labeled%widths(block + 1), block + 1), &
                weights(:, :labeled%widths(block + 1), block + 1), &
                rows(:, :labeled%widths(block), block), &
                weights(:, :labeled%widths(block), block), &
                stiffness%lower(block)%values, labeled%lower(block)%values)
        end do
        call identity_blocks(labeled%widths, mass_scale, labeled_mass, info)
        if (info /= cas3d_phase_transform_ok) return
        quotient%physical_unknowns = sum(stiffness%widths)
        quotient%envelope_unknowns = sum(labeled%widths)
        quotient%quotient_rank = (h1_dofs + l2_dofs) &
            * reached_rows(map)
        quotient%nullity = quotient%envelope_unknowns - quotient%quotient_rank
        quotient%peak_block_width = maxval(labeled%widths)
        info = cas3d_phase_transform_invalid
        if (quotient%quotient_rank /= quotient%physical_unknowns) return
        call move_alloc(labeled%widths, stiffness%widths)
        call move_alloc(labeled%diagonal, stiffness%diagonal)
        call move_alloc(labeled%lower, stiffness%lower)
        mass = labeled_mass
        info = cas3d_phase_transform_ok
    contains
        pure integer function physical_normal_width(index) result(count)
            integer, intent(in) :: index

            count = 0
            if (index <= h1_dofs) count = map%physical_mode_count
        end function physical_normal_width

        pure integer function physical_width(index) result(count)
            integer, intent(in) :: index

            count = physical_normal_width(index) + map%physical_mode_count
        end function physical_width
    end subroutine apply_cas3d_phase_envelope_block_congruence

    ! Rank of the per-component map: the number of physical modes reached by
    ! a label with a nonzero weight. Each pair of labels spans its two
    ! modes and a single-support label spans its mode, so this integer count
    ! is the exact rank; valid_map requires every mode to be reached.
    pure integer function reached_rows(map) result(count)
        type(cas3d_phase_envelope_map_t), intent(in) :: map
        integer :: row

        count = 0
        do row = 1, map%physical_mode_count
            if (any(map%physical_row == row .and. map%normal_weight /= 0.0_dp) &
                .and. any(map%physical_row == row &
                .and. map%eta_weight /= 0.0_dp)) count = count + 1
        end do
    end function reached_rows

    pure subroutine block_congruence(left_rows, left_weights, right_rows, &
            right_weights, source, target)
        integer, intent(in) :: left_rows(:, :), right_rows(:, :)
        real(dp), intent(in) :: left_weights(:, :), right_weights(:, :)
        real(dp), intent(in) :: source(:, :)
        real(dp), intent(out) :: target(:, :)
        integer :: first, left, right, second

        ! target = C_left^T source C_right with at most two supports per
        ! column of C.
        target = 0.0_dp
        do second = 1, size(target, 2)
            do first = 1, size(target, 1)
                do right = 1, 2
                    if (right_rows(right, second) == 0) cycle
                    do left = 1, 2
                        if (left_rows(left, first) == 0) cycle
                        target(first, second) = target(first, second) &
                            + left_weights(left, first) &
                            * right_weights(right, second) &
                            * source(left_rows(left, first), &
                            right_rows(right, second))
                    end do
                end do
            end do
        end do
    end subroutine block_congruence

    subroutine identity_blocks(widths, scale, blocks, info)
        integer, intent(in) :: widths(:)
        real(dp), intent(in) :: scale
        type(variable_block_tridiagonal_t), intent(out) :: blocks
        integer, intent(out) :: info
        integer :: allocation_status, block, column

        info = cas3d_phase_transform_allocation
        allocate (blocks%widths, source=widths, stat=allocation_status)
        if (allocation_status /= 0) return
        allocate (blocks%diagonal(size(widths)), &
            blocks%lower(max(0, size(widths) - 1)), stat=allocation_status)
        if (allocation_status /= 0) return
        do block = 1, size(widths)
            allocate (blocks%diagonal(block)%values(widths(block), &
                widths(block)), source=0.0_dp, stat=allocation_status)
            if (allocation_status /= 0) return
            do column = 1, widths(block)
                blocks%diagonal(block)%values(column, column) = scale
            end do
            if (block == size(widths)) cycle
            allocate (blocks%lower(block)%values(widths(block + 1), &
                widths(block)), source=0.0_dp, stat=allocation_status)
            if (allocation_status /= 0) return
        end do
        info = cas3d_phase_transform_ok
    end subroutine identity_blocks

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
