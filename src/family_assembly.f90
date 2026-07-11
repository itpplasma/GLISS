module family_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use block_tridiagonal, only: apply_block_tridiagonal, &
        block_factor_t, block_tridiagonal_t, factorize_shifted, &
        solve_factored
    use family_point_assembly, only: assemble_direct_surface_resolved, &
        assemble_transformed_surface_resolved, resolve_normal_stored_power
    use phase_assembly_policy, only: phase_assembly_direct, &
        phase_assembly_transformed
    use radial_space_policy, only: radial_space_config_t, radial_space_ok, &
        validate_radial_space
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
    public :: phase_assembly_direct
    public :: phase_assembly_transformed

    type, public :: surface_geometry_t
        real(dp), allocatable :: fields(:, :, :)
        real(dp), allocatable :: drive(:, :)
    end type surface_geometry_t

    type, public :: family_assembly_options_t
        integer :: field_periods = 1
        integer :: parity_class = 0
        integer :: phase_assembly = phase_assembly_transformed
        type(radial_space_config_t) :: radial_space
    end type family_assembly_options_t

    public :: assemble_family_blocks
    public :: assemble_family_stiffness
    public :: family_negative_count
    public :: iterate_block_eigenvalue
    public :: iterate_family_eigenvalue
    public :: lowest_family_eigenvalue

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*)
            real(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine dsyev
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dposv
    end interface

contains

    subroutine lowest_family_eigenvalue(geometry, mode_m, mode_n, &
            radial_step, lowest, info, options, normal_stored_power)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: lowest
        integer, intent(out) :: info
        type(family_assembly_options_t), intent(in), optional :: options
        real(dp), intent(in), optional :: normal_stored_power(:)
        real(dp), allocatable :: stiffness(:, :), eigenvalues(:), work(:)
        integer :: unknowns

        call assemble_family_stiffness(geometry, mode_m, mode_n, &
            radial_step, stiffness, info, options, normal_stored_power)
        if (info /= 0) return
        unknowns = size(stiffness, 1)
        allocate (eigenvalues(unknowns), work(8 * unknowns))
        call dsyev("N", "U", unknowns, stiffness, unknowns, eigenvalues, &
            work, size(work), info)
        if (info /= 0) return
        lowest = eigenvalues(1) / radial_step
    end subroutine lowest_family_eigenvalue

    subroutine assemble_family_blocks(geometry, mode_m, mode_n, &
            radial_step, blocks, info, options, normal_stored_power)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        type(block_tridiagonal_t), intent(out) :: blocks
        integer, intent(out) :: info
        type(family_assembly_options_t), intent(in), optional :: options
        real(dp), intent(in), optional :: normal_stored_power(:)
        real(dp), allocatable :: element(:, :)
        integer, allocatable :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), allocatable :: mode_stored_power(:), trial_stored_power(:)
        integer :: trials, intervals, nodes, i, periods, selector
        integer :: phase_assembly
        type(radial_space_config_t) :: radial_space

        call validate_family_inputs(geometry, mode_m, mode_n, radial_step, info)
        if (info /= 0) return
        call resolve_options(options, periods, selector, phase_assembly, &
            radial_space, info)
        if (info /= 0) return
        call resolve_normal_stored_power(normal_stored_power, size(mode_m), &
            mode_stored_power, info)
        if (info /= 0) then
            info = -2
            return
        end if
        call build_trial_tables(mode_m, mode_n, mode_stored_power, selector, &
            trial_m, trial_n, trial_parity, trial_stored_power)
        trials = size(trial_m)
        intervals = size(geometry)
        nodes = intervals - 1
        allocate (blocks%diag(trials, trials, nodes), source=0.0_dp)
        allocate (blocks%off(trials, trials, nodes - 1), source=0.0_dp)
        allocate (element(2 * trials, 2 * trials))
        do i = 1, intervals
            call condensed_element(geometry(i), trial_m, trial_n, &
                trial_parity, trial_stored_power, periods, phase_assembly, &
                radial_space, &
                (real(i, dp) - 0.5_dp) * radial_step, radial_step, &
                element, info)
            if (info /= 0) return
            if (i > 1) then
                blocks%diag(:, :, i - 1) = blocks%diag(:, :, i - 1) &
                    + element(1:trials, 1:trials)
            end if
            if (i <= nodes) then
                blocks%diag(:, :, i) = blocks%diag(:, :, i) &
                    + element(trials + 1:, trials + 1:)
            end if
            if (i > 1 .and. i <= nodes) then
                blocks%off(:, :, i - 1) = blocks%off(:, :, i - 1) &
                    + element(trials + 1:, 1:trials)
            end if
        end do
        info = 0
    end subroutine assemble_family_blocks

    subroutine family_negative_count(geometry, mode_m, mode_n, &
            radial_step, shift, count, info, options, normal_stored_power)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step, shift
        integer, intent(out) :: count
        integer, intent(out) :: info
        type(family_assembly_options_t), intent(in), optional :: options
        real(dp), intent(in), optional :: normal_stored_power(:)
        type(block_tridiagonal_t) :: blocks
        type(block_factor_t) :: factor

        count = -1
        call assemble_family_blocks(geometry, mode_m, mode_n, &
            radial_step, blocks, info, options, normal_stored_power)
        if (info /= 0) return
        call factorize_shifted(blocks, shift * radial_step, factor, info)
        if (info /= 0) return
        count = factor%negative_count
    end subroutine family_negative_count

    subroutine iterate_family_eigenvalue(geometry, mode_m, mode_n, &
            radial_step, shift, eigenvalue, info, options, normal_stored_power)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step, shift
        real(dp), intent(out) :: eigenvalue
        integer, intent(out) :: info
        type(family_assembly_options_t), intent(in), optional :: options
        real(dp), intent(in), optional :: normal_stored_power(:)
        type(block_tridiagonal_t) :: blocks

        call assemble_family_blocks(geometry, mode_m, mode_n, &
            radial_step, blocks, info, options, normal_stored_power)
        if (info /= 0) return
        call iterate_block_eigenvalue(blocks, radial_step, shift, &
            eigenvalue, info)
    end subroutine iterate_family_eigenvalue

    subroutine iterate_block_eigenvalue(blocks, radial_step, shift, &
            eigenvalue, info)
        type(block_tridiagonal_t), intent(in) :: blocks
        real(dp), intent(in) :: radial_step, shift
        real(dp), intent(out) :: eigenvalue
        integer, intent(out) :: info
        type(block_factor_t) :: factor
        real(dp), allocatable :: vector(:, :)
        real(dp) :: rayleigh, previous
        integer :: trials, nodes, iteration, t, j

        call factorize_shifted(blocks, shift * radial_step, factor, info)
        if (info /= 0) return
        trials = size(blocks%diag, 1)
        nodes = size(blocks%diag, 3)
        allocate (vector(trials, nodes))
        do j = 1, nodes
            do t = 1, trials
                vector(t, j) = 1.0_dp + 0.1_dp * real(t, dp) &
                    + 0.01_dp * real(j, dp)
            end do
        end do
        vector = vector / norm2(vector)
        previous = huge(1.0_dp)
        do iteration = 1, 500
            call solve_factored(blocks, factor, vector, info)
            if (info /= 0) return
            vector = vector / norm2(vector)
            rayleigh = sum(vector * apply_block_tridiagonal(blocks, &
                vector))
            if (abs(rayleigh - previous) <= 1.0e-13_dp &
                * max(1.0_dp, abs(rayleigh))) exit
            previous = rayleigh
        end do
        if (abs(rayleigh - previous) > 1.0e-11_dp &
            * max(1.0_dp, abs(rayleigh))) then
            info = -1
            return
        end if
        eigenvalue = rayleigh / radial_step
        info = 0
    end subroutine iterate_block_eigenvalue

    subroutine assemble_family_stiffness(geometry, mode_m, mode_n, &
            radial_step, stiffness, info, options, normal_stored_power)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        type(family_assembly_options_t), intent(in), optional :: options
        real(dp), intent(in), optional :: normal_stored_power(:)
        real(dp), allocatable :: element(:, :)
        integer, allocatable :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), allocatable :: mode_stored_power(:), trial_stored_power(:)
        integer :: trials, intervals, i, a, b, row, column, periods, selector
        integer :: phase_assembly
        type(radial_space_config_t) :: radial_space

        call validate_family_inputs(geometry, mode_m, mode_n, radial_step, info)
        if (info /= 0) return
        call resolve_options(options, periods, selector, phase_assembly, &
            radial_space, info)
        if (info /= 0) return
        call resolve_normal_stored_power(normal_stored_power, size(mode_m), &
            mode_stored_power, info)
        if (info /= 0) then
            info = -2
            return
        end if
        call build_trial_tables(mode_m, mode_n, mode_stored_power, selector, &
            trial_m, trial_n, trial_parity, trial_stored_power)
        trials = size(trial_m)
        intervals = size(geometry)
        allocate (stiffness(trials * (intervals - 1), &
            trials * (intervals - 1)), source=0.0_dp)
        allocate (element(2 * trials, 2 * trials))
        do i = 1, intervals
            call condensed_element(geometry(i), trial_m, trial_n, &
                trial_parity, trial_stored_power, periods, phase_assembly, &
                radial_space, &
                (real(i, dp) - 0.5_dp) * radial_step, radial_step, &
                element, info)
            if (info /= 0) return
            do b = 1, 2 * trials
                column = global_index(i, b, trials, intervals)
                if (column == 0) cycle
                do a = 1, 2 * trials
                    row = global_index(i, a, trials, intervals)
                    if (row == 0) cycle
                    stiffness(row, column) = stiffness(row, column) &
                        + element(a, b)
                end do
            end do
        end do
        info = 0
    end subroutine assemble_family_stiffness

    pure subroutine resolve_options(options, periods, selector, &
            phase_assembly, radial_space, info)
        type(family_assembly_options_t), intent(in), optional :: options
        type(radial_space_config_t), intent(out) :: radial_space
        integer, intent(out) :: periods, selector, phase_assembly, info

        periods = 1
        selector = 0
        phase_assembly = phase_assembly_transformed
        if (present(options)) then
            periods = options%field_periods
            selector = options%parity_class
            phase_assembly = options%phase_assembly
            radial_space = options%radial_space
        end if
        info = -2
        if (periods < 1) return
        if (selector < 0 .or. selector > 2) return
        if (phase_assembly /= phase_assembly_transformed .and. &
            phase_assembly /= phase_assembly_direct) return
        call validate_radial_space(radial_space, info)
        if (info /= radial_space_ok) then
            info = -2
            return
        end if
        info = 0
    end subroutine resolve_options

    pure subroutine validate_family_inputs(geometry, mode_m, mode_n, &
            radial_step, info)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        integer, intent(out) :: info

        info = -2
        if (size(geometry) < 2) return
        if (size(mode_m) < 1) return
        if (size(mode_n) /= size(mode_m)) return
        if (.not. ieee_is_finite(radial_step)) return
        if (radial_step <= 0.0_dp) return
        info = 0
    end subroutine validate_family_inputs

    pure subroutine build_trial_tables(mode_m, mode_n, normal_stored_power, &
            selector, trial_m, trial_n, trial_parity, trial_stored_power)
        integer, intent(in) :: mode_m(:), mode_n(:), selector
        real(dp), intent(in) :: normal_stored_power(:)
        integer, allocatable, intent(out) :: trial_m(:), trial_n(:)
        integer, allocatable, intent(out) :: trial_parity(:)
        real(dp), allocatable, intent(out) :: trial_stored_power(:)
        integer :: k, parity, t, trials

        trials = 2 * size(mode_m)
        if (selector /= 0) trials = size(mode_m)
        allocate (trial_m(trials), trial_n(trials), &
            trial_parity(trials), trial_stored_power(trials))
        t = 0
        do k = 1, size(mode_m)
            do parity = 1, 2
                if (selector /= 0) then
                    if (parity /= selector) cycle
                end if
                t = t + 1
                trial_m(t) = mode_m(k)
                trial_n(t) = mode_n(k)
                trial_parity(t) = parity
                trial_stored_power(t) = normal_stored_power(k)
            end do
        end do
    end subroutine build_trial_tables

    pure function global_index(interval, local, trials, intervals) &
            result(index)
        integer, intent(in) :: interval, local, trials, intervals
        integer :: index
        integer :: node, trial

        if (local <= trials) then
            node = interval - 1
            trial = local
        else
            node = interval
            trial = local - trials
        end if
        if (node == 0 .or. node == intervals) then
            index = 0
        else
            index = (node - 1) * trials + trial
        end if
    end function global_index

    subroutine condensed_element(surface, trial_m, trial_n, &
            trial_parity, normal_stored_power, field_periods, phase_assembly, &
            radial_space, radial_coordinate, radial_step, element, info)
        type(surface_geometry_t), intent(in) :: surface
        integer, intent(in) :: trial_m(:), trial_n(:), trial_parity(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: field_periods, phase_assembly
        type(radial_space_config_t), intent(in) :: radial_space
        real(dp), intent(in) :: radial_coordinate, radial_step
        real(dp), intent(out) :: element(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: full(:, :)
        integer :: trials, k

        trials = size(trial_m)
        allocate (full(3 * trials, 3 * trials), source=0.0_dp)
        if (phase_assembly == phase_assembly_direct) then
            call assemble_direct_surface_resolved(surface%fields, &
                surface%drive, trial_m, trial_n, trial_parity, &
                normal_stored_power, field_periods, radial_space, &
                radial_coordinate, radial_step, full, info)
        else
            call assemble_transformed_surface_resolved(surface%fields, &
                surface%drive, trial_m, trial_n, trial_parity, &
                normal_stored_power, field_periods, radial_space, &
                radial_coordinate, radial_step, full, info)
        end if
        if (info /= 0) return
        call condense_tangential(full, trials, element, info)
        do k = 1, size(element, 1)
            element(:, k) = element(:, k) * radial_step
        end do
    end subroutine condensed_element

    subroutine condense_tangential(full, trials, element, info)
        real(dp), intent(in) :: full(:, :)
        integer, intent(in) :: trials
        real(dp), intent(out) :: element(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: yy(:, :), yx(:, :)
        integer :: two_k

        two_k = 2 * trials
        yy = full(two_k + 1:, two_k + 1:)
        yx = full(two_k + 1:, 1:two_k)
        call dposv("U", trials, two_k, yy, trials, yx, trials, info)
        if (info /= 0) return
        element = full(1:two_k, 1:two_k) &
            - matmul(full(1:two_k, two_k + 1:), yx)
    end subroutine condense_tangential

end module family_assembly
