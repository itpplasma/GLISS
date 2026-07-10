module family_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use two_component_kernel, only: two_component_components
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    type, public :: surface_geometry_t
        real(dp), allocatable :: fields(:, :, :)
        real(dp), allocatable :: drive(:, :)
    end type surface_geometry_t

    public :: assemble_family_stiffness
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
            radial_step, lowest, info)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: lowest
        integer, intent(out) :: info
        real(dp), allocatable :: stiffness(:, :), eigenvalues(:), work(:)
        integer :: trials, unknowns

        trials = 2 * size(mode_m)
        unknowns = trials * (size(geometry) - 1)
        call assemble_family_stiffness(geometry, mode_m, mode_n, &
            radial_step, stiffness, info)
        if (info /= 0) return
        allocate (eigenvalues(unknowns), work(8 * unknowns))
        call dsyev("N", "U", unknowns, stiffness, unknowns, eigenvalues, &
            work, size(work), info)
        if (info /= 0) return
        lowest = eigenvalues(1) / radial_step
    end subroutine lowest_family_eigenvalue

    subroutine assemble_family_stiffness(geometry, mode_m, mode_n, &
            radial_step, stiffness, info)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: element(:, :)
        integer :: trials, intervals, i, a, b, row, column

        trials = 2 * size(mode_m)
        intervals = size(geometry)
        allocate (stiffness(trials * (intervals - 1), &
            trials * (intervals - 1)), source=0.0_dp)
        allocate (element(2 * trials, 2 * trials))
        do i = 1, intervals
            call condensed_element(geometry(i), mode_m, mode_n, &
                radial_step, element, info)
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

    subroutine condensed_element(surface, mode_m, mode_n, radial_step, &
            element, info)
        type(surface_geometry_t), intent(in) :: surface
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: element(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: full(:, :)
        integer :: trials, n_theta, n_zeta, j, l, k
        real(dp) :: weight

        trials = 2 * size(mode_m)
        n_theta = size(surface%fields, 1)
        n_zeta = size(surface%fields, 2)
        allocate (full(3 * trials, 3 * trials), source=0.0_dp)
        weight = 1.0_dp / real(n_theta * n_zeta, dp)
        do l = 1, n_zeta
            do j = 1, n_theta
                call accumulate_point(surface%fields(j, l, :), &
                    surface%drive(j, l), mode_m, mode_n, &
                    (real(j, dp) - 1.0_dp) / real(n_theta, dp), &
                    (real(l, dp) - 1.0_dp) / real(n_zeta, dp), &
                    radial_step, weight, full)
            end do
        end do
        call condense_tangential(full, trials, element, info)
        do k = 1, size(element, 1)
            element(:, k) = element(:, k) * radial_step
        end do
    end subroutine condensed_element

    subroutine accumulate_point(fields, drive, mode_m, mode_n, theta, &
            zeta, radial_step, weight, full)
        real(dp), intent(in) :: fields(:), drive, theta, zeta
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step, weight
        real(dp), intent(inout) :: full(:, :)
        real(dp) :: rows(4, 6 * size(mode_m))
        real(dp) :: phase, cosine, sine
        real(dp) :: value, dvalue, dother
        real(dp) :: c1_of(6), c2_of(6), c3_of(6)
        integer :: trials, k, parity, trial, entry_index
        real(dp) :: unit_inputs(6)

        do entry_index = 1, 6
            unit_inputs = 0.0_dp
            unit_inputs(entry_index) = 1.0_dp
            call two_component_components(fields(1), fields(2), &
                fields(3), fields(4), fields(5), fields(6), &
                fields(7), fields(8), fields(9), fields(10), &
                fields(11), fields(12), fields(13), &
                unit_inputs(1), unit_inputs(2), unit_inputs(3), &
                unit_inputs(4), unit_inputs(5), unit_inputs(6), &
                c1_of(entry_index), c2_of(entry_index), &
                c3_of(entry_index))
        end do
        trials = 2 * size(mode_m)
        rows = 0.0_dp
        do k = 1, size(mode_m)
            phase = two_pi * (real(mode_m(k), dp) * theta &
                - real(mode_n(k), dp) * zeta)
            cosine = cos(phase)
            sine = sin(phase)
            do parity = 1, 2
                if (parity == 1) then
                    value = cosine
                    dvalue = -sine
                    dother = cosine
                else
                    value = sine
                    dvalue = cosine
                    dother = -sine
                end if
                trial = 2 * (k - 1) + parity
                do entry_index = 1, 6
                    call add_linear(rows, trial, trials, entry_index, &
                        value, dvalue, dother, mode_m(k), mode_n(k), &
                        radial_step, c1_of(entry_index), &
                        c2_of(entry_index), c3_of(entry_index))
                end do
                rows(4, trial) = rows(4, trial) + 0.5_dp * value
                rows(4, trials + trial) = rows(4, trials + trial) &
                    + 0.5_dp * value
            end do
        end do
        call rank_updates(rows, drive, weight * abs(fields(7)), full)
    end subroutine accumulate_point

    subroutine add_linear(rows, trial, trials, entry_index, value, &
            dvalue, dother, m, n, radial_step, c1, c2, c3)
        real(dp), intent(inout) :: rows(:, :)
        integer, intent(in) :: trial, trials, entry_index, m, n
        real(dp), intent(in) :: value, dvalue, dother, radial_step
        real(dp), intent(in) :: c1, c2, c3

        select case (entry_index)
        case (1)
            call apply(rows, trial, trials, value * 0.5_dp, &
                value * 0.5_dp, 0.0_dp, c1, c2, c3)
        case (2)
            call apply(rows, trial, trials, -value / radial_step, &
                value / radial_step, 0.0_dp, c1, c2, c3)
        case (3)
            call apply(rows, trial, trials, &
                two_pi * real(m, dp) * dvalue * 0.5_dp, &
                two_pi * real(m, dp) * dvalue * 0.5_dp, 0.0_dp, c1, c2, &
                c3)
        case (4)
            call apply(rows, trial, trials, &
                -two_pi * real(n, dp) * dvalue * 0.5_dp, &
                -two_pi * real(n, dp) * dvalue * 0.5_dp, 0.0_dp, c1, c2, &
                c3)
        case (5)
            call apply(rows, trial, trials, 0.0_dp, 0.0_dp, &
                two_pi * real(m, dp) * dother, c1, c2, c3)
        case (6)
            call apply(rows, trial, trials, 0.0_dp, 0.0_dp, &
                -two_pi * real(n, dp) * dother, c1, c2, c3)
        end select
    end subroutine add_linear

    subroutine apply(rows, trial, trials, left_factor, right_factor, &
            tangential_factor, c1, c2, c3)
        real(dp), intent(inout) :: rows(:, :)
        integer, intent(in) :: trial, trials
        real(dp), intent(in) :: left_factor, right_factor
        real(dp), intent(in) :: tangential_factor, c1, c2, c3

        rows(1, trial) = rows(1, trial) + c1 * left_factor
        rows(1, trials + trial) = rows(1, trials + trial) &
            + c1 * right_factor
        rows(1, 2 * trials + trial) = rows(1, 2 * trials + trial) &
            + c1 * tangential_factor
        rows(2, trial) = rows(2, trial) + c2 * left_factor
        rows(2, trials + trial) = rows(2, trials + trial) &
            + c2 * right_factor
        rows(2, 2 * trials + trial) = rows(2, 2 * trials + trial) &
            + c2 * tangential_factor
        rows(3, trial) = rows(3, trial) + c3 * left_factor
        rows(3, trials + trial) = rows(3, trials + trial) &
            + c3 * right_factor
        rows(3, 2 * trials + trial) = rows(3, 2 * trials + trial) &
            + c3 * tangential_factor
    end subroutine apply

    subroutine rank_updates(rows, drive, weight, full)
        real(dp), intent(in) :: rows(:, :), drive, weight
        real(dp), intent(inout) :: full(:, :)
        integer :: a, b, component

        do b = 1, size(full, 2)
            do a = 1, size(full, 1)
                do component = 1, 3
                    full(a, b) = full(a, b) + weight &
                        * rows(component, a) * rows(component, b)
                end do
                full(a, b) = full(a, b) - weight * drive &
                    * rows(4, a) * rows(4, b)
            end do
        end do
    end subroutine rank_updates

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
