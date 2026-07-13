module starwall_fourier_coupling
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use dynamic_family_layout, only: dynamic_family_layout_t, &
        normal_global_index
    use fourier_phase_kind, only: phase_cosine, phase_sine, valid_phase_kind
    use trial_space_topology, only: trial_space_topology_t
    implicit none
    private

    integer, parameter, public :: starwall_fourier_ok = 0
    integer, parameter, public :: starwall_fourier_invalid_input = 1
    integer, parameter, public :: starwall_fourier_underresolved = 2
    integer, parameter, public :: starwall_fourier_nonsymmetric = 3
    integer, parameter, public :: starwall_fourier_invalid_layout = 4

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: add_starwall_fourier_stiffness
    public :: build_starwall_fourier_map

contains

    subroutine build_starwall_fourier_map(nu, nv, topology, map, info)
        integer, intent(in) :: nu, nv
        type(trial_space_topology_t), intent(in) :: topology
        real(dp), allocatable, intent(out) :: map(:, :)
        integer, intent(out) :: info
        real(dp) :: phase, u, v
        integer :: i, k, node, trial, trials

        info = starwall_fourier_invalid_input
        if (.not. valid_topology(topology)) return
        trials = size(topology%poloidal)
        if (duplicate_normal_trial(topology)) return
        if (.not. grid_size_is_representable(nu, nv)) return
        if (.not. grid_resolves_modes(nu, nv, topology)) then
            info = starwall_fourier_underresolved
            return
        end if

        allocate (map(nu * nv, trials), source=0.0_dp)
        do k = 1, nv
            v = real(k - 1, dp) / real(nv, dp)
            do i = 1, nu
                u = real(i - 1, dp) / real(nu, dp)
                node = i + nu * (k - 1)
                do trial = 1, trials
                    if (.not. topology%active(1, trial)) cycle
                    ! v spans the full torus, so physical n has no N_FP factor.
                    phase = two_pi * (real(topology%poloidal(trial), dp) * u &
                        - real(topology%toroidal(trial), dp) * v)
                    if (topology%normal_phase(trial) == phase_cosine) &
                        map(node, trial) = cos(phase)
                    if (topology%normal_phase(trial) == phase_sine) &
                        map(node, trial) = sin(phase)
                end do
            end do
        end do
        info = starwall_fourier_ok
    end subroutine build_starwall_fourier_map

    subroutine add_starwall_fourier_stiffness(nu, nv, nodal, topology, &
            layout, stiffness, info)
        integer, intent(in) :: nu, nv
        real(dp), intent(in) :: nodal(:, :)
        type(trial_space_topology_t), intent(in) :: topology
        type(dynamic_family_layout_t), intent(in) :: layout
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: fourier(:, :), map(:, :)
        integer :: matrix_info

        info = starwall_fourier_invalid_layout
        if (.not. valid_boundary_layout(layout, topology)) return
        if (any(shape(stiffness) /= layout%total_unknowns)) return
        call validate_symmetric_matrix(stiffness, matrix_info)
        if (matrix_info /= starwall_fourier_ok) then
            info = matrix_info
            return
        end if
        info = starwall_fourier_invalid_input
        if (.not. grid_size_is_representable(nu, nv)) return
        if (any(shape(nodal) /= nu * nv)) return
        call validate_symmetric_matrix(nodal, matrix_info)
        if (matrix_info /= starwall_fourier_ok) then
            info = matrix_info
            return
        end if

        call build_starwall_fourier_map(nu, nv, topology, map, info)
        if (info /= starwall_fourier_ok) return
        fourier = matmul(transpose(map), matmul(nodal, map))
        fourier = 0.5_dp * (fourier + transpose(fourier))
        if (.not. all(ieee_is_finite(fourier))) then
            info = starwall_fourier_invalid_input
            return
        end if

        call add_boundary_block(fourier, topology, layout, stiffness, info)
    end subroutine add_starwall_fourier_stiffness

    pure subroutine validate_symmetric_matrix(matrix, info)
        real(dp), intent(in) :: matrix(:, :)
        integer, intent(out) :: info
        real(dp) :: scale

        info = starwall_fourier_invalid_input
        if (size(matrix, 1) /= size(matrix, 2)) return
        if (.not. all(ieee_is_finite(matrix))) return
        scale = max(1.0_dp, maxval(abs(matrix)))
        if (maxval(abs(matrix - transpose(matrix))) &
            > 1024.0_dp * epsilon(1.0_dp) * scale) then
            info = starwall_fourier_nonsymmetric
            return
        end if
        info = starwall_fourier_ok
    end subroutine validate_symmetric_matrix

    subroutine add_boundary_block(fourier, topology, layout, stiffness, info)
        real(dp), intent(in) :: fourier(:, :)
        type(trial_space_topology_t), intent(in) :: topology
        type(dynamic_family_layout_t), intent(in) :: layout
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp) :: updated(size(fourier, 1), size(fourier, 2))
        integer :: a, b, global_a, global_b

        updated = 0.0_dp
        do b = 1, size(fourier, 2)
            if (.not. topology%active(1, b)) cycle
            global_b = normal_global_index(layout, layout%intervals, b)
            do a = 1, size(fourier, 1)
                if (.not. topology%active(1, a)) cycle
                global_a = normal_global_index(layout, layout%intervals, a)
                updated(a, b) = stiffness(global_a, global_b) + fourier(a, b)
            end do
        end do
        info = starwall_fourier_invalid_input
        if (.not. all(ieee_is_finite(updated))) return
        do b = 1, size(fourier, 2)
            if (.not. topology%active(1, b)) cycle
            global_b = normal_global_index(layout, layout%intervals, b)
            do a = 1, size(fourier, 1)
                if (.not. topology%active(1, a)) cycle
                global_a = normal_global_index(layout, layout%intervals, a)
                stiffness(global_a, global_b) = updated(a, b)
            end do
        end do
        info = starwall_fourier_ok
    end subroutine add_boundary_block

    pure function valid_topology(topology) result(valid)
        type(trial_space_topology_t), intent(in) :: topology
        logical :: valid
        integer :: trial, trials

        valid = .false.
        if (.not. allocated(topology%poloidal) &
            .or. .not. allocated(topology%toroidal) &
            .or. .not. allocated(topology%parity) &
            .or. .not. allocated(topology%normal_phase) &
            .or. .not. allocated(topology%active)) return
        trials = size(topology%poloidal)
        if (trials < 1) return
        if (size(topology%toroidal) /= trials &
            .or. size(topology%parity) /= trials &
            .or. size(topology%normal_phase) /= trials) return
        if (any(shape(topology%active) /= [3, trials])) return
        if (any(topology%poloidal < 0)) return
        if (any(topology%normal_phase /= topology%parity)) return
        do trial = 1, trials
            if (.not. valid_phase_kind(topology%normal_phase(trial))) return
        end do
        valid = .true.
    end function valid_topology

    pure function duplicate_normal_trial(topology) result(duplicate)
        type(trial_space_topology_t), intent(in) :: topology
        logical :: duplicate
        integer :: first, second

        duplicate = .false.
        do second = 2, size(topology%poloidal)
            if (.not. topology%active(1, second)) cycle
            do first = 1, second - 1
                if (.not. topology%active(1, first)) cycle
                duplicate = duplicate .or. &
                    topology%poloidal(first) == topology%poloidal(second) &
                    .and. topology%toroidal(first) == topology%toroidal(second) &
                    .and. topology%normal_phase(first) &
                    == topology%normal_phase(second)
            end do
        end do
    end function duplicate_normal_trial

    pure function grid_size_is_representable(nu, nv) result(valid)
        integer, intent(in) :: nu, nv
        logical :: valid
        integer(int64) :: nodes

        valid = .false.
        if (nu < 3 .or. nv < 3) return
        nodes = int(nu, int64) * int(nv, int64)
        valid = nodes <= int(huge(1), int64)
    end function grid_size_is_representable

    pure function grid_resolves_modes(nu, nv, topology) result(resolved)
        integer, intent(in) :: nu, nv
        type(trial_space_topology_t), intent(in) :: topology
        logical :: resolved
        integer :: trial
        integer(int64) :: m, n

        resolved = .false.
        do trial = 1, size(topology%poloidal)
            if (.not. topology%active(1, trial)) cycle
            m = int(topology%poloidal(trial), int64)
            n = int(topology%toroidal(trial), int64)
            if (int(nu, int64) <= 2_int64 * m) return
            if (int(nv, int64) <= 2_int64 * abs(n)) return
        end do
        resolved = .true.
    end function grid_resolves_modes

    pure function valid_boundary_layout(layout, topology) result(valid)
        type(dynamic_family_layout_t), intent(in) :: layout
        type(trial_space_topology_t), intent(in) :: topology
        logical :: valid
        integer :: trial

        valid = .false.
        if (.not. valid_topology(topology)) return
        if (.not. layout%outer_normal_retained) return
        if (layout%trials /= size(topology%poloidal)) return
        if (layout%intervals < 2 .or. layout%total_unknowns < 1) return
        if (.not. allocated(layout%active)) return
        if (any(shape(layout%active) /= shape(topology%active))) return
        if (any(layout%active .neqv. topology%active)) return
        do trial = 1, layout%trials
            if (.not. topology%active(1, trial)) cycle
            if (normal_global_index(layout, layout%intervals, trial) < 1) return
        end do
        valid = .true.
    end function valid_boundary_layout

end module starwall_fourier_coupling
