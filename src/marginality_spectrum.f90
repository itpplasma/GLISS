module marginality_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, &
        ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_problem_ok, &
        compatible_two_component_problem_t
    use dense_spectrum_support, only: dense_spectrum_ok, &
        refine_dense_eigenpair
    use field_profile_identities, only: compute_field_profile_identities, &
        field_profile_identities_ok, field_profile_identity_result_t
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated, &
        symmetric_eigensolver_ok
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    use variable_generalized_solver, only: variable_generalized_inertia, &
        variable_generalized_ok
    implicit none
    private

    integer, parameter, public :: marginality_spectrum_ok = 0
    integer, parameter, public :: marginality_spectrum_invalid = 1
    integer, parameter, public :: marginality_spectrum_compute_error = 2
    real(dp), parameter :: zero_floor = 1.0e-12_dp

    type, public :: marginality_spectrum_result_t
        logical :: has_eigenpair = .false.
        integer :: field_periods = 0
        integer :: mode_count = 0
        integer :: radial_surfaces = 0
        integer :: parity_class = 0
        integer :: degree = 0
        integer :: negative_count = 0
        real(dp) :: lowest_eigenvalue = 0.0_dp
        real(dp) :: certificate = 0.0_dp
        real(dp) :: eigenpair_residual = 0.0_dp
        real(dp) :: force_balance_residual = 0.0_dp
    end type marginality_spectrum_result_t

    public :: compute_marginality_spectrum
    public :: compute_phase_envelope_spectrum

contains

    subroutine compute_marginality_spectrum(equilibrium, mode_m, mode_n, &
            normal_stored_power, parity_class, degree, n_theta, n_zeta, &
            solve_eigenpair, result, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: normal_stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        logical, intent(in) :: solve_eigenpair
        type(marginality_spectrum_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        call compute_spectrum(equilibrium, mode_m, mode_n, &
            normal_stored_power, size(mode_m), parity_class, degree, &
            n_theta, n_zeta, solve_eigenpair, result, info, message)
    end subroutine compute_marginality_spectrum

    subroutine compute_phase_envelope_spectrum(equilibrium, base_m, base_n, &
            envelope_m, envelope_n, parity_class, degree, n_theta, n_zeta, &
            solve_eigenpair, result, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: base_m, base_n
        integer, intent(in) :: envelope_m(:), envelope_n(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        logical, intent(in) :: solve_eigenpair
        type(marginality_spectrum_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer, allocatable :: labeled_m(:), labeled_n(:)
        integer, allocatable :: mode_m(:), mode_n(:)
        real(dp), allocatable :: stored_power(:)
        integer :: allocation_status, labeled_count, mode

        call build_labeled_sidebands(base_m, base_n, &
            equilibrium%field_periods, envelope_m, envelope_n, labeled_m, &
            labeled_n, info, message)
        if (info /= marginality_spectrum_ok) return
        labeled_count = size(labeled_m)
        call unique_mode_table(labeled_m, labeled_n, mode_m, mode_n, info)
        if (info /= marginality_spectrum_ok) then
            message = "phase-envelope mode allocation failed"
            info = marginality_spectrum_compute_error
            return
        end if
        allocate (stored_power(size(mode_m)), stat=allocation_status)
        if (allocation_status /= 0) then
            message = "phase-envelope power allocation failed"
            info = marginality_spectrum_compute_error
            return
        end if
        do mode = 1, size(mode_m)
            stored_power(mode) = 0.0_dp
            if (mode_m(mode) > 0) stored_power(mode) = &
                1.0_dp - 0.5_dp * real(mode_m(mode), dp)
        end do
        call compute_spectrum(equilibrium, mode_m, mode_n, stored_power, &
            labeled_count, parity_class, degree, n_theta, n_zeta, &
            solve_eigenpair, result, info, message)
    end subroutine compute_phase_envelope_spectrum

    subroutine compute_spectrum(equilibrium, mode_m, mode_n, stored_power, &
            reported_mode_count, parity_class, degree, n_theta, n_zeta, &
            solve_eigenpair, result, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:), reported_mode_count
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        logical, intent(in) :: solve_eigenpair
        type(marginality_spectrum_result_t), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(compatible_two_component_problem_t) :: problem
        type(field_profile_identity_result_t) :: identities

        result = marginality_spectrum_result_t()
        call validate_inputs(equilibrium, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, n_zeta, info, message)
        if (info /= marginality_spectrum_ok) return
        call build_compatible_two_component_problem(equilibrium, mode_m, &
            mode_n, stored_power, parity_class, degree, n_theta, n_zeta, &
            problem, info)
        if (info /= compatible_problem_ok) then
            info = marginality_spectrum_compute_error
            message = "compatible FEEC marginality assembly failed"
            return
        end if
        call solve_problem(problem, solve_eigenpair, result, info, message)
        if (info /= marginality_spectrum_ok) return
        call compute_field_profile_identities(equilibrium, n_theta, n_zeta, &
            identities, info)
        if (info /= field_profile_identities_ok) then
            info = marginality_spectrum_compute_error
            message = "equilibrium force-balance diagnostic failed"
            return
        end if
        result%has_eigenpair = solve_eigenpair
        result%field_periods = equilibrium%field_periods
        result%mode_count = reported_mode_count
        result%radial_surfaces = size(equilibrium%s)
        result%parity_class = parity_class
        result%degree = degree
        result%force_balance_residual = &
            maxval(identities%general_force_balance_deviation)
        info = marginality_spectrum_ok
        message = ""
    end subroutine compute_spectrum

    subroutine solve_problem(problem, solve_eigenpair, result, info, message)
        type(compatible_two_component_problem_t), intent(in) :: problem
        logical, intent(in) :: solve_eigenpair
        type(marginality_spectrum_result_t), intent(inout) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(fixed_boundary_solver_controls_t) :: controls
        type(variable_block_tridiagonal_t) :: block_k, block_m
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        real(dp), allocatable :: stiffness(:, :), mass(:, :), vector(:)
        real(dp) :: eigenvalue, residual, resolution
        integer :: allocation_status, local_info

        call pack_problem(problem, block_k, block_m, info)
        if (info /= marginality_spectrum_ok) then
            message = "compatible FEEC matrix packing failed"
            return
        end if
        call variable_generalized_inertia(block_k, block_m, -zero_floor, &
            result%negative_count, local_info)
        if (local_info /= variable_generalized_ok) then
            info = marginality_spectrum_compute_error
            message = "zero-shift compatible FEEC inertia failed"
            return
        end if
        if (.not. solve_eigenpair) then
            result%lowest_eigenvalue = ieee_value(0.0_dp, ieee_quiet_nan)
            result%certificate = ieee_value(0.0_dp, ieee_quiet_nan)
            result%eigenpair_residual = ieee_value(0.0_dp, ieee_quiet_nan)
            info = marginality_spectrum_ok
            message = ""
            return
        end if
        allocate (stiffness(size(problem%stiffness, 1), &
            size(problem%stiffness, 2)), mass(size(problem%mass, 1), &
            size(problem%mass, 2)), stat=allocation_status)
        if (allocation_status /= 0) then
            info = marginality_spectrum_compute_error
            message = "compatible FEEC eigensolve allocation failed"
            return
        end if
        stiffness = problem%stiffness
        mass = problem%mass
        call solve_symmetric_generalized_allocated(stiffness, mass, &
            eigenvalues, eigenvectors, local_info, .true.)
        if (local_info /= symmetric_eigensolver_ok) then
            info = marginality_spectrum_compute_error
            message = "compatible FEEC dense eigensolve failed"
            return
        end if
        call refine_dense_eigenpair(block_k, block_m, controls, eigenvalues, &
            1, eigenvectors(:, 1), eigenvalue, vector, residual, resolution, &
            local_info)
        if (local_info /= dense_spectrum_ok) then
            info = marginality_spectrum_compute_error
            message = "compatible FEEC eigenpair refinement failed"
            return
        end if
        result%lowest_eigenvalue = eigenvalue
        result%eigenpair_residual = residual
        result%certificate = residual + resolution
        info = marginality_spectrum_ok
        message = ""
    end subroutine solve_problem

    subroutine pack_problem(problem, block_k, block_m, info)
        type(compatible_two_component_problem_t), intent(in) :: problem
        type(variable_block_tridiagonal_t), intent(out) :: block_k, block_m
        integer, intent(out) :: info
        integer :: width(1), local_info

        info = marginality_spectrum_compute_error
        width(1) = size(problem%stiffness, 1)
        call pack_variable_blocks(problem%stiffness, width, block_k, &
            local_info)
        if (local_info /= variable_block_ok) return
        call pack_variable_blocks(problem%mass, width, block_m, local_info)
        if (local_info /= variable_block_ok) return
        info = marginality_spectrum_ok
    end subroutine pack_problem

    subroutine validate_inputs(equilibrium, mode_m, mode_n, stored_power, &
            parity_class, degree, n_theta, n_zeta, info, message)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: parity_class, degree, n_theta, n_zeta
        integer, intent(out) :: info
        character(len=*), intent(out) :: message

        info = marginality_spectrum_invalid
        if (.not. valid_mode_table(mode_m, mode_n, stored_power)) then
            message = "mode table is invalid"
        else if (parity_class < 1 .or. parity_class > 2) then
            message = "parity class must be 1 or 2"
        else if (degree < 1 .or. degree > 4) then
            message = "FEEC degree must be between 1 and 4"
        else if (n_theta < 8 .or. n_zeta < 8) then
            message = "angular resolutions must be at least 8"
        else if (equilibrium%field_periods < 1) then
            message = "field periods must be positive"
        else if (size(equilibrium%s) < 2) then
            message = "equilibrium requires at least two radial surfaces"
        else if (angular_grid_aliases(equilibrium, mode_m, mode_n, n_theta, &
                n_zeta)) then
            message = "mode table aliases the angular quadrature"
        else
            info = marginality_spectrum_ok
            message = ""
        end if
    end subroutine validate_inputs

    function valid_mode_table(mode_m, mode_n, stored_power) result(valid)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: stored_power(:)
        logical :: valid
        integer :: first, second

        valid = .false.
        if (size(mode_m) < 1) return
        if (size(mode_n) /= size(mode_m)) return
        if (size(stored_power) /= size(mode_m)) return
        if (any(mode_m < 0)) return
        if (.not. all(ieee_is_finite(stored_power))) return
        do first = 1, size(mode_m)
            if (mode_m(first) == 0 .and. mode_n(first) < 0) return
            do second = 1, first - 1
                if (mode_m(first) == mode_m(second) &
                    .and. mode_n(first) == mode_n(second)) return
            end do
        end do
        valid = .true.
    end function valid_mode_table

    subroutine unique_mode_table(source_m, source_n, mode_m, mode_n, info)
        integer, intent(in) :: source_m(:), source_n(:)
        integer, allocatable, intent(out) :: mode_m(:), mode_n(:)
        integer, intent(out) :: info
        integer, allocatable :: work_m(:), work_n(:)
        integer :: allocation_status, count, output, source, previous
        logical :: duplicate

        info = marginality_spectrum_compute_error
        allocate (work_m(size(source_m)), work_n(size(source_n)), &
            stat=allocation_status)
        if (allocation_status /= 0) return
        count = 0
        do source = 1, size(source_m)
            duplicate = .false.
            do previous = 1, count
                if (source_m(source) == work_m(previous) &
                    .and. source_n(source) == work_n(previous)) then
                    duplicate = .true.
                    exit
                end if
            end do
            if (duplicate) cycle
            count = count + 1
            work_m(count) = source_m(source)
            work_n(count) = source_n(source)
        end do
        allocate (mode_m(count), mode_n(count), stat=allocation_status)
        if (allocation_status /= 0) return
        do output = 1, count
            mode_m(output) = work_m(output)
            mode_n(output) = work_n(output)
        end do
        info = marginality_spectrum_ok
    end subroutine unique_mode_table

    subroutine build_labeled_sidebands(base_m, base_n, field_periods, &
            envelope_m, envelope_n, mode_m, mode_n, info, message)
        integer, intent(in) :: base_m, base_n, field_periods
        integer, intent(in) :: envelope_m(:), envelope_n(:)
        integer, allocatable, intent(out) :: mode_m(:), mode_n(:)
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        integer(int64) :: first_m, first_n, second_m, second_n
        integer(int64) :: base_m64, base_n64, delta_m, delta_n, periods64
        integer :: allocation_status, envelope, first, labeled, second

        info = marginality_spectrum_invalid
        message = "phase-envelope table is invalid"
        if (field_periods < 1 .or. base_m < 0) return
        if (base_m == 0 .and. base_n < 0) return
        if (size(envelope_m) < 1) return
        if (size(envelope_n) /= size(envelope_m)) return
        if (envelope_m(1) /= 0 .or. envelope_n(1) /= 0) return
        if (size(envelope_m) > (huge(labeled) + 1_int64) / 2_int64) return
        do first = 1, size(envelope_m)
            if (envelope_m(first) < 0) return
            do second = 1, first - 1
                if (envelope_m(first) == envelope_m(second) &
                    .and. envelope_n(first) == envelope_n(second)) return
            end do
        end do
        allocate (mode_m(2 * size(envelope_m) - 1), &
            mode_n(2 * size(envelope_n) - 1), stat=allocation_status)
        if (allocation_status /= 0) then
            info = marginality_spectrum_compute_error
            message = "phase-envelope mode allocation failed"
            return
        end if
        mode_m(1) = base_m
        mode_n(1) = base_n
        base_m64 = int(base_m, int64)
        base_n64 = int(base_n, int64)
        periods64 = int(field_periods, int64)
        labeled = 1
        do envelope = 2, size(envelope_m)
            delta_m = int(envelope_m(envelope), int64)
            delta_n = periods64 * int(envelope_n(envelope), int64)
            first_m = base_m64 + delta_m
            first_n = base_n64 + delta_n
            second_m = base_m64 - delta_m
            second_n = base_n64 - delta_n
            call canonicalize_mode(first_m, first_n)
            call canonicalize_mode(second_m, second_n)
            if (.not. mode_fits_default_integer(first_m, first_n)) then
                message = "phase-envelope sideband exceeds integer range"
                return
            end if
            if (.not. mode_fits_default_integer(second_m, second_n)) then
                message = "phase-envelope sideband exceeds integer range"
                return
            end if
            mode_m(labeled + 1) = int(first_m)
            mode_n(labeled + 1) = int(first_n)
            mode_m(labeled + 2) = int(second_m)
            mode_n(labeled + 2) = int(second_n)
            labeled = labeled + 2
        end do
        info = marginality_spectrum_ok
        message = ""
    end subroutine build_labeled_sidebands

    pure subroutine canonicalize_mode(mode_m, mode_n)
        integer(int64), intent(inout) :: mode_m, mode_n

        if (mode_m < 0_int64) then
            mode_m = -mode_m
            mode_n = -mode_n
        else if (mode_m == 0_int64 .and. mode_n < 0_int64) then
            mode_n = -mode_n
        end if
    end subroutine canonicalize_mode

    pure function mode_fits_default_integer(mode_m, mode_n) result(fits)
        integer(int64), intent(in) :: mode_m, mode_n
        logical :: fits
        integer(int64), parameter :: maximum = int(huge(0), int64)
        integer(int64), parameter :: minimum = -maximum - 1_int64

        fits = mode_m >= 0_int64 .and. mode_m <= maximum &
            .and. mode_n >= minimum .and. mode_n <= maximum
    end function mode_fits_default_integer

    function angular_grid_aliases(equilibrium, mode_m, mode_n, n_theta, &
            n_zeta) result(aliases)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: mode_m(:), mode_n(:), n_theta, n_zeta
        logical :: aliases
        integer(int64) :: poloidal_bandwidth, toroidal_bandwidth

        aliases = .true.
        if (.not. allocated(equilibrium%poloidal_modes)) return
        if (.not. allocated(equilibrium%toroidal_modes)) return
        if (size(equilibrium%poloidal_modes) < 1) return
        if (size(equilibrium%toroidal_modes) < 1) return
        poloidal_bandwidth = 2_int64 * int(maxval(mode_m), int64) &
            + int(maxval(abs(equilibrium%poloidal_modes)), int64)
        toroidal_bandwidth = 2_int64 &
            * maxval(abs(int(mode_n, int64))) &
            + maxval(abs(int(equilibrium%toroidal_modes, int64)))
        aliases = poloidal_bandwidth >= int(n_theta, int64) &
            .or. toroidal_bandwidth >= int(n_zeta, int64)
    end function angular_grid_aliases

end module marginality_spectrum
