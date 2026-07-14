module cartesian_harmonic_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use axis_regular_harmonic_spline, only: &
        axis_regular_harmonic_allocation_error, axis_regular_harmonic_field_t, &
        axis_regular_harmonic_ok, evaluate_axis_regular_harmonics, &
        fit_axis_regular_harmonics
    use gvec_cas3d_reconstruction, only: reconstruct_harmonic_grid, &
        reconstruction_ok
    use gvec_cas3d_types, only: harmonic_pair_t
    use radial_cubic_spline, only: radial_cubic_spline_grid_t
    implicit none
    private

    integer, parameter, public :: cartesian_harmonic_ok = 0
    integer, parameter, public :: cartesian_harmonic_invalid = -1
    integer, parameter, public :: cartesian_harmonic_allocation_error = -2

    type, public :: cartesian_harmonic_spline_t
        integer, allocatable :: poloidal_modes(:), toroidal_modes(:)
        type(axis_regular_harmonic_field_t) :: coefficients
    end type cartesian_harmonic_spline_t

    type, public :: cartesian_jet_grid_t
        real(dp), allocatable :: value(:, :, :), radial(:, :, :)
        real(dp), allocatable :: poloidal(:, :, :), toroidal(:, :, :)
        real(dp), allocatable :: radial_radial(:, :, :)
        real(dp), allocatable :: radial_poloidal(:, :, :)
        real(dp), allocatable :: radial_toroidal(:, :, :)
        real(dp), allocatable :: poloidal_poloidal(:, :, :)
        real(dp), allocatable :: poloidal_toroidal(:, :, :)
        real(dp), allocatable :: toroidal_toroidal(:, :, :)
    end type cartesian_jet_grid_t

    public :: evaluate_cartesian_harmonic_spline
    public :: fit_cartesian_harmonic_spline

contains

    subroutine fit_cartesian_harmonic_spline(grid, poloidal_modes, &
            toroidal_modes, x, y, z, spline, info)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        type(harmonic_pair_t), intent(in) :: x, y, z
        type(cartesian_harmonic_spline_t), intent(out) :: spline
        integer, intent(out) :: info
        real(dp), allocatable :: samples(:, :)
        integer, allocatable :: coefficient_modes(:)
        integer :: allocation_status, coefficient_count, fit_info

        info = cartesian_harmonic_invalid
        if (.not. position_pairs_are_valid(grid, poloidal_modes, &
            toroidal_modes, x, y, z)) return
        if (.not. safe_coefficient_count(size(poloidal_modes), &
            size(toroidal_modes), coefficient_count)) return
        allocate (spline%poloidal_modes(size(poloidal_modes)), &
            spline%toroidal_modes(size(toroidal_modes)), &
            samples(size(grid%nodes), coefficient_count), &
            coefficient_modes(coefficient_count), stat=allocation_status)
        if (allocation_status /= 0) then
            info = cartesian_harmonic_allocation_error
            return
        end if
        spline%poloidal_modes = poloidal_modes
        spline%toroidal_modes = toroidal_modes
        call flatten_position_pairs(x, y, z, poloidal_modes, samples, &
            coefficient_modes)
        call fit_axis_regular_harmonics(grid, coefficient_modes, samples, &
            spline%coefficients, fit_info)
        if (fit_info == axis_regular_harmonic_allocation_error) then
            info = cartesian_harmonic_allocation_error
            return
        end if
        if (fit_info /= axis_regular_harmonic_ok) return
        info = cartesian_harmonic_ok
    end subroutine fit_cartesian_harmonic_spline

    subroutine evaluate_cartesian_harmonic_spline(grid, spline, coordinate, &
            theta, zeta_period, jet, info)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        type(cartesian_harmonic_spline_t), intent(in) :: spline
        real(dp), intent(in) :: coordinate, theta(:), zeta_period(:)
        type(cartesian_jet_grid_t), intent(out) :: jet
        integer, intent(out) :: info
        real(dp), allocatable :: coefficients(:), radial(:), second(:)
        integer :: allocation_status, coefficient_count, component, fit_info

        jet = cartesian_jet_grid_t()
        info = cartesian_harmonic_invalid
        if (.not. spline_is_valid(spline)) return
        if (size(theta) < 1 .or. size(zeta_period) < 1) return
        if (.not. all(ieee_is_finite(theta)) &
            .or. .not. all(ieee_is_finite(zeta_period))) return
        coefficient_count = size(spline%coefficients%poloidal_modes)
        allocate (coefficients(coefficient_count), radial(coefficient_count), &
            second(coefficient_count), stat=allocation_status)
        if (allocation_status /= 0) then
            info = cartesian_harmonic_allocation_error
            return
        end if
        call evaluate_axis_regular_harmonics(grid, spline%coefficients, &
            coordinate, coefficients, radial, second, fit_info)
        if (fit_info /= axis_regular_harmonic_ok) return
        call allocate_jet(jet, size(theta), size(zeta_period), &
            allocation_status)
        if (allocation_status /= 0) then
            jet = cartesian_jet_grid_t()
            info = cartesian_harmonic_allocation_error
            return
        end if
        do component = 1, 3
            call reconstruct_component(spline, component, coefficients, &
                radial, second, theta, zeta_period, jet, fit_info)
            if (fit_info /= reconstruction_ok) then
                jet = cartesian_jet_grid_t()
                return
            end if
        end do
        if (.not. jet_is_finite(jet)) then
            jet = cartesian_jet_grid_t()
            return
        end if
        info = cartesian_harmonic_ok
    end subroutine evaluate_cartesian_harmonic_spline

    subroutine reconstruct_component(spline, component, coefficients, &
            radial, second, theta, zeta_period, jet, info)
        type(cartesian_harmonic_spline_t), intent(in) :: spline
        integer, intent(in) :: component
        real(dp), intent(in) :: coefficients(:), radial(:), second(:)
        real(dp), intent(in) :: theta(:), zeta_period(:)
        type(cartesian_jet_grid_t), intent(inout) :: jet
        integer, intent(out) :: info
        type(harmonic_pair_t) :: value_pair, radial_pair, second_pair
        real(dp), allocatable :: value(:, :), dt(:, :), dz(:, :)
        real(dp), allocatable :: dtt(:, :), dtz(:, :), dzz(:, :)
        real(dp), allocatable :: ds(:, :), dst(:, :), dsz(:, :)
        real(dp), allocatable :: dss(:, :), discard_t(:, :), discard_z(:, :)

        call unflatten_component(spline, component, coefficients, value_pair)
        call unflatten_component(spline, component, radial, radial_pair)
        call unflatten_component(spline, component, second, second_pair)
        call reconstruct_harmonic_grid(value_pair, 1, &
            spline%poloidal_modes, spline%toroidal_modes, theta, zeta_period, &
            value, dt, dz, info, dtt, dtz, dzz)
        if (info /= reconstruction_ok) return
        call reconstruct_harmonic_grid(radial_pair, 1, &
            spline%poloidal_modes, spline%toroidal_modes, theta, zeta_period, &
            ds, dst, dsz, info)
        if (info /= reconstruction_ok) return
        call reconstruct_harmonic_grid(second_pair, 1, &
            spline%poloidal_modes, spline%toroidal_modes, theta, zeta_period, &
            dss, discard_t, discard_z, info)
        if (info /= reconstruction_ok) return
        jet%value(:, :, component) = value
        jet%radial(:, :, component) = ds
        jet%poloidal(:, :, component) = dt
        jet%toroidal(:, :, component) = dz
        jet%radial_radial(:, :, component) = dss
        jet%radial_poloidal(:, :, component) = dst
        jet%radial_toroidal(:, :, component) = dsz
        jet%poloidal_poloidal(:, :, component) = dtt
        jet%poloidal_toroidal(:, :, component) = dtz
        jet%toroidal_toroidal(:, :, component) = dzz
    end subroutine reconstruct_component

    subroutine flatten_position_pairs(x, y, z, poloidal_modes, samples, modes)
        type(harmonic_pair_t), intent(in) :: x, y, z
        integer, intent(in) :: poloidal_modes(:)
        real(dp), intent(out) :: samples(:, :)
        integer, intent(out) :: modes(:)
        integer :: column

        column = 0
        call flatten_pair(x, poloidal_modes, samples, modes, column)
        call flatten_pair(y, poloidal_modes, samples, modes, column)
        call flatten_pair(z, poloidal_modes, samples, modes, column)
    end subroutine flatten_position_pairs

    subroutine flatten_pair(pair, poloidal_modes, samples, modes, column)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: poloidal_modes(:)
        real(dp), intent(inout) :: samples(:, :)
        integer, intent(inout) :: modes(:)
        integer, intent(inout) :: column
        integer :: m, n

        do n = 1, size(pair%cosine, 3)
            do m = 1, size(poloidal_modes)
                column = column + 1
                samples(:, column) = pair%cosine(:, m, n)
                modes(column) = poloidal_modes(m)
                column = column + 1
                samples(:, column) = pair%sine(:, m, n)
                modes(column) = poloidal_modes(m)
            end do
        end do
    end subroutine flatten_pair

    subroutine unflatten_component(spline, component, values, pair)
        type(cartesian_harmonic_spline_t), intent(in) :: spline
        integer, intent(in) :: component
        real(dp), intent(in) :: values(:)
        type(harmonic_pair_t), intent(out) :: pair
        integer :: column, m, n, component_columns

        allocate (pair%cosine(1, size(spline%poloidal_modes), &
            size(spline%toroidal_modes)), pair%sine(1, &
            size(spline%poloidal_modes), size(spline%toroidal_modes)))
        component_columns = 2 * size(spline%poloidal_modes) &
            * size(spline%toroidal_modes)
        column = (component - 1) * component_columns
        do n = 1, size(spline%toroidal_modes)
            do m = 1, size(spline%poloidal_modes)
                column = column + 1
                pair%cosine(1, m, n) = values(column)
                column = column + 1
                pair%sine(1, m, n) = values(column)
            end do
        end do
    end subroutine unflatten_component

    subroutine allocate_jet(jet, n_theta, n_zeta, allocation_status)
        type(cartesian_jet_grid_t), intent(out) :: jet
        integer, intent(in) :: n_theta, n_zeta
        integer, intent(out) :: allocation_status

        allocate (jet%value(n_theta, n_zeta, 3), &
            jet%radial(n_theta, n_zeta, 3), &
            jet%poloidal(n_theta, n_zeta, 3), &
            jet%toroidal(n_theta, n_zeta, 3), &
            jet%radial_radial(n_theta, n_zeta, 3), &
            jet%radial_poloidal(n_theta, n_zeta, 3), &
            jet%radial_toroidal(n_theta, n_zeta, 3), &
            jet%poloidal_poloidal(n_theta, n_zeta, 3), &
            jet%poloidal_toroidal(n_theta, n_zeta, 3), &
            jet%toroidal_toroidal(n_theta, n_zeta, 3), &
            stat=allocation_status)
    end subroutine allocate_jet

    function position_pairs_are_valid(grid, poloidal_modes, toroidal_modes, &
            x, y, z) result(valid)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        type(harmonic_pair_t), intent(in) :: x, y, z
        logical :: valid
        integer :: expected(3)

        valid = allocated(grid%nodes) .and. size(poloidal_modes) > 0 &
            .and. size(toroidal_modes) > 0
        if (.not. valid) return
        expected = [size(grid%nodes), size(poloidal_modes), &
            size(toroidal_modes)]
        valid = pair_has_shape(x, expected) .and. pair_has_shape(y, expected) &
            .and. pair_has_shape(z, expected)
    end function position_pairs_are_valid

    function pair_has_shape(pair, expected) result(valid)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: expected(3)
        logical :: valid

        valid = allocated(pair%cosine) .and. allocated(pair%sine)
        if (.not. valid) return
        valid = all(shape(pair%cosine) == expected) &
            .and. all(shape(pair%sine) == expected) &
            .and. all(ieee_is_finite(pair%cosine)) &
            .and. all(ieee_is_finite(pair%sine))
    end function pair_has_shape

    function safe_coefficient_count(m_count, n_count, count) result(valid)
        integer, intent(in) :: m_count, n_count
        integer, intent(out) :: count
        logical :: valid
        integer :: maximum_count

        count = 0
        maximum_count = huge(count)
        valid = m_count > 0 .and. n_count > 0
        if (.not. valid) return
        if (m_count > maximum_count / 6) then
            valid = .false.
            return
        end if
        count = 6 * m_count
        if (n_count > maximum_count / count) then
            count = 0
            valid = .false.
            return
        end if
        count = count * n_count
    end function safe_coefficient_count

    function spline_is_valid(spline) result(valid)
        type(cartesian_harmonic_spline_t), intent(in) :: spline
        logical :: valid
        integer :: expected

        valid = allocated(spline%poloidal_modes) &
            .and. allocated(spline%toroidal_modes) &
            .and. allocated(spline%coefficients%poloidal_modes)
        if (.not. valid) return
        if (.not. safe_coefficient_count(size(spline%poloidal_modes), &
            size(spline%toroidal_modes), expected)) then
            valid = .false.
            return
        end if
        valid = size(spline%coefficients%poloidal_modes) == expected
    end function spline_is_valid

    function jet_is_finite(jet) result(finite)
        type(cartesian_jet_grid_t), intent(in) :: jet
        logical :: finite

        finite = all(ieee_is_finite(jet%value)) &
            .and. all(ieee_is_finite(jet%radial)) &
            .and. all(ieee_is_finite(jet%poloidal)) &
            .and. all(ieee_is_finite(jet%toroidal)) &
            .and. all(ieee_is_finite(jet%radial_radial)) &
            .and. all(ieee_is_finite(jet%radial_poloidal)) &
            .and. all(ieee_is_finite(jet%radial_toroidal)) &
            .and. all(ieee_is_finite(jet%poloidal_poloidal)) &
            .and. all(ieee_is_finite(jet%poloidal_toroidal)) &
            .and. all(ieee_is_finite(jet%toroidal_toroidal))
    end function jet_is_finite

end module cartesian_harmonic_spline
