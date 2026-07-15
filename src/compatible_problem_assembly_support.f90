module compatible_problem_assembly_support
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: compatible_support_ok = 0
    integer, parameter, public :: compatible_support_invalid = -1
    integer, parameter, public :: compatible_support_allocation = -2

    public :: apply_stored_power
    public :: apply_l2_stored_power
    public :: build_active_indices
    public :: build_uniform_breaks
    public :: evaluate_generalized_eigenpair
    public :: mode_table_is_unique
    public :: quadratic_form
    public :: replicate_indexed_values
    public :: scale_matrix
    public :: scale_tensor
    public :: scatter_matrix
    public :: sum_tensor
    public :: symmetrize_matrix
    public :: symmetrize_tensor

contains

    subroutine apply_l2_stored_power(coordinate, stored_power, l2, indices, &
            values, info)
        real(dp), intent(in) :: coordinate, stored_power(:), l2(:)
        integer, intent(in) :: indices(:)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: info
        real(dp) :: scale
        integer :: basis, trial

        info = compatible_support_invalid
        if (coordinate <= 0.0_dp) return
        if (size(stored_power) /= size(values, 2)) return
        if (size(indices) /= size(values, 1)) return
        if (any(indices < 1) .or. any(indices > size(l2))) return
        do trial = 1, size(stored_power)
            scale = coordinate**(-stored_power(trial))
            if (.not. ieee_is_finite(scale)) return
            do basis = 1, size(indices)
                values(basis, trial) = scale * l2(indices(basis))
            end do
        end do
        if (.not. all(ieee_is_finite(values))) return
        info = compatible_support_ok
    end subroutine apply_l2_stored_power

    subroutine apply_stored_power(coordinate, stored_power, h1, dh1, indices, &
            values, derivatives, info)
        real(dp), intent(in) :: coordinate, stored_power(:), h1(:), dh1(:)
        integer, intent(in) :: indices(:)
        real(dp), intent(out) :: values(:, :), derivatives(:, :)
        integer, intent(out) :: info
        real(dp) :: scale
        integer :: trial

        info = compatible_support_invalid
        if (coordinate <= 0.0_dp) return
        if (size(stored_power) /= size(values, 2)) return
        if (size(indices) /= size(values, 1)) return
        if (any(shape(derivatives) /= shape(values))) return
        if (size(dh1) /= size(h1)) return
        if (any(indices < 1) .or. any(indices > size(h1))) return
        do trial = 1, size(stored_power)
            scale = coordinate**(-stored_power(trial))
            if (.not. ieee_is_finite(scale)) return
            call set_stored_power_column(coordinate, stored_power(trial), &
                scale, h1, dh1, indices, values(:, trial), &
                derivatives(:, trial))
        end do
        if (.not. all(ieee_is_finite(values))) return
        if (.not. all(ieee_is_finite(derivatives))) return
        info = compatible_support_ok
    end subroutine apply_stored_power

    subroutine set_stored_power_column(coordinate, power, scale, h1, dh1, &
            indices, values, derivatives)
        real(dp), intent(in) :: coordinate, power, scale, h1(:), dh1(:)
        integer, intent(in) :: indices(:)
        real(dp), intent(out) :: values(:), derivatives(:)
        integer :: basis, source

        do basis = 1, size(indices)
            source = indices(basis)
            values(basis) = scale * h1(source)
            derivatives(basis) = scale &
                * (dh1(source) - power * h1(source) / coordinate)
        end do
    end subroutine set_stored_power_column

    subroutine build_active_indices(values, indices, info, derivatives)
        real(dp), intent(in) :: values(:)
        integer, allocatable, intent(out) :: indices(:)
        integer, intent(out) :: info
        real(dp), optional, intent(in) :: derivatives(:)
        integer :: allocation_status, count, index

        info = compatible_support_invalid
        if (present(derivatives)) then
            if (size(derivatives) /= size(values)) return
        end if
        count = 0
        do index = 1, size(values)
            if (is_active(index)) count = count + 1
        end do
        allocate (indices(count), stat=allocation_status)
        if (allocation_status /= 0) then
            info = compatible_support_allocation
            return
        end if
        count = 0
        do index = 1, size(values)
            if (.not. is_active(index)) cycle
            count = count + 1
            indices(count) = index
        end do
        info = compatible_support_ok
    contains
        function is_active(index) result(active)
            integer, intent(in) :: index
            logical :: active

            active = values(index) /= 0.0_dp
            if (present(derivatives)) then
                if (derivatives(index) /= 0.0_dp) active = .true.
            end if
        end function is_active
    end subroutine build_active_indices

    subroutine replicate_indexed_values(values, indices, repeated, info)
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: indices(:)
        real(dp), intent(out) :: repeated(:, :)
        integer, intent(out) :: info
        integer :: basis, trial

        info = compatible_support_invalid
        if (size(repeated, 1) /= size(indices)) return
        if (any(indices < 1) .or. any(indices > size(values))) return
        do trial = 1, size(repeated, 2)
            do basis = 1, size(indices)
                repeated(basis, trial) = values(indices(basis))
            end do
        end do
        info = compatible_support_ok
    end subroutine replicate_indexed_values

    subroutine build_uniform_breaks(intervals, breaks, info)
        integer, intent(in) :: intervals
        real(dp), intent(out) :: breaks(:)
        integer, intent(out) :: info
        integer :: boundary

        info = compatible_support_invalid
        if (intervals < 1 .or. size(breaks) /= intervals + 1) return
        do boundary = 0, intervals
            breaks(boundary + 1) = real(boundary, dp) / real(intervals, dp)
        end do
        info = compatible_support_ok
    end subroutine build_uniform_breaks

    subroutine evaluate_generalized_eigenpair(stiffness, mass, vector, &
            eigenvalue, kinetic, potential, residual, info)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :), vector(:)
        real(dp), intent(in) :: eigenvalue
        real(dp), intent(out) :: kinetic, potential, residual
        integer, intent(out) :: info
        real(dp) :: stiffness_value, mass_value, difference
        real(dp) :: stiffness_norm, mass_norm, difference_norm
        integer :: column, row

        info = compatible_support_invalid
        kinetic = 0.0_dp
        potential = 0.0_dp
        residual = huge(1.0_dp)
        if (.not. compatible_dimensions(stiffness, mass, vector)) return
        if (.not. ieee_is_finite(eigenvalue)) return
        stiffness_norm = 0.0_dp
        mass_norm = 0.0_dp
        difference_norm = 0.0_dp
        do row = 1, size(vector)
            stiffness_value = 0.0_dp
            mass_value = 0.0_dp
            do column = 1, size(vector)
                stiffness_value = stiffness_value &
                    + stiffness(row, column) * vector(column)
                mass_value = mass_value + mass(row, column) * vector(column)
            end do
            kinetic = kinetic + vector(row) * mass_value
            potential = potential + vector(row) * stiffness_value
            stiffness_norm = stiffness_norm + stiffness_value**2
            mass_norm = mass_norm + mass_value**2
            difference = stiffness_value - eigenvalue * mass_value
            difference_norm = difference_norm + difference**2
        end do
        residual = sqrt(difference_norm) / max(tiny(1.0_dp), &
            sqrt(stiffness_norm) + abs(eigenvalue) * sqrt(mass_norm))
        info = compatible_support_ok
    end subroutine evaluate_generalized_eigenpair

    subroutine quadratic_form(matrix, vector, value, info)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp), intent(out) :: value
        integer, intent(out) :: info
        integer :: column, row

        info = compatible_support_invalid
        value = 0.0_dp
        if (size(matrix, 1) /= size(vector)) return
        if (size(matrix, 2) /= size(vector)) return
        if (.not. all(ieee_is_finite(matrix))) return
        if (.not. all(ieee_is_finite(vector))) return
        do column = 1, size(vector)
            do row = 1, size(vector)
                value = value + vector(row) * matrix(row, column) &
                    * vector(column)
            end do
        end do
        info = compatible_support_ok
    end subroutine quadratic_form

    function compatible_dimensions(stiffness, mass, vector) result(valid)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :), vector(:)
        logical :: valid
        integer :: dimension

        dimension = size(vector)
        valid = dimension >= 1
        valid = valid .and. all(shape(stiffness) == [dimension, dimension])
        valid = valid .and. all(shape(mass) == [dimension, dimension])
        valid = valid .and. all(ieee_is_finite(stiffness))
        valid = valid .and. all(ieee_is_finite(mass))
        valid = valid .and. all(ieee_is_finite(vector))
    end function compatible_dimensions

    subroutine scatter_matrix(map, local, scale, global)
        integer, intent(in) :: map(:)
        real(dp), intent(in) :: local(:, :), scale
        real(dp), intent(inout) :: global(:, :)
        integer :: a, b

        do b = 1, size(map)
            if (map(b) == 0) cycle
            do a = 1, size(map)
                if (map(a) == 0) cycle
                global(map(a), map(b)) = global(map(a), map(b)) &
                    + scale * local(a, b)
            end do
        end do
    end subroutine scatter_matrix

    subroutine sum_tensor(tensor, matrix)
        real(dp), intent(in) :: tensor(:, :, :)
        real(dp), intent(out) :: matrix(:, :)
        integer :: a, b, term

        matrix = 0.0_dp
        do term = 1, size(tensor, 3)
            do b = 1, size(matrix, 2)
                do a = 1, size(matrix, 1)
                    matrix(a, b) = matrix(a, b) + tensor(a, b, term)
                end do
            end do
        end do
    end subroutine sum_tensor

    subroutine symmetrize_matrix(matrix)
        real(dp), intent(inout) :: matrix(:, :)
        real(dp) :: average
        integer :: a, b

        do b = 1, size(matrix, 2)
            do a = 1, b - 1
                average = 0.5_dp * (matrix(a, b) + matrix(b, a))
                matrix(a, b) = average
                matrix(b, a) = average
            end do
        end do
    end subroutine symmetrize_matrix

    subroutine symmetrize_tensor(tensor)
        real(dp), intent(inout) :: tensor(:, :, :)
        integer :: term

        do term = 1, size(tensor, 3)
            call symmetrize_matrix(tensor(:, :, term))
        end do
    end subroutine symmetrize_tensor

    subroutine scale_matrix(matrix, scale)
        real(dp), intent(inout) :: matrix(:, :)
        real(dp), intent(in) :: scale
        integer :: a, b

        do b = 1, size(matrix, 2)
            do a = 1, size(matrix, 1)
                matrix(a, b) = scale * matrix(a, b)
            end do
        end do
    end subroutine scale_matrix

    subroutine scale_tensor(tensor, scale)
        real(dp), intent(inout) :: tensor(:, :, :)
        real(dp), intent(in) :: scale
        integer :: a, b, term

        do term = 1, size(tensor, 3)
            do b = 1, size(tensor, 2)
                do a = 1, size(tensor, 1)
                    tensor(a, b, term) = scale * tensor(a, b, term)
                end do
            end do
        end do
    end subroutine scale_tensor

    pure function mode_table_is_unique(mode_m, mode_n) result(unique)
        integer, intent(in) :: mode_m(:), mode_n(:)
        logical :: unique
        integer :: first, second

        unique = .false.
        if (size(mode_n) /= size(mode_m)) return
        do first = 1, size(mode_m)
            do second = 1, first - 1
                if (mode_m(first) == mode_m(second) &
                    .and. mode_n(first) == mode_n(second)) return
            end do
        end do
        unique = .true.
    end function mode_table_is_unique

end module compatible_problem_assembly_support
