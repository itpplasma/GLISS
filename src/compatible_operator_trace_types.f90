module compatible_operator_trace_types
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    type, public :: compatible_radial_point_trace_t
        real(dp) :: coordinate = 0.0_dp
        real(dp) :: weight = 0.0_dp
        logical :: term_mask(4) = .false.
        logical :: assembles_mass = .false.
        integer, allocatable :: map(:)
        real(dp), allocatable :: fields(:, :, :), drive(:, :)
        real(dp), allocatable :: h1(:, :), dh1(:, :), l2(:, :)
        real(dp), allocatable :: stiffness_terms(:, :, :), mass(:, :)
    end type compatible_radial_point_trace_t

    type, public :: compatible_cell_trace_t
        integer :: cell = 0
        type(compatible_radial_point_trace_t), allocatable :: points(:)
    end type compatible_cell_trace_t

    public :: build_trace_radial_mass

contains

    subroutine build_trace_radial_mass(h1, l2, mass)
        real(dp), intent(in) :: h1(:, :), l2(:, :)
        real(dp), allocatable, intent(out) :: mass(:, :)
        real(dp) :: basis(size(h1, 1) * size(h1, 2) &
            + size(l2, 1) * size(l2, 2))
        integer :: a, b, basis_index, h1_columns, trial, trials

        trials = size(h1, 2)
        h1_columns = size(h1, 1) * trials
        do basis_index = 1, size(h1, 1)
            do trial = 1, trials
                basis((basis_index - 1) * trials + trial) = &
                    h1(basis_index, trial)
            end do
        end do
        do basis_index = 1, size(l2, 1)
            do trial = 1, trials
                basis(h1_columns + (basis_index - 1) * trials + trial) = &
                    l2(basis_index, trial)
            end do
        end do
        allocate (mass(size(basis), size(basis)), source=0.0_dp)
        do b = 1, size(basis)
            do a = 1, size(basis)
                if (modulo(a - 1, trials) /= modulo(b - 1, trials)) cycle
                if ((a <= h1_columns) .neqv. (b <= h1_columns)) cycle
                mass(a, b) = basis(a) * basis(b)
            end do
        end do
    end subroutine build_trace_radial_mass

end module compatible_operator_trace_types
