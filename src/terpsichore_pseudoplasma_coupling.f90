module terpsichore_pseudoplasma_coupling
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamic_family_layout, only: dynamic_family_layout_t, &
        normal_global_index
    use terpsichore_matrix_fixture, only: terpsichore_matrix_fixture_t
    use terpsichore_pseudoplasma_fixture, only: &
        terpsichore_pseudoplasma_fixture_is_valid, &
        terpsichore_pseudoplasma_fixture_t
    use terpsichore_pseudoplasma_stiffness, only: &
        assemble_terpsichore_pseudoplasma_stiffness, &
        pseudoplasma_stiffness_ok
    use vacuum_schur, only: eliminate_vacuum, vacuum_schur_ok
    implicit none
    private

    integer, parameter, public :: pseudoplasma_coupling_ok = 0
    integer, parameter, public :: pseudoplasma_coupling_invalid = -1
    integer, parameter, public :: pseudoplasma_coupling_not_spd = -2

    public :: add_terpsichore_pseudoplasma_schur

contains

    subroutine add_terpsichore_pseudoplasma_schur(plasma, vacuum, layout, &
            stiffness, effective, response, info)
        type(terpsichore_matrix_fixture_t), intent(in) :: plasma
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: vacuum
        type(dynamic_family_layout_t), intent(in) :: layout
        real(dp), intent(inout) :: stiffness(:, :)
        real(dp), allocatable, intent(out) :: effective(:, :), response(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: vacuum_stiffness(:, :)
        integer :: mode, modes, schur_info, target

        info = pseudoplasma_coupling_invalid
        if (.not. compatible_inputs(plasma, vacuum, layout, stiffness)) return
        modes = vacuum%modes
        call assemble_terpsichore_pseudoplasma_stiffness(vacuum, &
            vacuum_stiffness, schur_info)
        if (schur_info /= pseudoplasma_stiffness_ok) return
        call eliminate_vacuum(vacuum_stiffness(:modes, :modes), &
            vacuum_stiffness(modes + 1:, modes + 1:), &
            vacuum_stiffness(:modes, modes + 1:), effective, response, &
            schur_info)
        if (schur_info /= vacuum_schur_ok) then
            info = pseudoplasma_coupling_not_spd
            return
        end if
        do mode = 1, modes
            target = normal_global_index(layout, layout%intervals, mode)
            call add_interface_row(layout, mode, target, effective, stiffness)
        end do
        info = pseudoplasma_coupling_ok
    end subroutine add_terpsichore_pseudoplasma_schur

    subroutine add_interface_row(layout, mode, target, effective, stiffness)
        type(dynamic_family_layout_t), intent(in) :: layout
        integer, intent(in) :: mode, target
        real(dp), intent(in) :: effective(:, :)
        real(dp), intent(inout) :: stiffness(:, :)
        integer :: other, other_target

        do other = 1, layout%trials
            other_target = normal_global_index(layout, layout%intervals, other)
            stiffness(target, other_target) = stiffness(target, other_target) &
                + effective(mode, other)
        end do
    end subroutine add_interface_row

    pure function compatible_inputs(plasma, vacuum, layout, stiffness) &
            result(valid)
        type(terpsichore_matrix_fixture_t), intent(in) :: plasma
        type(terpsichore_pseudoplasma_fixture_t), intent(in) :: vacuum
        type(dynamic_family_layout_t), intent(in) :: layout
        real(dp), intent(in) :: stiffness(:, :)
        logical :: valid
        integer :: mode

        valid = terpsichore_pseudoplasma_fixture_is_valid(vacuum)
        if (.not. valid) return
        valid = plasma%intervals == vacuum%plasma_intervals &
            .and. plasma%modes == vacuum%modes
        if (.not. valid) return
        valid = allocated(plasma%mode_m)
        if (.not. valid) return
        valid = allocated(plasma%mode_n)
        if (.not. valid) return
        valid = size(plasma%mode_m) == plasma%modes &
            .and. size(plasma%mode_n) == plasma%modes
        if (.not. valid) return
        valid = all(plasma%mode_m == vacuum%mode_m) &
            .and. all(plasma%mode_n == vacuum%mode_n)
        if (.not. valid) return
        valid = layout%trials == plasma%modes &
            .and. layout%intervals == plasma%intervals &
            .and. layout%outer_normal_retained
        if (.not. valid) return
        valid = all(shape(stiffness) == layout%total_unknowns)
        if (.not. valid) return
        valid = all(ieee_is_finite(stiffness))
        if (.not. valid) return
        do mode = 1, plasma%modes
            valid = normal_global_index(layout, layout%intervals, mode) > 0
            if (.not. valid) return
        end do
    end function compatible_inputs

end module terpsichore_pseudoplasma_coupling
