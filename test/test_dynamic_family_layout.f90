program test_dynamic_family_layout
    use, intrinsic :: iso_fortran_env, only: error_unit
    use dynamic_family_layout, only: build_dynamic_family_layout, &
        dynamic_family_layout_t, dynamic_layout_ok, eta_global_index, &
        mu_global_index, normal_global_index
    implicit none

    type(dynamic_family_layout_t) :: layout
    integer :: info

    call build_dynamic_family_layout(3, 4, layout, info)
    call require(info == dynamic_layout_ok, "valid dynamic layout failed")
    call require(layout%normal_unknowns == 9, "normal count is wrong")
    call require(layout%eta_unknowns == 12, "eta count is wrong")
    call require(layout%mu_unknowns == 12, "mu count is wrong")
    call require(layout%total_unknowns == 33, "total count is wrong")
    call require(normal_global_index(layout, 0, 1) == 0, &
        "axis normal coefficient was retained")
    call require(normal_global_index(layout, 4, 1) == 0, &
        "edge normal coefficient was retained")
    call require(normal_global_index(layout, 1, 1) == 1, &
        "first interior normal index is wrong")
    call require(normal_global_index(layout, 3, 3) == 9, &
        "last interior normal index is wrong")
    call require(eta_global_index(layout, 1, 1) == 10, &
        "first eta index is wrong")
    call require(eta_global_index(layout, 4, 3) == 21, &
        "last eta index is wrong")
    call require(mu_global_index(layout, 1, 1) == 22, &
        "first mu index is wrong")
    call require(mu_global_index(layout, 4, 3) == 33, &
        "last mu index is wrong")
    call build_dynamic_family_layout(0, 4, layout, info)
    call require(info /= dynamic_layout_ok, "zero-trial layout was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_dynamic_family_layout
