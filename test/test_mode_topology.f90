program test_mode_topology
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use mode_topology, only: axis_form_function, &
        axis_form_function_slope, build_mode_family, family_count, &
        mode_family_t, modes_coupled, nonzero_family_count
    implicit none

    type(mode_family_t) :: family, other
    integer :: periods, index_a, index_b, i, j, m, info
    integer :: covered(-8:8)
    real(dp) :: s, step, central

    do periods = 1, 12
        if (family_count(periods) /= periods / 2 + 1) then
            call fail("family count does not include the zero family")
        end if
        if (nonzero_family_count(periods) /= periods / 2) then
            call fail("nonzero family count is not floor(NT/2)")
        end if
    end do

    do periods = 2, 8
        covered = 0
        do index_a = 1, nonzero_family_count(periods)
            call build_mode_family(periods, index_a, 2, 8, family, info)
            if (info /= 0) call fail("valid family was rejected")
            if (size(family%toroidal) == 0) call fail("family is empty")
            do i = 1, size(family%toroidal)
                do j = 1, size(family%toroidal)
                    if (.not. modes_coupled(family%toroidal(i), &
                        family%toroidal(j), periods)) then
                        call fail("family is not closed under coupling")
                    end if
                end do
                if (family%poloidal(i) > 0 .or. family%toroidal(i) >= 0) &
                    then
                    if (modulo(family%toroidal(i), periods) == 0) then
                        call fail("family contains the excluded residue")
                    end if
                end if
                covered(family%toroidal(i)) = &
                    covered(family%toroidal(i)) + 1
            end do
            do index_b = index_a + 1, nonzero_family_count(periods)
                call build_mode_family(periods, index_b, 2, 8, other, info)
                if (info /= 0) call fail("valid comparison family was rejected")
                do i = 1, size(family%toroidal)
                    do j = 1, size(other%toroidal)
                        if (modes_coupled(family%toroidal(i), &
                            other%toroidal(j), periods)) then
                            call fail("distinct families are coupled")
                        end if
                    end do
                end do
            end do
        end do
        do i = -8, 8
            if (modulo(i, periods) == 0) cycle
            if (2 * (modulo(i, periods)) == periods .and. covered(i) &
                == 0) cycle
            if (covered(i) == 0) call fail("families do not cover a mode")
        end do
    end do

    call check_zero_family()
    if (modes_coupled(0, 0, 0)) call fail("zero periods coupled modes")
    if (modes_coupled(0, 0, -1)) call fail("negative periods coupled modes")
    call build_mode_family(3, -1, 2, 2, family, info)
    if (info == 0) call fail("negative family index was accepted")
    call build_mode_family(3, 2, 2, 2, family, info)
    if (info == 0) call fail("family index above maximum was accepted")
    call build_mode_family(3, 1, -1, 2, family, info)
    if (info == 0) call fail("negative truncation was accepted")

    do m = 0, 6
        s = 0.4_dp
        if (abs(axis_form_function(m, s) - s**(0.5_dp * m) * (1 - s)) > &
            1.0e-15_dp) call fail("form function value is wrong")
        if (abs(axis_form_function(m, 1.0_dp)) > 1.0e-15_dp) then
            call fail("form function does not vanish at the boundary")
        end if
        step = 1.0e-6_dp
        central = (axis_form_function(m, s + step) - &
            axis_form_function(m, s - step)) / (2.0_dp * step)
        if (abs(axis_form_function_slope(m, s) - central) > 1.0e-8_dp) &
            then
            call fail("form function slope is wrong")
        end if
    end do
    write (*, "(a)") "PASS"

contains

    subroutine check_zero_family()
        type(mode_family_t) :: zero_family, nonzero_family
        integer :: i, j, local_info

        call build_mode_family(3, 0, 2, 3, zero_family, local_info)
        if (local_info /= 0) call fail("zero family was rejected")
        if (size(zero_family%poloidal) /= 8) &
            call fail("zero family has the wrong size")
        if (any(modulo(zero_family%toroidal, 3) /= 0)) &
            call fail("zero family contains another residue")
        if (count(zero_family%poloidal == 0 .and. &
            zero_family%toroidal == 0) /= 1) &
            call fail("zero harmonic is not represented once")
        if (any(zero_family%poloidal == 0 .and. &
            zero_family%toroidal < 0)) &
            call fail("zero family duplicates negative axis modes")
        call build_mode_family(3, 1, 2, 3, nonzero_family, local_info)
        if (local_info /= 0) call fail("nonzero family comparison failed")
        do i = 1, size(zero_family%toroidal)
            do j = 1, size(nonzero_family%toroidal)
                if (modes_coupled(zero_family%toroidal(i), &
                    nonzero_family%toroidal(j), 3)) &
                    call fail("zero and nonzero families are coupled")
            end do
        end do
        call build_mode_family(1, 0, 1, 1, zero_family, local_info)
        if (local_info /= 0) call fail("one-period zero family was rejected")
    end subroutine check_zero_family

    subroutine fail(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") message
        error stop 1
    end subroutine fail

end program test_mode_topology
