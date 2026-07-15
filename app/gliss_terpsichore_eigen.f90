program gliss_terpsichore_eigen
    use, intrinsic :: iso_fortran_env, only: error_unit
    use terpsichore_fixed_boundary_spectrum, only: &
        solve_terpsichore_fixed_boundary_file, &
        terpsichore_fixed_boundary_result_t, terpsichore_fixed_spectrum_ok
    use terpsichore_pseudoplasma_spectrum, only: &
        solve_terpsichore_pseudoplasma_files, &
        terpsichore_pseudoplasma_result_t, &
        terpsichore_pseudoplasma_spectrum_ok
    implicit none

    type(terpsichore_fixed_boundary_result_t) :: fixed_result
    type(terpsichore_pseudoplasma_result_t) :: free_result
    character(len=1024) :: matrix_file, vacuum_file
    character(len=128) :: message
    integer :: info, vacuum_intervals

    call parse_arguments(matrix_file, vacuum_intervals, vacuum_file)
    if (vacuum_intervals == 0) then
        call solve_terpsichore_fixed_boundary_file(matrix_file, fixed_result, &
            info, message)
        if (info /= terpsichore_fixed_spectrum_ok) call fail(message, 1)
        write (error_unit, "(a,i0)") "negative count at zero shift: ", &
            fixed_result%negative_count
        call write_fixed_result(fixed_result)
    else
        call solve_terpsichore_pseudoplasma_files(matrix_file, &
            vacuum_intervals, vacuum_file, free_result, info, message)
        if (info /= terpsichore_pseudoplasma_spectrum_ok) &
            call fail(message, 1)
        write (error_unit, "(a,i0)") "negative count at zero shift: ", &
            free_result%negative_count
        call write_free_result(free_result)
    end if

contains

    subroutine parse_arguments(matrix_path, vacuum_count, vacuum_path)
        character(len=*), intent(out) :: matrix_path, vacuum_path
        integer, intent(out) :: vacuum_count
        character(len=1024) :: argument
        integer :: count, status

        matrix_path = ""
        vacuum_path = ""
        count = command_argument_count()
        if (count < 2 .or. count > 3) call usage()
        call get_command_argument(1, matrix_path, status=status)
        if (status /= 0 .or. len_trim(matrix_path) == 0) call usage()
        call get_command_argument(2, argument, status=status)
        if (status /= 0) call fail("IVAC must be a nonnegative decimal integer", 2)
        if (.not. decimal_integer(argument)) &
            call fail("IVAC must be a nonnegative decimal integer", 2)
        read (argument, *, iostat=status) vacuum_count
        if (status /= 0) call fail("IVAC is outside the supported integer range", 2)
        if (vacuum_count < 0) &
            call fail("IVAC is outside the supported integer range", 2)
        if (vacuum_count == 0 .and. count /= 2) &
            call fail("IVAC=0 does not accept FORT.24", 2)
        if (vacuum_count > 0 .and. count /= 3) &
            call fail("IVAC>0 requires FORT.24", 2)
        if (vacuum_count == 0) return
        call get_command_argument(3, vacuum_path, status=status)
        if (status /= 0) call fail("IVAC>0 requires FORT.24", 2)
        if (len_trim(vacuum_path) == 0) call fail("IVAC>0 requires FORT.24", 2)
    end subroutine parse_arguments

    pure function decimal_integer(argument) result(valid)
        character(len=*), intent(in) :: argument
        logical :: valid
        integer :: i, digit

        valid = len_trim(argument) > 0
        if (.not. valid) return
        do i = 1, len_trim(argument)
            digit = iachar(argument(i:i))
            valid = digit >= iachar("0") .and. digit <= iachar("9")
            if (.not. valid) return
        end do
    end function decimal_integer

    subroutine write_fixed_result(result)
        type(terpsichore_fixed_boundary_result_t), intent(in) :: result

        write (*, "(a)") "unknowns,negative_count,eigenvalue,certificate," // &
            "residual,resolution"
        write (*, '(i0,",",i0,4(",",es24.16))') result%unknowns, &
            result%negative_count, result%eigenvalue, result%certificate, &
            result%residual, result%resolution
    end subroutine write_fixed_result

    subroutine write_free_result(result)
        type(terpsichore_pseudoplasma_result_t), intent(in) :: result

        write (*, "(a)") "unknowns,negative_count,eigenvalue,certificate," // &
            "residual,resolution,growth_rate,reference_eigenvalue," // &
            "reference_potential,computed_potential,reference_kinetic," // &
            "computed_kinetic,reference_residual,mode_overlap"
        write (*, '(i0,",",i0,12(",",es24.16))') result%unknowns, &
            result%negative_count, result%eigenvalue, result%certificate, &
            result%residual, result%resolution, result%growth_rate, &
            result%reference_eigenvalue, result%reference_potential, &
            result%computed_potential, result%reference_kinetic, &
            result%computed_kinetic, result%reference_residual, &
            result%mode_overlap
    end subroutine write_free_result

    subroutine usage()
        write (error_unit, "(a)") &
            "usage: gliss_terpsichore_eigen FORT.23 IVAC [FORT.24]"
        flush (error_unit)
        stop 2
    end subroutine usage

    subroutine fail(message, status)
        character(len=*), intent(in) :: message
        integer, intent(in) :: status

        write (error_unit, "(a)") trim(message)
        flush (error_unit)
        if (status == 2) stop 2
        stop 1
    end subroutine fail

end program gliss_terpsichore_eigen
