program gliss_axisymmetric_fkg
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use family_assembly, only: condensed_surface_coefficients, &
        surface_geometry_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: equilibrium_is_axisymmetric, &
        gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none

    integer, parameter :: n_theta = 64, n_zeta = 8
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(surface_geometry_t) :: surface
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    real(dp), allocatable :: f_matrix(:, :), g_matrix(:, :), k_matrix(:, :)
    integer, allocatable :: mode_m(:), mode_n(:), parity(:)
    character(len=1024) :: filename, token
    real(dp) :: step
    integer :: arguments, i, info, j, m, mode, poloidal_max
    integer :: toroidal_mode

    interface
        subroutine terminate_process(status) bind(C, name="exit")
            import c_int
            integer(c_int), value :: status
        end subroutine terminate_process
    end interface

    arguments = command_argument_count()
    if (arguments /= 3) call fail_usage("expected exactly three arguments")
    call read_argument(1, "EXPORT_FILE", filename)
    call read_integer_argument(2, "N", toroidal_mode)
    call read_integer_argument(3, "MMAX", poloidal_max)
    if (toroidal_mode <= 0) call fail_usage("N must be positive")
    if (poloidal_max < 1) call fail_usage("MMAX must be positive")

    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) call fail("equilibrium export could not be read")
    if (.not. equilibrium%has_chart_metric) &
        call fail("equilibrium export lacks g_st/g_sz chart metrics")
    if (equilibrium%field_periods /= 1) &
        call fail("axisymmetric comparison requires N_FP=1")
    if (.not. equilibrium_is_axisymmetric(equilibrium)) &
        call fail("equilibrium contains nonaxisymmetric harmonics")
    if (2 * poloidal_max + maxval(abs(equilibrium%poloidal_modes)) &
        >= n_theta) call fail_usage("MMAX aliases the fixed angular quadrature")

    allocate (mode_m(2 * poloidal_max + 1), mode_n(2 * poloidal_max + 1))
    mode_m(1) = 0
    mode_n(1) = toroidal_mode
    mode = 1
    do m = 1, poloidal_max
        mode = mode + 1
        mode_m(mode) = m
        mode_n(mode) = -toroidal_mode
        mode = mode + 1
        mode_m(mode) = m
        mode_n(mode) = toroidal_mode
    end do
    allocate (parity(size(mode_m)), source=1)
    call build_kernel_geometry(equilibrium, n_theta, n_zeta, fields, &
        drive, info)
    if (info /= mercier_ok) call fail("kernel geometry could not be built")
    step = 1.0_dp / real(size(equilibrium%s), dp)

    write (*, "(a)") "s,matrix,row_m,row_n,column_m,column_n,real"
    do i = 1, size(equilibrium%s)
        surface%fields = fields(:, :, :, i)
        surface%drive = drive(:, :, i)
        call condensed_surface_coefficients(surface, mode_m, mode_n, parity, &
            equilibrium%s(i), step, 1, f_matrix, k_matrix, g_matrix, info)
        if (info /= 0) call fail("surface coefficient extraction failed")
        do j = 1, size(mode_m)
            do m = 1, size(mode_m)
                call write_entry(equilibrium%s(i), "F", m, j, f_matrix(m, j))
                call write_entry(equilibrium%s(i), "K", m, j, k_matrix(m, j))
                call write_entry(equilibrium%s(i), "G", m, j, g_matrix(m, j))
            end do
        end do
    end do

contains

    subroutine write_entry(s, name, row, column, value)
        real(dp), intent(in) :: s, value
        character(len=1), intent(in) :: name
        integer, intent(in) :: row, column

        write (*, "(es24.16,2a,4(a,i0),a,es24.16)") s, ",", name, ",", &
            mode_m(row), ",", mode_n(row), ",", mode_m(column), ",", &
            mode_n(column), ",", value
    end subroutine write_entry

    subroutine fail_usage(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_axisymmetric_fkg: " // trim(message)
        write (error_unit, "(a)") &
            "usage: gliss_axisymmetric_fkg EXPORT_FILE N MMAX"
        call terminate_process(2_c_int)
    end subroutine fail_usage

    subroutine fail(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_axisymmetric_fkg: " // trim(message)
        call terminate_process(1_c_int)
    end subroutine fail

    subroutine read_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        integer :: status

        call get_command_argument(position, value, status=status)
        if (status /= 0) call fail_usage(trim(name) // " is too long")
        if (len_trim(value) == 0) call fail_usage(trim(name) // " is empty")
    end subroutine read_argument

    subroutine read_integer_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        integer, intent(out) :: value
        integer :: status

        call read_argument(position, name, token)
        read (token, *, iostat=status) value
        if (status /= 0) call fail_usage(trim(name) // " must be an integer")
    end subroutine read_integer_argument

end program gliss_axisymmetric_fkg
