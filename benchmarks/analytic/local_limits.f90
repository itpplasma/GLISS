program benchmark_local_limits
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use local_mode_model, only: assemble_local_mode
    use helical_cylinder_limit, only: helical_vertical_margin
    use symmetric_eigensolver, only: solve_three_component_modes
    implicit none
    real(dp) :: stiffness(3, 3), mass(3, 3), values(3), vectors(3, 3)
    real(dp) :: angle, fraction, wave_vector(3)
    integer :: i, info, unit
    character(len=1024) :: filename

    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) error stop 'Supply output CSV path'
    open (newunit=unit, file=trim(filename), status='replace')
    write (unit, '(a)') 'angle_rad,slow_s_minus2,alfven_s_minus2,fast_s_minus2'
    do i = 0, 40
        angle = acos(-1.0_dp)*real(i, dp)/80.0_dp
        wave_vector(1) = sin(angle)
        wave_vector(2) = 0.0_dp
        wave_vector(3) = cos(angle)
        call assemble_local_mode(wave_vector, &
            [0.0_dp, 0.0_dp, 1.0_dp], 1.0e5_dp, 2.0_dp, &
            5.0_dp/3.0_dp, 0.0_dp, [1.0_dp, 0.0_dp, 0.0_dp], &
            stiffness, mass)
        call solve_three_component_modes(stiffness, mass, values, &
            vectors, info)
        if (info /= 0) error stop 'Local eigensolve failed'
        write (unit, '(3(es25.16,a),es25.16)') angle, ',', values(1), &
            ',', values(2), ',', values(3)
    end do
    close (unit)
    call get_command_argument(2, filename)
    if (len_trim(filename) == 0) error stop 'Supply helical CSV path'
    open (newunit=unit, file=trim(filename), status='replace')
    write (unit, '(a)') 'external_iota_fraction,normalized_vertical_margin'
    do i = 0, 40
        fraction = real(i, dp)/40.0_dp
        write (unit, '(es25.16,a,es25.16)') fraction, ',', &
            helical_vertical_margin(2.0_dp, fraction)
    end do
    close (unit)
end program benchmark_local_limits
