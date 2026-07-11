program gliss_newcomb
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use newcomb_limit, only: cylinder_profiles_t, &
        lowest_artificial_stiffness_level
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: resolutions(6) = [25, 50, 100, 200, 400, 800]
    type(cylinder_profiles_t) :: profiles
    real(dp) :: lowest
    integer :: i, info

    profiles%length = 6.0_dp * pi
    profiles%b_axial = 1.0_dp
    profiles%b_linear = 0.3_dp
    profiles%b_cubic = 0.4_dp

    write (*, "(a)") "n_radial,m,n,artificial_stiffness_level"
    do i = 1, size(resolutions)
        call lowest_artificial_stiffness_level(profiles, 1, 1, 0.5_dp, &
            resolutions(i), lowest, info)
        if (info /= 0) error stop "m=1 Newcomb solve failed"
        write (*, "(i0, 2(',', i0), ',', es24.16)") resolutions(i), 1, &
            1, lowest
        call lowest_artificial_stiffness_level(profiles, 2, 1, 0.5_dp, &
            resolutions(i), lowest, info)
        if (info /= 0) error stop "m=2 Newcomb solve failed"
        write (*, "(i0, 2(',', i0), ',', es24.16)") resolutions(i), 2, &
            1, lowest
    end do
end program gliss_newcomb
