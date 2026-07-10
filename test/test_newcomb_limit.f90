program test_newcomb_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use newcomb_limit, only: cylinder_profiles_t, &
        lowest_eigenvalue_single_mode
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    type(cylinder_profiles_t) :: profiles
    real(dp) :: unstable_coarse, unstable_fine, stable_mode

    profiles%length = 6.0_dp * pi
    profiles%b_axial = 1.0_dp
    profiles%b_linear = 0.3_dp
    profiles%b_cubic = 0.4_dp

    call lowest_eigenvalue_single_mode(profiles, 1, 1, 0.5_dp, 100, &
        unstable_coarse)
    call lowest_eigenvalue_single_mode(profiles, 1, 1, 0.5_dp, 200, &
        unstable_fine)
    call lowest_eigenvalue_single_mode(profiles, 2, 1, 0.5_dp, 100, &
        stable_mode)

    if (unstable_coarse >= 0.0_dp) then
        write (error_unit, "(a, es12.4)") &
            "resonant Suydam-unstable mode is not unstable: ", &
            unstable_coarse
        error stop 1
    end if
    if (unstable_fine > unstable_coarse) then
        write (error_unit, "(a)") &
            "instability does not deepen under refinement"
        error stop 1
    end if
    if (stable_mode <= 0.0_dp) then
        write (error_unit, "(a, es12.4)") &
            "non-resonant mode is not stable: ", stable_mode
        error stop 1
    end if
    write (*, "(a)") "PASS"
end program test_newcomb_limit
