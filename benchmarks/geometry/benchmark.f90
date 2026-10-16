program geometry_timing
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use cartesian_primitive_geometry_original, only: old_build => &
        build_primitive_geometry_point, old_point_t => primitive_geometry_point_t
    use cartesian_primitive_geometry, only: new_build => &
        build_primitive_geometry_point, new_point_t => primitive_geometry_point_t
    implicit none
    integer, parameter :: points = 4096, passes = 2048, trials = 11
    real(dp) :: jets(27, points), flux(2, points)
    integer :: periods(points)
    real(dp) :: old_seconds(trials), new_seconds(trials)
    real(dp) :: old_checksum, new_checksum, maximum_error, minimum_frame_ratio
    real(dp) :: reference(22), actual(22), frame_scale, phase
    type(old_point_t) :: old_point
    type(new_point_t) :: new_point
    integer :: point, component, trial, info, negative_count
    integer(int64) :: clock_rate

    ! A changing, diagonally dominant, nonorthogonal frame and arbitrary
    ! symmetric second jets. Reflect all z components on half the samples.
    do point = 1, points
        phase = 0.0037_dp * real(point, dp)
        jets(1:9, point) = [1.2_dp, 0.2_dp, -0.1_dp, 0.3_dp, 1.1_dp, &
            0.25_dp, -0.15_dp, 0.4_dp, 1.3_dp]
        do component = 1, 9
            jets(component, point) = jets(component, point) + &
                0.08_dp * sin(phase * real(component + 1, dp))
        end do
        do component = 10, 27
            jets(component, point) = 0.4_dp * &
                sin(phase * real(component + 1, dp)) + &
                0.05_dp * cos(phase * real(component - 3, dp))
        end do
        if (mod(point, 2) == 0) jets(3:27:3, point) = -jets(3:27:3, point)
        flux(1, point) = 7.0_dp + 0.1_dp * sin(phase)
        flux(2, point) = -3.0_dp + 0.1_dp * cos(phase)
        periods(point) = 1 + mod(point, 7)
    end do

    maximum_error = 0.0_dp
    minimum_frame_ratio = huge(1.0_dp)
    negative_count = 0
    do point = 1, points
        call evaluate_old(point, old_point, info)
        if (info /= 0) error stop 'original rejected a regular fixture'
        call evaluate_new(point, new_point, info)
        if (info /= 0) error stop 'generated caller rejected a regular fixture'
        call flatten_old(old_point, reference)
        call flatten_new(new_point, actual)
        if (.not. all(ieee_is_finite(reference))) error stop 'nonfinite original'
        if (.not. all(ieee_is_finite(actual))) error stop 'nonfinite current'
        maximum_error = max(maximum_error, &
            maxval(abs(actual - reference) / max(1.0_dp, abs(actual), abs(reference))))
        frame_scale = sqrt(sum(jets(1:3, point)**2)) * &
            sqrt(sum(jets(4:6, point)**2)) * sqrt(sum(jets(7:9, point)**2))
        minimum_frame_ratio = min(minimum_frame_ratio, &
            abs(new_point%signed_jacobian) / frame_scale)
        if (new_point%signed_jacobian < 0.0_dp) negative_count = negative_count + 1
    end do
    if (maximum_error > 1.0e-13_dp) error stop 'old/new geometries disagree'
    if (negative_count /= points / 2) error stop 'both orientation signs required'
    write (*, '(a,es24.16)') 'maximum normalized output difference: ', maximum_error
    write (*, '(a,es24.16)') 'minimum normalized frame determinant: ', minimum_frame_ratio
    write (*, '(a,i0)') 'negative orientation samples: ', negative_count
    write (*, '(a,i0)') 'calls per timed block: ', points * passes
    call system_clock(count_rate=clock_rate)
    write (*, '(a,i0)') 'clock ticks per second: ', clock_rate

    call time_old(16, phase, old_checksum)
    call time_new(16, phase, new_checksum)
    do trial = 1, trials
        if (mod(trial, 2) == 1) then
            call time_old(passes, old_seconds(trial), old_checksum)
            call time_new(passes, new_seconds(trial), new_checksum)
        else
            call time_new(passes, new_seconds(trial), new_checksum)
            call time_old(passes, old_seconds(trial), old_checksum)
        end if
        if (abs(old_checksum - new_checksum) > 1.0e-12_dp * &
            max(1.0_dp, abs(old_checksum), abs(new_checksum))) &
            error stop 'timed checksums disagree'
        write (*, '(i3,2f14.8,2es24.16)') trial, old_seconds(trial), &
            new_seconds(trial), old_checksum, new_checksum
    end do
    write (*, '(a,f14.8)') 'original median seconds: ', median(old_seconds)
    write (*, '(a,f14.8)') 'current median seconds: ', median(new_seconds)
    write (*, '(a,f14.8)') 'current/original median ratio: ', &
        median(new_seconds) / median(old_seconds)

contains

    subroutine evaluate_old(index, geometry, status)
        integer, intent(in) :: index
        type(old_point_t), intent(out) :: geometry
        integer, intent(out) :: status
        call old_build(jets(1:3, index), jets(4:6, index), jets(7:9, index), &
            jets(10:12, index), jets(13:15, index), jets(16:18, index), &
            jets(19:21, index), jets(22:24, index), jets(25:27, index), &
            periods(index), flux(1, index), flux(2, index), geometry, status)
    end subroutine evaluate_old

    subroutine evaluate_new(index, geometry, status)
        integer, intent(in) :: index
        type(new_point_t), intent(out) :: geometry
        integer, intent(out) :: status
        call new_build(jets(1:3, index), jets(4:6, index), jets(7:9, index), &
            jets(10:12, index), jets(13:15, index), jets(16:18, index), &
            jets(19:21, index), jets(22:24, index), jets(25:27, index), &
            periods(index), flux(1, index), flux(2, index), geometry, status)
    end subroutine evaluate_new

    subroutine time_old(repeats, seconds, checksum)
        integer, intent(in) :: repeats
        real(dp), intent(out) :: seconds, checksum
        type(old_point_t) :: geometry
        integer :: pass, index, status
        integer(int64) :: start, finish
        checksum = 0.0_dp
        call system_clock(start)
        do pass = 1, repeats
            do index = 1, points
                call evaluate_old(index, geometry, status)
                if (status /= 0) error stop 'timed original failure'
                checksum = checksum + geometry%signed_jacobian + geometry%mod_b &
                    + geometry%metric(1, 1) + geometry%second_form(1, 2)
            end do
        end do
        call system_clock(finish)
        seconds = real(finish - start, dp) / real(clock_rate, dp)
    end subroutine time_old

    subroutine time_new(repeats, seconds, checksum)
        integer, intent(in) :: repeats
        real(dp), intent(out) :: seconds, checksum
        type(new_point_t) :: geometry
        integer :: pass, index, status
        integer(int64) :: start, finish
        checksum = 0.0_dp
        call system_clock(start)
        do pass = 1, repeats
            do index = 1, points
                call evaluate_new(index, geometry, status)
                if (status /= 0) error stop 'timed generated caller failure'
                checksum = checksum + geometry%signed_jacobian + geometry%mod_b &
                    + geometry%metric(1, 1) + geometry%second_form(1, 2)
            end do
        end do
        call system_clock(finish)
        seconds = real(finish - start, dp) / real(clock_rate, dp)
    end subroutine time_new

    subroutine flatten_old(geometry, values)
        type(old_point_t), intent(in) :: geometry
        real(dp), intent(out) :: values(22)
        values(1:9) = reshape(geometry%metric, [9])
        values(10:13) = [geometry%signed_jacobian, geometry%jacobian_s, &
            geometry%jacobian_theta, geometry%jacobian_zeta]
        values(14:15) = geometry%b_contravariant
        values(16:17) = geometry%b_covariant
        values(18) = geometry%mod_b
        values(19:22) = reshape(geometry%second_form, [4])
    end subroutine flatten_old

    subroutine flatten_new(geometry, values)
        type(new_point_t), intent(in) :: geometry
        real(dp), intent(out) :: values(22)
        values(1:9) = reshape(geometry%metric, [9])
        values(10:13) = [geometry%signed_jacobian, geometry%jacobian_s, &
            geometry%jacobian_theta, geometry%jacobian_zeta]
        values(14:15) = geometry%b_contravariant
        values(16:17) = geometry%b_covariant
        values(18) = geometry%mod_b
        values(19:22) = reshape(geometry%second_form, [4])
    end subroutine flatten_new

    function median(values) result(value)
        real(dp), intent(in) :: values(:)
        real(dp) :: ordered(size(values)), value, temporary
        integer :: i, j
        ordered = values
        do i = 2, size(ordered)
            do j = i, 2, -1
                if (ordered(j) >= ordered(j - 1)) exit
                temporary = ordered(j)
                ordered(j) = ordered(j - 1)
                ordered(j - 1) = temporary
            end do
        end do
        value = ordered((size(ordered) + 1) / 2)
    end function median
end program geometry_timing
