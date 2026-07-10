module eigenvalue_tracking
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use block_tridiagonal, only: block_factor_t, block_tridiagonal_t, &
        factorize_shifted
    use family_assembly, only: assemble_family_blocks, &
        iterate_block_eigenvalue, surface_geometry_t
    implicit none
    private

    public :: certified_lowest_eigenvalue

contains

    subroutine certified_lowest_eigenvalue(geometry, mode_m, mode_n, &
            radial_step, eigenvalue, certificate_width, info, &
            class_selector)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: eigenvalue, certificate_width
        integer, intent(out) :: info
        integer, intent(in), optional :: class_selector
        type(block_tridiagonal_t) :: blocks
        real(dp) :: lower, upper, mid, scale, rayleigh, width_floor
        real(dp) :: window
        integer :: attempt, counted, count_upper, iterate_info
        integer :: below, above

        call assemble_family_blocks(geometry, mode_m, mode_n, &
            radial_step, blocks, info, class_selector)
        if (info /= 0) return
        call gershgorin_bounds(blocks, radial_step, lower, upper)
        scale = max(1.0_dp, abs(lower), abs(upper))
        lower = lower - 1.0e-3_dp * scale
        upper = upper + 1.0e-3_dp * scale
        width_floor = 1.0e-9_dp * scale
        count_upper = size(blocks%diag, 1) * size(blocks%diag, 3)
        do attempt = 1, 200
            if (count_upper == 1 .or. upper - lower <= width_floor) then
                call iterate_block_eigenvalue(blocks, radial_step, &
                    upper, rayleigh, iterate_info)
                if (iterate_info == 0) then
                    window = 1.0e-9_dp * max(1.0_dp, abs(rayleigh))
                    call shifted_count(blocks, radial_step, &
                        rayleigh - window, below, info)
                    if (info /= 0) return
                    call shifted_count(blocks, radial_step, &
                        rayleigh + window, above, info)
                    if (info /= 0) return
                    if (below == 0 .and. above >= 1) then
                        eigenvalue = rayleigh
                        certificate_width = window
                        info = 0
                        return
                    end if
                    if (below >= 1) then
                        upper = rayleigh - window
                        count_upper = below
                        cycle
                    end if
                end if
            end if
            mid = 0.5_dp * (lower + upper)
            call shifted_count(blocks, radial_step, mid, counted, info)
            if (info /= 0) return
            if (counted == 0) then
                lower = mid
            else
                upper = mid
                count_upper = counted
            end if
        end do
        info = -1
    end subroutine certified_lowest_eigenvalue

    pure subroutine gershgorin_bounds(blocks, radial_step, lower, upper)
        type(block_tridiagonal_t), intent(in) :: blocks
        real(dp), intent(in) :: radial_step
        real(dp), intent(out) :: lower, upper
        real(dp) :: center, radius
        integer :: trials, nodes, i, r, c

        trials = size(blocks%diag, 1)
        nodes = size(blocks%diag, 3)
        lower = huge(1.0_dp)
        upper = -huge(1.0_dp)
        do i = 1, nodes
            do r = 1, trials
                center = blocks%diag(r, r, i)
                radius = 0.0_dp
                do c = 1, trials
                    if (c /= r) then
                        radius = radius + abs(blocks%diag(r, c, i))
                    end if
                    if (i > 1) then
                        radius = radius + abs(blocks%off(r, c, i - 1))
                    end if
                    if (i < nodes) then
                        radius = radius + abs(blocks%off(c, r, i))
                    end if
                end do
                lower = min(lower, center - radius)
                upper = max(upper, center + radius)
            end do
        end do
        lower = lower / radial_step
        upper = upper / radial_step
    end subroutine gershgorin_bounds

    subroutine shifted_count(blocks, radial_step, shift, counted, info)
        type(block_tridiagonal_t), intent(in) :: blocks
        real(dp), intent(in) :: radial_step, shift
        integer, intent(out) :: counted, info
        type(block_factor_t) :: factor
        real(dp) :: nudged

        call factorize_shifted(blocks, shift * radial_step, factor, &
            info)
        if (info /= 0) then
            nudged = shift + 1.0e-12_dp * max(1.0_dp, abs(shift))
            call factorize_shifted(blocks, nudged * radial_step, &
                factor, info)
            if (info /= 0) return
        end if
        counted = factor%negative_count
    end subroutine shifted_count

end module eigenvalue_tracking
