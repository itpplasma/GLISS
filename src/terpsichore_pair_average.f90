module terpsichore_pair_average
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: terpsichore_pair_ok = 0
    integer, parameter, public :: terpsichore_pair_invalid = -1

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: terpsichore_pair_averages

contains

    ! Pair-averaged trigonometric sums over the TERPSICHORE tensor
    ! angular grid, in the naive double-loop convention:
    !   normal_normal(a, b)   = s Sum_p field(p) sin_a(p) sin_b(p)
    !   normal_tangent(a, b)  = s Sum_p field(p) sin_a(p) cos_b(p)
    !   tangent_tangent(a, b) = s Sum_p field(p) cos_a(p) cos_b(p)
    ! with angle_x(p) = m_x theta_p - n_x zeta_p and s = 2/points.
    ! Product-to-sum turns each pair entry into two lookups of the
    ! separable field transforms on the difference and sum mode table,
    ! so the cost is O(points * table + modes^2) per field instead of
    ! O(points * modes^2).
    subroutine terpsichore_pair_averages(field, poloidal_points, &
            toroidal_points, field_periods, mode_m, mode_n, &
            normal_normal, normal_tangent, tangent_tangent, info)
        real(dp), intent(in) :: field(:)
        integer, intent(in) :: poloidal_points, toroidal_points
        integer, intent(in) :: field_periods
        integer, intent(in) :: mode_m(:), mode_n(:)
        real(dp), intent(out) :: normal_normal(:, :)
        real(dp), intent(out) :: normal_tangent(:, :)
        real(dp), intent(out) :: tangent_tangent(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: cos_table(:, :), sin_table(:, :)
        real(dp) :: cos_diff, sin_diff, cos_sum, sin_sum
        integer :: modes, m_span, n_span, a, b, dm, dn, sm, sn

        info = terpsichore_pair_invalid
        modes = size(mode_m)
        if (size(mode_n) /= modes .or. modes < 1) return
        if (poloidal_points < 1 .or. toroidal_points < 1) return
        if (field_periods < 1) return
        if (size(field) /= poloidal_points * toroidal_points) return
        if (.not. all(ieee_is_finite(field))) return
        if (any(shape(normal_normal) /= modes)) return
        if (any(shape(normal_tangent) /= modes)) return
        if (any(shape(tangent_tangent) /= modes)) return

        m_span = 2 * maxval(abs(mode_m))
        n_span = 2 * maxval(abs(mode_n))
        call build_field_transforms(field, poloidal_points, &
            toroidal_points, field_periods, m_span, n_span, cos_table, &
            sin_table)

        do b = 1, modes
            do a = 1, modes
                dm = mode_m(a) - mode_m(b)
                dn = mode_n(a) - mode_n(b)
                sm = mode_m(a) + mode_m(b)
                sn = mode_n(a) + mode_n(b)
                call lookup(cos_table, sin_table, n_span, dm, dn, &
                    cos_diff, sin_diff)
                call lookup(cos_table, sin_table, n_span, sm, sn, &
                    cos_sum, sin_sum)
                normal_normal(a, b) = 0.5_dp * (cos_diff - cos_sum)
                normal_tangent(a, b) = 0.5_dp * (sin_diff + sin_sum)
                tangent_tangent(a, b) = 0.5_dp * (cos_diff + cos_sum)
            end do
        end do
        info = terpsichore_pair_ok
    end subroutine terpsichore_pair_averages

    ! cos_table(dm, dn) = s Sum_p field(p) cos(dm theta_p - dn zeta_p)
    ! for dm >= 0; negative dm is served by the lookup parity flip.
    subroutine build_field_transforms(field, poloidal_points, &
            toroidal_points, field_periods, m_span, n_span, cos_table, &
            sin_table)
        real(dp), intent(in) :: field(:)
        integer, intent(in) :: poloidal_points, toroidal_points
        integer, intent(in) :: field_periods, m_span, n_span
        real(dp), allocatable, intent(out) :: cos_table(:, :)
        real(dp), allocatable, intent(out) :: sin_table(:, :)
        real(dp), allocatable :: stage_cos(:, :), stage_sin(:, :)
        real(dp) :: theta, zeta, angle, scale, cos_z, sin_z
        integer :: j, k, dm, dn, point

        allocate (stage_cos(0:m_span, toroidal_points))
        allocate (stage_sin(0:m_span, toroidal_points))
        stage_cos = 0.0_dp
        stage_sin = 0.0_dp
        do k = 1, toroidal_points
            do j = 1, poloidal_points
                point = j + (k - 1) * poloidal_points
                theta = two_pi * real(j - 1, dp) &
                    / real(poloidal_points, dp)
                do dm = 0, m_span
                    angle = real(dm, dp) * theta
                    stage_cos(dm, k) = stage_cos(dm, k) &
                        + field(point) * cos(angle)
                    stage_sin(dm, k) = stage_sin(dm, k) &
                        + field(point) * sin(angle)
                end do
            end do
        end do

        allocate (cos_table(0:m_span, -n_span:n_span))
        allocate (sin_table(0:m_span, -n_span:n_span))
        cos_table = 0.0_dp
        sin_table = 0.0_dp
        scale = 2.0_dp / real(poloidal_points * toroidal_points, dp)
        do dn = -n_span, n_span
            do k = 1, toroidal_points
                zeta = two_pi * real(k - 1, dp) &
                    / real(toroidal_points * field_periods, dp)
                angle = -real(dn, dp) * zeta
                cos_z = cos(angle)
                sin_z = sin(angle)
                do dm = 0, m_span
                    cos_table(dm, dn) = cos_table(dm, dn) &
                        + stage_cos(dm, k) * cos_z &
                        - stage_sin(dm, k) * sin_z
                    sin_table(dm, dn) = sin_table(dm, dn) &
                        + stage_sin(dm, k) * cos_z &
                        + stage_cos(dm, k) * sin_z
                end do
            end do
        end do
        cos_table = scale * cos_table
        sin_table = scale * sin_table
    end subroutine build_field_transforms

    pure subroutine lookup(cos_table, sin_table, n_span, dm, dn, &
            cos_value, sin_value)
        real(dp), intent(in) :: cos_table(0:, -n_span:), &
            sin_table(0:, -n_span:)
        integer, intent(in) :: n_span, dm, dn
        real(dp), intent(out) :: cos_value, sin_value

        if (dm >= 0) then
            cos_value = cos_table(dm, dn)
            sin_value = sin_table(dm, dn)
        else
            cos_value = cos_table(-dm, -dn)
            sin_value = -sin_table(-dm, -dn)
        end if
    end subroutine lookup

end module terpsichore_pair_average
