module newcomb_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use two_component_kernel, only: two_component_components
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    type, public :: cylinder_profiles_t
        real(dp) :: length = 0.0_dp
        real(dp) :: b_axial = 0.0_dp
        real(dp) :: b_linear = 0.0_dp
        real(dp) :: b_cubic = 0.0_dp
    end type cylinder_profiles_t

    public :: lowest_eigenvalue_single_mode

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*)
            real(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

contains

    subroutine lowest_eigenvalue_single_mode(profiles, poloidal_mode, &
            toroidal_mode, minor_radius, n_radial, lowest)
        type(cylinder_profiles_t), intent(in) :: profiles
        integer, intent(in) :: poloidal_mode, toroidal_mode, n_radial
        real(dp), intent(in) :: minor_radius
        real(dp), intent(out) :: lowest
        real(dp), allocatable :: stiffness(:, :), work(:)
        real(dp), allocatable :: eigenvalues(:)
        real(dp) :: element(2, 2), step
        integer :: i, info, unknowns

        unknowns = n_radial - 1
        allocate (stiffness(unknowns, unknowns), source=0.0_dp)
        step = minor_radius / real(n_radial, dp)
        do i = 1, n_radial
            call element_matrix(profiles, poloidal_mode, toroidal_mode, &
                (real(i, dp) - 0.5_dp) * step, step, element)
            if (i > 1) then
                stiffness(i - 1, i - 1) = stiffness(i - 1, i - 1) &
                    + element(1, 1)
            end if
            if (i > 1 .and. i < n_radial) then
                stiffness(i - 1, i) = stiffness(i - 1, i) + element(1, 2)
                stiffness(i, i - 1) = stiffness(i, i - 1) + element(2, 1)
            end if
            if (i < n_radial) then
                stiffness(i, i) = stiffness(i, i) + element(2, 2)
            end if
        end do

        allocate (eigenvalues(unknowns), work(8 * unknowns))
        call dsyev("N", "U", unknowns, stiffness, unknowns, eigenvalues, &
            work, size(work), info)
        if (info /= 0) then
            lowest = huge(1.0_dp)
            return
        end if
        lowest = eigenvalues(1) / step
    end subroutine lowest_eigenvalue_single_mode

    subroutine element_matrix(profiles, poloidal_mode, toroidal_mode, &
            radius_mid, step, element)
        type(cylinder_profiles_t), intent(in) :: profiles
        integer, intent(in) :: poloidal_mode, toroidal_mode
        real(dp), intent(in) :: radius_mid, step
        real(dp), intent(out) :: element(2, 2)
        real(dp) :: coeff_xx, coeff_xd, coeff_dd
        real(dp) :: value_weight(2), slope_weight(2)
        integer :: a, b

        call reduced_quadratic(profiles, poloidal_mode, toroidal_mode, &
            radius_mid, coeff_xx, coeff_xd, coeff_dd)
        value_weight = [0.5_dp, 0.5_dp]
        slope_weight = [-1.0_dp / step, 1.0_dp / step]
        do b = 1, 2
            do a = 1, 2
                element(a, b) = step * (coeff_xx * value_weight(a) &
                    * value_weight(b) + coeff_xd * (value_weight(a) &
                    * slope_weight(b) + slope_weight(a) &
                    * value_weight(b)) + coeff_dd * slope_weight(a) &
                    * slope_weight(b))
            end do
        end do
    end subroutine element_matrix

    subroutine reduced_quadratic(profiles, poloidal_mode, toroidal_mode, &
            radius_mid, coeff_xx, coeff_xd, coeff_dd)
        type(cylinder_profiles_t), intent(in) :: profiles
        integer, intent(in) :: poloidal_mode, toroidal_mode
        real(dp), intent(in) :: radius_mid
        real(dp), intent(out) :: coeff_xx, coeff_xd, coeff_dd
        real(dp) :: geometry(14), angle_m, angle_n
        real(dp) :: c1_x, c2_x, c2_y, c3_x, c3_d, c3_y
        real(dp) :: q_xx, q_xd, q_dd, q_xy, q_dy, q_yy

        call cylinder_geometry(profiles, radius_mid, geometry)
        angle_m = two_pi * real(poloidal_mode, dp)
        angle_n = two_pi * real(toroidal_mode, dp)

        call component_coefficients(geometry, -angle_m, angle_n, &
            c1_x, c2_x, c2_y, c3_x, c3_d, c3_y)

        q_xx = 0.5_dp * (c1_x**2 + c2_x**2 + c3_x**2 &
            - geometry(13) * geometry(14))
        q_xd = 0.5_dp * c3_x * c3_d
        q_dd = 0.5_dp * c3_d**2
        q_xy = 0.5_dp * (c2_x * c2_y + c3_x * c3_y)
        q_dy = 0.5_dp * c3_d * c3_y
        q_yy = 0.5_dp * (c2_y**2 + c3_y**2)

        coeff_xx = (q_xx - q_xy**2 / q_yy) * abs(geometry(7))
        coeff_xd = (q_xd - q_xy * q_dy / q_yy) * abs(geometry(7))
        coeff_dd = (q_dd - q_dy**2 / q_yy) * abs(geometry(7))
    end subroutine reduced_quadratic

    subroutine component_coefficients(geometry, d_theta, d_zeta, c1_x, &
            c2_x, c2_y, c3_x, c3_d, c3_y)
        real(dp), intent(in) :: geometry(14), d_theta, d_zeta
        real(dp), intent(out) :: c1_x, c2_x, c2_y, c3_x, c3_d, c3_y
        real(dp) :: c1, c2, c3

        call two_component_components(geometry(1), geometry(2), &
            geometry(3), geometry(4), geometry(5), geometry(6), &
            geometry(7), geometry(8), geometry(9), geometry(10), &
            geometry(11), geometry(12), 0.0_dp, 0.0_dp, 0.0_dp, &
            d_theta, d_zeta, 0.0_dp, 0.0_dp, c1, c2, c3)
        c1_x = c1
        call two_component_components(geometry(1), geometry(2), &
            geometry(3), geometry(4), geometry(5), geometry(6), &
            geometry(7), geometry(8), geometry(9), geometry(10), &
            geometry(11), geometry(12), 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, c1, c2, c3)
        c2_x = c2
        c3_x = c3
        call two_component_components(geometry(1), geometry(2), &
            geometry(3), geometry(4), geometry(5), geometry(6), &
            geometry(7), geometry(8), geometry(9), geometry(10), &
            geometry(11), geometry(12), 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, c1, c2, c3)
        c3_d = c3
        call two_component_components(geometry(1), geometry(2), &
            geometry(3), geometry(4), geometry(5), geometry(6), &
            geometry(7), geometry(8), geometry(9), geometry(10), &
            geometry(11), geometry(12), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, d_theta, d_zeta, c1, c2, c3)
        c2_y = c2
        c3_y = c3
    end subroutine component_coefficients

    subroutine cylinder_geometry(profiles, radius_mid, geometry)
        type(cylinder_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: radius_mid
        real(dp), intent(out) :: geometry(14)
        real(dp) :: r, b_theta, b_theta_slope_r, length

        length = profiles%length
        r = radius_mid
        b_theta = profiles%b_linear * r + profiles%b_cubic * r**3
        b_theta_slope_r = profiles%b_linear &
            + 3.0_dp * profiles%b_cubic * r**2
        geometry(1) = two_pi * r * profiles%b_axial
        geometry(2) = length * b_theta
        geometry(3) = two_pi * profiles%b_axial
        geometry(4) = length * b_theta_slope_r
        geometry(5) = length * profiles%b_axial
        geometry(6) = two_pi * r * b_theta
        geometry(7) = two_pi * length * r
        geometry(8) = sqrt(b_theta**2 + profiles%b_axial**2)
        geometry(9) = 1.0_dp
        geometry(10) = (b_theta_slope_r + b_theta / r) &
            * profiles%b_axial
        geometry(11) = -b_theta * (b_theta_slope_r + b_theta / r)
        geometry(12) = 0.0_dp
        geometry(13) = 2.0_dp * b_theta * (b_theta_slope_r + b_theta / r) &
            / r
        geometry(14) = 1.0_dp
    end subroutine cylinder_geometry

end module newcomb_limit
