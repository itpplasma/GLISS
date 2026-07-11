module newcomb_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use two_component_kernel, only: two_component_components
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
    integer, parameter, public :: newcomb_degenerate_reduction = 2
    integer, parameter, public :: newcomb_invalid_input = 1
    integer, parameter, public :: newcomb_ok = 0

    type, public :: cylinder_profiles_t
        real(dp) :: length = 0.0_dp
        real(dp) :: b_axial = 0.0_dp
        real(dp) :: b_linear = 0.0_dp
        real(dp) :: b_cubic = 0.0_dp
    end type cylinder_profiles_t

    public :: assemble_single_mode_stiffness
    public :: lowest_artificial_stiffness_level
    public :: newcomb_axis_regular_power
    public :: newcomb_quadratic_coefficients

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

    subroutine assemble_single_mode_stiffness(profiles, poloidal_mode, &
            toroidal_mode, minor_radius, n_radial, stiffness, info)
        type(cylinder_profiles_t), intent(in) :: profiles
        integer, intent(in) :: poloidal_mode, toroidal_mode, n_radial
        real(dp), intent(in) :: minor_radius
        real(dp), allocatable, intent(out) :: stiffness(:, :)
        integer, intent(out) :: info
        real(dp) :: element(2, 2), step
        integer :: first_node, i, left_index, right_index, unknowns

        if (.not. valid_inputs(profiles, minor_radius, n_radial)) then
            allocate (stiffness(0, 0))
            info = newcomb_invalid_input
            return
        end if
        first_node = 1
        if (newcomb_axis_regular_power(poloidal_mode) == 0) first_node = 0
        unknowns = n_radial - first_node
        allocate (stiffness(unknowns, unknowns), source=0.0_dp)
        step = minor_radius / real(n_radial, dp)
        do i = 1, n_radial
            call element_matrix(profiles, poloidal_mode, toroidal_mode, &
                (real(i, dp) - 0.5_dp) * step, step, element, info)
            if (info /= newcomb_ok) return
            left_index = i - first_node
            right_index = left_index + 1
            call add_element(stiffness, left_index, right_index, element)
        end do
        info = newcomb_ok
    end subroutine assemble_single_mode_stiffness

    subroutine lowest_artificial_stiffness_level(profiles, poloidal_mode, &
            toroidal_mode, minor_radius, n_radial, lowest, info)
        type(cylinder_profiles_t), intent(in) :: profiles
        integer, intent(in) :: poloidal_mode, toroidal_mode, n_radial
        real(dp), intent(in) :: minor_radius
        real(dp), intent(out) :: lowest
        integer, intent(out) :: info
        real(dp), allocatable :: stiffness(:, :), work(:)
        real(dp), allocatable :: eigenvalues(:)
        real(dp) :: step
        integer :: unknowns

        call assemble_single_mode_stiffness(profiles, poloidal_mode, &
            toroidal_mode, minor_radius, n_radial, stiffness, info)
        if (info /= newcomb_ok) then
            lowest = 0.0_dp
            return
        end if
        unknowns = size(stiffness, 1)
        step = minor_radius / real(n_radial, dp)
        allocate (eigenvalues(unknowns), work(8 * unknowns))
        call dsyev("N", "U", unknowns, stiffness, unknowns, eigenvalues, &
            work, size(work), info)
        if (info /= 0) then
            lowest = 0.0_dp
            info = newcomb_degenerate_reduction
            return
        end if
        lowest = eigenvalues(1) / step
        info = newcomb_ok
    end subroutine lowest_artificial_stiffness_level

    pure function newcomb_axis_regular_power(poloidal_mode) result(power)
        integer, intent(in) :: poloidal_mode
        integer :: power

        power = 1
        if (poloidal_mode /= 0) power = abs(poloidal_mode) - 1
    end function newcomb_axis_regular_power

    subroutine element_matrix(profiles, poloidal_mode, toroidal_mode, &
            radius_mid, step, element, info)
        type(cylinder_profiles_t), intent(in) :: profiles
        integer, intent(in) :: poloidal_mode, toroidal_mode
        real(dp), intent(in) :: radius_mid, step
        real(dp), intent(out) :: element(2, 2)
        integer, intent(out) :: info
        real(dp) :: coeff_xx, coeff_xd, coeff_dd
        real(dp) :: value_weight(2), slope_weight(2)
        integer :: a, b

        call newcomb_quadratic_coefficients(profiles, poloidal_mode, &
            toroidal_mode, radius_mid, coeff_xx, coeff_xd, coeff_dd, info)
        if (info /= newcomb_ok) return
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

    subroutine newcomb_quadratic_coefficients(profiles, poloidal_mode, &
            toroidal_mode, radius_mid, coeff_xx, coeff_xd, coeff_dd, info)
        type(cylinder_profiles_t), intent(in) :: profiles
        integer, intent(in) :: poloidal_mode, toroidal_mode
        real(dp), intent(in) :: radius_mid
        real(dp), intent(out) :: coeff_xx, coeff_xd, coeff_dd
        integer, intent(out) :: info
        real(dp) :: geometry(14), angle_m, angle_n
        real(dp) :: c1_x, c2_x, c2_y, c3_x, c3_d, c3_y
        real(dp) :: q_xx, q_xd, q_dd, q_xy, q_dy, q_yy

        coeff_xx = 0.0_dp
        coeff_xd = 0.0_dp
        coeff_dd = 0.0_dp
        if (.not. valid_inputs(profiles, radius_mid, 2)) then
            info = newcomb_invalid_input
            return
        end if
        call cylinder_geometry(profiles, radius_mid, geometry)
        if (.not. ieee_is_finite(geometry(8))) then
            info = newcomb_invalid_input
            return
        end if
        if (geometry(8) <= tiny(1.0_dp)) then
            info = newcomb_invalid_input
            return
        end if
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

        if (.not. ieee_is_finite(q_yy)) then
            info = newcomb_degenerate_reduction
            return
        end if
        if (q_yy <= tiny(1.0_dp)) then
            info = newcomb_degenerate_reduction
            return
        end if

        coeff_xx = (q_xx - q_xy**2 / q_yy) * abs(geometry(7))
        coeff_xd = (q_xd - q_xy * q_dy / q_yy) * abs(geometry(7))
        coeff_dd = (q_dd - q_dy**2 / q_yy) * abs(geometry(7))
        if (.not. ieee_is_finite(coeff_xx)) then
            info = newcomb_degenerate_reduction
            return
        end if
        if (.not. ieee_is_finite(coeff_xd)) then
            info = newcomb_degenerate_reduction
            return
        end if
        if (.not. ieee_is_finite(coeff_dd)) then
            info = newcomb_degenerate_reduction
            return
        end if
        info = newcomb_ok
    end subroutine newcomb_quadratic_coefficients

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

    subroutine add_element(stiffness, left_index, right_index, element)
        real(dp), intent(inout) :: stiffness(:, :)
        integer, intent(in) :: left_index, right_index
        real(dp), intent(in) :: element(2, 2)

        if (valid_index(left_index, size(stiffness, 1))) then
            stiffness(left_index, left_index) = &
                stiffness(left_index, left_index) + element(1, 1)
            if (valid_index(right_index, size(stiffness, 1))) then
                stiffness(left_index, right_index) = &
                    stiffness(left_index, right_index) + element(1, 2)
                stiffness(right_index, left_index) = &
                    stiffness(right_index, left_index) + element(2, 1)
            end if
        end if
        if (valid_index(right_index, size(stiffness, 1))) &
            stiffness(right_index, right_index) = &
            stiffness(right_index, right_index) + element(2, 2)
    end subroutine add_element

    pure function valid_index(index, extent) result(valid)
        integer, intent(in) :: index, extent
        logical :: valid

        valid = index >= 1 .and. index <= extent
    end function valid_index

    pure function valid_inputs(profiles, radius, n_radial) result(valid)
        type(cylinder_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: radius
        integer, intent(in) :: n_radial
        logical :: valid

        valid = .false.
        if (n_radial < 2) return
        if (.not. ieee_is_finite(radius)) return
        if (radius <= 0.0_dp) return
        if (.not. ieee_is_finite(profiles%length)) return
        if (profiles%length <= 0.0_dp) return
        if (.not. ieee_is_finite(profiles%b_axial)) return
        if (.not. ieee_is_finite(profiles%b_linear)) return
        if (.not. ieee_is_finite(profiles%b_cubic)) return
        valid = .true.
    end function valid_inputs

end module newcomb_limit
