program test_eigenvalue_tracking
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use eigenvalue_tracking, only: certified_lowest_eigenvalue
    use family_assembly, only: family_assembly_options_t, &
        lowest_family_eigenvalue, surface_geometry_t
    use newcomb_limit, only: cylinder_profiles_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: n_radial = 100, n_theta = 64, n_zeta = 32
    type(cylinder_profiles_t) :: profiles
    type(surface_geometry_t), allocatable :: geometry(:)
    real(dp) :: step
    integer :: i

    profiles%length = 6.0_dp * pi
    profiles%b_axial = 1.0_dp
    profiles%b_linear = 0.3_dp
    profiles%b_cubic = 0.4_dp
    step = 0.5_dp / real(n_radial, dp)

    allocate (geometry(n_radial))
    do i = 1, n_radial
        call cylinder_surface((real(i, dp) - 0.5_dp) * step, geometry(i))
    end do

    call check_tracks_dense(geometry, [1], [1], step, 0)
    call check_tracks_dense(geometry, [1, 2], [1, 1], step, 0)
    call check_tracks_dense(geometry, [2], [1], step, 0)
    call check_tracks_dense(geometry, [1], [1], step, 1)
    call check_tracks_dense(geometry, [1], [1], step, 2)
    call check_tracks_dense(geometry, [1], [1], step, 2, [0.25_dp])
    write (*, "(a)") "PASS"

contains

    subroutine check_tracks_dense(geometry, mode_m, mode_n, step, &
            selector, normal_stored_power)
        type(surface_geometry_t), intent(in) :: geometry(:)
        integer, intent(in) :: mode_m(:), mode_n(:), selector
        real(dp), intent(in) :: step
        real(dp), intent(in), optional :: normal_stored_power(:)
        real(dp) :: dense, certified, residual, width, tolerance
        integer :: info
        type(family_assembly_options_t) :: options

        options%parity_class = selector

        call lowest_family_eigenvalue(geometry, mode_m, mode_n, step, &
            dense, info, options, normal_stored_power)
        call require(info == 0, "dense reference solve failed")
        call certified_lowest_eigenvalue(geometry, mode_m, mode_n, &
            step, certified, width, info, options, normal_stored_power, &
            residual)
        call require(info == 0, "certified tracking failed")
        call require(residual <= 1.0e-10_dp, &
            "certified eigenpair residual is not tight")
        call require(width <= 1.0e-8_dp * max(1.0_dp, abs(certified)), &
            "certificate window is not tight")
        tolerance = width + 1.0e-7_dp * abs(dense)
        call require(abs(certified - dense) <= tolerance, &
            "certified eigenvalue disagrees with the dense solve")
    end subroutine check_tracks_dense

    subroutine cylinder_surface(radius, surface)
        real(dp), intent(in) :: radius
        type(surface_geometry_t), intent(out) :: surface
        real(dp) :: b_theta, b_theta_slope, fields(13), drive
        real(dp) :: length
        integer :: component

        length = profiles%length
        b_theta = profiles%b_linear * radius + profiles%b_cubic &
            * radius**3
        b_theta_slope = profiles%b_linear + 3.0_dp * profiles%b_cubic &
            * radius**2
        fields(1) = 2.0_dp * pi * radius * profiles%b_axial
        fields(2) = length * b_theta
        fields(3) = 2.0_dp * pi * profiles%b_axial
        fields(4) = length * b_theta_slope
        fields(5) = length * profiles%b_axial
        fields(6) = 2.0_dp * pi * radius * b_theta
        fields(7) = 2.0_dp * pi * length * radius
        fields(8) = sqrt(b_theta**2 + profiles%b_axial**2)
        fields(9) = 1.0_dp
        fields(10) = (b_theta_slope + b_theta / radius) &
            * profiles%b_axial
        fields(11) = -b_theta * (b_theta_slope + b_theta / radius)
        fields(12) = 0.0_dp
        fields(13) = 0.0_dp
        drive = 2.0_dp * b_theta * (b_theta_slope + b_theta / radius) &
            / radius
        allocate (surface%fields(n_theta, n_zeta, 13))
        allocate (surface%drive(n_theta, n_zeta))
        do component = 1, size(fields)
            surface%fields(:, :, component) = fields(component)
        end do
        surface%drive = drive
    end subroutine cylinder_surface

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_eigenvalue_tracking
