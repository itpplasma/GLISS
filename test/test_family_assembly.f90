program test_family_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use block_tridiagonal, only: block_tridiagonal_t
    use family_assembly, only: assemble_family_blocks, &
        assemble_family_stiffness, family_negative_count, &
        iterate_family_eigenvalue, lowest_family_eigenvalue, &
        surface_geometry_t
    use newcomb_limit, only: cylinder_profiles_t, &
        lowest_eigenvalue_single_mode
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: n_radial = 100, n_theta = 64, n_zeta = 32
    type(cylinder_profiles_t) :: profiles
    type(surface_geometry_t), allocatable :: geometry(:)
    real(dp) :: reference, family_value, pair_value, step
    real(dp) :: iterated
    integer :: i, info, negatives

    profiles%length = 6.0_dp * pi
    profiles%b_axial = 1.0_dp
    profiles%b_linear = 0.3_dp
    profiles%b_cubic = 0.4_dp
    step = 0.5_dp / real(n_radial, dp)

    allocate (geometry(n_radial))
    do i = 1, n_radial
        call cylinder_surface((real(i, dp) - 0.5_dp) * step, geometry(i))
    end do

    call lowest_eigenvalue_single_mode(profiles, 1, 1, 0.5_dp, &
        n_radial, reference)
    call lowest_family_eigenvalue(geometry, [1], [1], step, &
        family_value, info)
    call require(info == 0, "family assembly failed")
    call require(abs(family_value - reference) < 1.0e-6_dp * &
        abs(reference), &
        "single-mode family disagrees with the 1D assembly")

    call lowest_family_eigenvalue(geometry, [1, 2], [1, 1], step, &
        pair_value, info)
    call require(info == 0, "two-mode family assembly failed")
    call require(abs(pair_value - family_value) < 1.0e-10_dp * &
        abs(family_value), &
        "decoupled modes change the lowest eigenvalue")

    call check_block_assembly(geometry, step)

    call iterate_family_eigenvalue(geometry, [1], [1], step, &
        1.05_dp * family_value, iterated, info)
    call require(info == 0, "inverse iteration failed")
    ! The dual-parity pair is degenerate up to quadrature round-off;
    ! the iterate may land anywhere in that split.
    call require(abs(iterated - family_value) < 1.0e-7_dp * &
        abs(family_value), &
        "inverse iteration disagrees with the dense solve")

    call family_negative_count(geometry, [1], [1], step, 0.0_dp, &
        negatives, info)
    call require(info == 0, "inertia count failed")
    call require(negatives >= 1, &
        "the unstable mode is invisible to the inertia count")
    call family_negative_count(geometry, [1], [1], step, &
        2.0_dp * family_value, negatives, info)
    call require(info == 0, "shifted inertia count failed")
    call require(negatives == 0, &
        "the inertia count below the lowest eigenvalue is not zero")
    write (*, "(a)") "PASS"

contains

    subroutine check_block_assembly(geometry, step)
        type(surface_geometry_t), intent(in) :: geometry(:)
        real(dp), intent(in) :: step
        type(block_tridiagonal_t) :: blocks
        real(dp), allocatable :: stiffness(:, :), dense_block(:, :)
        real(dp) :: deviation
        integer :: trials, nodes, info, i

        call assemble_family_stiffness(geometry, [1], [1], step, &
            stiffness, info)
        call require(info == 0, "dense assembly failed")
        call assemble_family_blocks(geometry, [1], [1], step, blocks, &
            info)
        call require(info == 0, "block assembly failed")
        trials = size(blocks%diag, 1)
        nodes = size(blocks%diag, 3)
        deviation = 0.0_dp
        allocate (dense_block(trials, trials))
        do i = 1, nodes
            dense_block = stiffness(trials * (i - 1) + 1:trials * i, &
                trials * (i - 1) + 1:trials * i)
            deviation = max(deviation, maxval(abs(dense_block &
                - blocks%diag(:, :, i))))
            if (i < nodes) then
                dense_block = stiffness(trials * i + 1:trials * (i + 1), &
                    trials * (i - 1) + 1:trials * i)
                deviation = max(deviation, maxval(abs(dense_block &
                    - blocks%off(:, :, i))))
            end if
        end do
        call require(deviation == 0.0_dp, &
            "block assembly differs from the dense assembly")
    end subroutine check_block_assembly

    subroutine cylinder_surface(radius, surface)
        real(dp), intent(in) :: radius
        type(surface_geometry_t), intent(out) :: surface
        real(dp) :: b_theta, b_theta_slope, fields(13), drive
        real(dp) :: length

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
        surface%fields = spread(spread(fields, 1, n_theta), 2, n_zeta)
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

end program test_family_assembly
