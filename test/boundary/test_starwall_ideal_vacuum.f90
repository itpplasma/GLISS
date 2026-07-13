program test_starwall_ideal_vacuum
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use starwall_ideal_vacuum, only: assemble_starwall_ideal_vacuum, &
        starwall_degenerate_surface, starwall_invalid_input, starwall_ok, &
        starwall_surfaces_not_nested
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)

    call test_coaxial_cylinder_limit()
    call test_wall_ordering()
    call test_invalid_surfaces()
    write (*, "(a)") "PASS"

contains

    subroutine test_coaxial_cylinder_limit()
        integer, parameter :: nu = 8, nv = 30
        real(dp), parameter :: major_radius = 12.0_dp
        real(dp), parameter :: plasma_radius = 1.0_dp
        real(dp), parameter :: wall_radius = 1.5_dp
        real(dp), allocatable :: plasma(:, :, :), wall(:, :, :)
        real(dp), allocatable :: stiffness(:, :), response(:, :), displacement(:)
        real(dp) :: energy, expected, fp, length, relative_error
        integer :: harmonic, info

        call torus(major_radius, plasma_radius, nu, nv, plasma)
        call torus(major_radius, wall_radius, nu, nv, wall)
        length = 2.0_dp * pi * major_radius
        fp = length * plasma_radius
        call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
            response, info, wall)
        call require(info == starwall_ok, "valid nested tori were rejected")
        call require(all(shape(stiffness) == [nu * nv, nu * nv]), &
            "vacuum stiffness has the wrong shape")
        call require(size(response, 2) == nu * nv, &
            "vacuum response has the wrong displacement dimension")
        call require(maxval(abs(stiffness - transpose(stiffness))) < &
            1.0e-11_dp * maxval(abs(stiffness)), &
            "vacuum stiffness is not symmetric")

        allocate (displacement(nu * nv))
        displacement = 1.0_dp
        energy = 0.5_dp * dot_product(displacement, &
            matmul(stiffness, displacement))
        expected = pi * fp**2 / (length * log(wall_radius / plasma_radius))
        relative_error = abs(energy / expected - 1.0_dp)
        call require(relative_error < 0.04_dp, &
            "Lust-Martensen energy misses the cylindrical limit")

        harmonic = 2
        call cosine_mode(nu, nv, harmonic, displacement)
        energy = 0.5_dp * dot_product(displacement, &
            matmul(stiffness, displacement))
        expected = pi * real(harmonic, dp) * fp**2 / (2.0_dp * length) &
            * (wall_radius**(2 * harmonic) + plasma_radius**(2 * harmonic)) &
            / (wall_radius**(2 * harmonic) - plasma_radius**(2 * harmonic))
        relative_error = abs(energy / expected - 1.0_dp)
        call require(relative_error < 0.04_dp, &
            "periodic vacuum energy misses the cylindrical limit")
    end subroutine test_coaxial_cylinder_limit

    subroutine test_wall_ordering()
        integer, parameter :: nu = 6, nv = 8
        real(dp), allocatable :: plasma(:, :, :), near_wall(:, :, :)
        real(dp), allocatable :: far_wall(:, :, :), stiffness(:, :), response(:, :)
        real(dp), allocatable :: displacement(:)
        real(dp) :: far_energy, fp, near_energy, open_energy
        integer :: info

        call torus(8.0_dp, 1.0_dp, nu, nv, plasma)
        call torus(8.0_dp, 1.3_dp, nu, nv, near_wall)
        call torus(8.0_dp, 2.0_dp, nu, nv, far_wall)
        fp = 16.0_dp * pi
        allocate (displacement(nu * nv))
        call cosine_mode(nu, nv, 1, displacement)
        call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
            response, info, near_wall)
        call require(info == starwall_ok, "near-wall assembly failed")
        near_energy = 0.5_dp * dot_product(displacement, &
            matmul(stiffness, displacement))
        call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
            response, info, far_wall)
        call require(info == starwall_ok, "far-wall assembly failed")
        far_energy = 0.5_dp * dot_product(displacement, &
            matmul(stiffness, displacement))
        call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
            response, info)
        call require(info == starwall_ok, "open-vacuum assembly failed")
        open_energy = 0.5_dp * dot_product(displacement, &
            matmul(stiffness, displacement))
        call require(near_energy > far_energy .and. &
            far_energy > open_energy .and. open_energy > 0.0_dp, &
            "conducting-wall ordering is not stabilizing")
    end subroutine test_wall_ordering

    subroutine test_invalid_surfaces()
        real(dp), allocatable :: plasma(:, :, :), wall(:, :, :)
        real(dp), allocatable :: stiffness(:, :), response(:, :)
        integer :: info

        call torus(8.0_dp, 1.0_dp, 6, 8, plasma)
        call torus(8.0_dp, 0.5_dp, 6, 8, wall)
        call assemble_starwall_ideal_vacuum(plasma, 1.0_dp, 0.0_dp, stiffness, &
            response, info, wall)
        call require(info == starwall_surfaces_not_nested, &
            "an inner wall was accepted as an enclosing wall")
        plasma(1, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call assemble_starwall_ideal_vacuum(plasma, 1.0_dp, 0.0_dp, stiffness, &
            response, info)
        call require(info == starwall_invalid_input, &
            "a nonfinite plasma surface was accepted")
        plasma = 0.0_dp
        call assemble_starwall_ideal_vacuum(plasma, 1.0_dp, 0.0_dp, stiffness, &
            response, info)
        call require(info == starwall_degenerate_surface, &
            "a degenerate plasma surface was accepted")
    end subroutine test_invalid_surfaces

    subroutine torus(major_radius, minor_radius, nu, nv, surface)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: nu, nv
        real(dp), allocatable, intent(out) :: surface(:, :, :)
        real(dp) :: phi, radius, theta
        integer :: i, k

        allocate (surface(3, nu, nv))
        do k = 1, nv
            phi = 2.0_dp * pi * real(k - 1, dp) / real(nv, dp)
            do i = 1, nu
                theta = 2.0_dp * pi * real(i - 1, dp) / real(nu, dp)
                radius = major_radius + minor_radius * cos(theta)
                surface(:, i, k) = [radius * cos(phi), radius * sin(phi), &
                    minor_radius * sin(theta)]
            end do
        end do
    end subroutine torus

    subroutine cosine_mode(nu, nv, harmonic, values)
        integer, intent(in) :: nu, nv, harmonic
        real(dp), intent(out) :: values(nu * nv)
        integer :: i, k

        do k = 1, nv
            do i = 1, nu
                values(i + nu * (k - 1)) = cos(2.0_dp * pi &
                    * real(harmonic * (i - 1), dp) / real(nu, dp))
            end do
        end do
    end subroutine cosine_mode

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_starwall_ideal_vacuum
