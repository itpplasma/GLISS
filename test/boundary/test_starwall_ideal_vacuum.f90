program test_starwall_ideal_vacuum
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use starwall_cylinder_limit, only: build_circular_torus, &
        coaxial_cylinder_energy, cylinder_gate_ok, poloidal_cosine
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

        call build_circular_torus(major_radius, plasma_radius, nu, nv, plasma, info)
        call require(info == cylinder_gate_ok, "plasma fixture construction failed")
        call build_circular_torus(major_radius, wall_radius, nu, nv, wall, info)
        call require(info == cylinder_gate_ok, "wall fixture construction failed")
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
        energy = quadratic_energy(stiffness, displacement)
        expected = coaxial_cylinder_energy(fp, length, plasma_radius, &
            wall_radius, 0)
        relative_error = abs(energy / expected - 1.0_dp)
        call require(relative_error < 0.04_dp, &
            "Lust-Martensen energy misses the cylindrical limit")

        harmonic = 2
        call poloidal_cosine(nu, nv, harmonic, displacement, info)
        call require(info == cylinder_gate_ok, "mode construction failed")
        energy = quadratic_energy(stiffness, displacement)
        expected = coaxial_cylinder_energy(fp, length, plasma_radius, &
            wall_radius, harmonic)
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

        call build_circular_torus(8.0_dp, 1.0_dp, nu, nv, plasma, info)
        call require(info == cylinder_gate_ok, "plasma fixture construction failed")
        call build_circular_torus(8.0_dp, 1.3_dp, nu, nv, near_wall, info)
        call require(info == cylinder_gate_ok, "near-wall construction failed")
        call build_circular_torus(8.0_dp, 2.0_dp, nu, nv, far_wall, info)
        call require(info == cylinder_gate_ok, "far-wall construction failed")
        fp = 16.0_dp * pi
        allocate (displacement(nu * nv))
        call poloidal_cosine(nu, nv, 1, displacement, info)
        call require(info == cylinder_gate_ok, "mode construction failed")
        call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
            response, info, near_wall)
        call require(info == starwall_ok, "near-wall assembly failed")
        near_energy = quadratic_energy(stiffness, displacement)
        call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
            response, info, far_wall)
        call require(info == starwall_ok, "far-wall assembly failed")
        far_energy = quadratic_energy(stiffness, displacement)
        call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
            response, info)
        call require(info == starwall_ok, "open-vacuum assembly failed")
        open_energy = quadratic_energy(stiffness, displacement)
        call require(near_energy > far_energy .and. &
            far_energy > open_energy .and. open_energy > 0.0_dp, &
            "conducting-wall ordering is not stabilizing")
    end subroutine test_wall_ordering

    pure function quadratic_energy(matrix, vector) result(energy)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp) :: energy, image
        integer :: column, row

        energy = 0.0_dp
        do row = 1, size(matrix, 1)
            image = 0.0_dp
            do column = 1, size(matrix, 2)
                image = image + matrix(row, column) * vector(column)
            end do
            energy = energy + vector(row) * image
        end do
        energy = 0.5_dp * energy
    end function quadratic_energy

    subroutine test_invalid_surfaces()
        real(dp), allocatable :: plasma(:, :, :), wall(:, :, :)
        real(dp), allocatable :: stiffness(:, :), response(:, :)
        integer :: info

        call build_circular_torus(8.0_dp, 1.0_dp, 6, 8, plasma, info)
        call require(info == cylinder_gate_ok, "plasma fixture construction failed")
        call build_circular_torus(8.0_dp, 0.5_dp, 6, 8, wall, info)
        call require(info == cylinder_gate_ok, "inner-wall construction failed")
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

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program test_starwall_ideal_vacuum
