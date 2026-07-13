program gliss_starwall_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use starwall_cylinder_limit, only: build_circular_torus, &
        coaxial_cylinder_energy, cylinder_gate_ok, poloidal_cosine
    use starwall_ideal_vacuum, only: assemble_starwall_ideal_vacuum, starwall_ok
    implicit none

    integer, parameter :: nu = 8, nv = 30
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major_radius = 12.0_dp
    real(dp), parameter :: plasma_radius = 1.0_dp
    real(dp), allocatable :: plasma(:, :, :)
    real(dp) :: fp, length
    integer :: info

    if (command_argument_count() /= 0) call fail("no arguments are accepted")
    call build_circular_torus(major_radius, plasma_radius, nu, nv, plasma, info)
    if (info /= cylinder_gate_ok) call fail("plasma fixture construction failed")
    length = 2.0_dp * pi * major_radius
    fp = length * plasma_radius
    write (*, "(a)") "case,gate,wall_radius,harmonic,numerical_energy," // &
        "cylinder_reference_energy,relative_error,symmetry_defect"
    call report_case("finite_constant", "analytic", 1.5_dp, 0)
    call report_case("finite_harmonic", "analytic", 1.5_dp, 2)
    call report_case("near_wall", "ordering", 1.3_dp, 1)
    call report_case("far_wall", "ordering", 2.0_dp, 1)
    call report_case("open_vacuum", "ordering", 0.0_dp, 1)

contains

    subroutine report_case(name, gate, wall_radius, harmonic)
        character(len=*), intent(in) :: name, gate
        real(dp), intent(in) :: wall_radius
        integer, intent(in) :: harmonic
        real(dp), allocatable :: displacement(:), response(:, :)
        real(dp), allocatable :: stiffness(:, :), wall(:, :, :)
        real(dp) :: analytic, energy, error, symmetry
        integer :: status

        call poloidal_cosine(nu, nv, harmonic, displacement, status)
        if (status /= cylinder_gate_ok) call fail("mode construction failed")
        if (wall_radius > 0.0_dp) then
            call build_circular_torus(major_radius, wall_radius, nu, nv, &
                wall, status)
            if (status /= cylinder_gate_ok) call fail("wall construction failed")
            call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
                response, status, wall)
        else
            call assemble_starwall_ideal_vacuum(plasma, fp, 0.0_dp, stiffness, &
                response, status)
        end if
        if (status /= starwall_ok) call fail("vacuum assembly failed")
        energy = 0.5_dp * dot_product(displacement, &
            matmul(stiffness, displacement))
        analytic = coaxial_cylinder_energy(fp, length, plasma_radius, &
            wall_radius, harmonic)
        if (analytic <= 0.0_dp) call fail("analytic energy failed")
        error = abs(energy / analytic - 1.0_dp)
        symmetry = maxval(abs(stiffness - transpose(stiffness)))
        write (*, "(a,',',a,',',es24.16,',',i0,4(',',es24.16))") name, &
            gate, wall_radius, harmonic, energy, analytic, error, symmetry
    end subroutine report_case

    subroutine fail(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_starwall_diagnostic: " // message
        stop 2
    end subroutine fail

end program gliss_starwall_diagnostic
