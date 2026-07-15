module local_mode_model
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use physical_constants, only: vacuum_permeability
    implicit none
    private

    public :: vacuum_permeability

    public :: assemble_local_mode
    public :: benchmark_mode_energy
    public :: mode_energy
    public :: rayleigh_quotient

contains

    pure subroutine assemble_local_mode(wave_vector, magnetic_field, pressure, &
            density, adiabatic_index, pressure_curvature_drive, normal, &
            stiffness, mass)
        real(dp), intent(in) :: wave_vector(3), magnetic_field(3), pressure
        real(dp), intent(in) :: density, adiabatic_index
        real(dp), intent(in) :: pressure_curvature_drive, normal(3)
        real(dp), intent(out) :: stiffness(3, 3), mass(3, 3)
        real(dp) :: induction(3, 3)
        real(dp) :: parallel_wave_number
        integer :: i, j

        parallel_wave_number = dot_product(magnetic_field, wave_vector)
        induction = 0.0_dp
        do i = 1, 3
            induction(i, i) = parallel_wave_number
        end do
        do j = 1, 3
            do i = 1, 3
                induction(i, j) = induction(i, j) &
                    - magnetic_field(i) * wave_vector(j)
            end do
        end do

        stiffness = matmul(transpose(induction), induction) / &
            vacuum_permeability
        do j = 1, 3
            do i = 1, 3
                stiffness(i, j) = stiffness(i, j) &
                    + adiabatic_index * pressure * wave_vector(i) &
                    * wave_vector(j)
                stiffness(i, j) = stiffness(i, j) &
                    - pressure_curvature_drive * normal(i) * normal(j)
            end do
        end do

        mass = 0.0_dp
        do i = 1, 3
            mass(i, i) = density
        end do
    end subroutine assemble_local_mode

    pure function mode_energy(stiffness, displacement) result(energy)
        real(dp), intent(in) :: stiffness(3, 3), displacement(3)
        real(dp) :: energy, image(3)

        image = matmul(stiffness, displacement)
        energy = 0.5_dp * dot_product(displacement, image)
    end function mode_energy

    pure function rayleigh_quotient(stiffness, mass, displacement) &
            result(omega_squared)
        real(dp), intent(in) :: stiffness(3, 3), mass(3, 3)
        real(dp), intent(in) :: displacement(3)
        real(dp) :: omega_squared, stiffness_image(3), mass_image(3)

        stiffness_image = matmul(stiffness, displacement)
        mass_image = matmul(mass, displacement)
        omega_squared = dot_product(displacement, stiffness_image) &
            / dot_product(displacement, mass_image)
    end function rayleigh_quotient

    pure function benchmark_mode_energy(drive) result(energy) &
            bind(c, name="gvstab_benchmark_mode_energy")
        real(c_double), intent(in), value :: drive
        real(c_double) :: energy
        real(dp) :: stiffness(3, 3), mass(3, 3)
        real(dp), parameter :: wave_vector(3) = [0.0_dp, 0.0_dp, 1.0_dp]
        real(dp), parameter :: magnetic_field(3) = [0.0_dp, 0.0_dp, 1.0e-3_dp]
        real(dp), parameter :: normal(3) = [1.0_dp, 0.0_dp, 0.0_dp]
        real(dp), parameter :: displacement(3) = [1.0_dp, 0.0_dp, 0.0_dp]

        call assemble_local_mode(wave_vector, magnetic_field, 1.0_dp, &
            2.0_dp, 5.0_dp / 3.0_dp, drive, normal, stiffness, mass)
        energy = mode_energy(stiffness, displacement)
    end function benchmark_mode_energy

end module local_mode_model
