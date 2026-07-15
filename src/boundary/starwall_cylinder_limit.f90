module starwall_cylinder_limit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: cylinder_gate_ok = 0
    integer, parameter, public :: cylinder_gate_invalid_input = 1
    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: build_circular_torus, coaxial_cylinder_energy, poloidal_cosine

contains

    subroutine build_circular_torus(major_radius, minor_radius, nu, nv, &
            surface, info)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: nu, nv
        real(dp), allocatable, intent(out) :: surface(:, :, :)
        integer, intent(out) :: info
        real(dp) :: phi, radius, theta
        integer :: i, k

        info = cylinder_gate_invalid_input
        if (.not. ieee_is_finite(major_radius) .or. &
            .not. ieee_is_finite(minor_radius)) return
        if (minor_radius <= 0.0_dp .or. major_radius <= minor_radius) return
        if (nu < 3 .or. nv < 3) return
        allocate (surface(3, nu, nv))
        do k = 1, nv
            phi = 2.0_dp * pi * real(k - 1, dp) / real(nv, dp)
            do i = 1, nu
                theta = 2.0_dp * pi * real(i - 1, dp) / real(nu, dp)
                radius = major_radius + minor_radius * cos(theta)
                surface(1, i, k) = radius * cos(phi)
                surface(2, i, k) = radius * sin(phi)
                surface(3, i, k) = minor_radius * sin(theta)
            end do
        end do
        info = cylinder_gate_ok
    end subroutine build_circular_torus

    subroutine poloidal_cosine(nu, nv, harmonic, values, info)
        integer, intent(in) :: nu, nv, harmonic
        real(dp), allocatable, intent(out) :: values(:)
        integer, intent(out) :: info
        integer :: i, k

        info = cylinder_gate_invalid_input
        if (nu < 3 .or. nv < 3 .or. harmonic < 0) return
        allocate (values(nu * nv))
        do k = 1, nv
            do i = 1, nu
                values(i + nu * (k - 1)) = cos(2.0_dp * pi &
                    * real(harmonic * (i - 1), dp) / real(nu, dp))
            end do
        end do
        info = cylinder_gate_ok
    end subroutine poloidal_cosine

    function coaxial_cylinder_energy(fp, length, plasma_radius, wall_radius, &
            harmonic) result(energy)
        real(dp), intent(in) :: fp, length, plasma_radius, wall_radius
        integer, intent(in) :: harmonic
        real(dp) :: energy, ratio

        energy = -1.0_dp
        if (.not. ieee_is_finite(fp)) return
        if (.not. ieee_is_finite(length)) return
        if (.not. ieee_is_finite(plasma_radius)) return
        if (.not. ieee_is_finite(wall_radius)) return
        if (length <= 0.0_dp .or. plasma_radius <= 0.0_dp) return
        if (harmonic < 0 .or. wall_radius < 0.0_dp) return
        if (wall_radius > 0.0_dp .and. wall_radius <= plasma_radius) return
        if (harmonic == 0) then
            if (wall_radius == 0.0_dp) return
            energy = pi * fp**2 / (length * log(wall_radius / plasma_radius))
            return
        end if
        if (wall_radius > 0.0_dp) then
            ratio = (wall_radius**(2 * harmonic) &
                + plasma_radius**(2 * harmonic)) &
                / (wall_radius**(2 * harmonic) &
                - plasma_radius**(2 * harmonic))
        else
            ratio = 1.0_dp
        end if
        energy = pi * real(harmonic, dp) * fp**2 * ratio / (2.0_dp * length)
    end function coaxial_cylinder_energy

end module starwall_cylinder_limit
