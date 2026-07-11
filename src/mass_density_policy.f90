module mass_density_policy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: mass_density_ok = 0
    integer, parameter, public :: mass_density_invalid = -1
    integer, parameter, public :: mass_density_outside = -2

    type, public :: mass_density_profile_t
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: kilograms_per_cubic_metre(:)
    end type mass_density_profile_t

    public :: evaluate_mass_density
    public :: validate_mass_density_profile

contains

    pure subroutine validate_mass_density_profile(profile, info)
        type(mass_density_profile_t), intent(in) :: profile
        integer, intent(out) :: info

        info = mass_density_invalid
        if (.not. allocated(profile%s)) return
        if (.not. allocated(profile%kilograms_per_cubic_metre)) return
        if (size(profile%s) < 2) return
        if (size(profile%kilograms_per_cubic_metre) /= size(profile%s)) return
        if (.not. all(ieee_is_finite(profile%s))) return
        if (.not. all(ieee_is_finite( &
            profile%kilograms_per_cubic_metre))) return
        if (profile%s(1) < 0.0_dp .or. profile%s(size(profile%s)) > 1.0_dp) &
            return
        if (any(profile%s(2:) <= profile%s(:size(profile%s) - 1))) return
        if (any(profile%kilograms_per_cubic_metre <= 0.0_dp)) return
        info = mass_density_ok
    end subroutine validate_mass_density_profile

    pure subroutine evaluate_mass_density(profile, radial_coordinate, &
            density_kg_m3, info)
        type(mass_density_profile_t), intent(in) :: profile
        real(dp), intent(in) :: radial_coordinate
        real(dp), intent(out) :: density_kg_m3
        integer, intent(out) :: info
        real(dp) :: fraction
        integer :: interval

        density_kg_m3 = 0.0_dp
        call validate_mass_density_profile(profile, info)
        if (info /= mass_density_ok) return
        info = mass_density_outside
        if (.not. ieee_is_finite(radial_coordinate)) return
        if (radial_coordinate < profile%s(1)) return
        if (radial_coordinate > profile%s(size(profile%s))) return
        if (radial_coordinate == profile%s(size(profile%s))) then
            density_kg_m3 = &
                profile%kilograms_per_cubic_metre(size(profile%s))
            info = mass_density_ok
            return
        end if
        do interval = 1, size(profile%s) - 1
            if (radial_coordinate < profile%s(interval + 1)) exit
        end do
        fraction = (radial_coordinate - profile%s(interval)) &
            / (profile%s(interval + 1) - profile%s(interval))
        density_kg_m3 = (1.0_dp - fraction) &
            * profile%kilograms_per_cubic_metre(interval) &
            + fraction * profile%kilograms_per_cubic_metre(interval + 1)
        info = mass_density_ok
    end subroutine evaluate_mass_density

end module mass_density_policy
