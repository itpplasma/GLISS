module radial_space_policy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: radial_space_ok = 0
    integer, parameter, public :: radial_space_invalid = -1
    integer, parameter, public :: form_identity = 1
    integer, parameter, public :: form_s_power_edge = 2
    integer, parameter, public :: constraint_eliminate = 1

    type, public :: radial_space_config_t
        integer :: normal_degree = 1
        integer :: tangential_degree = 0
        integer :: form_policy = form_identity
        integer :: axis_constraint = constraint_eliminate
        integer :: edge_constraint = constraint_eliminate
    end type radial_space_config_t

    public :: evaluate_normal_basis
    public :: validate_radial_space

contains

    pure subroutine validate_radial_space(config, info)
        type(radial_space_config_t), intent(in) :: config
        integer, intent(out) :: info

        info = radial_space_invalid
        if (config%normal_degree /= 1) return
        if (config%tangential_degree /= 0) return
        if (config%form_policy /= form_identity .and. &
            config%form_policy /= form_s_power_edge) return
        if (config%axis_constraint /= constraint_eliminate) return
        if (config%edge_constraint /= constraint_eliminate) return
        info = radial_space_ok
    end subroutine validate_radial_space

    pure subroutine evaluate_normal_basis(config, poloidal_mode, s, &
            radial_step, local_coordinate, values, derivatives, info)
        type(radial_space_config_t), intent(in) :: config
        integer, intent(in) :: poloidal_mode
        real(dp), intent(in) :: s, radial_step, local_coordinate
        real(dp), intent(out) :: values(2), derivatives(2)
        integer, intent(out) :: info
        real(dp) :: hats(2), slopes(2), form, form_slope, half_m

        call validate_radial_space(config, info)
        if (info /= radial_space_ok) return
        info = radial_space_invalid
        if (poloidal_mode < 0) return
        if (s < 0.0_dp .or. s > 1.0_dp) return
        if (radial_step <= 0.0_dp) return
        if (local_coordinate < 0.0_dp .or. local_coordinate > 1.0_dp) return
        hats = [1.0_dp - local_coordinate, local_coordinate]
        slopes = [-1.0_dp / radial_step, 1.0_dp / radial_step]
        form = 1.0_dp
        form_slope = 0.0_dp
        if (config%form_policy == form_s_power_edge) then
            half_m = 0.5_dp * real(poloidal_mode, dp)
            form = s**half_m * (1.0_dp - s)
            if (poloidal_mode == 0) then
                form_slope = -1.0_dp
            else
                if (s == 0.0_dp .and. poloidal_mode == 1) return
                form_slope = half_m * s**(half_m - 1.0_dp) &
                    * (1.0_dp - s) - s**half_m
            end if
        end if
        values = form * hats
        derivatives = form_slope * hats + form * slopes
        info = radial_space_ok
    end subroutine evaluate_normal_basis

end module radial_space_policy
