module primitive_geometry_grid
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cartesian_harmonic_spline, only: cartesian_jet_grid_t
    use cartesian_primitive_geometry, only: build_primitive_geometry_point, &
        primitive_geometry_ok, primitive_geometry_point_t
    implicit none
    private

    integer, parameter, public :: primitive_geometry_grid_ok = 0
    integer, parameter, public :: primitive_geometry_grid_invalid = -1
    integer, parameter, public :: primitive_geometry_grid_allocation_error = -2

    type, public :: primitive_geometry_grid_t
        real(dp), allocatable :: metric(:, :, :, :)
        real(dp), allocatable :: signed_jacobian(:, :)
        real(dp), allocatable :: jacobian_s(:, :), jacobian_theta(:, :)
        real(dp), allocatable :: jacobian_zeta(:, :)
        real(dp), allocatable :: b_contravariant(:, :, :)
        real(dp), allocatable :: b_covariant(:, :, :), mod_b(:, :)
        real(dp), allocatable :: second_form(:, :, :, :)
    end type primitive_geometry_grid_t

    public :: build_primitive_geometry_grid

contains

    subroutine build_primitive_geometry_grid(jet, field_periods, &
            toroidal_flux_slope, poloidal_flux_slope, geometry, info)
        type(cartesian_jet_grid_t), intent(in) :: jet
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: toroidal_flux_slope, poloidal_flux_slope
        type(primitive_geometry_grid_t), intent(out) :: geometry
        integer, intent(out) :: info
        type(primitive_geometry_point_t) :: point
        integer :: allocation_status, j, k, point_info

        geometry = primitive_geometry_grid_t()
        info = primitive_geometry_grid_invalid
        if (field_periods < 1) return
        if (.not. ieee_is_finite(toroidal_flux_slope) &
            .or. .not. ieee_is_finite(poloidal_flux_slope)) return
        if (.not. jet_is_valid(jet)) return
        call allocate_geometry(geometry, size(jet%value, 1), &
            size(jet%value, 2), allocation_status)
        if (allocation_status /= 0) then
            geometry = primitive_geometry_grid_t()
            info = primitive_geometry_grid_allocation_error
            return
        end if
        do k = 1, size(jet%value, 2)
            do j = 1, size(jet%value, 1)
                call build_primitive_geometry_point(jet%radial(j, k, :), &
                    jet%poloidal(j, k, :), jet%toroidal(j, k, :), &
                    jet%radial_radial(j, k, :), &
                    jet%radial_poloidal(j, k, :), &
                    jet%radial_toroidal(j, k, :), &
                    jet%poloidal_poloidal(j, k, :), &
                    jet%poloidal_toroidal(j, k, :), &
                    jet%toroidal_toroidal(j, k, :), field_periods, &
                    toroidal_flux_slope, poloidal_flux_slope, point, &
                    point_info)
                if (point_info /= primitive_geometry_ok) then
                    geometry = primitive_geometry_grid_t()
                    return
                end if
                call store_point(point, j, k, geometry)
            end do
        end do
        info = primitive_geometry_grid_ok
    end subroutine build_primitive_geometry_grid

    subroutine allocate_geometry(geometry, n_theta, n_zeta, status)
        type(primitive_geometry_grid_t), intent(out) :: geometry
        integer, intent(in) :: n_theta, n_zeta
        integer, intent(out) :: status

        allocate (geometry%metric(n_theta, n_zeta, 3, 3), &
            geometry%signed_jacobian(n_theta, n_zeta), &
            geometry%jacobian_s(n_theta, n_zeta), &
            geometry%jacobian_theta(n_theta, n_zeta), &
            geometry%jacobian_zeta(n_theta, n_zeta), &
            geometry%b_contravariant(n_theta, n_zeta, 2), &
            geometry%b_covariant(n_theta, n_zeta, 2), &
            geometry%mod_b(n_theta, n_zeta), &
            geometry%second_form(n_theta, n_zeta, 2, 2), stat=status)
    end subroutine allocate_geometry

    subroutine store_point(point, j, k, geometry)
        type(primitive_geometry_point_t), intent(in) :: point
        integer, intent(in) :: j, k
        type(primitive_geometry_grid_t), intent(inout) :: geometry

        geometry%metric(j, k, :, :) = point%metric
        geometry%signed_jacobian(j, k) = point%signed_jacobian
        geometry%jacobian_s(j, k) = point%jacobian_s
        geometry%jacobian_theta(j, k) = point%jacobian_theta
        geometry%jacobian_zeta(j, k) = point%jacobian_zeta
        geometry%b_contravariant(j, k, :) = point%b_contravariant
        geometry%b_covariant(j, k, :) = point%b_covariant
        geometry%mod_b(j, k) = point%mod_b
        geometry%second_form(j, k, :, :) = point%second_form
    end subroutine store_point

    function jet_is_valid(jet) result(valid)
        type(cartesian_jet_grid_t), intent(in) :: jet
        logical :: valid
        integer :: expected(3)

        valid = allocated(jet%value)
        if (.not. valid) return
        expected = shape(jet%value)
        valid = expected(1) > 0 .and. expected(2) > 0 .and. expected(3) == 3
        if (.not. valid) return
        valid = array_is_valid(jet%radial, expected) &
            .and. array_is_valid(jet%poloidal, expected) &
            .and. array_is_valid(jet%toroidal, expected) &
            .and. array_is_valid(jet%radial_radial, expected) &
            .and. array_is_valid(jet%radial_poloidal, expected) &
            .and. array_is_valid(jet%radial_toroidal, expected) &
            .and. array_is_valid(jet%poloidal_poloidal, expected) &
            .and. array_is_valid(jet%poloidal_toroidal, expected) &
            .and. array_is_valid(jet%toroidal_toroidal, expected) &
            .and. all(ieee_is_finite(jet%value))
    end function jet_is_valid

    function array_is_valid(array, expected) result(valid)
        real(dp), allocatable, intent(in) :: array(:, :, :)
        integer, intent(in) :: expected(3)
        logical :: valid

        valid = allocated(array)
        if (.not. valid) return
        valid = all(shape(array) == expected) .and. all(ieee_is_finite(array))
    end function array_is_valid

end module primitive_geometry_grid
