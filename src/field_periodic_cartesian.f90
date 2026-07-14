module field_periodic_cartesian
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cartesian_harmonic_spline, only: cartesian_jet_grid_t
    implicit none
    private

    integer, parameter, public :: field_periodic_cartesian_ok = 0
    integer, parameter, public :: field_periodic_cartesian_invalid = -1

    public :: convert_field_periodic_jet

contains

    subroutine convert_field_periodic_jet(zeta_period, field_periods, &
            winding, jet, info)
        real(dp), intent(in) :: zeta_period(:)
        integer, intent(in) :: field_periods, winding
        type(cartesian_jet_grid_t), intent(inout) :: jet
        integer, intent(out) :: info
        real(dp) :: angle, rate
        integer :: j, k

        info = field_periodic_cartesian_invalid
        if (field_periods < 1 .or. size(zeta_period) < 1) return
        if (.not. all(ieee_is_finite(zeta_period))) return
        if (.not. jet_has_shape(jet, size(zeta_period))) return
        rate = 2.0_dp * acos(-1.0_dp) * real(winding, dp) &
            / real(field_periods, dp)
        do k = 1, size(zeta_period)
            angle = rate * zeta_period(k)
            do j = 1, size(jet%value, 1)
                call convert_point(jet, j, k, angle, rate)
            end do
        end do
        if (.not. jet_is_finite(jet)) return
        info = field_periodic_cartesian_ok
    end subroutine convert_field_periodic_jet

    subroutine convert_point(jet, j, k, angle, rate)
        type(cartesian_jet_grid_t), intent(inout) :: jet
        integer, intent(in) :: j, k
        real(dp), intent(in) :: angle, rate
        real(dp) :: value(3), radial(3), poloidal(3), toroidal(3)
        real(dp) :: radial_radial(3), radial_poloidal(3)
        real(dp) :: radial_toroidal(3), poloidal_poloidal(3)
        real(dp) :: poloidal_toroidal(3), toroidal_toroidal(3)

        value = jet%value(j, k, :)
        radial = jet%radial(j, k, :)
        poloidal = jet%poloidal(j, k, :)
        toroidal = jet%toroidal(j, k, :)
        radial_radial = jet%radial_radial(j, k, :)
        radial_poloidal = jet%radial_poloidal(j, k, :)
        radial_toroidal = jet%radial_toroidal(j, k, :)
        poloidal_poloidal = jet%poloidal_poloidal(j, k, :)
        poloidal_toroidal = jet%poloidal_toroidal(j, k, :)
        toroidal_toroidal = jet%toroidal_toroidal(j, k, :)
        jet%value(j, k, :) = rotate(value, angle)
        jet%radial(j, k, :) = rotate(radial, angle)
        jet%poloidal(j, k, :) = rotate(poloidal, angle)
        jet%toroidal(j, k, :) = rotate(toroidal &
            + rate * generator(value), angle)
        jet%radial_radial(j, k, :) = rotate(radial_radial, angle)
        jet%radial_poloidal(j, k, :) = rotate(radial_poloidal, angle)
        jet%radial_toroidal(j, k, :) = rotate( &
            radial_toroidal &
            + rate * generator(radial), angle)
        jet%poloidal_poloidal(j, k, :) = rotate(poloidal_poloidal, angle)
        jet%poloidal_toroidal(j, k, :) = rotate( &
            poloidal_toroidal &
            + rate * generator(poloidal), angle)
        jet%toroidal_toroidal(j, k, :) = rotate( &
            toroidal_toroidal &
            + 2.0_dp * rate * generator(toroidal) &
            + rate**2 * generator_squared(value), angle)
    end subroutine convert_point

    pure function rotate(vector, angle) result(rotated)
        real(dp), intent(in) :: vector(3), angle
        real(dp) :: rotated(3)
        real(dp) :: cosine, sine

        cosine = cos(angle)
        sine = sin(angle)
        rotated = [cosine * vector(1) - sine * vector(2), &
            sine * vector(1) + cosine * vector(2), vector(3)]
    end function rotate

    pure function generator(vector) result(generated)
        real(dp), intent(in) :: vector(3)
        real(dp) :: generated(3)

        generated = [-vector(2), vector(1), 0.0_dp]
    end function generator

    pure function generator_squared(vector) result(generated)
        real(dp), intent(in) :: vector(3)
        real(dp) :: generated(3)

        generated = [-vector(1), -vector(2), 0.0_dp]
    end function generator_squared

    function jet_has_shape(jet, n_zeta) result(valid)
        type(cartesian_jet_grid_t), intent(in) :: jet
        integer, intent(in) :: n_zeta
        logical :: valid
        integer :: expected(3)

        valid = jet_is_allocated(jet)
        if (.not. valid) return
        expected = [size(jet%value, 1), n_zeta, 3]
        valid = all(shape(jet%value) == expected) &
            .and. all(shape(jet%radial) == expected) &
            .and. all(shape(jet%poloidal) == expected) &
            .and. all(shape(jet%toroidal) == expected) &
            .and. all(shape(jet%radial_radial) == expected) &
            .and. all(shape(jet%radial_poloidal) == expected) &
            .and. all(shape(jet%radial_toroidal) == expected) &
            .and. all(shape(jet%poloidal_poloidal) == expected) &
            .and. all(shape(jet%poloidal_toroidal) == expected) &
            .and. all(shape(jet%toroidal_toroidal) == expected)
        if (valid) valid = jet_is_finite(jet)
    end function jet_has_shape

    function jet_is_allocated(jet) result(allocated_all)
        type(cartesian_jet_grid_t), intent(in) :: jet
        logical :: allocated_all

        allocated_all = allocated(jet%value) .and. allocated(jet%radial) &
            .and. allocated(jet%poloidal) .and. allocated(jet%toroidal) &
            .and. allocated(jet%radial_radial) &
            .and. allocated(jet%radial_poloidal) &
            .and. allocated(jet%radial_toroidal) &
            .and. allocated(jet%poloidal_poloidal) &
            .and. allocated(jet%poloidal_toroidal) &
            .and. allocated(jet%toroidal_toroidal)
    end function jet_is_allocated

    function jet_is_finite(jet) result(finite)
        type(cartesian_jet_grid_t), intent(in) :: jet
        logical :: finite

        finite = all(ieee_is_finite(jet%value)) &
            .and. all(ieee_is_finite(jet%radial)) &
            .and. all(ieee_is_finite(jet%poloidal)) &
            .and. all(ieee_is_finite(jet%toroidal)) &
            .and. all(ieee_is_finite(jet%radial_radial)) &
            .and. all(ieee_is_finite(jet%radial_poloidal)) &
            .and. all(ieee_is_finite(jet%radial_toroidal)) &
            .and. all(ieee_is_finite(jet%poloidal_poloidal)) &
            .and. all(ieee_is_finite(jet%poloidal_toroidal)) &
            .and. all(ieee_is_finite(jet%toroidal_toroidal))
    end function jet_is_finite

end module field_periodic_cartesian
