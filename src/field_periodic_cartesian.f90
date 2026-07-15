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
        real(dp) :: work(3)

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
        call rotate_store(value, angle, jet%value, j, k)
        call rotate_store(radial, angle, jet%radial, j, k)
        call rotate_store(poloidal, angle, jet%poloidal, j, k)
        work(1) = toroidal(1) - rate * value(2)
        work(2) = toroidal(2) + rate * value(1)
        work(3) = toroidal(3)
        call rotate_store(work, angle, jet%toroidal, j, k)
        call rotate_store(radial_radial, angle, jet%radial_radial, j, k)
        call rotate_store(radial_poloidal, angle, jet%radial_poloidal, j, k)
        work(1) = radial_toroidal(1) - rate * radial(2)
        work(2) = radial_toroidal(2) + rate * radial(1)
        work(3) = radial_toroidal(3)
        call rotate_store(work, angle, jet%radial_toroidal, j, k)
        call rotate_store(poloidal_poloidal, angle, &
            jet%poloidal_poloidal, j, k)
        work(1) = poloidal_toroidal(1) - rate * poloidal(2)
        work(2) = poloidal_toroidal(2) + rate * poloidal(1)
        work(3) = poloidal_toroidal(3)
        call rotate_store(work, angle, jet%poloidal_toroidal, j, k)
        work(1) = toroidal_toroidal(1) - 2.0_dp * rate * toroidal(2) &
            - rate**2 * value(1)
        work(2) = toroidal_toroidal(2) + 2.0_dp * rate * toroidal(1) &
            - rate**2 * value(2)
        work(3) = toroidal_toroidal(3)
        call rotate_store(work, angle, jet%toroidal_toroidal, j, k)
    end subroutine convert_point

    pure subroutine rotate_store(vector, angle, field, j, k)
        real(dp), intent(in) :: vector(3), angle
        real(dp), intent(inout) :: field(:, :, :)
        integer, intent(in) :: j, k
        real(dp) :: cosine, sine

        cosine = cos(angle)
        sine = sin(angle)
        field(j, k, 1) = cosine * vector(1) - sine * vector(2)
        field(j, k, 2) = sine * vector(1) + cosine * vector(2)
        field(j, k, 3) = vector(3)
    end subroutine rotate_store

    function jet_has_shape(jet, n_zeta) result(valid)
        type(cartesian_jet_grid_t), intent(in) :: jet
        integer, intent(in) :: n_zeta
        logical :: valid
        integer :: expected(3)

        valid = jet_is_allocated(jet)
        if (.not. valid) return
        expected(1) = size(jet%value, 1)
        expected(2) = n_zeta
        expected(3) = 3
        valid = tensor_has_shape(jet%value, expected) &
            .and. tensor_has_shape(jet%radial, expected) &
            .and. tensor_has_shape(jet%poloidal, expected) &
            .and. tensor_has_shape(jet%toroidal, expected) &
            .and. tensor_has_shape(jet%radial_radial, expected) &
            .and. tensor_has_shape(jet%radial_poloidal, expected) &
            .and. tensor_has_shape(jet%radial_toroidal, expected) &
            .and. tensor_has_shape(jet%poloidal_poloidal, expected) &
            .and. tensor_has_shape(jet%poloidal_toroidal, expected) &
            .and. tensor_has_shape(jet%toroidal_toroidal, expected)
        if (valid) valid = jet_is_finite(jet)
    end function jet_has_shape

    pure function tensor_has_shape(tensor, expected) result(valid)
        real(dp), intent(in) :: tensor(:, :, :)
        integer, intent(in) :: expected(3)
        logical :: valid

        valid = size(tensor, 1) == expected(1) &
            .and. size(tensor, 2) == expected(2) &
            .and. size(tensor, 3) == expected(3)
    end function tensor_has_shape

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
