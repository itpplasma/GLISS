module fixed_boundary_input_validation
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: valid_fixed_boundary_inputs

contains

    function valid_fixed_boundary_inputs(adiabatic_index, density_kg_m3, &
            zero_floor, mode_m, mode_n, radial_quadrature) result(valid)
        real(dp), intent(in) :: adiabatic_index, density_kg_m3, zero_floor
        integer, intent(in) :: mode_m(:), mode_n(:), radial_quadrature
        logical :: valid
        integer :: first, second

        valid = .false.
        if (.not. ieee_is_finite(adiabatic_index)) return
        if (.not. ieee_is_finite(density_kg_m3)) return
        if (.not. ieee_is_finite(zero_floor)) return
        if (adiabatic_index < 0.0_dp) return
        if (density_kg_m3 <= 0.0_dp .or. zero_floor <= 0.0_dp) return
        if (zero_floor > 0.125_dp * huge(zero_floor)) return
        if (size(mode_m) < 1 .or. size(mode_m) /= size(mode_n)) return
        if (radial_quadrature /= 1) return
        do first = 1, size(mode_m)
            if (mode_m(first) < 0) return
            if (mode_m(first) == 0 .and. mode_n(first) < 0) return
            do second = 1, first - 1
                if (mode_m(first) == mode_m(second) &
                    .and. mode_n(first) == mode_n(second)) return
            end do
        end do
        valid = .true.
    end function valid_fixed_boundary_inputs

end module fixed_boundary_input_validation
