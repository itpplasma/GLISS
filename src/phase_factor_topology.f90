module phase_factor_topology
    use fourier_phase_kind, only: phase_cosine, phase_sine
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: phase_factor_ok = 0
    integer, parameter, public :: phase_factor_invalid = -1
    public :: phase_cosine
    public :: phase_sine
    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    type, public :: phase_envelope_table_t
        integer :: field_periods = 0
        integer :: family_index = 0
        integer :: base_poloidal = 0
        integer, allocatable :: envelope_poloidal(:)
        integer, allocatable :: envelope_toroidal(:)
        integer, allocatable :: base_sign(:)
    end type phase_envelope_table_t

    public :: build_phase_envelope_table
    public :: evaluate_phase_envelope
    public :: phase_product_average
    public :: phase_product_coefficients

contains

    pure subroutine build_phase_envelope_table(field_periods, family_index, &
            base_poloidal, mode_m, mode_n, table, info)
        integer, intent(in) :: field_periods, family_index, base_poloidal
        integer, intent(in) :: mode_m(:), mode_n(:)
        type(phase_envelope_table_t), intent(out) :: table
        integer, intent(out) :: info
        integer :: mode, residue, sign

        info = phase_factor_invalid
        if (field_periods < 1) return
        if (family_index < 0 .or. family_index > field_periods / 2) return
        if (size(mode_m) == 0 .or. size(mode_m) /= size(mode_n)) return
        do mode = 1, size(mode_m)
            residue = modulo(mode_n(mode), field_periods)
            if (residue /= family_index .and. &
                residue /= field_periods - family_index) return
        end do
        table%field_periods = field_periods
        table%family_index = family_index
        table%base_poloidal = base_poloidal
        allocate (table%envelope_poloidal(size(mode_m)), &
            table%envelope_toroidal(size(mode_m)), &
            table%base_sign(size(mode_m)))
        do mode = 1, size(mode_m)
            residue = modulo(mode_n(mode), field_periods)
            if (residue == family_index) then
                sign = 1
                table%envelope_poloidal(mode) = mode_m(mode) &
                    - base_poloidal
                table%envelope_toroidal(mode) = &
                    (mode_n(mode) - family_index) / field_periods
            else if (residue == field_periods - family_index) then
                sign = -1
                table%envelope_poloidal(mode) = mode_m(mode) &
                    + base_poloidal
                table%envelope_toroidal(mode) = &
                    (mode_n(mode) + family_index) / field_periods
            end if
            table%base_sign(mode) = sign
        end do
        info = phase_factor_ok
    end subroutine build_phase_envelope_table

    pure subroutine evaluate_phase_envelope(table, mode, theta, zeta, &
            cosine, sine, info)
        type(phase_envelope_table_t), intent(in) :: table
        integer, intent(in) :: mode
        real(dp), intent(in) :: theta, zeta
        real(dp), intent(out) :: cosine, sine
        integer, intent(out) :: info
        real(dp) :: base_phase, envelope_phase
        real(dp) :: base_cosine, base_sine, envelope_cosine, envelope_sine
        real(dp) :: sign

        info = phase_factor_invalid
        cosine = ieee_value(cosine, ieee_quiet_nan)
        sine = ieee_value(sine, ieee_quiet_nan)
        if (table%field_periods < 1) return
        if (.not. allocated(table%envelope_poloidal)) return
        if (.not. allocated(table%envelope_toroidal)) return
        if (.not. allocated(table%base_sign)) return
        if (size(table%envelope_poloidal) /= &
            size(table%envelope_toroidal)) return
        if (size(table%envelope_poloidal) /= size(table%base_sign)) return
        if (mode < 1 .or. mode > size(table%envelope_poloidal)) return
        if (abs(table%base_sign(mode)) /= 1) return
        if (.not. ieee_is_finite(theta)) return
        if (.not. ieee_is_finite(zeta)) return

        base_phase = two_pi * (real(table%base_poloidal, dp) * theta &
            - real(table%family_index, dp) * zeta &
            / real(table%field_periods, dp))
        envelope_phase = two_pi &
            * (real(table%envelope_poloidal(mode), dp) * theta &
            - real(table%envelope_toroidal(mode), dp) * zeta)
        base_cosine = cos(base_phase)
        base_sine = sin(base_phase)
        envelope_cosine = cos(envelope_phase)
        envelope_sine = sin(envelope_phase)
        sign = real(table%base_sign(mode), dp)
        cosine = envelope_cosine * base_cosine &
            - sign * envelope_sine * base_sine
        sine = envelope_sine * base_cosine &
            + sign * envelope_cosine * base_sine
        info = phase_factor_ok
    end subroutine evaluate_phase_envelope

    pure function phase_product_average(first_kind, second_kind, &
            first_phase, second_phase, first_n, second_n, field_periods) &
            result(average)
        integer, intent(in) :: first_kind, second_kind
        real(dp), intent(in) :: first_phase, second_phase
        integer, intent(in) :: first_n, second_n, field_periods
        real(dp) :: average
        real(dp) :: products(2, 2)

        if (field_periods < 1) then
            average = ieee_value(average, ieee_quiet_nan)
            return
        end if
        if (first_kind < phase_cosine .or. first_kind > phase_sine) then
            average = ieee_value(average, ieee_quiet_nan)
            return
        end if
        if (second_kind < phase_cosine .or. second_kind > phase_sine) then
            average = ieee_value(average, ieee_quiet_nan)
            return
        end if
        call phase_product_coefficients(cos(first_phase), sin(first_phase), &
            cos(second_phase), sin(second_phase), first_n, second_n, &
            field_periods, products)
        average = products(first_kind, second_kind)
    end function phase_product_average

    pure subroutine phase_product_coefficients(first_cosine, first_sine, &
            second_cosine, second_sine, first_n, second_n, field_periods, &
            products)
        real(dp), intent(in) :: first_cosine, first_sine
        real(dp), intent(in) :: second_cosine, second_sine
        integer, intent(in) :: first_n, second_n, field_periods
        real(dp), intent(out) :: products(2, 2)
        real(dp) :: difference_cosine, difference_sine
        real(dp) :: sum_cosine, sum_sine

        if (field_periods < 1) then
            products = ieee_value(products, ieee_quiet_nan)
            return
        end if
        products = 0.0_dp
        difference_cosine = first_cosine * second_cosine &
            + first_sine * second_sine
        difference_sine = first_sine * second_cosine &
            - first_cosine * second_sine
        sum_cosine = first_cosine * second_cosine &
            - first_sine * second_sine
        sum_sine = first_sine * second_cosine &
            + first_cosine * second_sine
        if (modulo(first_n - second_n, field_periods) == 0) then
            products(phase_cosine, phase_cosine) = difference_cosine
            products(phase_cosine, phase_sine) = -difference_sine
            products(phase_sine, phase_cosine) = difference_sine
            products(phase_sine, phase_sine) = difference_cosine
        end if
        if (modulo(first_n + second_n, field_periods) == 0) then
            products(phase_cosine, phase_cosine) = &
                products(phase_cosine, phase_cosine) + sum_cosine
            products(phase_cosine, phase_sine) = &
                products(phase_cosine, phase_sine) + sum_sine
            products(phase_sine, phase_cosine) = &
                products(phase_sine, phase_cosine) + sum_sine
            products(phase_sine, phase_sine) = &
                products(phase_sine, phase_sine) - sum_cosine
        end if
        products = 0.5_dp * products
    end subroutine phase_product_coefficients

end module phase_factor_topology
