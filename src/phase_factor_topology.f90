module phase_factor_topology
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: phase_factor_ok = 0
    integer, parameter, public :: phase_factor_invalid = -1
    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    type, public :: phase_factor_table_t
        integer :: field_periods = 0
        integer :: family_index = 0
        integer :: base_poloidal = 0
        integer, allocatable :: envelope_poloidal(:)
        integer, allocatable :: envelope_toroidal(:)
        integer, allocatable :: base_sign(:)
    end type phase_factor_table_t

    public :: build_phase_factor_table
    public :: evaluate_phase_factor

contains

    pure subroutine build_phase_factor_table(field_periods, family_index, &
            base_poloidal, mode_m, mode_n, table, info)
        integer, intent(in) :: field_periods, family_index, base_poloidal
        integer, intent(in) :: mode_m(:), mode_n(:)
        type(phase_factor_table_t), intent(out) :: table
        integer, intent(out) :: info
        integer :: mode, residue, sign

        info = phase_factor_invalid
        if (field_periods < 2) return
        if (family_index < 1 .or. family_index > field_periods / 2) return
        if (2 * family_index == field_periods) return
        if (size(mode_m) == 0 .or. size(mode_m) /= size(mode_n)) return
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
            else
                return
            end if
            table%base_sign(mode) = sign
        end do
        info = phase_factor_ok
    end subroutine build_phase_factor_table

    pure subroutine evaluate_phase_factor(table, mode, theta, zeta, &
            cosine, sine)
        type(phase_factor_table_t), intent(in) :: table
        integer, intent(in) :: mode
        real(dp), intent(in) :: theta, zeta
        real(dp), intent(out) :: cosine, sine
        real(dp) :: base_phase, envelope_phase
        real(dp) :: base_cosine, base_sine, envelope_cosine, envelope_sine
        real(dp) :: sign

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
    end subroutine evaluate_phase_factor

end module phase_factor_topology
