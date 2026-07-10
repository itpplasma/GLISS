module mode_topology
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    type, public :: mode_family_t
        integer :: field_periods = 0
        integer :: family_index = 0
        integer, allocatable :: poloidal(:)
        integer, allocatable :: toroidal(:)
    end type mode_family_t

    public :: build_mode_family
    public :: family_count
    public :: modes_coupled
    public :: axis_form_function
    public :: axis_form_function_slope

contains

    pure function family_count(field_periods) result(count)
        integer, intent(in) :: field_periods
        integer :: count

        count = field_periods / 2
    end function family_count

    pure function modes_coupled(first, second, field_periods) &
            result(coupled)
        integer, intent(in) :: first, second, field_periods
        logical :: coupled

        coupled = mod(abs(first - second), field_periods) == 0 .or. &
            mod(abs(first + second), field_periods) == 0
    end function modes_coupled

    pure subroutine build_mode_family(field_periods, family_index, &
            poloidal_max, toroidal_max, family, info)
        integer, intent(in) :: field_periods, family_index
        integer, intent(in) :: poloidal_max, toroidal_max
        type(mode_family_t), intent(out) :: family
        integer, intent(out) :: info
        integer :: m, n, count

        info = -1
        if (field_periods < 2) return
        if (family_index < 1 .or. &
            family_index > family_count(field_periods)) return
        if (poloidal_max < 0 .or. toroidal_max < 0) return
        family%field_periods = field_periods
        family%family_index = family_index
        count = 0
        do m = 0, poloidal_max
            do n = -toroidal_max, toroidal_max
                if (.not. family_member(m, n, field_periods, &
                    family_index)) cycle
                count = count + 1
            end do
        end do
        allocate (family%poloidal(count), family%toroidal(count))
        count = 0
        do m = 0, poloidal_max
            do n = -toroidal_max, toroidal_max
                if (.not. family_member(m, n, field_periods, &
                    family_index)) cycle
                count = count + 1
                family%poloidal(count) = m
                family%toroidal(count) = n
            end do
        end do
        info = 0
    end subroutine build_mode_family

    pure function family_member(m, n, field_periods, family_index) &
            result(member)
        integer, intent(in) :: m, n, field_periods, family_index
        logical :: member
        integer :: residue

        residue = residue_of(n, field_periods)
        member = residue == modulo(family_index, field_periods) .or. &
            residue == modulo(-family_index, field_periods)
        if (m == 0 .and. n < 0) member = .false.
    end function family_member

    pure function residue_of(n, field_periods) result(residue)
        integer, intent(in) :: n, field_periods
        integer :: residue

        residue = modulo(n, field_periods)
    end function residue_of

    pure function axis_form_function(poloidal_mode, s) result(value)
        integer, intent(in) :: poloidal_mode
        real(dp), intent(in) :: s
        real(dp) :: value

        value = s**(0.5_dp * real(poloidal_mode, dp)) * (1.0_dp - s)
    end function axis_form_function

    pure function axis_form_function_slope(poloidal_mode, s) result(value)
        integer, intent(in) :: poloidal_mode
        real(dp), intent(in) :: s
        real(dp) :: value
        real(dp) :: half_m

        half_m = 0.5_dp * real(poloidal_mode, dp)
        if (poloidal_mode == 0) then
            value = -1.0_dp
        else
            value = half_m * s**(half_m - 1.0_dp) * (1.0_dp - s) &
                - s**half_m
        end if
    end function axis_form_function_slope

end module mode_topology
