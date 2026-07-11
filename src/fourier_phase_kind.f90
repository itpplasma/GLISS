module fourier_phase_kind
    implicit none
    private

    integer, parameter, public :: phase_cosine = 1
    integer, parameter, public :: phase_sine = 2

    public :: opposite_phase_kind
    public :: valid_phase_kind

contains

    pure function opposite_phase_kind(phase) result(opposite)
        integer, intent(in) :: phase
        integer :: opposite

        opposite = 0
        if (phase == phase_cosine) opposite = phase_sine
        if (phase == phase_sine) opposite = phase_cosine
    end function opposite_phase_kind

    pure function valid_phase_kind(phase) result(valid)
        integer, intent(in) :: phase
        logical :: valid

        valid = phase == phase_cosine .or. phase == phase_sine
    end function valid_phase_kind

end module fourier_phase_kind
