module trial_space_topology
    use fourier_phase_kind, only: &
        opposite_phase_kind, trial_phase_cosine => phase_cosine, &
        trial_phase_sine => phase_sine, valid_phase_kind
    implicit none
    private

    integer, parameter, public :: trial_topology_ok = 0
    integer, parameter, public :: trial_topology_invalid = -1
    integer, parameter, public :: trial_component_normal = 1
    integer, parameter, public :: trial_component_eta = 2
    integer, parameter, public :: trial_component_mu = 3
    public :: trial_phase_cosine
    public :: trial_phase_sine

    type, public :: trial_space_topology_t
        integer, allocatable :: poloidal(:)
        integer, allocatable :: toroidal(:)
        integer, allocatable :: parity(:)
        integer, allocatable :: normal_phase(:)
        integer, allocatable :: tangential_phase(:)
        logical, allocatable :: active(:, :)
    end type trial_space_topology_t

    public :: build_trial_space_topology

contains

    pure subroutine build_trial_space_topology(poloidal, toroidal, parity, &
            topology, info)
        integer, intent(in) :: poloidal(:), toroidal(:), parity(:)
        type(trial_space_topology_t), intent(out) :: topology
        integer, intent(out) :: info
        integer :: trial, trials

        info = trial_topology_invalid
        trials = size(poloidal)
        if (trials < 1) return
        if (size(toroidal) /= trials .or. size(parity) /= trials) return
        if (any(poloidal < 0)) return
        do trial = 1, trials
            if (.not. valid_phase_kind(parity(trial))) return
        end do
        allocate (topology%poloidal, source=poloidal)
        allocate (topology%toroidal, source=toroidal)
        allocate (topology%parity, source=parity)
        allocate (topology%normal_phase, source=parity)
        allocate (topology%tangential_phase(trials))
        allocate (topology%active(3, trials), source=.true.)
        do trial = 1, trials
            topology%tangential_phase(trial) = &
                opposite_phase_kind(parity(trial))
            if (poloidal(trial) /= 0 .or. toroidal(trial) /= 0) cycle
            topology%active(trial_component_normal, trial) = &
                topology%normal_phase(trial) /= trial_phase_sine
            topology%active(trial_component_eta, trial) = &
                topology%tangential_phase(trial) /= trial_phase_sine
            topology%active(trial_component_mu, trial) = &
                topology%tangential_phase(trial) /= trial_phase_sine
        end do
        info = trial_topology_ok
    end subroutine build_trial_space_topology

end module trial_space_topology
