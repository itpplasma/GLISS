program test_trial_space_topology
    use, intrinsic :: iso_fortran_env, only: error_unit
    use trial_space_topology, only: build_trial_space_topology, &
        trial_component_eta, trial_component_mu, trial_component_normal, &
        trial_phase_cosine, trial_phase_sine, trial_space_topology_t, &
        trial_topology_ok
    implicit none

    type(trial_space_topology_t) :: topology
    integer :: info

    call build_trial_space_topology([0, 0, 1, 0], [0, 0, 0, 3], &
        [1, 2, 1, 2], topology, info)
    call require(info == trial_topology_ok, "valid trial topology failed")
    call require(all(topology%normal_phase == [trial_phase_cosine, &
        trial_phase_sine, trial_phase_cosine, trial_phase_sine]), &
        "normal phase topology is wrong")
    call require(all(topology%tangential_phase == [trial_phase_sine, &
        trial_phase_cosine, trial_phase_sine, trial_phase_cosine]), &
        "tangential phase topology is wrong")
    call require(topology%active(trial_component_normal, 1), &
        "even zero normal trial is inactive")
    call require(.not. topology%active(trial_component_eta, 1) .and. &
        .not. topology%active(trial_component_mu, 1), &
        "identically zero even tangential trials are active")
    call require(.not. topology%active(trial_component_normal, 2), &
        "identically zero odd normal trial is active")
    call require(topology%active(trial_component_eta, 2) .and. &
        topology%active(trial_component_mu, 2), &
        "odd zero tangential trials are inactive")
    call require(all(topology%active(:, 3:4)), &
        "nonconstant harmonic lost a physical component")
    call build_trial_space_topology([0], [0, 1], [1], topology, info)
    call require(info /= trial_topology_ok, "mismatched trials were accepted")
    call build_trial_space_topology([-1], [0], [1], topology, info)
    call require(info /= trial_topology_ok, "negative poloidal mode was accepted")
    call build_trial_space_topology([0], [0], [3], topology, info)
    call require(info /= trial_topology_ok, "invalid parity was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_trial_space_topology
