program test_terpsichore_topology
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use terpsichore_topology, only: convert_terpsichore_mask, &
        parity_xi_cosine, parity_xi_sine, terpsichore_mode_mask_t, &
        terpsichore_mode_selection_t, terpsichore_topology_config_t, &
        terpsichore_topology_invalid, terpsichore_topology_ok
    implicit none

    type(terpsichore_topology_config_t) :: config
    type(terpsichore_mode_mask_t) :: mask
    type(terpsichore_mode_selection_t) :: selection
    integer :: info

    config%equilibrium_periods = 3
    config%stability_periods = 3
    config%parfac = 0.0_dp
    config%qn = 2.0_dp
    mask%poloidal_min = 0
    mask%toroidal_min = -1
    allocate (mask%selected(3, 3), source=.false.)
    mask%selected(2:3, 1) = .true.
    mask%selected(1:2, 3) = .true.
    call convert_terpsichore_mask(config, mask, selection, info)
    call require(info == terpsichore_topology_ok, &
        "valid ragged mask was rejected")
    call require(selection%field_periods == 3, &
        "equilibrium period count was not preserved")
    call require(selection%parity_class == parity_xi_sine, &
        "PARFAC zero selected the wrong parity")
    call require(all(selection%poloidal == [1, 2, 0, 1]), &
        "ragged mask poloidal order is wrong")
    call require(all(selection%toroidal == [-1, -1, 1, 1]), &
        "ragged mask toroidal order is wrong")
    call require(all(selection%stored_variable_power == [2.0_dp, 0.0_dp, &
        0.0_dp, 2.0_dp]), "QN was not restricted to shifted m=1")

    config%equilibrium_periods = 6
    config%stability_periods = 3
    config%poloidal_shift = 1
    config%parfac = 0.5_dp
    config%qn = 0.25_dp
    mask%selected = .false.
    mask%selected(1, 1) = .true.
    mask%selected(3, 3) = .true.
    call convert_terpsichore_mask(config, mask, selection, info)
    call require(info == terpsichore_topology_ok, &
        "shifted mask was rejected")
    call require(selection%parity_class == parity_xi_cosine, &
        "PARFAC one-half selected the wrong parity")
    call require(all(selection%poloidal == [1, 3]), &
        "MPINIT conversion is wrong")
    call require(all(selection%toroidal == [-2, 2]), &
        "NSTA conversion is wrong")
    call require(all(selection%stored_variable_power == [0.25_dp, 0.0_dp]), &
        "shifted QN selection is wrong")

    config%stability_periods = 4
    call convert_terpsichore_mask(config, mask, selection, info)
    call require(info == terpsichore_topology_invalid, &
        "non-divisor NSTA was accepted")
    config%stability_periods = 3
    config%parfac = 0.25_dp
    call convert_terpsichore_mask(config, mask, selection, info)
    call require(info == terpsichore_topology_invalid, &
        "unsupported PARFAC was accepted")
    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_topology
