program test_terpsichore_eigen_diagnostics
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use terpsichore_eigen_diagnostics, only: &
        compute_terpsichore_eigen_diagnostics, &
        terpsichore_eigen_diagnostics_invalid, &
        terpsichore_eigen_diagnostics_ok, terpsichore_eigen_diagnostics_t
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_tridiagonal_t
    implicit none

    type(variable_block_tridiagonal_t) :: stiffness, mass
    type(terpsichore_eigen_diagnostics_t) :: diagnostics
    real(dp) :: eigenvector(2), reference(2)
    integer :: info

    call pack_variable_blocks(reshape([-2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], &
        [2, 2]), [2], stiffness, info)
    call pack_variable_blocks(reshape([2.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], &
        [2, 2]), [2], mass, info)
    eigenvector = [1.0_dp / sqrt(2.0_dp), 0.0_dp]
    reference = [1.0_dp, 0.0_dp]
    call compute_terpsichore_eigen_diagnostics(stiffness, mass, -1.0_dp, &
        eigenvector, reference, -2.0_dp, 2.0_dp, 4.0_dp, diagnostics, info)
    call require(info == terpsichore_eigen_diagnostics_ok, &
        "valid eigen diagnostics failed")
    call require(abs(diagnostics%growth_rate - 0.5_dp) < 1.0e-14_dp, &
        "growth-rate normalization is wrong")
    call require(diagnostics%reference_quotient == -1.0_dp &
        .and. diagnostics%computed_potential == -2.0_dp &
        .and. diagnostics%computed_kinetic == 2.0_dp, &
        "energy diagnostics are wrong")
    call require(diagnostics%reference_residual < 1.0e-14_dp &
        .and. abs(diagnostics%mode_overlap - 1.0_dp) < 1.0e-14_dp, &
        "eigenvector diagnostics are wrong")
    call compute_terpsichore_eigen_diagnostics(stiffness, mass, -1.0_dp, &
        eigenvector, reference, -2.0_dp, 0.0_dp, 4.0_dp, diagnostics, info)
    call require(info == terpsichore_eigen_diagnostics_invalid, &
        "nonpositive reference kinetic energy was accepted")
    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_terpsichore_eigen_diagnostics
