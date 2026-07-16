program test_cas3d_coefficient_mass
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cas3d_coefficient_mass, only: &
        cas3d2mn_envelope_mass_scale, cas3d2_physical_mass_scale
    implicit none

    real(dp) :: envelope_mass(3, 3), inverse_transform(3, 3)
    real(dp) :: inverse_transpose(3, 3), work(3, 3)
    real(dp) :: expected_physical(3, 3), identity3(3, 3)
    real(dp) :: physical_mass(3, 3), scale

    scale = cas3d2_physical_mass_scale(48, 36, 24, 1.0_dp)
    call require(abs(scale - 1.0_dp / 41472.0_dp) < 1.0e-20_dp, &
        "CAS3D physical-sideband coefficient scale differs")
    identity3 = identity(3)
    envelope_mass = cas3d2mn_envelope_mass_scale(48, 36, 24, 1.0_dp) &
        * identity3
    inverse_transform = 0.0_dp
    inverse_transform(1, 1) = 1.0_dp
    inverse_transform(2:3, 2:3) = reshape( &
        [1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], [2, 2])
    work = matmul(envelope_mass, inverse_transform)
    inverse_transpose = transpose(inverse_transform)
    physical_mass = matmul(inverse_transpose, work)
    expected_physical = 0.0_dp
    expected_physical(1, 1) = 0.5_dp
    expected_physical(2, 2) = 1.0_dp
    expected_physical(3, 3) = 1.0_dp
    call require(maxval(abs(physical_mass - scale * expected_physical)) &
        == 0.0_dp, &
        "carrier-plus-pair envelope congruence differs")
    call require(cas3d2_physical_mass_scale(48, 36, 24, 10.0_dp) &
        == 1000.0_dp * scale, "reference-length cubic scaling differs")
    call require(cas3d2_physical_mass_scale(0, 36, 24, 1.0_dp) == 0.0_dp, &
        "zero radial interval count was accepted")
    call require(cas3d2_physical_mass_scale(48, 0, 24, 1.0_dp) == 0.0_dp, &
        "zero poloidal point count was accepted")
    call require(cas3d2_physical_mass_scale(48, 36, 0, 1.0_dp) == 0.0_dp, &
        "zero toroidal point count was accepted")
    call require(cas3d2_physical_mass_scale(48, 36, 24, 0.0_dp) == 0.0_dp, &
        "zero reference length was accepted")
    write (*, "(a)") "PASS"

contains

    pure function identity(size_value) result(matrix)
        integer, intent(in) :: size_value
        real(dp) :: matrix(size_value, size_value)
        integer :: index

        matrix = 0.0_dp
        do index = 1, size_value
            matrix(index, index) = 1.0_dp
        end do
    end function identity

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_cas3d_coefficient_mass
