program test_field_profile_identities
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use field_profile_identities, only: compute_field_profile_identities, &
        field_profile_identities_ok, field_profile_identity_result_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use cylinder_fixture, only: create_cylinder_fixture
    implicit none

    character(len=*), parameter :: fixture = "field_profile_cylinder.nc"
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(field_profile_identity_result_t) :: result
    integer :: info
    real(dp) :: coarse_force_residual

    call create_cylinder_fixture(fixture)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "cylinder fixture was rejected")
    call compute_field_profile_identities(equilibrium, 32, 16, result, info)
    call require(info == field_profile_identities_ok, &
        "field-profile identity computation failed")
    call require(maxval(result%toroidal_flux_deviation) < 1.0e-10_dp, &
        "toroidal flux derivative disagrees with the field")
    call require(maxval(result%poloidal_flux_deviation) < 1.0e-10_dp, &
        "poloidal flux derivative disagrees with the field")
    call require(maxval(result%covariant_theta_deviation) < 1.0e-10_dp, &
        "exported B_theta average disagrees with the field")
    call require(maxval(result%covariant_zeta_deviation) < 1.0e-10_dp, &
        "exported B_zeta average disagrees with the field")
    call require(maxval(result%iota_flux_deviation) < 1.0e-10_dp, &
        "iota disagrees with the exported flux derivatives")
    call require(maxval(result%ampere_theta_deviation) < 1.0e-10_dp, &
        "toroidal current profile violates Ampere's law")
    call require(maxval(result%ampere_zeta_deviation) < 1.0e-10_dp, &
        "poloidal current profile violates Ampere's law")
    call require(maxval(result%exported_jacobian_deviation) < 1.0e-10_dp, &
        "exported profiles violate the Boozer Jacobian identity")
    coarse_force_residual = maxval(result%general_force_balance_deviation)
    call require(coarse_force_residual < 1.0e-4_dp, &
        "exported profiles violate force balance")

    equilibrium%toroidal_flux = equilibrium%toroidal_flux &
        + 1.0e-3_dp * equilibrium%s
    call compute_field_profile_identities(equilibrium, 32, 16, result, info)
    call require(maxval(result%toroidal_flux_deviation) > 1.0e-6_dp, &
        "diagnostic ignored a corrupted toroidal-flux profile")

    open (unit=13, file=fixture, status="old")
    close (13, status="delete")
    call create_cylinder_fixture(fixture, surfaces=65)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "refined cylinder fixture was rejected")
    call compute_field_profile_identities(equilibrium, 32, 16, result, info)
    call require(maxval(result%general_force_balance_deviation) &
        < coarse_force_residual / 3.5_dp, &
        "force-balance residual did not converge at second order")
    open (unit=13, file=fixture, status="old")
    close (13, status="delete")
    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_field_profile_identities
