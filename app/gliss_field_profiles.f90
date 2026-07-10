program gliss_field_profiles
    use, intrinsic :: iso_fortran_env, only: error_unit
    use field_profile_identities, only: compute_field_profile_identities, &
        field_profile_identities_ok, field_profile_identity_result_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(field_profile_identity_result_t) :: result
    character(len=1024) :: filename
    integer :: info, i

    if (command_argument_count() /= 1) then
        write (error_unit, "(a)") "usage: gliss_field_profiles EXPORT_FILE"
        error stop 2
    end if
    call get_command_argument(1, filename)
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) then
        write (error_unit, "(a, i0)") "reader error ", info
        error stop 1
    end if
    call compute_field_profile_identities(equilibrium, 64, 64, result, info)
    if (info /= field_profile_identities_ok) then
        write (error_unit, "(a, i0)") "field-profile error ", info
        error stop 1
    end if

    write (*, "(a)") "s,toroidal_flux_deviation," // &
        "poloidal_flux_deviation,covariant_theta_deviation," // &
        "covariant_zeta_deviation,iota_flux_deviation," // &
        "ampere_theta_deviation,ampere_zeta_deviation," // &
        "exported_jacobian_deviation,general_force_balance_deviation"
    do i = 1, size(result%s)
        write (*, "(es24.16, 9(',', es24.16))") result%s(i), &
            result%toroidal_flux_deviation(i), &
            result%poloidal_flux_deviation(i), &
            result%covariant_theta_deviation(i), &
            result%covariant_zeta_deviation(i), &
            result%iota_flux_deviation(i), &
            result%ampere_theta_deviation(i), &
            result%ampere_zeta_deviation(i), &
            result%exported_jacobian_deviation(i), &
            result%general_force_balance_deviation(i)
    end do
end program gliss_field_profiles
