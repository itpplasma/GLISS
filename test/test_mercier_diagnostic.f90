program test_mercier_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: compute_mercier, mercier_ok, &
        mercier_result_t
    use netcdf_c_api, only: nc_close_file, nc_create_netcdf4, &
        nc_def_dimension, nc_def_scalar, nc_def_variable, nc_double, &
        nc_end_definitions, nc_int64, nc_noerr, nc_put_global_text, &
        nc_put_integer, nc_put_real
    implicit none

    integer, parameter :: ns = 33, nm = 2, nn = 3
    integer, parameter :: field_count = 13, profile_count = 6
    character(len=11), parameter :: field_names(field_count) = &
        [character(len=11) :: "mod_B", "xhat", "yhat", "zhat", "Jac", &
        "g_tt", "g_tz", "g_zz", "II_tt", "II_tz", "II_zz", &
        "B_contra_t", "B_contra_z"]
    character(len=11), parameter :: profile_names(profile_count) = &
        [character(len=11) :: "p", "B_theta_avg", "B_zeta_avg", &
        "Phi", "chi", "iota"]
    character(len=*), parameter :: fixture = "mercier_cylinder.nc"

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: mu0 = 4.0e-7_dp * pi
    real(dp), parameter :: minor_radius = 0.5_dp
    real(dp), parameter :: period_length = 6.0_dp * pi
    real(dp), parameter :: axis_radius = 3.0_dp
    real(dp), parameter :: b_axial = 1.0_dp
    real(dp), parameter :: b_linear = 0.3_dp
    real(dp), parameter :: b_cubic = 0.4_dp

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(mercier_result_t) :: result
    integer :: info, i
    real(dp) :: s_value, ratio, expected_ratio, expected_shear

    call create_cylinder_fixture(fixture)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "cylinder fixture was rejected")
    call compute_mercier(equilibrium, 32, 16, result, info)
    call require(info == mercier_ok, "mercier computation failed")

    do i = 1, ns
        s_value = equilibrium%s(i)
        expected_shear = shear_term_expected()
        call require(abs(result%d_shear(i) - expected_shear) < &
            1.0e-8_dp * expected_shear, "shear term is wrong")
        ratio = result%d_mercier(i) / result%d_shear(i)
        expected_ratio = suydam_ratio_expected(radius_of(s_value))
        call require(abs(ratio - expected_ratio) < 1.0e-3_dp * &
            max(1.0_dp, abs(expected_ratio)), &
            "Mercier does not match the Suydam ratio")
        call require(abs(result%d_geodesic(i)) < 1.0e-6_dp * &
            (abs(result%d_well(i)) + expected_shear), &
            "geodesic term does not vanish")
        call require(result%iota_deviation(i) < 1.0e-10_dp, &
            "iota convention check failed")
        call require(result%boozer_deviation(i) < 1.0e-10_dp, &
            "covariant components are not flux functions")
        call require(result%force_balance_residual(i) < 1.0e-3_dp, &
            "force balance residual is too large")
        call require(result%jacobian_identity_deviation(i) < 1.0e-10_dp, &
            "Boozer Jacobian identity is violated")
    end do
    call require(any(result%d_mercier < 0.0_dp), &
        "expected Mercier-unstable surfaces were not found")

    open (unit=13, file=fixture, status="old")
    close (13, status="delete")
    write (*, "(a)") "PASS"

contains

    pure function radius_of(s_value) result(radius)
        real(dp), intent(in) :: s_value
        real(dp) :: radius

        radius = minor_radius * sqrt(s_value)
    end function radius_of

    pure function b_poloidal(radius) result(value)
        real(dp), intent(in) :: radius
        real(dp) :: value

        value = b_linear * radius + b_cubic * radius**3
    end function b_poloidal

    pure function pressure_of(radius) result(value)
        real(dp), intent(in) :: radius
        real(dp) :: value

        value = (integral_at(minor_radius) - integral_at(radius)) / mu0 &
            + 1.0e2_dp
    end function pressure_of

    pure function integral_at(radius) result(value)
        real(dp), intent(in) :: radius
        real(dp) :: value

        value = b_linear**2 * radius**2 &
            + 1.5_dp * b_linear * b_cubic * radius**4 &
            + (2.0_dp / 3.0_dp) * b_cubic**2 * radius**6
    end function integral_at

    pure function iota_of(s_value) result(value)
        real(dp), intent(in) :: s_value
        real(dp) :: value

        value = period_length * (b_linear + b_cubic * minor_radius**2 &
            * s_value) / (2.0_dp * pi * b_axial)
    end function iota_of

    pure function shear_term_expected() result(value)
        real(dp) :: value
        real(dp) :: iota_slope, psi_slope

        iota_slope = period_length * b_cubic * minor_radius**2 &
            / (2.0_dp * pi * b_axial)
        psi_slope = -minor_radius**2 * b_axial / 2.0_dp
        value = (iota_slope / psi_slope)**2 / (16.0_dp * pi**2)
    end function shear_term_expected

    pure function suydam_ratio_expected(radius) result(value)
        real(dp), intent(in) :: radius
        real(dp) :: value
        real(dp) :: pressure_slope, shear_log

        pressure_slope = -b_poloidal(radius) &
            * (2.0_dp * b_linear * radius + 4.0_dp * b_cubic * radius**3) &
            / (mu0 * radius)
        shear_log = 2.0_dp * b_cubic * radius &
            / (b_linear + b_cubic * radius**2)
        value = 1.0_dp + 8.0_dp * mu0 * pressure_slope &
            / (radius * b_axial**2 * shear_log**2)
    end function suydam_ratio_expected

    subroutine create_cylinder_fixture(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid, dim_s, dim_m, dim_n
        integer :: id_nfp, id_beta, id_winding, id_m, id_n, id_rho
        integer :: id_profiles(profile_count)
        integer :: id_cos(field_count), id_sin(field_count)
        real(dp) :: rho_values(ns), s_values(ns), radius
        real(dp) :: profile_values(ns)
        real(dp) :: cosine(nn, nm, ns), sine(nn, nm, ns)
        integer :: field, profile, i

        call require_netcdf(nc_create_netcdf4(filename, ncid))
        call require_netcdf(nc_def_dimension(ncid, "s", ns, dim_s))
        call require_netcdf(nc_def_dimension(ncid, "m", nm, dim_m))
        call require_netcdf(nc_def_dimension(ncid, "n", nn, dim_n))
        call require_netcdf(nc_def_scalar(ncid, "N_FP", nc_int64, id_nfp))
        call require_netcdf(nc_def_scalar(ncid, "beta_avg", nc_double, &
            id_beta))
        call require_netcdf(nc_def_scalar(ncid, "winding", nc_int64, &
            id_winding))
        call require_netcdf(nc_def_variable(ncid, "m", nc_int64, [dim_m], &
            id_m))
        call require_netcdf(nc_def_variable(ncid, "n", nc_int64, [dim_n], &
            id_n))
        call require_netcdf(nc_def_variable(ncid, "rho", nc_double, &
            [dim_s], id_rho))
        do profile = 1, profile_count
            call require_netcdf(nc_def_variable(ncid, &
                trim(profile_names(profile)), nc_double, [dim_s], &
                id_profiles(profile)))
        end do
        do field = 1, field_count
            call require_netcdf(nc_def_variable(ncid, &
                trim(field_names(field)) // "_mnc", nc_double, &
                [dim_s, dim_m, dim_n], id_cos(field)))
            call require_netcdf(nc_def_variable(ncid, &
                trim(field_names(field)) // "_mns", nc_double, &
                [dim_s, dim_m, dim_n], id_sin(field)))
        end do
        call require_netcdf(nc_put_global_text(ncid, &
            "stellarator_symmetry", "False"))
        call require_netcdf(nc_end_definitions(ncid))

        do i = 1, ns
            s_values(i) = (real(i, dp) - 0.5_dp) / real(ns, dp)
        end do
        rho_values = sqrt(s_values)
        call require_netcdf(nc_put_integer(ncid, id_nfp, 1))
        call require_netcdf(nc_put_real(ncid, id_beta, 0.0_dp))
        call require_netcdf(nc_put_integer(ncid, id_winding, 1))
        call require_netcdf(nc_put_integer(ncid, id_m, [0, 1]))
        call require_netcdf(nc_put_integer(ncid, id_n, [0, 1, -1]))
        call require_netcdf(nc_put_real(ncid, id_rho, rho_values))

        do profile = 1, profile_count
            do i = 1, ns
                radius = radius_of(s_values(i))
                select case (trim(profile_names(profile)))
                case ("p")
                    profile_values(i) = pressure_of(radius)
                case ("B_theta_avg")
                    profile_values(i) = 2.0_dp * pi * radius &
                        * b_poloidal(radius)
                case ("B_zeta_avg")
                    profile_values(i) = period_length * b_axial
                case ("Phi")
                    profile_values(i) = -s_values(i) * pi &
                        * minor_radius**2 * b_axial
                case ("chi")
                    profile_values(i) = -minor_radius**2 * b_axial &
                        * iota_of(s_values(i)) * s_values(i) / 2.0_dp
                case ("iota")
                    profile_values(i) = iota_of(s_values(i))
                end select
            end do
            call require_netcdf(nc_put_real(ncid, id_profiles(profile), &
                profile_values))
        end do

        do field = 1, field_count
            cosine = 0.0_dp
            sine = 0.0_dp
            do i = 1, ns
                radius = radius_of((real(i, dp) - 0.5_dp) / real(ns, dp))
                select case (trim(field_names(field)))
                case ("mod_B")
                    cosine(1, 1, i) = sqrt(b_poloidal(radius)**2 &
                        + b_axial**2)
                case ("Jac")
                    cosine(1, 1, i) = -pi * minor_radius**2 * period_length
                case ("g_tt")
                    cosine(1, 1, i) = (2.0_dp * pi * radius)**2
                case ("g_zz")
                    cosine(1, 1, i) = period_length**2
                case ("B_contra_t")
                    cosine(1, 1, i) = b_poloidal(radius) &
                        / (2.0_dp * pi * radius)
                case ("B_contra_z")
                    cosine(1, 1, i) = b_axial / period_length
                case ("xhat")
                    cosine(3, 1, i) = axis_radius
                    cosine(2, 2, i) = radius / 2.0_dp
                    cosine(3, 2, i) = radius / 2.0_dp
                case ("yhat")
                    sine(3, 1, i) = axis_radius
                    sine(3, 2, i) = radius / 2.0_dp
                    sine(2, 2, i) = -radius / 2.0_dp
                case ("zhat")
                    sine(1, 2, i) = radius
                end select
            end do
            call require_netcdf(nc_put_real(ncid, id_cos(field), cosine))
            call require_netcdf(nc_put_real(ncid, id_sin(field), sine))
        end do
        call require_netcdf(nc_close_file(ncid))
    end subroutine create_cylinder_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

    subroutine require_netcdf(status)
        integer, intent(in) :: status

        if (status /= nc_noerr) then
            write (error_unit, "(a, i0)") "netcdf error ", status
            error stop 1
        end if
    end subroutine require_netcdf

end program test_mercier_diagnostic
