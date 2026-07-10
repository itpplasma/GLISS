module cylinder_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use netcdf_c_api, only: nc_close_file, nc_create_netcdf4, &
        nc_def_dimension, nc_def_scalar, nc_def_variable, nc_double, &
        nc_end_definitions, nc_int64, nc_noerr, nc_put_global_text, &
        nc_put_integer, nc_put_real
    implicit none
    private

    integer, parameter, public :: fixture_ns = 33
    integer, parameter :: nm = 2, nn = 3
    integer, parameter :: field_count = 13, profile_count = 6
    character(len=11), parameter :: field_names(field_count) = &
        [character(len=11) :: "mod_B", "xhat", "yhat", "zhat", "Jac", &
        "g_tt", "g_tz", "g_zz", "II_tt", "II_tz", "II_zz", &
        "B_contra_t", "B_contra_z"]
    character(len=11), parameter :: profile_names(profile_count) = &
        [character(len=11) :: "p", "B_theta_avg", "B_zeta_avg", &
        "Phi", "chi", "iota"]

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: fixture_mu0 = 4.0e-7_dp * pi
    real(dp), parameter, public :: fixture_minor_radius = 0.5_dp
    real(dp), parameter, public :: fixture_period_length = 6.0_dp * pi
    real(dp), parameter :: mu0 = fixture_mu0
    real(dp), parameter :: minor_radius = fixture_minor_radius
    real(dp), parameter :: period_length = fixture_period_length
    real(dp), parameter :: axis_radius = 3.0_dp
    real(dp), parameter, public :: fixture_b_axial = 1.0_dp
    real(dp), parameter, public :: fixture_b_linear = 0.3_dp
    real(dp), parameter, public :: fixture_b_cubic = 0.4_dp
    real(dp), parameter :: b_axial = fixture_b_axial
    real(dp), parameter :: b_linear = fixture_b_linear
    real(dp), parameter :: b_cubic = fixture_b_cubic

    public :: create_cylinder_fixture
    public :: radius_of, b_poloidal, pressure_of, iota_of

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

    subroutine create_cylinder_fixture(filename, chart_shift, surfaces)
        character(len=*), intent(in) :: filename
        real(dp), intent(in), optional :: chart_shift
        integer, intent(in), optional :: surfaces
        integer :: ncid, dim_s, dim_m, dim_n
        integer :: id_nfp, id_beta, id_winding, id_m, id_n, id_rho
        integer :: id_profiles(profile_count)
        integer :: id_cos(field_count), id_sin(field_count)
        integer :: id_gst_cos, id_gst_sin, id_gsz_cos, id_gsz_sin
        real(dp), allocatable :: rho_values(:), s_values(:)
        real(dp), allocatable :: profile_values(:)
        real(dp), allocatable :: cosine(:, :, :), sine(:, :, :)
        real(dp) :: shift, angle, rotated_cos, rotated_sin, radius
        logical :: shifted
        integer :: ns, field, profile, i, j

        ns = fixture_ns
        if (present(surfaces)) ns = surfaces
        allocate (rho_values(ns), s_values(ns), profile_values(ns))
        allocate (cosine(nn, nm, ns), sine(nn, nm, ns))

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
        shifted = present(chart_shift)
        shift = 0.0_dp
        if (shifted) shift = chart_shift
        do field = 1, field_count
            call require_netcdf(nc_def_variable(ncid, &
                trim(field_names(field)) // "_mnc", nc_double, &
                [dim_s, dim_m, dim_n], id_cos(field)))
            call require_netcdf(nc_def_variable(ncid, &
                trim(field_names(field)) // "_mns", nc_double, &
                [dim_s, dim_m, dim_n], id_sin(field)))
        end do
        if (shifted) then
            call require_netcdf(nc_def_variable(ncid, "g_st_mnc", &
                nc_double, [dim_s, dim_m, dim_n], id_gst_cos))
            call require_netcdf(nc_def_variable(ncid, "g_st_mns", &
                nc_double, [dim_s, dim_m, dim_n], id_gst_sin))
            call require_netcdf(nc_def_variable(ncid, "g_sz_mnc", &
                nc_double, [dim_s, dim_m, dim_n], id_gsz_cos))
            call require_netcdf(nc_def_variable(ncid, "g_sz_mns", &
                nc_double, [dim_s, dim_m, dim_n], id_gsz_sin))
        end if
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
            select case (trim(field_names(field)))
            case ("xhat", "yhat", "zhat")
                if (shifted) then
                    do i = 1, ns
                        angle = 2.0_dp * pi * shift * s_values(i)
                        do j = 1, nn
                            rotated_cos = cosine(j, 2, i) * cos(angle) &
                                - sine(j, 2, i) * sin(angle)
                            rotated_sin = cosine(j, 2, i) * sin(angle) &
                                + sine(j, 2, i) * cos(angle)
                            cosine(j, 2, i) = rotated_cos
                            sine(j, 2, i) = rotated_sin
                        end do
                    end do
                end if
            end select
            call require_netcdf(nc_put_real(ncid, id_cos(field), cosine))
            call require_netcdf(nc_put_real(ncid, id_sin(field), sine))
        end do
        if (shifted) then
            cosine = 0.0_dp
            sine = 0.0_dp
            do i = 1, ns
                radius = radius_of(s_values(i))
                cosine(1, 1, i) = -shift * (2.0_dp * pi * radius)**2
            end do
            call require_netcdf(nc_put_real(ncid, id_gst_cos, cosine))
            call require_netcdf(nc_put_real(ncid, id_gst_sin, sine))
            cosine = 0.0_dp
            call require_netcdf(nc_put_real(ncid, id_gsz_cos, cosine))
            call require_netcdf(nc_put_real(ncid, id_gsz_sin, sine))
        end if
        call require_netcdf(nc_close_file(ncid))
    end subroutine create_cylinder_fixture

    subroutine require_netcdf(status)
        integer, intent(in) :: status

        if (status /= nc_noerr) then
            write (error_unit, "(a, i0)") "netcdf error ", status
            error stop 1
        end if
    end subroutine require_netcdf

end module cylinder_fixture
