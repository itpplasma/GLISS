program gvec_cas3d_compare_pointwise
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_netcdf, only: read_dimension, read_integer_scalar, &
        read_real_tensor, read_real_vector, reader_ok
    use gvec_cas3d_reader, only: read_gvec_cas3d_file
    use gvec_cas3d_reconstruction, only: &
        periodic_sixth_order_derivatives, reconstruct_harmonic_grid, &
        reconstruction_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t
    use netcdf_c_api, only: nc_close_file, nc_noerr, nc_open_read
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    character(len=1024) :: fourier_file, pointwise_file
    real(dp), allocatable :: pointwise_rho(:), pointwise_s(:)
    real(dp), allocatable :: theta(:), zeta_period(:)
    integer :: argument_status, field_periods, info, ncid, winding
    integer :: radial_count, theta_count, zeta_count

    call read_arguments(fourier_file, pointwise_file)
    call read_gvec_cas3d_file(trim(fourier_file), equilibrium, info)
    call require(info == reader_ok, "Fourier export could not be read")
    info = nc_open_read(trim(pointwise_file), ncid)
    call require(info == nc_noerr, "pointwise export could not be opened")
    call read_pointwise_grid(ncid, radial_count, theta_count, zeta_count, &
        field_periods, winding, pointwise_rho, pointwise_s, theta, zeta_period)
    call require(field_periods == equilibrium%field_periods, &
        "field-period counts differ")
    call require(winding == equilibrium%winding, "winding values differ")
    call require(radial_count == size(equilibrium%s), &
        "radial surface counts differ")
    call require(maxval(abs(pointwise_rho - equilibrium%rho)) < 1.0e-12_dp, &
        "rho coordinates differ")
    call require(maxval(abs(pointwise_s - equilibrium%s)) < 1.0e-12_dp, &
        "s coordinates differ")
    call require_normalized_grid(theta, "theta")
    call require_normalized_grid(zeta_period, "zeta")

    write (*, "(a)") "field,max_absolute_error,max_relative_error," // &
        "max_absolute_theta_derivative_error," // &
        "max_relative_theta_derivative_error," // &
        "max_absolute_zeta_derivative_error," // &
        "max_relative_zeta_derivative_error"
    call compare_field(ncid, "mod_B", equilibrium%mod_b, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "xhat", equilibrium%xhat, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "yhat", equilibrium%yhat, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "zhat", equilibrium%zhat, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "Jac", equilibrium%jacobian, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "g_tt", equilibrium%g_tt, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "g_tz", equilibrium%g_tz, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "g_zz", equilibrium%g_zz, equilibrium, &
        theta, zeta_period)
    call compare_field(ncid, "II_tt", equilibrium%second_form_tt, &
        equilibrium, theta, zeta_period)
    call compare_field(ncid, "II_tz", equilibrium%second_form_tz, &
        equilibrium, theta, zeta_period)
    call compare_field(ncid, "II_zz", equilibrium%second_form_zz, &
        equilibrium, theta, zeta_period)
    call compare_field(ncid, "B_contra_t", &
        equilibrium%b_contravariant_theta, equilibrium, theta, zeta_period)
    call compare_field(ncid, "B_contra_z", &
        equilibrium%b_contravariant_zeta, equilibrium, theta, zeta_period)
    info = nc_close_file(ncid)
    call require(info == nc_noerr, "pointwise export could not be closed")

contains

    subroutine read_arguments(fourier_filename, pointwise_filename)
        character(len=*), intent(out) :: fourier_filename, pointwise_filename

        call get_command_argument(1, fourier_filename, status=argument_status)
        call require(argument_status == 0, &
            "usage: gvec_cas3d_compare_pointwise FOURIER POINTWISE")
        call require(len_trim(fourier_filename) > 0, &
            "usage: gvec_cas3d_compare_pointwise FOURIER POINTWISE")
        call get_command_argument(2, pointwise_filename, status=argument_status)
        call require(argument_status == 0, &
            "usage: gvec_cas3d_compare_pointwise FOURIER POINTWISE")
        call require(len_trim(pointwise_filename) > 0, &
            "usage: gvec_cas3d_compare_pointwise FOURIER POINTWISE")
    end subroutine read_arguments

    subroutine read_pointwise_grid(file_id, ns, ntheta, nzeta, nfp, winding, &
            rho, s, theta_coordinates, zeta_coordinates)
        integer, intent(in) :: file_id
        integer, intent(out) :: ns, ntheta, nzeta, nfp, winding
        real(dp), allocatable, intent(out) :: rho(:), s(:)
        real(dp), allocatable, intent(out) :: theta_coordinates(:)
        real(dp), allocatable, intent(out) :: zeta_coordinates(:)

        call read_dimension(file_id, "s", ns, info)
        call require(info == reader_ok, "pointwise s dimension is invalid")
        call read_dimension(file_id, "theta", ntheta, info)
        call require(info == reader_ok, "pointwise theta dimension is invalid")
        call read_dimension(file_id, "zeta", nzeta, info)
        call require(info == reader_ok, "pointwise zeta dimension is invalid")
        call read_integer_scalar(file_id, "N_FP", nfp, info)
        call require(info == reader_ok, "pointwise N_FP is invalid")
        call read_integer_scalar(file_id, "winding", winding, info)
        call require(info == reader_ok, "pointwise winding is invalid")
        call read_real_vector(file_id, "rho", "s", ns, rho, info)
        call require(info == reader_ok, "pointwise rho is invalid")
        call read_real_vector(file_id, "s", "s", ns, s, info)
        call require(info == reader_ok, "pointwise s is invalid")
        call read_real_vector(file_id, "theta", "theta", ntheta, &
            theta_coordinates, info)
        call require(info == reader_ok, "pointwise theta is invalid")
        call read_real_vector(file_id, "zeta", "zeta", nzeta, &
            zeta_coordinates, info)
        call require(info == reader_ok, "pointwise zeta is invalid")
    end subroutine read_pointwise_grid

    subroutine require_normalized_grid(coordinates, name)
        real(dp), intent(in) :: coordinates(:)
        character(len=*), intent(in) :: name
        real(dp) :: expected_spacing
        integer :: point

        call require(size(coordinates) >= 2, trim(name) // " grid is too short")
        expected_spacing = 1.0_dp / real(size(coordinates), dp)
        do point = 1, size(coordinates)
            call require(abs(coordinates(point) - &
                real(point - 1, dp) * expected_spacing) < 1.0e-12_dp, &
                trim(name) // " grid is not periodic and uniform")
        end do
    end subroutine require_normalized_grid

    subroutine compare_field(file_id, name, pair, state, theta_coordinates, &
            zeta_coordinates)
        integer, intent(in) :: file_id
        character(len=*), intent(in) :: name
        type(harmonic_pair_t), intent(in) :: pair
        type(gvec_cas3d_equilibrium_t), intent(in) :: state
        real(dp), intent(in) :: theta_coordinates(:), zeta_coordinates(:)
        real(dp), allocatable :: pointwise(:, :, :)
        real(dp) :: maximum_absolute_error, maximum_relative_error
        real(dp) :: maximum_theta_error, maximum_theta_reference
        real(dp) :: maximum_zeta_error, maximum_zeta_reference
        real(dp) :: reference_scale
        integer :: surface

        call read_real_tensor(file_id, name, &
            [character(len=5) :: "s", "theta", "zeta"], &
            [size(state%s), size(theta_coordinates), size(zeta_coordinates)], &
            pointwise, info)
        call require(info == reader_ok, trim(name) // " tensor is invalid")
        maximum_absolute_error = 0.0_dp
        maximum_theta_error = 0.0_dp
        maximum_theta_reference = 0.0_dp
        maximum_zeta_error = 0.0_dp
        maximum_zeta_reference = 0.0_dp
        do surface = 1, size(state%s)
            call compare_surface(pair, state, surface, &
                pointwise(surface, :, :), theta_coordinates, zeta_coordinates, &
                maximum_absolute_error, maximum_theta_error, &
                maximum_theta_reference, maximum_zeta_error, &
                maximum_zeta_reference)
        end do
        reference_scale = maxval(abs(pointwise))
        call require(reference_scale > 0.0_dp, &
            trim(name) // " has zero reference scale")
        call require(maximum_theta_reference > 0.0_dp, &
            trim(name) // " has zero theta-derivative reference scale")
        call require(maximum_zeta_reference > 0.0_dp, &
            trim(name) // " has zero zeta-derivative reference scale")
        maximum_relative_error = maximum_absolute_error / reference_scale
        write (*, "(a,6(',',es24.16e3))") trim(name), &
            maximum_absolute_error, maximum_relative_error, &
            maximum_theta_error, maximum_theta_error / &
            maximum_theta_reference, maximum_zeta_error, &
            maximum_zeta_error / maximum_zeta_reference
    end subroutine compare_field

    subroutine compare_surface(pair, state, surface, pointwise, &
            theta_coordinates, zeta_coordinates, maximum_value_error, &
            maximum_theta_error, maximum_theta_reference, &
            maximum_zeta_error, maximum_zeta_reference)
        type(harmonic_pair_t), intent(in) :: pair
        type(gvec_cas3d_equilibrium_t), intent(in) :: state
        integer, intent(in) :: surface
        real(dp), intent(in) :: pointwise(:, :)
        real(dp), intent(in) :: theta_coordinates(:), zeta_coordinates(:)
        real(dp), intent(inout) :: maximum_value_error, maximum_theta_error
        real(dp), intent(inout) :: maximum_theta_reference
        real(dp), intent(inout) :: maximum_zeta_error, maximum_zeta_reference
        real(dp), allocatable :: reconstructed(:, :), derivative_theta(:, :)
        real(dp), allocatable :: derivative_zeta(:, :)
        real(dp), allocatable :: reference_theta(:, :), reference_zeta(:, :)

        call reconstruct_harmonic_grid(pair, surface, state%poloidal_modes, &
            state%toroidal_modes, theta_coordinates, zeta_coordinates, &
            reconstructed, derivative_theta, derivative_zeta, info)
        call require(info == reconstruction_ok, "reconstruction failed")
        call periodic_sixth_order_derivatives(pointwise, &
            1.0_dp / real(size(theta_coordinates), dp), &
            1.0_dp / real(size(zeta_coordinates), dp), reference_theta, &
            reference_zeta, info)
        call require(info == reconstruction_ok, &
            "pointwise derivative failed")
        maximum_value_error = max(maximum_value_error, &
            maxval(abs(reconstructed - pointwise)))
        maximum_theta_error = max(maximum_theta_error, &
            maxval(abs(derivative_theta - reference_theta)))
        maximum_theta_reference = max(maximum_theta_reference, &
            maxval(abs(reference_theta)))
        maximum_zeta_error = max(maximum_zeta_error, &
            maxval(abs(derivative_zeta - reference_zeta)))
        maximum_zeta_reference = max(maximum_zeta_reference, &
            maxval(abs(reference_zeta)))
    end subroutine compare_surface

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: " // message
            error stop 1
        end if
    end subroutine require

end program gvec_cas3d_compare_pointwise
