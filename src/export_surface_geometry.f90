module export_surface_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use gvec_cas3d_reconstruction, only: project_harmonic_grid, &
        reconstruct_harmonic_grid, reconstruction_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_half
    use nonuniform_derivative, only: first_derivative_nonuniform
    implicit none
    private

    integer, parameter, public :: mercier_ok = 0
    integer, parameter, public :: mercier_invalid_input = 1
    integer, parameter, public :: mercier_reconstruction_error = 2

    real(dp), parameter, public :: two_pi = 2.0_dp * acos(-1.0_dp)
    real(dp), parameter, public :: mu0 = 2.0_dp * two_pi * 1.0e-7_dp
    integer(int64), parameter :: metric_points_per_harmonic = 16_int64

    type, public :: surface_data_t
        real(dp), allocatable :: jacobian(:, :), g_tt(:, :), g_tz(:, :)
        real(dp), allocatable :: g_zz(:, :), b_theta(:, :), b_zeta(:, :)
        real(dp), allocatable :: g_st(:, :), g_sz(:, :)
        real(dp), allocatable :: mod_b(:, :)
        real(dp), allocatable :: area_element(:, :)
    end type surface_data_t

    type, public :: surface_profiles_t
        real(dp) :: flux_slope, poloidal_slope
        real(dp) :: flux_curvature, poloidal_curvature
        real(dp) :: covariant_theta, covariant_zeta
        real(dp) :: covariant_theta_slope, covariant_zeta_slope
        real(dp) :: pressure_slope
    end type surface_profiles_t

    public :: build_angular_grids
    public :: build_kernel_geometry
    public :: build_surface_kernel_fields
    public :: beta_from_positions
    public :: differentiate_pair
    public :: grid_mean
    public :: load_surface
    public :: solve_beta_derivatives
    public :: solve_beta_derivatives_modes
    public :: surface_derivatives
    public :: surface_values
    public :: validate_tangential_metric

contains

    subroutine build_kernel_geometry(equilibrium, n_theta, n_zeta, &
            fields, drive, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        real(dp), allocatable, intent(out) :: fields(:, :, :, :)
        real(dp), allocatable, intent(out) :: drive(:, :, :)
        integer, intent(out) :: info
        type(surface_data_t) :: surface
        type(harmonic_pair_t) :: jacobian_slope
        type(surface_profiles_t) :: profiles
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: covariant_theta(:), covariant_zeta(:)
        real(dp), allocatable :: flux_slope(:), poloidal_slope(:)
        real(dp), allocatable :: covariant_theta_slope(:)
        real(dp), allocatable :: covariant_zeta_slope(:)
        real(dp), allocatable :: pressure_slope(:)
        real(dp), allocatable :: flux_curvature(:), poloidal_curvature(:)
        integer :: ns, i

        info = mercier_invalid_input
        if (n_theta < 8 .or. n_zeta < 8) return
        if (equilibrium%radial_grid /= radial_grid_half) return
        ns = size(equilibrium%s)
        if (ns < 5) return
        call validate_tangential_metric(equilibrium, n_theta, n_zeta, info)
        if (info /= mercier_ok) return

        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        call differentiate_pair(equilibrium%s, equilibrium%jacobian, &
            jacobian_slope)
        allocate (covariant_theta(ns), covariant_zeta(ns))
        allocate (flux_slope(ns), poloidal_slope(ns))
        allocate (covariant_theta_slope(ns), covariant_zeta_slope(ns))
        allocate (pressure_slope(ns), flux_curvature(ns))
        allocate (poloidal_curvature(ns))
        allocate (fields(n_theta, n_zeta, 13, ns))
        allocate (drive(n_theta, n_zeta, ns))

        do i = 1, ns
            call load_surface(equilibrium, i, theta, zeta, surface, info)
            if (info /= mercier_ok) return
            covariant_theta(i) = grid_two_product_mean(surface%g_tt, &
                surface%b_theta, surface%g_tz, surface%b_zeta)
            covariant_zeta(i) = grid_two_product_mean(surface%g_tz, &
                surface%b_theta, surface%g_zz, surface%b_zeta)
            flux_slope(i) = grid_product_mean(surface%jacobian, &
                surface%b_zeta)
            poloidal_slope(i) = grid_product_mean(surface%jacobian, &
                surface%b_theta)
        end do
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%pressure, pressure_slope)
        call first_derivative_nonuniform(equilibrium%s, covariant_theta, &
            covariant_theta_slope)
        call first_derivative_nonuniform(equilibrium%s, covariant_zeta, &
            covariant_zeta_slope)
        call first_derivative_nonuniform(equilibrium%s, flux_slope, &
            flux_curvature)
        call first_derivative_nonuniform(equilibrium%s, poloidal_slope, &
            poloidal_curvature)

        do i = 1, ns
            call load_surface(equilibrium, i, theta, zeta, surface, info)
            if (info /= mercier_ok) return
            profiles = surface_profiles_t(flux_slope(i), &
                poloidal_slope(i), flux_curvature(i), &
                poloidal_curvature(i), covariant_theta(i), &
                covariant_zeta(i), covariant_theta_slope(i), &
                covariant_zeta_slope(i), pressure_slope(i))
            call fill_surface_fields(equilibrium, surface, profiles, &
                jacobian_slope, i, theta, zeta, fields(:, :, :, i), &
                drive(:, :, i), info)
            if (info /= mercier_ok) return
        end do
        info = mercier_ok
    end subroutine build_kernel_geometry

    subroutine validate_tangential_metric(equilibrium, n_theta, n_zeta, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        integer, intent(out) :: info
        real(dp), allocatable :: determinant(:, :), g_tt(:, :)
        real(dp), allocatable :: g_tz(:, :), g_zz(:, :), theta(:), zeta(:)
        integer :: check_theta, check_zeta, surface

        info = mercier_invalid_input
        call metric_validation_grid(equilibrium, n_theta, n_zeta, &
            check_theta, check_zeta, info)
        if (info /= mercier_ok) return
        call build_angular_grids(check_theta, check_zeta, theta, zeta)
        do surface = 1, size(equilibrium%s)
            call surface_values(equilibrium%g_tt, surface, equilibrium, &
                theta, zeta, g_tt, info)
            if (info /= mercier_ok) return
            call surface_values(equilibrium%g_tz, surface, equilibrium, &
                theta, zeta, g_tz, info)
            if (info /= mercier_ok) return
            call surface_values(equilibrium%g_zz, surface, equilibrium, &
                theta, zeta, g_zz, info)
            if (info /= mercier_ok) return
            determinant = g_tt * g_zz - g_tz**2
            if (.not. tangential_metric_is_positive( &
                g_tt, g_zz, determinant)) then
                info = mercier_invalid_input
                return
            end if
        end do
        info = mercier_ok
    end subroutine validate_tangential_metric

    subroutine metric_validation_grid(equilibrium, n_theta, n_zeta, &
            check_theta, check_zeta, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        integer, intent(out) :: check_theta, check_zeta, info
        integer(int64) :: maximum_grid_size, poloidal, toroidal

        info = mercier_invalid_input
        if (n_theta < 8 .or. n_zeta < 8) return
        if (.not. allocated(equilibrium%poloidal_modes)) return
        if (.not. allocated(equilibrium%toroidal_modes)) return
        if (size(equilibrium%poloidal_modes) < 1) return
        if (size(equilibrium%toroidal_modes) < 1) return
        poloidal = maxval(abs(int(equilibrium%poloidal_modes, int64)))
        toroidal = maxval(abs(int(equilibrium%toroidal_modes, int64)))
        maximum_grid_size = int(huge(check_theta), int64)
        if (poloidal > maximum_grid_size &
            / int(metric_points_per_harmonic, int64) - 1) return
        if (toroidal > maximum_grid_size &
            / int(metric_points_per_harmonic, int64) - 1) return
        check_theta = max(n_theta, int(metric_points_per_harmonic &
            * (poloidal + 1)))
        check_zeta = max(n_zeta, int(metric_points_per_harmonic &
            * (toroidal + 1)))
        info = mercier_ok
    end subroutine metric_validation_grid

    pure function tangential_metric_is_positive(g_tt, g_zz, determinant) &
            result(positive)
        real(dp), intent(in) :: g_tt(:, :), g_zz(:, :), determinant(:, :)
        logical :: positive

        positive = all(ieee_is_finite(g_tt))
        if (.not. positive) return
        positive = all(ieee_is_finite(g_zz))
        if (.not. positive) return
        positive = all(ieee_is_finite(determinant))
        if (.not. positive) return
        positive = all(g_tt > 0.0_dp) .and. all(g_zz > 0.0_dp) &
            .and. all(determinant > 0.0_dp)
    end function tangential_metric_is_positive

    function kernel_surface_inputs_are_valid(poloidal_modes, toroidal_modes, &
            surface, jacobian_slope, theta, zeta, fields, drive) result(valid)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: jacobian_slope(:, :), theta(:), zeta(:)
        real(dp), intent(in) :: fields(:, :, :), drive(:, :)
        logical :: valid
        integer :: expected(2)

        valid = size(poloidal_modes) > 0 .and. size(toroidal_modes) > 0
        if (.not. valid) return
        valid = size(theta) > 0 .and. size(zeta) > 0
        if (.not. valid) return
        valid = all(ieee_is_finite(theta)) .and. all(ieee_is_finite(zeta))
        if (.not. valid) return
        expected(1) = size(theta)
        expected(2) = size(zeta)
        valid = matrix_has_shape(jacobian_slope, expected) &
            .and. matrix_has_shape(drive, expected) &
            .and. size(fields, 1) == expected(1) &
            .and. size(fields, 2) == expected(2) &
            .and. size(fields, 3) == 13
        if (.not. valid) return
        valid = allocated(surface%jacobian) .and. allocated(surface%g_tt) &
            .and. allocated(surface%g_tz) .and. allocated(surface%g_zz) &
            .and. allocated(surface%b_theta) .and. allocated(surface%b_zeta) &
            .and. allocated(surface%g_st) .and. allocated(surface%g_sz) &
            .and. allocated(surface%mod_b)
        if (.not. valid) return
        valid = matrix_has_shape(surface%jacobian, expected) &
            .and. matrix_has_shape(surface%g_tt, expected) &
            .and. matrix_has_shape(surface%g_tz, expected) &
            .and. matrix_has_shape(surface%g_zz, expected) &
            .and. matrix_has_shape(surface%b_theta, expected) &
            .and. matrix_has_shape(surface%b_zeta, expected) &
            .and. matrix_has_shape(surface%g_st, expected) &
            .and. matrix_has_shape(surface%g_sz, expected) &
            .and. matrix_has_shape(surface%mod_b, expected)
        if (.not. valid) return
        valid = all(ieee_is_finite(jacobian_slope)) &
            .and. all(ieee_is_finite(surface%jacobian)) &
            .and. all(ieee_is_finite(surface%g_tt)) &
            .and. all(ieee_is_finite(surface%g_tz)) &
            .and. all(ieee_is_finite(surface%g_zz)) &
            .and. all(ieee_is_finite(surface%b_theta)) &
            .and. all(ieee_is_finite(surface%b_zeta)) &
            .and. all(ieee_is_finite(surface%g_st)) &
            .and. all(ieee_is_finite(surface%g_sz)) &
            .and. all(ieee_is_finite(surface%mod_b))
        if (.not. valid) return
        valid = all(surface%jacobian /= 0.0_dp) &
            .and. all(surface%mod_b > 0.0_dp)
    end function kernel_surface_inputs_are_valid

    subroutine fill_surface_fields(equilibrium, surface, profiles, &
            jacobian_slope, i, theta, zeta, fields, drive, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        type(surface_profiles_t), intent(in) :: profiles
        type(harmonic_pair_t), intent(in) :: jacobian_slope
        integer, intent(in) :: i
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(out) :: fields(:, :, :)
        real(dp), intent(out) :: drive(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: jac_slope_grid(:, :)
        real(dp), allocatable :: discard_a(:, :), discard_b(:, :)
        integer :: rec_info

        call reconstruct_harmonic_grid(jacobian_slope, i, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, &
            theta, zeta, jac_slope_grid, discard_a, discard_b, rec_info)
        if (rec_info /= reconstruction_ok) then
            info = mercier_reconstruction_error
            return
        end if
        call build_surface_kernel_fields(equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, equilibrium%has_chart_metric, &
            surface, profiles, jac_slope_grid, theta, zeta, fields, drive, &
            info)
    end subroutine fill_surface_fields

    subroutine build_surface_kernel_fields(poloidal_modes, toroidal_modes, &
            has_chart_metric, surface, profiles, jacobian_slope, theta, zeta, &
            fields, drive, info)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        logical, intent(in) :: has_chart_metric
        type(surface_data_t), intent(in) :: surface
        type(surface_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: jacobian_slope(:, :), theta(:), zeta(:)
        real(dp), intent(out) :: fields(:, :, :), drive(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: beta_values(:, :), beta_theta(:, :)
        real(dp), allocatable :: beta_zeta(:, :), grad_s2(:, :)

        info = mercier_invalid_input
        if (.not. kernel_surface_inputs_are_valid(poloidal_modes, &
            toroidal_modes, surface, jacobian_slope, theta, zeta, fields, &
            drive)) return
        if (.not. surface_profiles_are_finite(profiles)) return
        call solve_beta_derivatives_modes(poloidal_modes, toroidal_modes, &
            surface, theta, zeta, profiles%covariant_theta_slope, &
            profiles%covariant_zeta_slope, profiles%pressure_slope, &
            profiles%poloidal_slope, profiles%flux_slope, beta_values, &
            beta_theta, beta_zeta, info=info)
        if (info /= mercier_ok) return
        fields(:, :, 1) = profiles%flux_slope
        fields(:, :, 2) = profiles%poloidal_slope
        fields(:, :, 3) = profiles%flux_curvature
        fields(:, :, 4) = profiles%poloidal_curvature
        fields(:, :, 5) = surface%g_tz * surface%b_theta &
            + surface%g_zz * surface%b_zeta
        fields(:, :, 6) = surface%g_tt * surface%b_theta &
            + surface%g_tz * surface%b_zeta
        fields(:, :, 7) = surface%jacobian
        fields(:, :, 8) = surface%mod_b
        grad_s2 = (surface%g_tt * surface%g_zz - surface%g_tz**2) &
            / surface%jacobian**2
        if (.not. all(ieee_is_finite(grad_s2)) &
            .or. any(grad_s2 <= 0.0_dp)) then
            info = mercier_invalid_input
            return
        end if
        fields(:, :, 9) = grad_s2
        fields(:, :, 10) = ((beta_zeta &
            - profiles%covariant_zeta_slope) * fields(:, :, 6) &
            + (profiles%covariant_theta_slope - beta_theta) &
            * fields(:, :, 5)) / surface%jacobian
        fields(:, :, 11) = mu0 * profiles%pressure_slope
        if (has_chart_metric) then
            fields(:, :, 12) = ((surface%g_tz * surface%b_theta &
                + surface%g_zz * surface%b_zeta) * surface%g_st &
                - (surface%g_tt * surface%b_theta &
                + surface%g_tz * surface%b_zeta) * surface%g_sz) &
                / (surface%jacobian * surface%mod_b)
            fields(:, :, 13) = surface%b_theta * surface%g_st &
                + surface%b_zeta * surface%g_sz
        else
            fields(:, :, 12) = 0.0_dp
            fields(:, :, 13) = beta_values
        end if
        ! The current-curvature group enters with a minus sign in the
        ! export chart; pinned against the geometric drive in
        ! derivations/drive_machinery_identity.wl.
        drive = (fields(:, :, 10)**2 &
            + (mu0 * profiles%pressure_slope)**2 * fields(:, :, 9)) &
            / (surface%mod_b**2 * fields(:, :, 9)) &
            - (profiles%flux_curvature * profiles%covariant_zeta_slope &
            + profiles%poloidal_curvature &
            * profiles%covariant_theta_slope &
            - profiles%flux_curvature * beta_zeta &
            - profiles%poloidal_curvature * beta_theta) &
            / surface%jacobian &
            - mu0 * profiles%pressure_slope * jacobian_slope &
            / surface%jacobian
        info = mercier_ok
        if (has_chart_metric) then
            call add_drive_chart_term(poloidal_modes, toroidal_modes, &
                surface, theta, zeta, &
                profiles%covariant_theta_slope, &
                profiles%covariant_zeta_slope, beta_theta, beta_zeta, &
                profiles%poloidal_slope, profiles%flux_slope, drive, &
                info)
        end if
    end subroutine build_surface_kernel_fields

    pure function surface_profiles_are_finite(profiles) result(finite)
        type(surface_profiles_t), intent(in) :: profiles
        logical :: finite

        finite = ieee_is_finite(profiles%flux_slope) &
            .and. ieee_is_finite(profiles%poloidal_slope) &
            .and. ieee_is_finite(profiles%flux_curvature) &
            .and. ieee_is_finite(profiles%poloidal_curvature) &
            .and. ieee_is_finite(profiles%covariant_theta) &
            .and. ieee_is_finite(profiles%covariant_zeta) &
            .and. ieee_is_finite(profiles%covariant_theta_slope) &
            .and. ieee_is_finite(profiles%covariant_zeta_slope) &
            .and. ieee_is_finite(profiles%pressure_slope)
    end function surface_profiles_are_finite

    subroutine add_drive_chart_term(poloidal_modes, toroidal_modes, surface, &
            theta, zeta, &
            covariant_theta_slope, covariant_zeta_slope, beta_theta, &
            beta_zeta, poloidal_flux_slope, toroidal_flux_slope, drive, &
            info)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: beta_theta(:, :), beta_zeta(:, :)
        real(dp), intent(in) :: poloidal_flux_slope, toroidal_flux_slope
        real(dp), intent(inout) :: drive(:, :)
        integer, intent(out) :: info
        type(harmonic_pair_t) :: operand_pair
        real(dp), allocatable :: operand(:, :), grad_s2(:, :)
        real(dp), allocatable :: upper_ts(:, :), upper_zs(:, :)
        real(dp), allocatable :: operand_values(:, :)
        real(dp), allocatable :: operand_theta(:, :), operand_zeta(:, :)
        real(dp) :: operand_cosine(size(poloidal_modes), size(toroidal_modes))
        real(dp) :: operand_sine(size(poloidal_modes), size(toroidal_modes))
        integer :: rec_info

        upper_ts = (surface%g_sz * surface%g_tz &
            - surface%g_st * surface%g_zz) / surface%jacobian**2
        upper_zs = (surface%g_st * surface%g_tz &
            - surface%g_sz * surface%g_tt) / surface%jacobian**2
        grad_s2 = (surface%g_tt * surface%g_zz - surface%g_tz**2) &
            / surface%jacobian**2
        if (.not. all(ieee_is_finite(grad_s2)) &
            .or. any(grad_s2 <= 0.0_dp)) then
            info = mercier_invalid_input
            return
        end if
        operand = ((covariant_theta_slope - beta_theta) * upper_ts &
            - (beta_zeta - covariant_zeta_slope) * upper_zs) / grad_s2
        call project_harmonic_grid(operand, poloidal_modes, toroidal_modes, &
            theta, zeta, operand_cosine, &
            operand_sine)
        allocate (operand_pair%cosine(1, size(operand_cosine, 1), &
            size(operand_cosine, 2)))
        allocate (operand_pair%sine(1, size(operand_sine, 1), &
            size(operand_sine, 2)))
        operand_pair%cosine(1, :, :) = operand_cosine
        operand_pair%sine(1, :, :) = operand_sine
        call reconstruct_harmonic_grid(operand_pair, 1, &
            poloidal_modes, toroidal_modes, &
            theta, zeta, operand_values, operand_theta, operand_zeta, &
            rec_info)
        if (rec_info /= reconstruction_ok) then
            info = mercier_reconstruction_error
            return
        end if
        drive = drive + (poloidal_flux_slope * operand_theta &
            + toroidal_flux_slope * operand_zeta) / surface%jacobian
        info = mercier_ok
    end subroutine add_drive_chart_term

    subroutine load_surface(equilibrium, i, theta, zeta, surface, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: i
        real(dp), intent(in) :: theta(:), zeta(:)
        type(surface_data_t), intent(out) :: surface
        integer, intent(out) :: info

        call surface_values(equilibrium%jacobian, i, equilibrium, theta, &
            zeta, surface%jacobian, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_tt, i, equilibrium, theta, zeta, &
            surface%g_tt, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_tz, i, equilibrium, theta, zeta, &
            surface%g_tz, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_zz, i, equilibrium, theta, zeta, &
            surface%g_zz, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%b_contravariant_theta, i, &
            equilibrium, theta, zeta, surface%b_theta, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%b_contravariant_zeta, i, &
            equilibrium, theta, zeta, surface%b_zeta, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%mod_b, i, equilibrium, theta, zeta, &
            surface%mod_b, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_st, i, equilibrium, theta, zeta, &
            surface%g_st, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_sz, i, equilibrium, theta, zeta, &
            surface%g_sz, info)
        if (info /= mercier_ok) return

        surface%area_element = sqrt(max(surface%g_tt * surface%g_zz &
            - surface%g_tz**2, 0.0_dp))
        info = mercier_ok
    end subroutine load_surface

    subroutine surface_values(pair, i, equilibrium, theta, zeta, values, &
            info)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: i
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: discard_theta(:, :), discard_zeta(:, :)
        integer :: rec_info

        call reconstruct_harmonic_grid(pair, i, equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, theta, zeta, values, discard_theta, &
            discard_zeta, rec_info)
        info = merge(mercier_ok, mercier_reconstruction_error, &
            rec_info == reconstruction_ok)
    end subroutine surface_values

    subroutine solve_beta_derivatives(equilibrium, surface, theta, zeta, &
            covariant_theta_slope, covariant_zeta_slope, pressure_slope, &
            poloidal_flux_slope, toroidal_flux_slope, beta_values, &
            beta_theta, beta_zeta, beta_harmonics)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: pressure_slope, poloidal_flux_slope
        real(dp), intent(in) :: toroidal_flux_slope
        real(dp), allocatable, intent(out) :: beta_values(:, :)
        real(dp), allocatable, intent(out) :: beta_theta(:, :)
        real(dp), allocatable, intent(out) :: beta_zeta(:, :)
        type(harmonic_pair_t), intent(out), optional :: beta_harmonics
        type(harmonic_pair_t) :: beta_pair

        call solve_beta_derivatives_modes(equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, surface, theta, zeta, &
            covariant_theta_slope, covariant_zeta_slope, pressure_slope, &
            poloidal_flux_slope, toroidal_flux_slope, beta_values, &
            beta_theta, beta_zeta, beta_pair)
        if (present(beta_harmonics)) beta_harmonics = beta_pair
    end subroutine solve_beta_derivatives

    subroutine solve_beta_derivatives_modes(poloidal_modes, toroidal_modes, &
            surface, theta, zeta, covariant_theta_slope, &
            covariant_zeta_slope, pressure_slope, poloidal_flux_slope, &
            toroidal_flux_slope, beta_values, beta_theta, beta_zeta, &
            beta_harmonics, info)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: pressure_slope, poloidal_flux_slope
        real(dp), intent(in) :: toroidal_flux_slope
        real(dp), allocatable, intent(out) :: beta_values(:, :)
        real(dp), allocatable, intent(out) :: beta_theta(:, :)
        real(dp), allocatable, intent(out) :: beta_zeta(:, :)
        type(harmonic_pair_t), intent(out), optional :: beta_harmonics
        integer, intent(out), optional :: info
        type(harmonic_pair_t) :: beta_pair
        real(dp), allocatable :: rhs(:, :)
        real(dp) :: rhs_cosine(size(poloidal_modes), size(toroidal_modes))
        real(dp) :: rhs_sine(size(poloidal_modes), size(toroidal_modes))
        real(dp) :: denominator, scale
        integer :: allocation_status, mode_m, mode_n, rec_info

        if (present(info)) info = mercier_invalid_input
        if (size(poloidal_modes) < 1 .or. size(toroidal_modes) < 1) return
        if (size(theta) < 1 .or. size(zeta) < 1) return
        if (.not. all(ieee_is_finite(theta)) &
            .or. .not. all(ieee_is_finite(zeta))) return
        rhs = surface%jacobian * (mu0 * pressure_slope &
            + covariant_zeta_slope * surface%b_zeta &
            + covariant_theta_slope * surface%b_theta)
        call project_harmonic_grid(rhs, poloidal_modes, toroidal_modes, theta, &
            zeta, rhs_cosine, rhs_sine)
        allocate (beta_pair%cosine(1, size(rhs_cosine, 1), &
            size(rhs_cosine, 2)), beta_pair%sine(1, size(rhs_sine, 1), &
            size(rhs_sine, 2)), stat=allocation_status)
        if (allocation_status /= 0) return
        scale = abs(toroidal_flux_slope) + abs(poloidal_flux_slope)
        do mode_n = 1, size(toroidal_modes)
            do mode_m = 1, size(poloidal_modes)
                denominator = two_pi * (real( &
                    poloidal_modes(mode_m), dp) &
                    * poloidal_flux_slope - real( &
                    toroidal_modes(mode_n), dp) &
                    * toroidal_flux_slope)
                if (abs(denominator) < 1.0e-10_dp * scale) then
                    beta_pair%cosine(1, mode_m, mode_n) = 0.0_dp
                    beta_pair%sine(1, mode_m, mode_n) = 0.0_dp
                else
                    beta_pair%sine(1, mode_m, mode_n) = &
                        rhs_cosine(mode_m, mode_n) / denominator
                    beta_pair%cosine(1, mode_m, mode_n) = &
                        -rhs_sine(mode_m, mode_n) / denominator
                end if
            end do
        end do
        call reconstruct_harmonic_grid(beta_pair, 1, &
            poloidal_modes, toroidal_modes, theta, zeta, beta_values, &
            beta_theta, beta_zeta, rec_info)
        if (rec_info /= reconstruction_ok) return
        if (present(beta_harmonics)) beta_harmonics = beta_pair
        if (present(info)) info = mercier_ok
    end subroutine solve_beta_derivatives_modes

    pure subroutine differentiate_pair(s, pair, slope_pair)
        real(dp), intent(in) :: s(:)
        type(harmonic_pair_t), intent(in) :: pair
        type(harmonic_pair_t), intent(out) :: slope_pair
        integer :: m, n

        allocate (slope_pair%cosine, mold=pair%cosine)
        allocate (slope_pair%sine, mold=pair%sine)
        do n = 1, size(pair%cosine, 3)
            do m = 1, size(pair%cosine, 2)
                call first_derivative_nonuniform(s, pair%cosine(:, m, n), &
                    slope_pair%cosine(:, m, n))
                call first_derivative_nonuniform(s, pair%sine(:, m, n), &
                    slope_pair%sine(:, m, n))
            end do
        end do
    end subroutine differentiate_pair

    subroutine beta_from_positions(equilibrium, xhat_s, yhat_s, zhat_s, &
            i, theta, zeta, surface, beta_positions, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(harmonic_pair_t), intent(in) :: xhat_s, yhat_s, zhat_s
        integer, intent(in) :: i
        real(dp), intent(in) :: theta(:), zeta(:)
        type(surface_data_t), intent(in) :: surface
        real(dp), allocatable, intent(out) :: beta_positions(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: es_x(:, :), es_y(:, :), es_z(:, :)
        real(dp), allocatable :: eu_x(:, :), eu_y(:, :), eu_z(:, :)
        real(dp), allocatable :: ev_x(:, :), ev_y(:, :), ev_z(:, :)
        real(dp), allocatable :: discard_a(:, :), discard_b(:, :)

        call surface_values(xhat_s, i, equilibrium, theta, zeta, es_x, &
            info)
        if (info /= mercier_ok) return
        call surface_values(yhat_s, i, equilibrium, theta, zeta, es_y, &
            info)
        if (info /= mercier_ok) return
        call surface_values(zhat_s, i, equilibrium, theta, zeta, es_z, &
            info)
        if (info /= mercier_ok) return
        call surface_derivatives(equilibrium%xhat, i, equilibrium, theta, &
            zeta, discard_a, eu_x, ev_x, info)
        if (info /= mercier_ok) return
        call surface_derivatives(equilibrium%yhat, i, equilibrium, theta, &
            zeta, discard_a, eu_y, ev_y, info)
        if (info /= mercier_ok) return
        call surface_derivatives(equilibrium%zhat, i, equilibrium, theta, &
            zeta, discard_b, eu_z, ev_z, info)
        if (info /= mercier_ok) return
        beta_positions = (es_x * eu_x + es_y * eu_y + es_z * eu_z) &
            * surface%b_theta &
            + (es_x * ev_x + es_y * ev_y + es_z * ev_z) * surface%b_zeta
    end subroutine beta_from_positions

    subroutine surface_derivatives(pair, i, equilibrium, theta, zeta, &
            values, derivative_theta, derivative_zeta, info)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: i
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        real(dp), allocatable, intent(out) :: derivative_theta(:, :)
        real(dp), allocatable, intent(out) :: derivative_zeta(:, :)
        integer, intent(out) :: info
        integer :: rec_info

        call reconstruct_harmonic_grid(pair, i, equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, theta, zeta, values, &
            derivative_theta, derivative_zeta, rec_info)
        info = merge(mercier_ok, mercier_reconstruction_error, &
            rec_info == reconstruction_ok)
    end subroutine surface_derivatives

    pure function grid_mean(values) result(mean)
        real(dp), intent(in) :: values(:, :)
        real(dp) :: mean

        mean = sum(values) / real(size(values), dp)
    end function grid_mean

    pure function grid_product_mean(first, second) result(mean)
        real(dp), intent(in) :: first(:, :), second(:, :)
        real(dp) :: mean
        integer :: column, row

        mean = 0.0_dp
        do column = 1, size(first, 2)
            do row = 1, size(first, 1)
                mean = mean + first(row, column) * second(row, column)
            end do
        end do
        mean = mean / real(size(first), dp)
    end function grid_product_mean

    pure function grid_two_product_mean(first_a, first_b, second_a, second_b) &
            result(mean)
        real(dp), intent(in) :: first_a(:, :), first_b(:, :)
        real(dp), intent(in) :: second_a(:, :), second_b(:, :)
        real(dp) :: mean
        integer :: column, row

        mean = 0.0_dp
        do column = 1, size(first_a, 2)
            do row = 1, size(first_a, 1)
                mean = mean + first_a(row, column) * first_b(row, column) &
                    + second_a(row, column) * second_b(row, column)
            end do
        end do
        mean = mean / real(size(first_a), dp)
    end function grid_two_product_mean

    pure function matrix_has_shape(matrix, expected) result(valid)
        real(dp), intent(in) :: matrix(:, :)
        integer, intent(in) :: expected(2)
        logical :: valid

        valid = size(matrix, 1) == expected(1) &
            .and. size(matrix, 2) == expected(2)
    end function matrix_has_shape

    pure subroutine build_angular_grids(n_theta, n_zeta, theta, zeta)
        integer, intent(in) :: n_theta, n_zeta
        real(dp), allocatable, intent(out) :: theta(:), zeta(:)
        integer :: j

        allocate (theta(n_theta), zeta(n_zeta))
        do j = 1, n_theta
            theta(j) = real(j - 1, dp) / real(n_theta, dp)
        end do
        do j = 1, n_zeta
            zeta(j) = real(j - 1, dp) / real(n_zeta, dp)
        end do
    end subroutine build_angular_grids

end module export_surface_geometry
