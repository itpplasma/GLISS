program gliss_primitive_geometry_trace
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use export_surface_geometry, only: build_angular_grids, differentiate_pair
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_reconstruction, only: reconstruct_harmonic_grid, &
        reconstruction_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    use primitive_equilibrium_spline, only: evaluate_primitive_equilibrium, &
        fit_primitive_equilibrium, primitive_equilibrium_ok, &
        primitive_equilibrium_spline_t
    use primitive_geometry_grid, only: primitive_geometry_grid_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
    use two_component_kernel, only: two_component_components
    implicit none

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(primitive_equilibrium_spline_t) :: spline
    type(primitive_geometry_grid_t) :: geometry
    type(harmonic_pair_t) :: jacobian_s_pair
    character(len=1024) :: filename
    real(dp), allocatable :: theta(:), zeta(:)
    real(dp), allocatable :: jacobian(:, :), jacobian_s(:, :)
    real(dp), allocatable :: g_tt(:, :), g_tz(:, :), g_zz(:, :)
    real(dp), allocatable :: g_st(:, :), g_sz(:, :), mod_b(:, :)
    real(dp), allocatable :: b_theta(:, :), b_zeta(:, :)
    real(dp), allocatable :: second_tt(:, :), second_tz(:, :)
    real(dp), allocatable :: second_zz(:, :)
    real(dp), allocatable :: discard_theta(:, :), discard_zeta(:, :)
    real(dp), allocatable :: covariant_theta(:, :), covariant_zeta(:, :)
    real(dp), allocatable :: legacy_fields(:, :, :, :), legacy_drive(:, :, :)
    real(dp), allocatable :: primitive_fields(:, :, :), primitive_drive(:, :)
    real(dp), allocatable :: geometric_drive(:, :)
    real(dp) :: pressure, pressure_slope
    integer :: info, n_theta, n_zeta, surface

    call read_arguments(filename, n_theta, n_zeta)
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    call require(info == reader_ok, "export could not be read")
    call require(equilibrium%has_chart_metric, &
        "export does not contain radial chart metrics")
    call fit_primitive_equilibrium(equilibrium, spline, info)
    call require(info == primitive_equilibrium_ok, &
        "primitive equilibrium fit failed")
    call build_angular_grids(n_theta, n_zeta, theta, zeta)
    call build_kernel_geometry(equilibrium, n_theta, n_zeta, legacy_fields, &
        legacy_drive, info)
    call require(info == mercier_ok, "legacy kernel geometry failed")
    call differentiate_pair(equilibrium%s, equilibrium%jacobian, &
        jacobian_s_pair)

    write (*, "(a)") "s,field,relative_l2,normalized_max,scaled_max," // &
        "absolute_max,reference_rms,reference_max,comparison_scale"
    do surface = 1, size(equilibrium%s)
        call evaluate_primitive_equilibrium(spline, equilibrium%s(surface), &
            theta, zeta, geometry, pressure, pressure_slope, info)
        call require(info == primitive_equilibrium_ok, &
            "primitive equilibrium evaluation failed")
        call reconstruct_reference_fields(surface)
        call evaluate_primitive_kernel_surface(spline, &
            equilibrium%s(surface), theta, zeta, primitive_fields, &
            primitive_drive, info, geometric_drive=geometric_drive)
        call require(info == primitive_kernel_ok, &
            "primitive kernel geometry failed")
        covariant_theta = g_tt * b_theta + g_tz * b_zeta
        covariant_zeta = g_tz * b_theta + g_zz * b_zeta
        call write_comparison("Jac", geometry%signed_jacobian, jacobian, &
            equilibrium%s(surface))
        call write_comparison("Jac_s", geometry%jacobian_s, jacobian_s, &
            equilibrium%s(surface))
        call write_comparison("g_tt", geometry%metric(:, :, 2, 2), g_tt, &
            equilibrium%s(surface))
        call write_comparison("g_tz", geometry%metric(:, :, 2, 3), g_tz, &
            equilibrium%s(surface))
        call write_comparison("g_zz", geometry%metric(:, :, 3, 3), g_zz, &
            equilibrium%s(surface))
        call write_comparison("g_st", geometry%metric(:, :, 1, 2), g_st, &
            equilibrium%s(surface))
        call write_comparison("g_sz", geometry%metric(:, :, 1, 3), g_sz, &
            equilibrium%s(surface))
        call write_comparison("mod_B", geometry%mod_b, mod_b, &
            equilibrium%s(surface))
        call write_comparison("B_contra_t", geometry%b_contravariant(:, :, 1), &
            b_theta, equilibrium%s(surface))
        call write_comparison("B_contra_z", geometry%b_contravariant(:, :, 2), &
            b_zeta, equilibrium%s(surface))
        call write_comparison("B_covariant_t", geometry%b_covariant(:, :, 1), &
            covariant_theta, equilibrium%s(surface), &
            maxval(abs(covariant_zeta)))
        call write_comparison("B_covariant_z", geometry%b_covariant(:, :, 2), &
            covariant_zeta, equilibrium%s(surface))
        call write_comparison("II_tt", geometry%second_form(:, :, 1, 1), &
            second_tt, equilibrium%s(surface))
        call write_comparison("II_tz", geometry%second_form(:, :, 1, 2), &
            second_tz, equilibrium%s(surface))
        call write_comparison("II_zz", geometry%second_form(:, :, 2, 2), &
            second_zz, equilibrium%s(surface))
        call write_kernel_comparisons(surface)
        call write_comparison("kernel_drive_geometric", geometric_drive, &
            primitive_drive, equilibrium%s(surface))
        call write_component_identity(equilibrium%s(surface))
    end do

contains

    subroutine read_arguments(path, poloidal_points, toroidal_points)
        character(len=*), intent(out) :: path
        integer, intent(out) :: poloidal_points, toroidal_points
        character(len=64) :: value
        integer :: argument_status, parse_status

        call require(command_argument_count() == 1 &
            .or. command_argument_count() == 3, &
            "usage: gliss_primitive_geometry_trace EXPORT [NTHETA NZETA]")
        call get_command_argument(1, path, status=argument_status)
        call require(argument_status == 0 .and. len_trim(path) > 0, &
            "usage: gliss_primitive_geometry_trace EXPORT [NTHETA NZETA]")
        poloidal_points = 64
        toroidal_points = 64
        if (command_argument_count() == 1) return
        call get_command_argument(2, value, status=argument_status)
        call require(argument_status == 0, "NTHETA could not be read")
        read (value, *, iostat=parse_status) poloidal_points
        call require(parse_status == 0 .and. poloidal_points >= 8, &
            "NTHETA must be an integer of at least 8")
        call get_command_argument(3, value, status=argument_status)
        call require(argument_status == 0, "NZETA could not be read")
        read (value, *, iostat=parse_status) toroidal_points
        call require(parse_status == 0 .and. toroidal_points >= 8, &
            "NZETA must be an integer of at least 8")
    end subroutine read_arguments

    subroutine reconstruct_reference_fields(radial_surface)
        integer, intent(in) :: radial_surface

        call reconstruct_pair(equilibrium%jacobian, radial_surface, jacobian)
        call reconstruct_pair(jacobian_s_pair, radial_surface, jacobian_s)
        call reconstruct_pair(equilibrium%g_tt, radial_surface, g_tt)
        call reconstruct_pair(equilibrium%g_tz, radial_surface, g_tz)
        call reconstruct_pair(equilibrium%g_zz, radial_surface, g_zz)
        call reconstruct_pair(equilibrium%g_st, radial_surface, g_st)
        call reconstruct_pair(equilibrium%g_sz, radial_surface, g_sz)
        call reconstruct_pair(equilibrium%mod_b, radial_surface, mod_b)
        call reconstruct_pair(equilibrium%b_contravariant_theta, &
            radial_surface, b_theta)
        call reconstruct_pair(equilibrium%b_contravariant_zeta, &
            radial_surface, b_zeta)
        call reconstruct_pair(equilibrium%second_form_tt, radial_surface, &
            second_tt)
        call reconstruct_pair(equilibrium%second_form_tz, radial_surface, &
            second_tz)
        call reconstruct_pair(equilibrium%second_form_zz, radial_surface, &
            second_zz)
    end subroutine reconstruct_reference_fields

    subroutine write_kernel_comparisons(radial_surface)
        integer, intent(in) :: radial_surface
        character(len=16), parameter :: names(13) = [character(len=16) :: &
            "flux_t_slope", "flux_p_slope", "flux_t_curve", &
            "flux_p_curve", "current_i", "current_j", "signed_sqrtg", &
            "mod_b", "grad_s2", "j_dot_b", "pressure_slope", &
            "sigma_tilde", "beta_tilde"]
        integer :: field

        do field = 1, size(names)
            call write_comparison("kernel_" // trim(names(field)), &
                primitive_fields(:, :, field), &
                legacy_fields(:, :, field, radial_surface), &
                equilibrium%s(radial_surface))
        end do
        call write_comparison("kernel_drive", primitive_drive, &
            legacy_drive(:, :, radial_surface), equilibrium%s(radial_surface))
    end subroutine write_kernel_comparisons

    subroutine write_component_identity(s)
        real(dp), intent(in) :: s
        real(dp), allocatable :: direct(:, :, :), kernel(:, :, :)
        real(dp) :: r_contra(3), r_cov(3), values(6)
        integer :: j, k

        values = [0.37_dp, -0.21_dp, 0.16_dp, -0.11_dp, 0.29_dp, -0.24_dp]
        allocate (direct(size(theta), size(zeta), 3), &
            kernel(size(theta), size(zeta), 3))
        do k = 1, size(zeta)
            do j = 1, size(theta)
                call two_component_components( &
                    primitive_fields(j, k, 1), primitive_fields(j, k, 2), &
                    primitive_fields(j, k, 3), primitive_fields(j, k, 4), &
                    primitive_fields(j, k, 5), primitive_fields(j, k, 6), &
                    primitive_fields(j, k, 7), primitive_fields(j, k, 8), &
                    primitive_fields(j, k, 9), primitive_fields(j, k, 10), &
                    primitive_fields(j, k, 11), primitive_fields(j, k, 12), &
                    primitive_fields(j, k, 13), values(1), values(2), &
                    values(3), values(4), values(5), values(6), &
                    kernel(j, k, 1), kernel(j, k, 2), kernel(j, k, 3))
                r_contra(1) = (primitive_fields(j, k, 2) * values(3) &
                    + primitive_fields(j, k, 1) * values(4)) &
                    / primitive_fields(j, k, 7)
                r_contra(2) = (values(6) &
                    - primitive_fields(j, k, 4) * values(1) &
                    - primitive_fields(j, k, 2) * values(2)) &
                    / primitive_fields(j, k, 7)
                r_contra(3) = (-values(5) &
                    - primitive_fields(j, k, 3) * values(1) &
                    - primitive_fields(j, k, 1) * values(2)) &
                    / primitive_fields(j, k, 7)
                r_cov = matmul(geometry%metric(j, k, :, :), r_contra)
                direct(j, k, 1) = r_contra(1) &
                    / sqrt(primitive_fields(j, k, 9))
                direct(j, k, 2) = (primitive_fields(j, k, 6) * r_cov(3) &
                    - primitive_fields(j, k, 5) * r_cov(2)) &
                    / (primitive_fields(j, k, 7) &
                    * sqrt(primitive_fields(j, k, 9)) &
                    * primitive_fields(j, k, 8)) &
                    - primitive_fields(j, k, 10) * values(1) &
                    / (sqrt(primitive_fields(j, k, 9)) &
                    * primitive_fields(j, k, 8))
                direct(j, k, 3) = (primitive_fields(j, k, 6) * r_contra(2) &
                    + primitive_fields(j, k, 5) * r_contra(3) &
                    - primitive_fields(j, k, 11) * values(1)) &
                    / primitive_fields(j, k, 8)
            end do
        end do
        call write_comparison("kernel_C1_geometric", kernel(:, :, 1), &
            direct(:, :, 1), s)
        call write_comparison("kernel_C2_geometric", kernel(:, :, 2), &
            direct(:, :, 2), s)
        call write_comparison("kernel_C3_geometric", kernel(:, :, 3), &
            direct(:, :, 3), s)
    end subroutine write_component_identity

    subroutine reconstruct_pair(pair, radial_surface, values)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: radial_surface
        real(dp), allocatable, intent(out) :: values(:, :)

        call reconstruct_harmonic_grid(pair, radial_surface, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, theta, &
            zeta, values, discard_theta, discard_zeta, info)
        call require(info == reconstruction_ok, &
            "exported harmonic reconstruction failed")
    end subroutine reconstruct_pair

    subroutine write_comparison(name, actual, reference, s, external_scale)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :), reference(:, :), s
        real(dp), intent(in), optional :: external_scale
        real(dp) :: absolute_max, comparison_scale, difference_norm
        real(dp) :: normalized_max, reference_max, reference_norm
        real(dp) :: relative_l2, scaled_max

        call require(all(shape(actual) == shape(reference)), &
            "trace comparison shapes differ")
        call require(all(ieee_is_finite(actual)) &
            .and. all(ieee_is_finite(reference)), &
            "trace comparison contains nonfinite values")
        reference_norm = sqrt(sum(reference**2) / real(size(reference), dp))
        reference_max = maxval(abs(reference))
        difference_norm = sqrt(sum((actual - reference)**2) &
            / real(size(reference), dp))
        absolute_max = maxval(abs(actual - reference))
        relative_l2 = difference_norm / max(reference_norm, tiny(1.0_dp))
        normalized_max = absolute_max / max(reference_max, tiny(1.0_dp))
        comparison_scale = reference_max
        if (present(external_scale)) comparison_scale = external_scale
        comparison_scale = max(comparison_scale, tiny(1.0_dp))
        scaled_max = absolute_max / comparison_scale
        write (*, "(es24.16e3,',',a,7(',',es24.16e3))") s, trim(name), &
            relative_l2, normalized_max, scaled_max, absolute_max, &
            reference_norm, reference_max, comparison_scale
    end subroutine write_comparison

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "gliss_primitive_geometry_trace: " // message
        error stop 2
    end subroutine require

end program gliss_primitive_geometry_trace
