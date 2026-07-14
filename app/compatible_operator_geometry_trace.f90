module compatible_operator_geometry_trace
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_operator_trace_types, only: compatible_cell_trace_t
    use export_surface_geometry, only: build_angular_grids
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use primitive_equilibrium_spline, only: evaluate_primitive_equilibrium, &
        fit_primitive_equilibrium, primitive_equilibrium_ok, &
        primitive_equilibrium_spline_t
    use primitive_geometry_grid, only: primitive_geometry_grid_t
    implicit none
    private

    integer, parameter, public :: operator_geometry_trace_ok = 0
    integer, parameter, public :: operator_geometry_trace_invalid = -1
    integer, parameter, public :: operator_geometry_trace_evaluation_error = -2

    public :: write_compatible_operator_geometry

contains

    subroutine write_compatible_operator_geometry(equilibrium, traces, &
            n_theta, n_zeta, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(compatible_cell_trace_t), intent(in) :: traces(:)
        integer, intent(in) :: n_theta, n_zeta
        integer, intent(out) :: info
        type(primitive_equilibrium_spline_t) :: spline
        type(primitive_geometry_grid_t) :: geometry
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp) :: pressure, pressure_slope
        integer :: cell_index, local_info, point

        info = operator_geometry_trace_invalid
        if (n_theta < 1 .or. n_zeta < 1 .or. size(traces) < 1) return
        call fit_primitive_equilibrium(equilibrium, spline, local_info)
        if (local_info /= primitive_equilibrium_ok) then
            info = operator_geometry_trace_evaluation_error
            return
        end if
        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        do cell_index = 1, size(traces)
            if (.not. allocated(traces(cell_index)%points)) return
            do point = 1, size(traces(cell_index)%points)
                call evaluate_primitive_equilibrium(spline, &
                    traces(cell_index)%points(point)%coordinate, theta, zeta, &
                    geometry, pressure, pressure_slope, local_info)
                if (local_info /= primitive_equilibrium_ok) then
                    info = operator_geometry_trace_evaluation_error
                    return
                end if
                if (.not. geometry_matches_kernel(geometry, &
                    traces(cell_index)%points(point)%fields)) return
                call write_geometry(traces(cell_index)%cell, point, geometry)
            end do
        end do
        info = operator_geometry_trace_ok
    end subroutine write_compatible_operator_geometry

    function geometry_matches_kernel(geometry, fields) result(matches)
        type(primitive_geometry_grid_t), intent(in) :: geometry
        real(dp), intent(in) :: fields(:, :, :)
        logical :: matches
        real(dp) :: inverse(6), scale
        integer :: j, k

        matches = size(fields, 1) == size(geometry%signed_jacobian, 1) &
            .and. size(fields, 2) == size(geometry%signed_jacobian, 2) &
            .and. size(fields, 3) >= 9
        if (.not. matches) return
        do k = 1, size(fields, 2)
            do j = 1, size(fields, 1)
                call invert_symmetric_metric( &
                    geometry%metric(j, k, 1, 1), &
                    geometry%metric(j, k, 1, 2), &
                    geometry%metric(j, k, 1, 3), &
                    geometry%metric(j, k, 2, 2), &
                    geometry%metric(j, k, 2, 3), &
                    geometry%metric(j, k, 3, 3), inverse, matches)
                if (.not. matches) return
                scale = max(abs(fields(j, k, 7)), tiny(1.0_dp))
                matches = abs(geometry%signed_jacobian(j, k) &
                    - fields(j, k, 7)) <= 64.0_dp * epsilon(1.0_dp) * scale
                if (.not. matches) return
                scale = max(abs(fields(j, k, 8)), tiny(1.0_dp))
                matches = abs(geometry%mod_b(j, k) - fields(j, k, 8)) &
                    <= 64.0_dp * epsilon(1.0_dp) * scale
                if (.not. matches) return
                scale = max(abs(fields(j, k, 9)), tiny(1.0_dp))
                matches = abs(inverse(1) - fields(j, k, 9)) &
                    <= 256.0_dp * epsilon(1.0_dp) * scale
                if (.not. matches) return
            end do
        end do
    end function geometry_matches_kernel

    subroutine write_geometry(cell, point, geometry)
        integer, intent(in) :: cell, point
        type(primitive_geometry_grid_t), intent(in) :: geometry
        real(dp) :: inverse(6)
        integer :: j, k
        logical :: valid

        do k = 1, size(geometry%signed_jacobian, 2)
            do j = 1, size(geometry%signed_jacobian, 1)
                call invert_symmetric_metric( &
                    geometry%metric(j, k, 1, 1), &
                    geometry%metric(j, k, 1, 2), &
                    geometry%metric(j, k, 1, 3), &
                    geometry%metric(j, k, 2, 2), &
                    geometry%metric(j, k, 2, 3), &
                    geometry%metric(j, k, 3, 3), inverse, valid)
                if (.not. valid) return
                write (*, "(a,4(i0,','),24(es24.16e3,:,','))") &
                    "GEOMETRY,", cell, point, j, k, &
                    geometry%metric(j, k, 1, 1), &
                    geometry%metric(j, k, 1, 2), &
                    geometry%metric(j, k, 1, 3), &
                    geometry%metric(j, k, 2, 2), &
                    geometry%metric(j, k, 2, 3), &
                    geometry%metric(j, k, 3, 3), &
                    geometry%metric_radial(j, k, 1, 1), &
                    geometry%metric_radial(j, k, 1, 2), &
                    geometry%metric_radial(j, k, 1, 3), &
                    geometry%metric_radial(j, k, 2, 2), &
                    geometry%metric_radial(j, k, 2, 3), &
                    geometry%metric_radial(j, k, 3, 3), inverse, &
                    geometry%signed_jacobian(j, k), &
                    geometry%b_contravariant(j, k, 1), &
                    geometry%b_contravariant(j, k, 2), &
                    geometry%b_covariant(j, k, 1), &
                    geometry%b_covariant(j, k, 2), geometry%mod_b(j, k)
            end do
        end do
    end subroutine write_geometry

    pure subroutine invert_symmetric_metric(g_ss, g_st, g_sz, g_tt, g_tz, &
            g_zz, inverse, valid)
        real(dp), intent(in) :: g_ss, g_st, g_sz, g_tt, g_tz, g_zz
        real(dp), intent(out) :: inverse(6)
        logical, intent(out) :: valid
        real(dp) :: determinant

        determinant = g_ss * (g_tt * g_zz - g_tz * g_tz) &
            - g_st * (g_st * g_zz - g_sz * g_tz) &
            + g_sz * (g_st * g_tz - g_sz * g_tt)
        valid = ieee_is_finite(determinant) .and. determinant > 0.0_dp
        inverse = 0.0_dp
        if (.not. valid) return
        inverse(1) = (g_tt * g_zz - g_tz * g_tz) / determinant
        inverse(2) = (g_sz * g_tz - g_st * g_zz) / determinant
        inverse(3) = (g_st * g_tz - g_sz * g_tt) / determinant
        inverse(4) = (g_ss * g_zz - g_sz * g_sz) / determinant
        inverse(5) = (g_st * g_sz - g_ss * g_tz) / determinant
        inverse(6) = (g_ss * g_tt - g_st * g_st) / determinant
        valid = all(ieee_is_finite(inverse))
    end subroutine invert_symmetric_metric

end module compatible_operator_geometry_trace
