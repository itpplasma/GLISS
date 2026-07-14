program test_export_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: b_poloidal, create_cylinder_fixture, &
        radius_of
    use family_assembly, only: lowest_family_eigenvalue, &
        surface_geometry_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use export_surface_geometry, only: build_angular_grids, &
        build_surface_kernel_fields, surface_data_t, surface_profiles_t
    use mercier_diagnostic, only: build_kernel_geometry, &
        mercier_invalid_input, mercier_ok
    implicit none

    character(len=*), parameter :: fixture = "export_assembly.nc"
    character(len=*), parameter :: shifted_fixture = &
        "export_assembly_shifted.nc"
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: chart_shift = 0.05_dp
    type(gvec_cas3d_equilibrium_t) :: equilibrium, hidden_invalid_metric
    type(gvec_cas3d_equilibrium_t) :: invalid_metric
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    real(dp) :: resonant, stable, shifted_resonant, shifted_stable
    real(dp) :: radius, beta_expected, beta_error
    real(dp) :: coarse_gap, fine_gap
    integer :: info, i, ns

    call check_local_covariant_fields()
    call create_cylinder_fixture(fixture)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "fixture was rejected")
    call require(.not. equilibrium%has_chart_metric, &
        "plain fixture must not announce chart metric data")
    call solve_modes(equilibrium, resonant, stable)

    call require(resonant < 0.0_dp, &
        "resonant Suydam mode is not unstable through the export chain")
    call require(stable > 0.0_dp, &
        "non-resonant mode is not stable through the export chain")

    call create_cylinder_fixture(shifted_fixture, chart_shift)
    call read_gvec_cas3d_file(shifted_fixture, equilibrium, info)
    call require(info == reader_ok, "shifted fixture was rejected")
    call require(equilibrium%has_chart_metric, &
        "shifted fixture must announce chart metric data")

    ns = size(equilibrium%s)
    call build_kernel_geometry(equilibrium, 32, 16, fields, drive, info)
    call require(info == mercier_ok, "shifted geometry build failed")
    beta_error = 0.0_dp
    do i = 1, ns
        radius = radius_of(equilibrium%s(i))
        beta_expected = -chart_shift * 2.0_dp * pi * radius &
            * b_poloidal(radius)
        beta_error = max(beta_error, maxval(abs(fields(:, :, 13, i) &
            - beta_expected)))
    end do
    call require(beta_error < 1.0e-12_dp, &
        "chart beta-tilde does not match the physical covariant B_s")
    invalid_metric = equilibrium
    invalid_metric%g_zz%cosine = 0.0_dp
    invalid_metric%g_zz%sine = 0.0_dp
    call build_kernel_geometry(invalid_metric, 32, 16, fields, drive, info)
    call require(info == mercier_invalid_input, &
        "nonpositive tangential metric was accepted")
    hidden_invalid_metric = equilibrium
    hidden_invalid_metric%poloidal_modes(2) = 4
    hidden_invalid_metric%g_zz%sine(:, 2, 1) = &
        2.0_dp * hidden_invalid_metric%g_zz%cosine(:, 1, 1)
    call build_kernel_geometry(hidden_invalid_metric, 8, 8, fields, drive, &
        info)
    call require(info == mercier_invalid_input, &
        "between-node negative metric was accepted")

    call solve_modes(equilibrium, shifted_resonant, shifted_stable)
    call require(shifted_resonant < 0.0_dp, &
        "resonant Suydam mode is not unstable in the shifted chart")
    call require(shifted_stable > 0.0_dp, &
        "non-resonant mode is not stable in the shifted chart")

    coarse_gap = abs(shifted_resonant - resonant) / abs(resonant)
    write (error_unit, "(a, es13.5)") "coarse gauge gap ", coarse_gap
    call require(coarse_gap < 1.0e-3_dp, &
        "shifted-chart eigenvalue deviates beyond discretization error")
    call measure_gap(2 * size(equilibrium%s), fine_gap)
    call require(fine_gap < coarse_gap / 2.5_dp, &
        "gauge gap does not vanish under radial refinement")

    open (unit=13, file=fixture, status="old")
    close (13, status="delete")
    open (unit=13, file=shifted_fixture, status="old")
    close (13, status="delete")
    write (*, "(a)") "PASS"

contains

    subroutine check_local_covariant_fields()
        integer, parameter :: n_theta = 8, n_zeta = 8
        integer, parameter :: poloidal_modes(2) = [0, 1]
        integer, parameter :: toroidal_modes(1) = [0]
        type(surface_data_t) :: surface
        type(surface_profiles_t) :: profiles
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp) :: jacobian_slope(n_theta, n_zeta)
        real(dp) :: fields(n_theta, n_zeta, 13)
        real(dp) :: drive(n_theta, n_zeta)
        real(dp) :: expected_theta(n_theta, n_zeta)
        real(dp) :: expected_zeta(n_theta, n_zeta)
        real(dp) :: expected_j_dot_b(n_theta, n_zeta)
        integer :: build_info, j

        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        allocate (surface%jacobian(n_theta, n_zeta), &
            surface%g_tt(n_theta, n_zeta), surface%g_tz(n_theta, n_zeta), &
            surface%g_zz(n_theta, n_zeta), &
            surface%b_theta(n_theta, n_zeta), &
            surface%b_zeta(n_theta, n_zeta), &
            surface%g_st(n_theta, n_zeta), surface%g_sz(n_theta, n_zeta), &
            surface%mod_b(n_theta, n_zeta), &
            surface%area_element(n_theta, n_zeta))
        do j = 1, n_zeta
            surface%g_tt(:, j) = 2.0_dp + 0.2_dp &
                * cos(2.0_dp * pi * theta)
            surface%g_tz(:, j) = 0.1_dp * sin(2.0_dp * pi * theta)
            surface%g_zz(:, j) = 3.0_dp
            surface%b_theta(:, j) = 0.4_dp
            surface%b_zeta(:, j) = -0.7_dp
        end do
        surface%jacobian = 1.8_dp
        surface%g_st = 0.0_dp
        surface%g_sz = 0.0_dp
        expected_theta = surface%g_tt * surface%b_theta &
            + surface%g_tz * surface%b_zeta
        expected_zeta = surface%g_tz * surface%b_theta &
            + surface%g_zz * surface%b_zeta
        surface%mod_b = sqrt(surface%b_theta * expected_theta &
            + surface%b_zeta * expected_zeta)
        surface%area_element = 1.0_dp
        profiles = surface_profiles_t(-1.26_dp, 0.72_dp, 0.0_dp, 0.0_dp, &
            sum(expected_theta) / real(n_theta * n_zeta, dp), &
            sum(expected_zeta) / real(n_theta * n_zeta, dp), &
            0.3_dp, -0.2_dp, -0.26_dp / (4.0_dp * pi * 1.0e-7_dp))
        expected_j_dot_b = (0.2_dp * expected_theta &
            + 0.3_dp * expected_zeta) / surface%jacobian
        jacobian_slope = 0.0_dp

        call build_surface_kernel_fields(poloidal_modes, toroidal_modes, &
            .true., surface, profiles, jacobian_slope, theta, zeta, fields, &
            drive, build_info)
        call require(build_info == mercier_ok, &
            "local covariant-field fixture was rejected")
        call require(maxval(abs(fields(:, :, 5) - expected_zeta)) &
            < 1.0e-14_dp, "kernel discarded local covariant B_zeta")
        call require(maxval(abs(fields(:, :, 6) - expected_theta)) &
            < 1.0e-14_dp, "kernel discarded local covariant B_theta")
        call require(maxval(abs(fields(:, :, 1) * fields(:, :, 5) &
            + fields(:, :, 2) * fields(:, :, 6) &
            - surface%mod_b**2 * surface%jacobian)) < 1.0e-14_dp, &
            "kernel covariant and contravariant fields are inconsistent")
        call require(maxval(abs(fields(:, :, 10) - expected_j_dot_b)) &
            < 1.0e-14_dp, "kernel j dot B discarded local covariant field")
    end subroutine check_local_covariant_fields

    subroutine measure_gap(surfaces, gap)
        integer, intent(in) :: surfaces
        real(dp), intent(out) :: gap
        character(len=*), parameter :: plain = "export_gap_plain.nc"
        character(len=*), parameter :: rotated = "export_gap_shift.nc"
        type(gvec_cas3d_equilibrium_t) :: fine
        real(dp) :: base_resonant, base_stable
        real(dp) :: turn_resonant, turn_stable
        integer :: read_info

        call create_cylinder_fixture(plain, surfaces=surfaces)
        call read_gvec_cas3d_file(plain, fine, read_info)
        call require(read_info == reader_ok, "fine fixture was rejected")
        call solve_modes(fine, base_resonant, base_stable)
        call create_cylinder_fixture(rotated, chart_shift, surfaces)
        call read_gvec_cas3d_file(rotated, fine, read_info)
        call require(read_info == reader_ok, &
            "fine shifted fixture was rejected")
        call solve_modes(fine, turn_resonant, turn_stable)
        gap = abs(turn_resonant - base_resonant) / abs(base_resonant)
        open (unit=14, file=plain, status="old")
        close (14, status="delete")
        open (unit=14, file=rotated, status="old")
        close (14, status="delete")
    end subroutine measure_gap

    subroutine solve_modes(equilibrium, resonant, stable)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(out) :: resonant, stable
        type(surface_geometry_t), allocatable :: geometry(:)
        real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
        real(dp) :: step
        integer :: info, i, ns

        call build_kernel_geometry(equilibrium, 32, 16, fields, drive, &
            info)
        call require(info == mercier_ok, "geometry build failed")
        ns = size(equilibrium%s)
        allocate (geometry(ns))
        do i = 1, ns
            geometry(i)%fields = fields(:, :, :, i)
            geometry(i)%drive = drive(:, :, i)
        end do
        step = 1.0_dp / real(ns, dp)
        call lowest_family_eigenvalue(geometry, [1], [1], step, &
            resonant, info)
        call require(info == 0, "resonant assembly failed")
        call lowest_family_eigenvalue(geometry, [2], [1], step, stable, &
            info)
        call require(info == 0, "stable assembly failed")
    end subroutine solve_modes

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_export_assembly
