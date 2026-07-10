program test_export_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: b_poloidal, create_cylinder_fixture, &
        radius_of
    use family_assembly, only: lowest_family_eigenvalue, &
        surface_geometry_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none

    character(len=*), parameter :: fixture = "export_assembly.nc"
    character(len=*), parameter :: shifted_fixture = &
        "export_assembly_shifted.nc"
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: chart_shift = 0.05_dp
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    real(dp) :: resonant, stable, shifted_resonant, shifted_stable
    real(dp) :: radius, beta_expected, beta_error
    real(dp) :: coarse_gap, fine_gap
    integer :: info, i, ns

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
