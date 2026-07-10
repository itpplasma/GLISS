program test_export_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    use family_assembly, only: lowest_family_eigenvalue, &
        surface_geometry_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none

    character(len=*), parameter :: fixture = "export_assembly.nc"
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(surface_geometry_t), allocatable :: geometry(:)
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    real(dp) :: resonant, stable, step
    integer :: info, i, ns

    call create_cylinder_fixture(fixture)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call require(info == reader_ok, "fixture was rejected")
    call build_kernel_geometry(equilibrium, 32, 16, fields, drive, info)
    call require(info == mercier_ok, "geometry build failed")

    ns = size(equilibrium%s)
    allocate (geometry(ns))
    do i = 1, ns
        geometry(i)%fields = fields(:, :, :, i)
        geometry(i)%drive = drive(:, :, i)
    end do
    step = 1.0_dp / real(ns, dp)

    call lowest_family_eigenvalue(geometry, [1], [1], step, resonant, &
        info)
    call require(info == 0, "resonant assembly failed")
    call lowest_family_eigenvalue(geometry, [2], [1], step, stable, info)
    call require(info == 0, "stable assembly failed")

    call require(resonant < 0.0_dp, &
        "resonant Suydam mode is not unstable through the export chain")
    call require(stable > 0.0_dp, &
        "non-resonant mode is not stable through the export chain")

    open (unit=13, file=fixture, status="old")
    close (13, status="delete")
    write (*, "(a)") "PASS"

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_export_assembly
