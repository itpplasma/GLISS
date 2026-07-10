program test_suydam_marginal
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    use family_assembly, only: family_negative_count, surface_geometry_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    implicit none

    ! Suydam threshold of the pressure-scaled fixture family at the
    ! iota = 1 surface (derivations/suydam_marginal_threshold.wl).
    ! With the criterion violated the indicial roots are complex and
    ! every resonant mode is unstable, so the continuum stability
    ! boundary sits at kappa* for each m; discrete boundaries approach
    ! it from above and slowly (the unstable layer narrows
    ! exponentially near the threshold), and the interval-constant
    ! tangential space keeps the discrete operator on the stable side
    ! of the independent pointwise-eliminated reduction.
    real(dp), parameter :: kappa_star = 0.31454704733350229_dp
    real(dp), parameter :: threshold = -1.0e-8_dp
    integer, parameter :: n_theta = 32, n_zeta = 32
    integer, parameter :: resolutions(3) = [96, 192, 384]
    character(len=*), parameter :: fixture = "suydam_marginal.nc"
    real(dp) :: boundary(3, 2)
    integer :: modes(2), level, which
    logical :: ok

    modes = [4, 6]
    do which = 1, 2
        do level = 1, 3
            call bisect_marginal(modes(which), resolutions(level), &
                boundary(level, which), ok)
            call require(ok, "no stability transition bracketed")
        end do
        write (error_unit, "(a, i0, a, 3f10.6)") "m=", modes(which), &
            " boundaries ", boundary(:, which)
    end do

    do which = 1, 2
        call require(boundary(1, which) > boundary(2, which) .and. &
            boundary(2, which) > boundary(3, which), &
            "the discrete boundary does not decrease under refinement")
        call require(boundary(3, which) > kappa_star, &
            "the discrete boundary undershoots the Suydam threshold")
        call require(boundary(3, which) - kappa_star < 0.75_dp * &
            (boundary(1, which) - kappa_star), &
            "the boundary does not contract toward the threshold")
    end do
    write (*, "(a)") "PASS"

contains

    subroutine bisect_marginal(mode, surfaces, marginal, bracketed)
        integer, intent(in) :: mode, surfaces
        real(dp), intent(out) :: marginal
        logical, intent(out) :: bracketed
        real(dp) :: low, high, middle
        integer :: iteration

        low = 0.10_dp
        high = 1.00_dp
        bracketed = .not. unstable(low, mode, surfaces) .and. &
            unstable(high, mode, surfaces)
        if (.not. bracketed) return
        do iteration = 1, 16
            middle = 0.5_dp * (low + high)
            if (unstable(middle, mode, surfaces)) then
                high = middle
            else
                low = middle
            end if
        end do
        marginal = 0.5_dp * (low + high)
    end subroutine bisect_marginal

    function unstable(fraction, mode, surfaces) result(is_unstable)
        real(dp), intent(in) :: fraction
        integer, intent(in) :: mode, surfaces
        logical :: is_unstable
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
        type(surface_geometry_t), allocatable :: geometry(:)
        real(dp) :: step
        integer :: info, negatives, i, ns

        call create_cylinder_fixture(fixture, surfaces=surfaces, &
            pressure_scale=fraction)
        call read_gvec_cas3d_file(fixture, equilibrium, info)
        call require(info == reader_ok, "fixture was rejected")
        call build_kernel_geometry(equilibrium, n_theta, n_zeta, fields, &
            drive, info)
        call require(info == mercier_ok, "geometry build failed")
        ns = size(equilibrium%s)
        allocate (geometry(ns))
        do i = 1, ns
            geometry(i)%fields = fields(:, :, :, i)
            geometry(i)%drive = drive(:, :, i)
        end do
        step = 1.0_dp / real(ns, dp)
        call family_negative_count(geometry, [mode], [mode], step, &
            threshold, negatives, info)
        call require(info == 0, "inertia count failed")
        is_unstable = negatives > 0
        open (unit=17, file=fixture, status="old")
        close (17, status="delete")
    end function unstable

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_suydam_marginal
