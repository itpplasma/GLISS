module plasma_vacuum_boundary
    ! Vacuum energy of the physical free-boundary problem on the plasma edge.
    !
    ! On the edge s = 1 of the full-torus chart (s, u, v), u = theta and
    ! v = (zeta_period + period) / N_FP, the perturbed normal flux is
    !     Q . grad(s) sqrt(g) du dv = (fp d/du + ft d/dv) xi^s du dv,
    ! with fp = -chi'(1) and ft = -Phi'(1) in the sign convention of
    ! cartesian_primitive_geometry. Each normal trial cos or sin(2 pi (m u -
    ! n v)) prescribes this flux; scalar_potential_vacuum returns the form E
    ! of integral_V |B|^2 dV over the vacuum region (units mu0 = 1), so the
    ! edge block of the stiffness K = 2 delta W is E / mu0. Tangential
    ! continuity of B leaves no surface-current term.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cartesian_harmonic_spline, only: cartesian_harmonic_ok, &
        cartesian_jet_grid_t, evaluate_cartesian_harmonic_spline
    use field_periodic_cartesian, only: convert_field_periodic_jet, &
        field_periodic_cartesian_ok
    use fourier_phase_kind, only: phase_cosine, phase_sine
    use physical_constants, only: vacuum_permeability
    use primitive_equilibrium_spline, only: primitive_equilibrium_spline_t
    use radial_cubic_spline, only: evaluate_radial_cubic_spline_field, &
        radial_cubic_spline_ok
    use scalar_potential_vacuum, only: assemble_exterior_vacuum, &
        edge_triangle_centroids, vacuum_bie_not_nested, vacuum_bie_ok
    use trial_space_topology, only: trial_component_normal, &
        trial_space_topology_t
    implicit none
    private

    integer, parameter, public :: plasma_vacuum_ok = 0
    integer, parameter, public :: plasma_vacuum_invalid = -1
    integer, parameter, public :: plasma_vacuum_geometry_error = -2
    integer, parameter, public :: plasma_vacuum_underresolved = -3
    integer, parameter, public :: plasma_vacuum_wall_not_nested = -4
    integer, parameter, public :: plasma_vacuum_singular = -5

    integer, parameter, public :: vacuum_wall_none = 0
    integer, parameter, public :: vacuum_wall_conformal = 1
    integer, parameter, public :: vacuum_wall_surface = 2

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    type, public :: plasma_vacuum_model_t
        ! Nodes of the full-torus edge mesh in u and v.
        integer :: nu = 0
        integer :: nv = 0
        integer :: wall_kind = vacuum_wall_none
        ! Normal offset of a conformal wall, in metres.
        real(dp) :: wall_distance = 0.0_dp
        ! Explicit wall nodes (3, wall_nu, wall_nv), Cartesian metres.
        real(dp), allocatable :: wall(:, :, :)
    end type plasma_vacuum_model_t

    public :: build_conformal_wall, build_plasma_boundary_surface
    public :: build_vacuum_edge_block, valid_plasma_vacuum_model

contains

    pure logical function valid_plasma_vacuum_model(model) result(valid)
        type(plasma_vacuum_model_t), intent(in) :: model

        valid = .false.
        if (model%nu < 3 .or. model%nv < 3) return
        if (model%nu > huge(1) / model%nv) return
        select case (model%wall_kind)
        case (vacuum_wall_none)
        case (vacuum_wall_conformal)
            if (.not. ieee_is_finite(model%wall_distance)) return
            if (model%wall_distance <= 0.0_dp) return
        case (vacuum_wall_surface)
            if (.not. allocated(model%wall)) return
            if (size(model%wall, 1) /= 3 .or. size(model%wall, 2) < 3 &
                .or. size(model%wall, 3) < 3) return
            if (.not. all(ieee_is_finite(model%wall))) return
        case default
            return
        end select
        valid = .true.
    end function valid_plasma_vacuum_model

    subroutine build_plasma_boundary_surface(spline, nu, nv, surface, &
            normal, fp, ft, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        integer, intent(in) :: nu, nv
        ! Nodes (3, nu, nv) and outward unit normals on the plasma edge.
        real(dp), allocatable, intent(out) :: surface(:, :, :), normal(:, :, :)
        real(dp), intent(out) :: fp, ft
        integer, intent(out) :: info
        type(cartesian_jet_grid_t) :: jet
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp) :: area(3), length, orientation, profiles(3), seconds(3)
        real(dp) :: tangent_u(3), tangent_v(3)
        real(dp) :: slopes(3)
        integer :: i, k, local_info

        info = plasma_vacuum_invalid
        fp = 0.0_dp
        ft = 0.0_dp
        if (spline%field_periods < 1 .or. nu < 3 .or. nv < 3) return
        allocate (theta(nu), zeta(nv))
        do i = 1, nu
            theta(i) = real(i - 1, dp) / real(nu, dp)
        end do
        ! The period-frame angle spans the N_FP periods of the full torus.
        do k = 1, nv
            zeta(k) = real(spline%field_periods * (k - 1), dp) / real(nv, dp)
        end do
        info = plasma_vacuum_geometry_error
        call evaluate_cartesian_harmonic_spline(spline%radial_grid, &
            spline%position, 1.0_dp, theta, zeta, jet, local_info)
        if (local_info /= cartesian_harmonic_ok) return
        call convert_field_periodic_jet(zeta, spline%field_periods, &
            spline%winding, jet, local_info)
        if (local_info /= field_periodic_cartesian_ok) return
        call evaluate_radial_cubic_spline_field(spline%radial_grid, &
            spline%profiles, 1.0_dp, profiles, slopes, seconds, local_info)
        if (local_info /= radial_cubic_spline_ok) return
        allocate (surface(3, nu, nv), normal(3, nu, nv))
        orientation = 0.0_dp
        do k = 1, nv
            do i = 1, nu
                surface(:, i, k) = jet%value(i, k, :)
                tangent_u = jet%poloidal(i, k, :)
                tangent_v = jet%toroidal(i, k, :)
                area = cross(tangent_u, tangent_v)
                length = norm2(area)
                if (.not. ieee_is_finite(length) .or. length <= tiny(1.0_dp)) &
                    return
                normal(:, i, k) = area / length
                ! Divergence theorem: sum r . dA is three times the volume.
                orientation = orientation + dot_product(surface(:, i, k), area)
            end do
        end do
        if (.not. ieee_is_finite(orientation) .or. orientation == 0.0_dp) return
        if (orientation < 0.0_dp) normal = -normal
        ! A frame whose rotation adds to the toroidal harmonics can trace the
        ! torus more than once; such a surface bounds no vacuum region.
        if (toroidal_turns(surface) /= 1) return
        fp = -slopes(2)
        ft = -slopes(1)
        if (.not. ieee_is_finite(fp) .or. .not. ieee_is_finite(ft)) return
        if (fp == 0.0_dp .and. ft == 0.0_dp) return
        info = plasma_vacuum_ok
    end subroutine build_plasma_boundary_surface

    pure subroutine build_conformal_wall(surface, normal, distance, wall, info)
        real(dp), intent(in) :: surface(:, :, :), normal(:, :, :), distance
        real(dp), allocatable, intent(out) :: wall(:, :, :)
        integer, intent(out) :: info

        info = plasma_vacuum_invalid
        if (any(shape(surface) /= shape(normal))) return
        if (.not. ieee_is_finite(distance) .or. distance <= 0.0_dp) return
        allocate (wall(3, size(surface, 2), size(surface, 3)))
        wall = surface + distance * normal
        info = plasma_vacuum_ok
    end subroutine build_conformal_wall

    subroutine build_vacuum_edge_block(spline, model, topology, block, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(plasma_vacuum_model_t), intent(in) :: model
        type(trial_space_topology_t), intent(in) :: topology
        ! Trial-by-trial vacuum stiffness E / mu0 of the edge normal
        ! coefficients; inactive normal trials have zero rows and columns.
        real(dp), allocatable, intent(out) :: block(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: flux(:, :), normal(:, :, :)
        real(dp), allocatable :: surface(:, :, :), wall(:, :, :)
        real(dp) :: fp, ft
        integer :: local_info

        info = plasma_vacuum_invalid
        if (.not. valid_plasma_vacuum_model(model)) return
        if (.not. valid_normal_trials(topology)) return
        if (.not. resolves_trials(model%nu, model%nv, topology)) then
            info = plasma_vacuum_underresolved
            return
        end if
        call build_plasma_boundary_surface(spline, model%nu, model%nv, &
            surface, normal, fp, ft, info)
        if (info /= plasma_vacuum_ok) return
        call trial_edge_flux(model%nu, model%nv, fp, ft, topology, flux)
        select case (model%wall_kind)
        case (vacuum_wall_conformal)
            call build_conformal_wall(surface, normal, model%wall_distance, &
                wall, info)
            if (info /= plasma_vacuum_ok) return
            call assemble_exterior_vacuum(surface, flux, block, local_info, &
                wall)
        case (vacuum_wall_surface)
            call assemble_exterior_vacuum(surface, flux, block, local_info, &
                model%wall)
        case default
            call assemble_exterior_vacuum(surface, flux, block, local_info)
        end select
        if (local_info == vacuum_bie_not_nested) then
            info = plasma_vacuum_wall_not_nested
            return
        else if (local_info /= vacuum_bie_ok) then
            info = plasma_vacuum_singular
            return
        end if
        block = block / vacuum_permeability
        info = plasma_vacuum_ok
    end subroutine build_vacuum_edge_block

    ! B.n dA / (du dv) = (fp d/du + ft d/dv) of each normal trial at the
    ! centroids of the edge triangles.
    subroutine trial_edge_flux(nu, nv, fp, ft, topology, flux)
        integer, intent(in) :: nu, nv
        real(dp), intent(in) :: fp, ft
        type(trial_space_topology_t), intent(in) :: topology
        real(dp), allocatable, intent(out) :: flux(:, :)
        real(dp), allocatable :: centroid(:, :)
        real(dp) :: phase, rate
        integer :: t, trial

        allocate (centroid(2, 2 * nu * nv), &
            flux(2 * nu * nv, size(topology%poloidal)))
        call edge_triangle_centroids(nu, nv, centroid)
        flux = 0.0_dp
        do trial = 1, size(topology%poloidal)
            if (.not. topology%active(trial_component_normal, trial)) cycle
            ! v spans the full torus, so physical n has no N_FP factor.
            rate = two_pi * (fp * real(topology%poloidal(trial), dp) &
                - ft * real(topology%toroidal(trial), dp))
            do t = 1, size(centroid, 2)
                phase = two_pi * (real(topology%poloidal(trial), dp) &
                    * centroid(1, t) - real(topology%toroidal(trial), dp) &
                    * centroid(2, t))
                if (topology%normal_phase(trial) == phase_cosine) &
                    flux(t, trial) = -rate * sin(phase)
                if (topology%normal_phase(trial) == phase_sine) &
                    flux(t, trial) = rate * cos(phase)
            end do
        end do
    end subroutine trial_edge_flux

    pure logical function valid_normal_trials(topology) result(valid)
        type(trial_space_topology_t), intent(in) :: topology
        integer :: trial

        valid = .false.
        if (.not. allocated(topology%poloidal) &
            .or. .not. allocated(topology%toroidal) &
            .or. .not. allocated(topology%normal_phase) &
            .or. .not. allocated(topology%active)) return
        if (size(topology%poloidal) < 1) return
        do trial = 1, size(topology%poloidal)
            if (.not. topology%active(trial_component_normal, trial)) cycle
            if (topology%normal_phase(trial) /= phase_cosine &
                .and. topology%normal_phase(trial) /= phase_sine) return
        end do
        valid = .true.
    end function valid_normal_trials

    ! The edge mesh must sample every normal trial above its Nyquist rate.
    pure logical function resolves_trials(nu, nv, topology) result(resolved)
        integer, intent(in) :: nu, nv
        type(trial_space_topology_t), intent(in) :: topology
        integer :: trial

        resolved = .false.
        do trial = 1, size(topology%poloidal)
            if (.not. topology%active(trial_component_normal, trial)) cycle
            if (nu <= 2 * abs(topology%poloidal(trial))) return
            if (nv <= 2 * abs(topology%toroidal(trial))) return
        end do
        resolved = .true.
    end function resolves_trials

    pure integer function toroidal_turns(surface) result(turns)
        real(dp), intent(in) :: surface(:, :, :)
        real(dp) :: angle, centre(2), previous, total
        integer :: k, nv

        nv = size(surface, 3)
        total = 0.0_dp
        previous = 0.0_dp
        do k = 1, nv + 1
            centre = sum(surface(1:2, :, modulo(k - 1, nv) + 1), dim=2)
            angle = atan2(centre(2), centre(1))
            if (k > 1) total = total + modulo(angle - previous &
                + 0.5_dp * two_pi, two_pi) - 0.5_dp * two_pi
            previous = angle
        end do
        turns = abs(nint(total / two_pi))
    end function toroidal_turns

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross

end module plasma_vacuum_boundary
