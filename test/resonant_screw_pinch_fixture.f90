! Exact cylindrical surface data with global production radial FEEC assembly.
! This is not a toroidal equilibrium import fixture.  Physical radius is
! r=a*sqrt(s), theta and zeta are turns, and physical axial period is L.
module resonant_screw_pinch_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use compatible_problem_assembly_support, only: apply_stored_power, &
        apply_tangential_axis_weight, build_active_indices, &
        replicate_indexed_values, scatter_matrix
    use compatible_radial_quadrature, only: accurate_nodes, accurate_weights, &
        axis_quadrature_points, build_axis_quadrature, build_constraint_quadrature
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    use radial_feec_complex, only: build_radial_feec_complex, &
        evaluate_radial_feec_complex, radial_feec_complex_t
    implicit none
    private

    real(dp), parameter, public :: pinch_radius = 0.7_dp
    real(dp), parameter, public :: pinch_blinear = 0.8_dp
    real(dp), parameter, public :: pinch_bz = 1.2_dp
    real(dp), parameter, public :: pinch_pressure_edge = 17000.0_dp
    real(dp), parameter, public :: pinch_density = 1.7e-7_dp
    real(dp), parameter, public :: pinch_lambda_scale = pinch_bz**2 &
        / (vacuum_permeability * pinch_density * pinch_radius**2)
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter, public :: pinch_mode = 2

    public :: assemble_resonant_screw_pinch

contains

    subroutine assemble_resonant_screw_pinch(cells, degree, parity, gamma, &
            stiffness, mass, normal_dofs, tangential_dofs)
        integer, intent(in) :: cells, degree, parity
        real(dp), intent(in) :: gamma
        real(dp), allocatable, intent(out) :: stiffness(:, :), mass(:, :)
        integer, intent(out) :: normal_dofs, tangential_dofs
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: breaks(:), nodes(:), weights(:)
        real(dp) :: s, weight
        integer :: index, status, dimension, pass, cell, point

        allocate (breaks(cells + 1))
        do index = 1, cells + 1
            breaks(index) = real(index - 1, dp) / real(cells, dp)
        end do
        call build_radial_feec_complex(breaks, degree, .true., .true., &
            complex, status)
        if (status /= 0) error stop 'resonant radial complex failed'
        normal_dofs = complex%h1_dofs
        tangential_dofs = complex%l2_dofs
        dimension = normal_dofs + 2 * tangential_dofs
        allocate (stiffness(dimension, dimension), source=0.0_dp)
        allocate (mass(dimension, dimension), source=0.0_dp)
        ! Match the current global production schedule: all stiffness terms
        ! and mass on the axis cell use its Gauss rule in sqrt(s). Only
        ! interior cells split accurate and reduced constraint terms.
        call build_axis_quadrature(breaks(2), &
            axis_quadrature_points(degree, pinch_mode), nodes, weights, status)
        if (status /= 0) error stop 'resonant axis quadrature failed'
        do point = 1, size(nodes)
            call accumulate(complex, nodes(point), weights(point), parity, &
                gamma, 0, stiffness, mass)
        end do
        do pass = 1, 2
            if (pass == 1) then
                nodes = accurate_nodes
                weights = accurate_weights
            else
                call build_constraint_quadrature(degree, nodes, weights, status)
                if (status /= 0) error stop 'resonant constraint rule failed'
            end if
            do cell = 2, cells
                do point = 1, size(nodes)
                    s = (real(cell, dp) - 0.5_dp &
                        + 0.5_dp * nodes(point)) / real(cells, dp)
                    weight = 0.5_dp * weights(point) / real(cells, dp)
                    call accumulate(complex, s, weight, parity, gamma, pass, &
                        stiffness, mass)
                end do
            end do
        end do
    end subroutine assemble_resonant_screw_pinch

    subroutine accumulate(complex, s, weight, parity, gamma, pass, stiffness, mass)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: s, weight, gamma
        integer, intent(in) :: parity, pass
        real(dp), intent(inout) :: stiffness(:, :), mass(:, :)
        real(dp) :: fields(8, 8, 13), zeros(8, 8), drive(8, 8), gamma_p(8, 8)
        real(dp), allocatable :: h1(:), dh1(:), l2(:), h(:, :), dh(:, :)
        real(dp), allocatable :: eta(:, :), nu(:, :)
        real(dp), allocatable :: local_k(:, :), local_m(:, :), terms(:, :, :)
        integer, allocatable :: hi(:), li(:), map(:)
        real(dp) :: length, b_squared, power(1)
        integer :: status, nh, nl, dimension, term
        integer :: modes(1), toroidal(1), kinds(1)

        call evaluate_radial_feec_complex(complex, s, h1, dh1, l2, status)
        if (status /= 0) error stop 'resonant basis evaluation failed'
        call build_active_indices(h1, hi, status, dh1)
        if (status /= 0) error stop 'resonant normal indices failed'
        call build_active_indices(l2, li, status)
        if (status /= 0) error stop 'resonant tangential indices failed'
        nh = size(hi)
        nl = size(li)
        dimension = nh + 2 * nl
        allocate (h(nh, 1), dh(nh, 1), eta(nl, 1), nu(nl, 1), map(dimension))
        power = 1.0_dp - 0.5_dp * real(pinch_mode, dp)
        call apply_stored_power(s, power, h1, dh1, hi, h, dh, status)
        if (status /= 0) error stop 'resonant normal axis factor failed'
        call apply_tangential_axis_weight(s, power, l2, li, eta, status)
        if (status /= 0) error stop 'resonant eta axis factor failed'
        call replicate_indexed_values(l2, li, nu, status)
        if (status /= 0) error stop 'resonant nu basis failed'
        map(:nh) = hi
        map(nh + 1:nh + nl) = complex%h1_dofs + li
        map(nh + nl + 1:) = complex%h1_dofs + complex%l2_dofs + li
        modes = pinch_mode
        toroidal = pinch_mode
        kinds = parity
        length = 2.0_dp * pi * pinch_bz / pinch_blinear
        b_squared = pinch_bz**2 + pinch_blinear**2 * pinch_radius**2 * s
        fields = 0.0_dp
        fields(:, :, 1) = pi * pinch_radius**2 * pinch_bz
        ! L=2*pi*Z/A makes P=T exactly; n=m is global resonance.
        fields(:, :, 2) = fields(:, :, 1)
        fields(:, :, 5) = length * pinch_bz
        fields(:, :, 6) = 2.0_dp * pi * pinch_blinear * pinch_radius**2 * s
        fields(:, :, 7) = pi * pinch_radius**2 * length
        fields(:, :, 8) = sqrt(b_squared)
        fields(:, :, 9) = 4.0_dp * s / pinch_radius**2
        fields(:, :, 10) = 2.0_dp * pinch_blinear * pinch_bz
        fields(:, :, 11) = -pinch_blinear**2 * pinch_radius**2
        drive = pinch_blinear**2 * pinch_radius**2 / s
        zeros = 0.0_dp
        gamma_p = gamma * (pinch_pressure_edge + pinch_blinear**2 &
            * pinch_radius**2 * (1.0_dp - s) / vacuum_permeability)
        allocate (local_k(dimension, dimension), source=0.0_dp)
        allocate (local_m(dimension, dimension), source=0.0_dp)
        allocate (terms(dimension, dimension, 5), source=0.0_dp)
        call assemble_compatible_compressible_stiffness_surface(fields, drive, &
            zeros, zeros, zeros, gamma_p, modes, toroidal, kinds, 1, &
            h, dh, eta, nu, weight, phase_assembly_transformed, local_k, &
            status, terms)
        if (status /= 0) error stop 'resonant production stiffness failed'
        local_k = 0.0_dp
        do term = 1, 5
            if (pass == 0 .or. &
                ((term == 3 .or. term == 5) .eqv. (pass == 2))) &
                local_k = local_k + terms(:, :, term)
        end do
        call scatter_matrix(map, local_k, 1.0_dp, stiffness)
        if (pass /= 2) then
            call assemble_compatible_physical_mass_surface(fields, pinch_density, &
                modes, toroidal, kinds, 1, h, eta, nu, weight, &
                phase_assembly_transformed, local_m, status)
            if (status /= 0) error stop 'resonant production mass failed'
            call scatter_matrix(map, local_m, 1.0_dp, mass)
        end if
    end subroutine accumulate

end module resonant_screw_pinch_fixture
