module exact_cylinder
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use radial_feec_complex, only: radial_feec_complex_t, &
        build_radial_feec_complex, evaluate_radial_feec_complex
    use compatible_radial_quadrature, only: accurate_nodes, accurate_weights, &
        build_constraint_quadrature
    use compatible_problem_assembly_support, only: build_active_indices, &
        apply_stored_power, replicate_indexed_values, scatter_matrix
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use phase_assembly_policy, only: phase_assembly_transformed
    implicit none
    private
    public :: assemble_exact_cylinder
contains
    subroutine assemble_exact_cylinder(ns, degree, parity, stiffness, mass)
        integer, intent(in) :: ns, degree, parity
        real(dp), allocatable, intent(out) :: stiffness(:, :), mass(:, :)
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: breaks(:), nodes(:), weights(:)
        real(dp) :: coordinate, weight
        integer :: i, cell, point, pass, dimension, info

        allocate (breaks(ns + 1))
        do i = 1, ns + 1
            breaks(i) = real(i - 1, dp)/real(ns, dp)
        end do
        call build_radial_feec_complex(breaks, degree, .true., .true., &
            complex, info)
        if (info /= 0) error stop 'Radial complex failed'
        dimension = complex%h1_dofs + 2*complex%l2_dofs
        allocate (stiffness(dimension, dimension), mass(dimension, dimension))
        stiffness = 0.0_dp
        mass = 0.0_dp
        do pass = 1, 2
            if (pass == 1) then
                nodes = accurate_nodes
                weights = accurate_weights
            else
                call build_constraint_quadrature(degree, nodes, weights, info)
                if (info /= 0) error stop 'Constraint quadrature failed'
            end if
            do cell = 1, ns
                do point = 1, size(nodes)
                    coordinate = (real(cell, dp) - 0.5_dp &
                        + 0.5_dp*nodes(point))/real(ns, dp)
                    weight = 0.5_dp*weights(point)/real(ns, dp)
                    call accumulate(complex, coordinate, weight, parity, &
                        pass, stiffness, mass)
                end do
            end do
        end do
    end subroutine assemble_exact_cylinder

    subroutine accumulate(complex, s, weight, parity, pass, stiffness, mass)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: s, weight
        integer, intent(in) :: parity, pass
        real(dp), intent(inout) :: stiffness(:, :), mass(:, :)
        real(dp), parameter :: pi = acos(-1.0_dp), radius = 0.5_dp
        real(dp), parameter :: length = 6.0_dp*pi
        real(dp) :: fields(8, 4, 13), zeros(8, 4), gamma_p(8, 4)
        real(dp), allocatable :: h1(:), dh1(:), l2(:), h(:, :), dh(:, :), l(:, :)
        real(dp), allocatable :: local_k(:, :), local_m(:, :), terms(:, :, :)
        integer, allocatable :: hi(:), li(:), map(:)
        integer :: info, nh, nl, dim, term, parity_table(1)

        call evaluate_radial_feec_complex(complex, s, h1, dh1, l2, info)
        if (info /= 0) error stop 'Radial evaluation failed'
        call build_active_indices(h1, hi, info, dh1)
        if (info /= 0) error stop 'H1 indices failed'
        call build_active_indices(l2, li, info)
        if (info /= 0) error stop 'L2 indices failed'
        parity_table(1) = parity
        nh = size(hi)
        nl = size(li)
        dim = nh + 2*nl
        allocate (h(nh, 1), dh(nh, 1), l(nl, 1), map(dim))
        call apply_stored_power(s, [-0.5_dp], h1, dh1, hi, h, dh, info)
        if (info /= 0) error stop 'Axis factor failed'
        call replicate_indexed_values(l2, li, l, info)
        if (info /= 0) error stop 'L2 values failed'
        map(:nh) = hi
        map(nh + 1:nh + nl) = complex%h1_dofs + li
        map(nh + nl + 1:) = complex%h1_dofs + complex%l2_dofs + li
        ! Exact right-handed (s=r²/a², theta/2pi, z/L) cylinder geometry.
        ! B=1 T, zero current, constant pressure100 Pa, physical rho=2.
        fields = 0.0_dp
        fields(:, :, 1) = pi*radius**2
        fields(:, :, 5) = length
        fields(:, :, 7) = pi*radius**2*length
        fields(:, :, 8) = 1.0_dp
        fields(:, :, 9) = 4.0_dp*s/radius**2
        zeros = 0.0_dp
        gamma_p = (5.0_dp/3.0_dp)*100.0_dp
        allocate (local_k(dim, dim), local_m(dim, dim), terms(dim, dim, 5))
        local_k = 0.0_dp
        local_m = 0.0_dp
        terms = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, &
            zeros, zeros, zeros, zeros, gamma_p, [3], [1], parity_table, 1, &
            h, dh, l, weight, phase_assembly_transformed, local_k, info, terms)
        if (info /= 0) error stop 'Production surface stiffness failed'
        local_k = 0.0_dp
        do term = 1, 5
            if ((term == 3 .or. term == 5) .eqv. (pass == 2)) &
                local_k = local_k + terms(:, :, term)
        end do
        call scatter_matrix(map, local_k, 1.0_dp, stiffness)
        if (pass == 1) then
            call assemble_compatible_physical_mass_surface(fields, 2.0_dp, &
                [3], [1], parity_table, 1, h, l, weight, &
                phase_assembly_transformed, local_m, info)
            if (info /= 0) error stop 'Production surface mass failed'
            call scatter_matrix(map, local_m, 1.0_dp, mass)
        end if
    end subroutine accumulate
end module exact_cylinder
