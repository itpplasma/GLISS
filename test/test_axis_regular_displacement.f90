! Exact-representation test for the axis-regular compatible FEEC spaces.
!
! In a straight cylinder (s=r^2/a^2, theta/2pi, z/L) the regular m=1
! displacement xi=(1-s) xi0 (cos kz, sin kz, 0) has
!   xi^s = (2/a) s^(1/2) (1-s) cos(psi),  eta = -(a Bz/2) s^(-1/2) (1-s) sin(psi)
! with psi=2 pi (theta - z/L).  Its weighted H1 and L2 coefficients are the
! polynomials (2/a) s (1-s) and -(a Bz/2)(1-s).  An unweighted L2 basis cannot
! represent the s^(-1/2) tangential behaviour, which made the |m|=1 discrete
! space non-conforming.  Independently of the discretization, the exact
! Rayleigh quotient in a uniform axial field is
!   [k^2 Bz^2/mu0 V/3 + (Bz^2/mu0 + gamma p) pi L] / (rho V/3),  V=pi a^2 L,
! and the kinetic energy rho V/3 does not depend on an added poloidal field.
! The latter exercises the regular third unknown nu=mu-(FP'/FT') sqrtg eta.
program test_axis_regular_displacement
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use radial_feec_complex, only: radial_feec_complex_t, &
        build_radial_feec_complex, evaluate_radial_feec_complex
    use compatible_radial_quadrature, only: accurate_nodes, accurate_weights, &
        build_constraint_quadrature
    use compatible_problem_assembly_support, only: build_active_indices, &
        apply_stored_power, apply_tangential_axis_weight, &
        replicate_indexed_values, scatter_matrix
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), radius = 0.5_dp
    real(dp), parameter :: length = 6.0_dp * pi, axial_field = 1.0_dp
    real(dp), parameter :: density = 2.0_dp, gamma_pressure = 500.0_dp / 3.0_dp
    real(dp), parameter :: power(1) = [0.5_dp]
    integer :: degree, surfaces

    do degree = 2, 4
        do surfaces = 3, 6, 3
            call check_case(surfaces, degree)
        end do
    end do
    write (*, '(a)') 'axis-regular displacement: all checks passed'

contains

    subroutine check_case(surfaces, degree)
        integer, intent(in) :: surfaces, degree
        type(radial_feec_complex_t) :: complex
        real(dp), allocatable :: breaks(:), stiffness(:, :), mass(:, :)
        real(dp), allocatable :: helical_mass(:, :), vector(:)
        real(dp) :: volume, wavenumber, expected, quotient, kinetic
        real(dp) :: helical_kinetic, exact_kinetic
        real(dp), allocatable :: image(:)
        integer :: i, info

        allocate (breaks(surfaces + 1))
        do i = 1, surfaces + 1
            breaks(i) = real(i - 1, dp) / real(surfaces, dp)
        end do
        call build_radial_feec_complex(breaks, degree, .true., .true., &
            complex, info)
        if (info /= 0) error stop 'radial complex failed'
        call assemble(complex, surfaces, degree, 0.0_dp, stiffness, mass)
        call assemble(complex, surfaces, degree, 0.3_dp, stiffness, &
            helical_mass, mass_only=.true.)
        call regular_displacement(complex, vector)

        volume = pi * radius**2 * length
        wavenumber = 2.0_dp * pi / length
        exact_kinetic = density * volume / 3.0_dp
        expected = (wavenumber**2 * axial_field**2 / vacuum_permeability &
            * volume / 3.0_dp + (axial_field**2 / vacuum_permeability &
            + gamma_pressure) * pi * length) / exact_kinetic
        image = matmul(mass, vector)
        kinetic = dot_product(vector, image)
        image = matmul(helical_mass, vector)
        helical_kinetic = dot_product(vector, image)
        image = matmul(stiffness, vector)
        quotient = dot_product(vector, image) / kinetic
        write (*, '(a,i0,a,i0,3(a,es11.3))') 'degree ', degree, &
            ' surfaces ', surfaces, &
            ' quotient error ', abs(quotient / expected - 1.0_dp), &
            ' kinetic error ', abs(kinetic / exact_kinetic - 1.0_dp), &
            ' helical kinetic error ', &
            abs(helical_kinetic / exact_kinetic - 1.0_dp)
        if (abs(kinetic / exact_kinetic - 1.0_dp) > 1.0e-12_dp) &
            error stop 'kinetic energy of the regular displacement is wrong'
        if (abs(quotient / expected - 1.0_dp) > 1.0e-11_dp) &
            error stop 'Rayleigh quotient of the regular displacement is wrong'
        if (abs(helical_kinetic / exact_kinetic - 1.0_dp) > 1.0e-12_dp) &
            error stop 'parallel unknown is not regular with poloidal field'
    end subroutine check_case

    ! Least-squares coefficients of the exact displacement in the weighted
    ! spaces; a nonzero residual means the space cannot represent it.
    subroutine regular_displacement(complex, vector)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), allocatable, intent(out) :: vector(:)
        integer, parameter :: samples = 97
        real(dp) :: normal(samples, complex%h1_dofs), normal_rhs(samples)
        real(dp) :: tangential(samples, complex%l2_dofs), tangential_rhs(samples)
        real(dp), allocatable :: h1(:), dh1(:), l2(:)
        real(dp) :: s, scale
        integer :: sample, info

        do sample = 1, samples
            s = (real(sample, dp) - 0.5_dp) / real(samples, dp)
            call evaluate_radial_feec_complex(complex, s, h1, dh1, l2, info)
            if (info /= 0) error stop 'radial evaluation failed'
            scale = s**(-power(1))
            normal(sample, :) = scale * h1
            tangential(sample, :) = scale * l2
            normal_rhs(sample) = 2.0_dp / radius * sqrt(s) * (1.0_dp - s)
            tangential_rhs(sample) = -0.5_dp * radius * axial_field &
                * (1.0_dp - s) / sqrt(s)
        end do
        allocate (vector(complex%h1_dofs + 2 * complex%l2_dofs), source=0.0_dp)
        call least_squares(normal, normal_rhs, vector(:complex%h1_dofs))
        call least_squares(tangential, tangential_rhs, &
            vector(complex%h1_dofs + 1:complex%h1_dofs + complex%l2_dofs))
    end subroutine regular_displacement

    subroutine least_squares(matrix, rhs, solution)
        real(dp), intent(in) :: matrix(:, :), rhs(:)
        real(dp), intent(out) :: solution(:)
        real(dp) :: work_matrix(size(matrix, 1), size(matrix, 2))
        real(dp) :: work_rhs(size(rhs), 1), work(64 * size(rhs))
        real(dp) :: residual(size(rhs))
        integer :: info, row, column

        interface
            subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, &
                    info)
                import dp
                character, intent(in) :: trans
                integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
                real(dp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
                integer, intent(out) :: info
            end subroutine dgels
        end interface

        work_matrix = matrix
        work_rhs(:, 1) = rhs
        call dgels('N', size(matrix, 1), size(matrix, 2), 1, work_matrix, &
            size(matrix, 1), work_rhs, size(rhs), work, size(work), info)
        if (info /= 0) error stop 'least-squares fit failed'
        solution = work_rhs(:size(matrix, 2), 1)
        residual = 0.0_dp
        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                residual(row) = residual(row) + matrix(row, column) * solution(column)
            end do
        end do
        do row = 1, size(rhs)
            residual(row) = residual(row) - rhs(row)
        end do
        if (maxval(abs(residual)) > 1.0e-11_dp * maxval(abs(rhs))) &
            error stop 'regular displacement is not in the discrete space'
    end subroutine least_squares

    subroutine assemble(complex, surfaces, degree, poloidal_field, &
            stiffness, mass, mass_only)
        type(radial_feec_complex_t), intent(in) :: complex
        integer, intent(in) :: surfaces, degree
        real(dp), intent(in) :: poloidal_field
        real(dp), allocatable, intent(inout) :: stiffness(:, :)
        real(dp), allocatable, intent(out) :: mass(:, :)
        logical, optional, intent(in) :: mass_only
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: coordinate, weight
        integer :: cell, point, pass, dimension, info
        logical :: with_stiffness

        with_stiffness = .true.
        if (present(mass_only)) with_stiffness = .not. mass_only
        dimension = complex%h1_dofs + 2 * complex%l2_dofs
        if (with_stiffness) then
            if (allocated(stiffness)) deallocate (stiffness)
            allocate (stiffness(dimension, dimension), source=0.0_dp)
        end if
        allocate (mass(dimension, dimension), source=0.0_dp)
        do pass = 1, 2
            if (pass == 1) then
                nodes = accurate_nodes
                weights = accurate_weights
            else
                call build_constraint_quadrature(degree, nodes, weights, info)
                if (info /= 0) error stop 'constraint quadrature failed'
            end if
            do cell = 1, surfaces
                do point = 1, size(nodes)
                    coordinate = (real(cell, dp) - 0.5_dp &
                        + 0.5_dp * nodes(point)) / real(surfaces, dp)
                    weight = 0.5_dp * weights(point) / real(surfaces, dp)
                    call accumulate(complex, coordinate, weight, pass, &
                        poloidal_field, with_stiffness, stiffness, mass)
                end do
            end do
        end do
    end subroutine assemble

    subroutine accumulate(complex, s, weight, pass, poloidal_field, &
            with_stiffness, stiffness, mass)
        type(radial_feec_complex_t), intent(in) :: complex
        real(dp), intent(in) :: s, weight, poloidal_field
        integer, intent(in) :: pass
        logical, intent(in) :: with_stiffness
        real(dp), intent(inout) :: stiffness(:, :), mass(:, :)
        real(dp) :: fields(8, 4, 13), zeros(8, 4), gamma_p(8, 4)
        real(dp), allocatable :: h1(:), dh1(:), l2(:), h(:, :), dh(:, :)
        real(dp), allocatable :: e(:, :), l(:, :), local_k(:, :), local_m(:, :)
        real(dp), allocatable :: terms(:, :, :)
        integer, allocatable :: hi(:), li(:), map(:)
        integer :: info, nh, nl, dim, term

        call evaluate_radial_feec_complex(complex, s, h1, dh1, l2, info)
        if (info /= 0) error stop 'radial evaluation failed'
        call build_active_indices(h1, hi, info, dh1)
        if (info /= 0) error stop 'H1 indices failed'
        call build_active_indices(l2, li, info)
        if (info /= 0) error stop 'L2 indices failed'
        nh = size(hi)
        nl = size(li)
        dim = nh + 2 * nl
        allocate (h(nh, 1), dh(nh, 1), e(nl, 1), l(nl, 1), map(dim))
        call apply_stored_power(s, power, h1, dh1, hi, h, dh, info)
        if (info /= 0) error stop 'normal axis factor failed'
        call apply_tangential_axis_weight(s, power, l2, li, e, info)
        if (info /= 0) error stop 'tangential axis factor failed'
        call replicate_indexed_values(l2, li, l, info)
        if (info /= 0) error stop 'L2 values failed'
        map(:nh) = hi
        map(nh + 1:nh + nl) = complex%h1_dofs + li
        map(nh + nl + 1:) = complex%h1_dofs + complex%l2_dofs + li
        ! Right-handed straight cylinder with Bz and B_theta=poloidal_field*r/a.
        fields = 0.0_dp
        fields(:, :, 1) = pi * radius**2 * axial_field
        fields(:, :, 2) = 0.5_dp * radius * length * poloidal_field
        fields(:, :, 5) = length * axial_field
        fields(:, :, 6) = 2.0_dp * pi * poloidal_field * radius * s
        fields(:, :, 7) = pi * radius**2 * length
        fields(:, :, 8) = sqrt(axial_field**2 + poloidal_field**2 * s)
        fields(:, :, 9) = 4.0_dp * s / radius**2
        zeros = 0.0_dp
        gamma_p = gamma_pressure
        allocate (local_k(dim, dim), local_m(dim, dim), terms(dim, dim, 5))
        local_k = 0.0_dp
        local_m = 0.0_dp
        terms = 0.0_dp
        if (with_stiffness) then
            call assemble_compatible_compressible_stiffness_surface(fields, &
                zeros, zeros, zeros, zeros, gamma_p, [1], [1], [1], 1, h, dh, &
                e, l, weight, phase_assembly_transformed, local_k, info, terms)
            if (info /= 0) error stop 'surface stiffness failed'
            local_k = 0.0_dp
            do term = 1, 5
                if ((term == 3 .or. term == 5) .eqv. (pass == 2)) &
                    local_k = local_k + terms(:, :, term)
            end do
            call scatter_matrix(map, local_k, 1.0_dp, stiffness)
        end if
        if (pass == 1) then
            call assemble_compatible_physical_mass_surface(fields, density, &
                [1], [1], [1], 1, h, e, l, weight, &
                phase_assembly_transformed, local_m, info)
            if (info /= 0) error stop 'surface mass failed'
            call scatter_matrix(map, local_m, 1.0_dp, mass)
        end if
    end subroutine accumulate

end program test_axis_regular_displacement
