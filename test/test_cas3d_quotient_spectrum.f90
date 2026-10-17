! Independent coefficient-space oracle: product-to-sum matrices below are
! specified directly, never obtained from the production envelope transform.
! Dense matrices are confined to these small tests. Orthonormalizing C^T
! removes the exact label null space before the dense reference eigensolve.
program test_cas3d_quotient_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use cas3d_phase_envelope_transform, only: &
        apply_cas3d_phase_envelope_block_congruence, &
        build_cas3d_phase_envelope_map, cas3d_phase_envelope_map_t, &
        cas3d_phase_transform_invalid, cas3d_phase_transform_ok, cas3d_quotient_t
    use compatible_block_storage, only: build_compatible_block_indices, &
                                        compatible_block_ok
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_problem_ok, &
        compatible_two_component_problem_t
    use cylinder_fixture, only: create_cylinder_fixture
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t
    use marginality_spectrum, only: compute_phase_envelope_spectrum, &
        marginality_normalization_cas3d2mn, marginality_spectrum_compute_error, &
        marginality_spectrum_invalid, marginality_spectrum_ok, &
                     marginality_spectrum_result_t, solve_compatible_marginality_problem
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated, &
                                     symmetric_eigensolver_ok
    use variable_block_tridiagonal, only: &
        pack_permuted_variable_blocks, pack_variable_blocks, variable_block_ok, &
        variable_block_tridiagonal_t
    implicit none

    character(len=*), parameter :: filename = 'cas3d_quotient_cylinder.nc'
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    integer :: case_id, parity, status

    call create_cylinder_fixture(filename, chart_shift=0.0_dp, surfaces=9)
    call read_gvec_cas3d_file(filename, equilibrium, status)
    call require(status == reader_ok, 'cylinder read failed')
    do case_id = 1, 3
        do parity = 1, 2
            call check_equilibrium(case_id, parity)
            call check_near_zero(case_id, parity)
        end do
    end do
    call check_equilibrium(4, 1)
    call check_zero_stiffness()
    call check_physical_normal_only()
    call check_empty_normal_component()
    call check_position_bandwidth_admission()
    open (newunit=status, file=filename, status='old')
    close (status, status='delete')
    write (*, '(a)') 'PASS: coefficient quotient eigenpairs and floor scaling'

contains

    subroutine tables(which, parity, pm, pn, lm, ln, orientation, em, en, cn, ce)
        integer, intent(in) :: which, parity
        integer, allocatable, intent(out) :: pm(:), pn(:), lm(:), ln(:)
        integer, allocatable, intent(out) :: orientation(:), em(:), en(:)
        real(dp), allocatable, intent(out) :: cn(:, :), ce(:, :)

        select case (which)
        case (1) ! Three distinct sidebands of a nonzero carrier.
            pm = [3, 4, 2]; pn = [2, 2, 2]
            lm = pm; ln = pn; orientation = [1, 1, 1]
            em = [0, 1]; en = [0, 0]
            allocate (cn(3, 3), source=0.0_dp)
            cn(1, 1) = 1.0_dp
            cn(2, 2:3) = [0.5_dp, -0.5_dp]
            cn(3, 2:3) = [0.5_dp, 0.5_dp]
            ce = cn
            ce(:, 3) = -cn(:, 3)
        case (2) ! Opposite envelope labels reach the same sidebands.
            pm = [1, 1, 1]; pn = [1, 0, 2]
            lm = [1, 1, 1, 1, 1]; ln = [1, 0, 2, 2, 0]
            orientation = [1, 1, 1, 1, 1]
            em = [0, 0, 0]; en = [0, -1, 1]
            allocate (cn(3, 5), source=0.0_dp)
            cn(1, 1) = 1.0_dp
            cn(2, 2:5) = [0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp]
            cn(3, 2:5) = [0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]
            ce = cn
            ce(:, 3) = -cn(:, 3)
            ce(:, 5) = -cn(:, 5)
        case (3) ! The +/- sidebands coincide after canonicalization.
            pm = [0, 1]; pn = [0, 0]
            lm = [0, 1, 1]; ln = [0, 0, 0]
            orientation = [1, 1, -1]
            em = [0, 1]; en = [0, 0]
            allocate (cn(2, 3), ce(2, 3), source=0.0_dp)
            cn(1, 1) = 1.0_dp; ce(1, 1) = 1.0_dp
            if (parity == 1) then
                cn(2, 2) = 1.0_dp; ce(2, 3) = 1.0_dp
            else
                cn(2, 3) = -1.0_dp; ce(2, 2) = 1.0_dp
            end if
        case (4) ! The zero carrier alone has one inactive component.
            pm = [0]; pn = [0]; lm = [0]; ln = [0]
            orientation = [1]; em = [0]; en = [0]
            allocate (cn(1, 1), ce(1, 1), source=1.0_dp)
        end select
    end subroutine tables

    subroutine full_map(cn, ce, h1, l2, c)
        real(dp), intent(in) :: cn(:, :), ce(:, :)
        integer, intent(in) :: h1, l2
        real(dp), allocatable, intent(out) :: c(:, :)
        integer :: b, p, e

        p = size(cn, 1); e = size(cn, 2)
        allocate (c((h1 + l2)*p, (h1 + l2)*e), source=0.0_dp)
        do b = 0, h1 + l2 - 1
            if (b < h1) then
                c(b*p + 1:(b + 1)*p, b*e + 1:(b + 1)*e) = cn
            else
                c(b*p + 1:(b + 1)*p, b*e + 1:(b + 1)*e) = ce
            end if
        end do
    end subroutine full_map

    subroutine remove_zero_phase_rows(pm, pn, parity, h1, l2, c)
        integer, intent(in) :: pm(:), pn(:), parity, h1, l2
        real(dp), allocatable, intent(inout) :: c(:, :)
        real(dp), allocatable :: active_map(:, :)
        logical, allocatable :: active(:)
        integer :: basis, mode, column

        allocate (active(size(c, 1)), source=.true.)
        do basis = 0, h1 + l2 - 1
            do mode = 1, size(pm)
                if (pm(mode) /= 0 .or. pn(mode) /= 0) cycle
                if ((basis < h1 .and. parity == 2) &
                    .or. (basis >= h1 .and. parity == 1)) &
                    active(basis*size(pm) + mode) = .false.
            end do
        end do
        allocate (active_map(count(active), size(c, 2)))
        do column = 1, size(c, 2)
            active_map(:, column) = pack(c(:, column), active)
        end do
        call move_alloc(active_map, c)
    end subroutine remove_zero_phase_rows

    subroutine dense_oracle(k, c, scale, values, floor)
        real(dp), intent(in) :: k(:, :), c(:, :), scale
        real(dp), allocatable, intent(out) :: values(:)
        real(dp), intent(out) :: floor
        real(dp), allocatable :: labeled(:, :), q(:, :), v(:), cq(:, :)
        real(dp), allocatable :: image(:, :)
        real(dp), allocatable :: reduced(:, :), mass(:, :), vectors(:, :)
        real(dp) :: norm
        integer :: i, j, rank, local_info

        call require(size(k, 1) == size(c, 1) &
                     .and. size(k, 2) == size(c, 1), 'oracle map dimensions differ')
        call require(all(ieee_is_finite(k)), 'physical oracle matrix nonfinite')
        image = matmul(k, c)
        labeled = matmul(transpose(c), image)
        floor = 1024.0_dp*epsilon(1.0_dp) &
                *maxval(sum(abs(labeled), dim=2))/scale
        allocate (q(size(c, 2), size(c, 1)), source=0.0_dp)
        rank = 0
        do i = 1, size(c, 1)
            v = c(i, :)
            do j = 1, rank
                v = v - dot_product(q(:, j), v)*q(:, j)
            end do
            norm = sqrt(dot_product(v, v))
            if (norm < 1.0e-12_dp) cycle
            rank = rank + 1
            q(:, rank) = v/norm
        end do
        cq = matmul(c, q(:, :rank))
        image = matmul(k, cq)
        reduced = matmul(transpose(cq), image)
        ! BLAS contractions can choose different accumulation orders for
        ! opposite triangles. Preserve the exact symmetric reference form.
        do j = 1, rank
            do i = j + 1, rank
                norm = 0.5_dp*(reduced(i, j) + reduced(j, i))
                reduced(i, j) = norm
                reduced(j, i) = norm
            end do
        end do
        allocate (mass(rank, rank), source=0.0_dp)
        do i = 1, rank
            mass(i, i) = scale
        end do
        call solve_symmetric_generalized_allocated(reduced, mass, values, &
                                                   vectors, local_info, .false.)
        if (local_info /= symmetric_eigensolver_ok) &
            write (error_unit, '(a,2i5)') 'dense quotient status/rank: ', &
            local_info, rank
        call require(local_info == symmetric_eigensolver_ok, &
                     'independent dense quotient solve failed')
    end subroutine dense_oracle

    subroutine check_result(result, values, floor)
        type(marginality_spectrum_result_t), intent(in) :: result
        real(dp), intent(in) :: values(:), floor
        real(dp) :: expected, tolerance
        integer :: first

        if (abs(result%zero_floor - floor) > 1.0e-12_dp*floor) &
            write (error_unit, '(a,2es24.15)') 'reported/coefficient floors: ', &
            result%zero_floor, floor
        call require(abs(result%zero_floor - floor) <= 1.0e-12_dp*floor, &
                     'reported floor is not the coefficient-pencil floor')
        call require(result%negative_count == count(values < -floor), &
                     'coefficient negative count differs from dense quotient')
        first = 1
        do while (first <= size(values))
            if (abs(values(first)) > floor) exit
            first = first + 1
        end do
        call require(first <= size(values), 'oracle has no resolved mode')
        expected = values(first)
        tolerance = max(4.0_dp*result%certificate, 2.0e-10_dp &
                        *abs(expected), 4.0_dp*floor)
        call require(abs(result%lowest_eigenvalue - expected) <= tolerance, &
                     'lowest coefficient eigenpair differs from dense quotient')
    end subroutine check_result

    subroutine check_equilibrium(which, parity)
        integer, intent(in) :: which, parity
        type(compatible_two_component_problem_t) :: physical
        type(marginality_spectrum_result_t) :: result, inertia
        integer, allocatable :: pm(:), pn(:), lm(:), ln(:), orientation(:)
        integer, allocatable :: em(:), en(:)
        real(dp), allocatable :: cn(:, :), ce(:, :), c(:, :), values(:)
        real(dp), allocatable :: powers(:)
        real(dp) :: scale, floor, length, baseline
        integer :: local_info, variant, nt, nz
        character(len=160) :: message

        call tables(which, parity, pm, pn, lm, ln, orientation, em, en, cn, ce)
        allocate (powers(size(pm)), source=0.0_dp)
        call build_compatible_two_component_problem(equilibrium, pm, pn, &
                                      powers, parity, 1, 13, 13, physical, local_info, &
                                                radial_cells=3, axis_conforming=.false.)
        call require(local_info == compatible_problem_ok, &
                     'small physical assembly failed')
        call full_map(cn, ce, physical%h1_dofs, physical%l2_dofs, c)
        call remove_zero_phase_rows(pm, pn, parity, physical%h1_dofs, &
                                    physical%l2_dofs, c)
        do variant = 1, 4
            nt = 7; nz = 11; length = 1.0_dp
            if (variant == 2) length = 2.0_dp
            if (variant >= 3) nt = 14
            if (variant == 3) nz = 22
            ! Independent transcription of the documented coefficient norm.
            scale = 0.5_dp*length**3/real(3*nt*nz, dp)
            call dense_oracle(physical%stiffness, c, scale, values, floor)
            call compute_phase_envelope_spectrum(equilibrium, pm(1), pn(1), &
                                em, en, parity, 1, 13, 13, .true., result, local_info, &
                     message, normalization_policy=marginality_normalization_cas3d2mn, &
                                        coefficient_n_theta=nt, coefficient_n_zeta=nz, &
                                                reference_length=length, radial_cells=3)
            if (local_info /= marginality_spectrum_ok) &
                write (error_unit, '(3i3,2a)') which, parity, variant, ': ', &
                trim(message)
            call require(local_info == marginality_spectrum_ok, &
                         'production coefficient solve failed')
            call check_result(result, values, floor)
            call require(result%quotient_rank == size(c, 1), &
                         'quotient rank did not exclude inactive sin(0) rows')
            call require(result%labeled_nullity == size(c, 2) - size(c, 1), &
                         'label nullity differs from independent active map')
            if (variant == 1) baseline = result%lowest_eigenvalue*scale
            call require(abs(result%lowest_eigenvalue*scale - baseline) &
                         <= 1.0e-9_dp*abs(baseline), &
                         'reference-length or coefficient-grid scaling differs')
            call compute_phase_envelope_spectrum(equilibrium, pm(1), pn(1), &
                              em, en, parity, 1, 13, 13, .false., inertia, local_info, &
                     message, normalization_policy=marginality_normalization_cas3d2mn, &
                                        coefficient_n_theta=nt, coefficient_n_zeta=nz, &
                                                reference_length=length, radial_cells=3)
            call require(local_info == marginality_spectrum_ok, &
                         'coefficient inertia-only solve failed')
            call require(inertia%negative_count == result%negative_count &
                         .and. inertia%zero_floor == result%zero_floor, &
                         'inertia and eigenpair paths use different floor semantics')
        end do
    end subroutine check_equilibrium

    subroutine check_near_zero(which, parity)
        integer, intent(in) :: which, parity
        type(compatible_two_component_problem_t) :: problem
        type(cas3d_phase_envelope_map_t) :: map
        type(cas3d_quotient_t) :: quotient
        type(marginality_spectrum_result_t) :: result
        integer, allocatable :: pm(:), pn(:), lm(:), ln(:), orientation(:)
        integer, allocatable :: em(:), en(:), widths(:), permutation(:)
        real(dp), allocatable :: cn(:, :), ce(:, :), c(:, :), k(:, :)
        real(dp), allocatable :: mass(:, :), values(:), vector(:)
        real(dp), allocatable :: component_vector(:), labeled(:, :), image(:)
        real(dp), allocatable :: product(:, :)
        real(dp) :: floor, scale, residual
        integer :: i, local_info, variant, n
        character(len=160) :: message

        call tables(which, parity, pm, pn, lm, ln, orientation, em, en, cn, ce)
        call full_map(cn, ce, 1, 2, c)
        n = 3*size(pm)
        allocate (k(n, n), mass(n, n), source=0.0_dp)
        do i = 1, n
            k(i, i) = real(i, dp)
            mass(i, i) = 1.0_dp
        end do
        k(1, 1) = -1.0e-14_dp ! Inside the marginal floor band.
        k(2, 2) = -1.0e-8_dp ! Resolved physical instability near zero.
        k(3, 3) = 1.0e-8_dp
        call build_cas3d_phase_envelope_map(lm, ln, orientation, pm, pn, &
                                            parity, map, local_info)
        call require(local_info == cas3d_phase_transform_ok, 'map failed')
        do variant = 1, 4
            scale = 10.0_dp**(2*modulo(variant - 1, 2) - 1)
            k(2, 2) = -1.0e-8_dp
            if (variant > 2) k(2, 2) = 1.0e-8_dp
            call dense_oracle(k, c, scale, values, floor)
            ! All unknowns in one block for packing, then the expected
            ! degree-one physical block layout [normal+eta, eta].
            call build_compatible_block_indices(1, 2, size(pm), size(pm), &
                                                1, widths, problem%sparse_block_index, &
                                                problem%sparse_local_index, local_info)
            call require(local_info == compatible_block_ok, 'layout failed')
            if (allocated(permutation)) deallocate (permutation)
            allocate (permutation(n))
            do i = 1, n
                permutation(sum(widths(:problem%sparse_block_index(i) - 1)) &
                            + problem%sparse_local_index(i)) = i
            end do
            call pack_permuted_variable_blocks(k, permutation, widths, &
                                               problem%sparse_stiffness, local_info)
            call require(local_info == variable_block_ok, 'K packing failed')
            call pack_permuted_variable_blocks(mass, permutation, widths, &
                                               problem%sparse_mass, local_info)
            call require(local_info == variable_block_ok, 'M packing failed')
            call apply_cas3d_phase_envelope_block_congruence(map, 1, 2, &
                                 scale, problem%sparse_stiffness, problem%sparse_mass, &
                                                             quotient, local_info)
            call require(local_info == cas3d_phase_transform_ok, &
                         'synthetic congruence failed')
            call build_compatible_block_indices(1, 2, size(lm), size(lm), &
                                                1, widths, problem%sparse_block_index, &
                                                problem%sparse_local_index, local_info)
            problem%normal_unknowns = size(lm)
            problem%eta_unknowns = 2*size(lm)
            problem%has_sparse_storage = .true.
            call solve_compatible_marginality_problem(problem, .true., &
                                        result, local_info, message, eigenvector=vector)
            call require(local_info == marginality_spectrum_ok, &
                         'synthetic near-zero solve failed')
            call check_result(result, values, floor)
            call require(ieee_is_finite(result%certificate), &
                         'near-zero eigenpair certificate is nonfinite')
            call require(result%certificate < 1.0e-10_dp/scale, &
                         'near-zero eigenpair certificate is too wide')
            component_vector = vector
            do i = 1, size(vector)
                component_vector(i) = vector(sum(widths( &
                                                 :problem%sparse_block_index(i) - 1)) &
                                             + problem%sparse_local_index(i))
            end do
            call require(abs(scale*dot_product(vector, vector) - 1.0_dp) &
                         < 1.0e-12_dp, 'coefficient vector normalization differs')
            product = matmul(k, c)
            labeled = matmul(transpose(c), product)
            image = matmul(labeled, component_vector)
            image = image - result%lowest_eigenvalue * scale * component_vector
            residual = sqrt(dot_product(image, image)) &
                       /(scale*sqrt(dot_product(component_vector, component_vector)))
            if (residual > max(2.0_dp*result%eigenpair_residual, &
                               4.0_dp*floor)) &
                write (error_unit, '(a,3i3,4es16.7)') 'residual case: ', &
                which, parity, variant, residual, result%eigenpair_residual, &
                floor, result%lowest_eigenvalue
            call require(residual <= max(2.0_dp*result%eigenpair_residual, &
                              4.0_dp*floor), 'independent coefficient residual differs')
        end do
    end subroutine check_near_zero

    subroutine check_zero_stiffness()
        type(compatible_two_component_problem_t) :: physical
        type(marginality_spectrum_result_t) :: result
        character(len=160) :: message
        integer :: local_info, normalization

        call build_compatible_two_component_problem(equilibrium, [0], [0], &
            [0.0_dp], 2, 1, 13, 13, physical, local_info, radial_cells=3)
        call require(local_info == compatible_problem_ok, &
            'zero-carrier eta-only assembly failed')
        call require(all(physical%stiffness == 0.0_dp), &
            'axisymmetric eta(0,0) symmetry mode is not exactly marginal')
        do normalization = 1, 2
            call compute_phase_envelope_spectrum(equilibrium, 0, 0, [0], [0], &
                2, 1, 13, 13, .false., result, local_info, message, &
                normalization_policy=normalization, radial_cells=3)
            call require(local_info == marginality_spectrum_ok, &
                'exact-zero pencil inertia failed')
            call require(result%negative_count == 0, &
                'exact-zero pencil has negative inertia')
            if (normalization == marginality_normalization_cas3d2mn) &
                call require(result%quotient_rank == 3 &
                    .and. result%labeled_nullity == 2, &
                    'eta-only quotient dimensions differ')
            call compute_phase_envelope_spectrum(equilibrium, 0, 0, [0], [0], &
                2, 1, 13, 13, .true., result, local_info, message, &
                normalization_policy=normalization, radial_cells=3)
            call require(local_info == marginality_spectrum_compute_error &
                .and. index(message, 'no resolved eigenvalue') > 0, &
                'exact-zero pencil returned a spurious resolved eigenpair')
        end do
        call build_compatible_two_component_problem(equilibrium, [0], [0], &
            [0.0_dp], 2, 1, 13, 13, physical, local_info, radial_cells=3, &
            sparse_storage=.true.)
        call require(local_info == compatible_problem_ok, &
            'zero-pencil mass-control assembly failed')
        physical%sparse_mass%diagonal(1)%values(1, 1) = -1.0_dp
        call solve_compatible_marginality_problem(physical, .false., result, &
            local_info, message)
        call require(local_info == marginality_spectrum_compute_error, &
            'zero stiffness accepted indefinite mass')
        physical%sparse_mass%diagonal(1)%values(1, 1) = 0.0_dp
        call solve_compatible_marginality_problem(physical, .false., result, &
            local_info, message)
        call require(local_info == marginality_spectrum_compute_error, &
            'zero stiffness accepted singular mass')
    end subroutine check_zero_stiffness

    subroutine check_physical_normal_only()
        type(compatible_two_component_problem_t) :: physical
        type(marginality_spectrum_result_t) :: result
        real(dp), allocatable :: values(:), vectors(:, :)
        character(len=160) :: message
        integer :: local_info

        call build_compatible_two_component_problem(equilibrium, [0], [0], &
            [0.0_dp], 1, 1, 13, 13, physical, local_info, radial_cells=3, &
            density_kg_m3=1.0_dp)
        call require(local_info == compatible_problem_ok, &
            'physical normal-only assembly failed')
        call require(physical%normal_unknowns == 2 &
            .and. physical%eta_unknowns == 0, 'normal-only topology differs')
        call solve_symmetric_generalized_allocated(physical%stiffness, &
            physical%mass, values, vectors, local_info)
        call require(local_info == symmetric_eigensolver_ok, &
            'physical normal-only dense oracle failed')
        call compute_phase_envelope_spectrum(equilibrium, 0, 0, [0], [0], &
            1, 1, 13, 13, .true., result, local_info, message, radial_cells=3)
        call require(local_info == marginality_spectrum_ok, &
            'physical normal-only sparse solve failed')
        call require(abs(result%lowest_eigenvalue - values(1)) &
            <= max(4.0_dp * result%certificate, 1.0e-9_dp * abs(values(1))), &
            'physical normal-only eigenpair differs from dense oracle')
    end subroutine check_physical_normal_only

    subroutine check_empty_normal_component()
        type(cas3d_phase_envelope_map_t) :: map
        type(cas3d_quotient_t) :: quotient
        type(variable_block_tridiagonal_t) :: stiffness, mass
        real(dp) :: diagonal(4, 4)
        integer :: i, local_info

        ! A sine normal component may be entirely absent while the cosine
        ! eta component retains physical unknowns. No eigenpair is needed.
        call build_cas3d_phase_envelope_map([0, 1, 1], [0, 0, 0], &
            [1, 1, -1], [0, 1], [0, 0], 2, map, local_info)
        call require(local_info == cas3d_phase_transform_ok, &
            'normal-empty map failed')
        diagonal = 0.0_dp
        do i = 1, 4
            diagonal(i, i) = real(i, dp)
        end do
        call pack_variable_blocks(diagonal, [2, 2], stiffness, local_info)
        call require(local_info == variable_block_ok, 'empty-normal K failed')
        call pack_variable_blocks(diagonal, [2, 2], mass, local_info)
        call require(local_info == variable_block_ok, 'empty-normal M failed')
        call apply_cas3d_phase_envelope_block_congruence(map, 1, 2, 1.0_dp, &
            stiffness, mass, quotient, local_info, normal_active=[.false.])
        call require(local_info == cas3d_phase_transform_invalid, &
            'wrong-size active mask accepted')
        call apply_cas3d_phase_envelope_block_congruence(map, 1, 2, 1.0_dp, &
            stiffness, mass, quotient, local_info, &
            normal_active=[.false., .false.])
        call require(local_info == cas3d_phase_transform_ok, &
            'normal-empty coefficient congruence failed')
        call require(quotient%quotient_rank == 4 .and. quotient%nullity == 5, &
            'normal-empty quotient dimensions differ')
        call require(all(stiffness%diagonal(1)%values(:3, :) == 0.0_dp) &
            .and. all(stiffness%lower(1)%values(:, :3) == 0.0_dp), &
            'inactive normal labels acquired stiffness')
    end subroutine check_empty_normal_component

    subroutine check_position_bandwidth_admission()
        type(gvec_cas3d_equilibrium_t) :: padded
        type(marginality_spectrum_result_t) :: result
        character(len=160) :: message, expected
        integer :: bandwidth, minimum, normalization, symmetry, local_info

        do bandwidth = 14, 18, 4
            padded = equilibrium
            ! Zero padding leaves the cylinder unchanged. The MDE admission
            ! contract uses the declared position table, including these rows.
            deallocate (padded%poloidal_modes, padded%toroidal_modes)
            allocate (padded%poloidal_modes(3), padded%toroidal_modes(5))
            padded%poloidal_modes(:2) = equilibrium%poloidal_modes
            padded%poloidal_modes(3) = bandwidth
            padded%toroidal_modes(:3) = equilibrium%toroidal_modes
            padded%toroidal_modes(4) = bandwidth
            padded%toroidal_modes(5) = -bandwidth
            call pad_harmonics(padded%xhat)
            call pad_harmonics(padded%yhat)
            call pad_harmonics(padded%zhat)
            call pad_harmonics(padded%mod_b)
            call pad_harmonics(padded%jacobian)
            call pad_harmonics(padded%g_tt)
            call pad_harmonics(padded%g_tz)
            call pad_harmonics(padded%g_zz)
            call pad_harmonics(padded%g_st)
            call pad_harmonics(padded%g_sz)
            call pad_harmonics(padded%second_form_tt)
            call pad_harmonics(padded%second_form_tz)
            call pad_harmonics(padded%second_form_zz)
            call pad_harmonics(padded%b_contravariant_theta)
            call pad_harmonics(padded%b_contravariant_zeta)
            minimum = 6 * bandwidth + 1
            write (expected, '(a,i0,a,i0)') &
                'cubic-Jacobian MDE requires angular_theta >= ', minimum, &
                ' and angular_zeta >= ', minimum
            do symmetry = 1, 2
                do normalization = 1, 2
                    call compute_phase_envelope_spectrum(padded, 3, 2, &
                        [0], [0], symmetry, 1, 96, 64, .true., result, &
                        local_info, message, normalization_policy=normalization, &
                        radial_cells=3)
                    if (local_info /= marginality_spectrum_invalid) &
                        write (error_unit, '(a,i0,2a)') &
                        'angular admission status ', local_info, ': ', trim(message)
                    call require(local_info == marginality_spectrum_invalid, &
                        'undersampled position table was not an input error')
                    call require(message == expected, &
                        'angular admission did not report actual minimum counts')
                    call compute_phase_envelope_spectrum(padded, 3, 2, &
                        [0], [0], symmetry, 1, minimum, minimum, .true., result, &
                        local_info, message, normalization_policy=normalization, &
                        radial_cells=3)
                    if (local_info /= marginality_spectrum_ok) &
                        write (error_unit, '(a,4(i0,1x),2a)') &
                        'admitted status/band/parity/norm ', local_info, bandwidth, &
                        symmetry, normalization, ': ', trim(message)
                    call require(local_info == marginality_spectrum_ok, &
                        'admitted position table failed at minimum angular grid')
                    call require(ieee_is_finite(result%lowest_eigenvalue) &
                        .and. ieee_is_finite(result%certificate), &
                        'admitted angular grid returned a nonfinite eigenpair')
                end do
            end do
            call compute_phase_envelope_spectrum(padded, 3, 2, [0], [0], &
                1, 1, minimum - 1, minimum, .false., result, local_info, &
                message, radial_cells=3)
            call require(local_info == marginality_spectrum_invalid &
                .and. message == expected, 'poloidal Nyquist boundary accepted')
            call compute_phase_envelope_spectrum(padded, 3, 2, [0], [0], &
                1, 1, minimum, minimum - 1, .false., result, local_info, &
                message, radial_cells=3)
            call require(local_info == marginality_spectrum_invalid &
                .and. message == expected, 'toroidal Nyquist boundary accepted')
        end do
    end subroutine check_position_bandwidth_admission

    subroutine pad_harmonics(pair)
        type(harmonic_pair_t), intent(inout) :: pair
        real(dp), allocatable :: coefficients(:, :, :)
        integer :: nn, nm, ns

        ns = size(pair%cosine, 1)
        nm = size(pair%cosine, 2)
        nn = size(pair%cosine, 3)
        allocate (coefficients(ns, nm + 1, nn + 2), source=0.0_dp)
        coefficients(:, :nm, :nn) = pair%cosine
        call move_alloc(coefficients, pair%cosine)
        allocate (coefficients(ns, nm + 1, nn + 2), source=0.0_dp)
        coefficients(:, :nm, :nn) = pair%sine
        call move_alloc(coefficients, pair%sine)
    end subroutine pad_harmonics

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, '(a)') 'FAIL: '//message
        error stop 1
    end subroutine require

end program test_cas3d_quotient_spectrum
