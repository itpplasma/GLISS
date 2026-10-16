program test_pressure_stiffness_tangent
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_dynamic, omp_get_max_threads, &
        omp_get_num_threads, omp_set_dynamic, omp_set_num_threads
    use compatible_three_component_problem, only: &
        build_compatible_three_component_problem, &
        build_compatible_pressure_stiffness_tangent, &
        build_compatible_pressure_trace_gradient, &
        compatible_three_component_problem_t, compatible_three_component_ok, &
        compatible_three_component_invalid
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_pressure_stiffness_tangent, &
        contract_compatible_pressure_stiffness
    use fixed_boundary_energy, only: fixed_boundary_energy_terms_t
    use fixed_boundary_spectrum, only: build_fixed_boundary_problem, &
        diagnose_fixed_boundary_energy, fixed_boundary_ok, &
        fixed_boundary_pressure_trace_jvp, fixed_boundary_pressure_trace_vjp, &
        fixed_boundary_problem_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use variable_block_tridiagonal, only: variable_block_to_dense, variable_block_ok
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    real(dp), allocatable :: direction(:)
    character(len=1024) :: directory
    integer :: info, parity, mesh, storage

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) directory = 'test/data'
    call read_gvec_cas3d_file(trim(directory) // '/solovev_q1.045.nc', &
        equilibrium, info)
    call require(info == reader_ok, 'pressure fixture read failed')
    direction = 0.03_dp * equilibrium%pressure * (1.0_dp + equilibrium%s)
    call check_parallel_invariance()
    do parity = 0, 2
        do mesh = 1, 2
            do storage = 0, 1
                call check_direction(parity, 6 * mesh, storage == 1)
            end do
        end do
    end do
    call check_duality()
    call check_nonuniform_mesh()
    call check_period_contraction()
    call check_invalid()
    print *, 'PASS'

contains

    subroutine check_parallel_invariance()
        real(dp), allocatable :: serial_k(:, :), serial_m(:, :), serial_dk(:, :)
        real(dp), allocatable :: parallel_k(:, :), parallel_m(:, :)
        real(dp), allocatable :: parallel_dk(:, :), tangent_m(:, :)
        real(dp), allocatable :: vectors(:, :), serial_gradient(:), gradient(:)
        real(dp) :: serial_jvp, jvp
        logical :: previous_dynamic
        integer :: previous_threads, actual_threads, row, repeat, status

        ! Concurrent geometry and matrix BLAS calls must preserve the same
        ! primal operator and pressure actions as one assembly worker.
        previous_threads = omp_get_max_threads()
        previous_dynamic = omp_get_dynamic()
        call omp_set_dynamic(.false.)
        call omp_set_num_threads(1)
        call assemble(equilibrium, 0, 6, .false., serial_k, serial_m)
        call assemble(equilibrium, 0, 6, .false., serial_dk, tangent_m, direction)
        allocate (vectors(size(serial_k, 1), 2))
        do row = 1, size(vectors, 1)
            vectors(row, 1) = sin(0.37_dp * real(row, dp))
            vectors(row, 2) = cos(0.29_dp * real(row, dp))
        end do
        serial_jvp = quadratic(vectors(:, 1), serial_dk) &
            + quadratic(vectors(:, 2), serial_dk)
        call build_compatible_pressure_trace_gradient(equilibrium, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], &
            0, 2, 64, 8, vectors, serial_gradient, status, radial_cells=6)
        call require(status == compatible_three_component_ok, &
            'serial pressure reverse assembly failed')
        call omp_set_num_threads(2)
        actual_threads = 0
        !$omp parallel reduction(max:actual_threads)
        actual_threads = omp_get_num_threads()
        !$omp end parallel
        call require(actual_threads == 2, 'pressure thread oracle requires two workers')
        do repeat = 1, 2
            call assemble(equilibrium, 0, 6, .false., parallel_k, parallel_m)
            call assemble(equilibrium, 0, 6, .false., parallel_dk, tangent_m, direction)
            call require(maxval(abs(parallel_k - serial_k)) &
                < 2.0e-12_dp * max(1.0_dp, maxval(abs(serial_k))), &
                'pressure primal matrix changed with assembly thread count')
            call require(maxval(abs(parallel_m - serial_m)) &
                < 2.0e-12_dp * max(1.0_dp, maxval(abs(serial_m))), &
                'pressure mass matrix changed with assembly thread count')
            call require(maxval(abs(tangent_m - serial_m)) &
                < 2.0e-12_dp * max(1.0_dp, maxval(abs(serial_m))), &
                'pressure tangent mass changed with assembly thread count')
            call require(maxval(abs(parallel_dk - serial_dk)) &
                < 2.0e-12_dp * max(1.0_dp, maxval(abs(serial_dk))), &
                'pressure tangent matrix changed with assembly thread count')
            jvp = quadratic(vectors(:, 1), parallel_dk) &
                + quadratic(vectors(:, 2), parallel_dk)
            call require(abs(jvp - serial_jvp) &
                < 2.0e-12_dp * max(1.0_dp, abs(serial_jvp)), &
                'pressure JVP changed with assembly thread count')
            call build_compatible_pressure_trace_gradient(equilibrium, 1.4_dp, &
                1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], &
                0, 2, 64, 8, vectors, gradient, status, radial_cells=6)
            call require(status == compatible_three_component_ok, &
                'parallel pressure reverse assembly failed')
            call require(maxval(abs(gradient - serial_gradient)) &
                < 2.0e-12_dp * max(1.0_dp, maxval(abs(serial_gradient))), &
                'pressure VJP changed with assembly thread count')
        end do
        call omp_set_num_threads(previous_threads)
        call omp_set_dynamic(previous_dynamic)
    end subroutine check_parallel_invariance

    subroutine assemble(data, parity, cells, sparse, matrix, mass, tangent, breaks)
        type(gvec_cas3d_equilibrium_t), intent(in) :: data
        integer, intent(in) :: parity, cells
        logical, intent(in) :: sparse
        real(dp), allocatable, intent(out) :: matrix(:, :), mass(:, :)
        real(dp), optional, intent(in) :: tangent(:)
        real(dp), optional, intent(in) :: breaks(:)
        type(compatible_three_component_problem_t) :: problem
        integer :: status

        if (present(tangent)) then
            call build_compatible_pressure_stiffness_tangent(data, 1.4_dp, 1.0_dp, &
                [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], parity, 2, 64, 8, &
                tangent, problem, status, sparse_storage=sparse, radial_cells=cells, &
                radial_breaks=breaks)
        else
            call build_compatible_three_component_problem(data, 1.4_dp, 1.0_dp, &
                [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], parity, 2, 64, 8, &
                problem, status, sparse_storage=sparse, radial_cells=cells, &
                radial_breaks=breaks)
        end if
        call require(status == compatible_three_component_ok, &
            'pressure assembly failed')
        if (sparse) then
            call variable_block_to_dense(problem%sparse_stiffness, matrix, status)
            call require(status == variable_block_ok, &
                'sparse pressure derivative failed')
            call variable_block_to_dense(problem%sparse_mass, mass, status)
            call require(status == variable_block_ok, 'sparse pressure mass failed')
        else
            call move_alloc(problem%stiffness, matrix)
            call move_alloc(problem%mass, mass)
        end if
    end subroutine assemble

    subroutine check_direction(parity, cells, sparse)
        integer, intent(in) :: parity, cells
        logical, intent(in) :: sparse
        real(dp), parameter :: steps(3) = [1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp]
        type(gvec_cas3d_equilibrium_t) :: plus, minus
        real(dp), allocatable :: derivative(:, :), mass(:, :), kp(:, :), km(:, :)
        real(dp), allocatable :: mp(:, :), mm(:, :)
        real(dp) :: error, scale, step
        integer :: index

        call assemble(equilibrium, parity, cells, sparse, derivative, mass, direction)
        scale = maxval(abs(derivative))
        call require(scale > 0.0_dp, 'pressure tangent is identically zero')
        plus = equilibrium
        minus = equilibrium
        do index = 1, size(steps)
            step = steps(index)
            plus%pressure = equilibrium%pressure + step * direction
            minus%pressure = equilibrium%pressure - step * direction
            call assemble(plus, parity, cells, sparse, kp, mp)
            call assemble(minus, parity, cells, sparse, km, mm)
            error = maxval(abs((kp - km) / (2.0_dp * step) - derivative)) / scale
            print '(a,i1,a,i2,a,l1,a,es10.2)', 'parity=', parity, ' cells=', cells, &
                ' sparse=', sparse, ' relative FD=', error
            call require(error < 2.0e-6_dp, 'full pressure K tangent disagrees with FD')
            call require(maxval(abs(mp - mass)) == 0.0_dp, 'pressure changed mass')
            call require(maxval(abs(mm - mass)) == 0.0_dp, 'pressure changed mass')
        end do
    end subroutine check_direction

    subroutine check_duality()
        real(dp), allocatable :: derivative(:, :), mass(:, :), vector(:), basis(:)
        real(dp), allocatable :: gradient(:), reverse_gradient(:), vectors(:, :)
        real(dp) :: jvp, vjp
        integer :: index

        call assemble(equilibrium, 1, 6, .false., derivative, mass, direction)
        allocate (vector(size(derivative, 1)), basis(size(direction)), &
            gradient(size(direction)))
        do index = 1, size(vector)
            vector(index) = sin(0.37_dp * real(index, dp))
        end do
        jvp = quadratic(vector, derivative)
        do index = 1, size(direction)
            basis = 0.0_dp
            basis(index) = 1.0_dp
            call assemble(equilibrium, 1, 6, .false., derivative, mass, basis)
            gradient(index) = quadratic(vector, derivative)
        end do
        vjp = dot_product(gradient, direction)
        call require(abs(jvp - vjp) < 1.0e-10_dp * max(abs(jvp), abs(vjp)), &
            'pressure sample JVP/VJP duality failed')
        allocate (vectors(size(vector), 1))
        vectors(:, 1) = vector
        call build_compatible_pressure_trace_gradient(equilibrium, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], &
            1, 2, 64, 8, vectors, reverse_gradient, info, radial_cells=6)
        call require(info == compatible_three_component_ok, &
            'contracted pressure reverse assembly failed')
        call require(maxval(abs(reverse_gradient - gradient)) &
            < 1.0e-10_dp * maxval(abs(gradient)), &
            'reverse gradient differs from independent cardinal assemblies')
        call require(abs(dot_product(reverse_gradient, direction) - jvp) &
            < 1.0e-10_dp * abs(jvp), 'contracted reverse/JVP duality failed')
    end subroutine check_duality

    subroutine check_nonuniform_mesh()
        real(dp), parameter :: breaks(7) = &
            [0.0_dp, 0.02_dp, 0.07_dp, 0.18_dp, 0.37_dp, 0.64_dp, 1.0_dp]
        real(dp), parameter :: steps(3) = [1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp]
        type(gvec_cas3d_equilibrium_t) :: plus, minus
        type(fixed_boundary_problem_t) :: fixed
        type(fixed_boundary_energy_terms_t) :: energy
        real(dp), allocatable :: stiffness(:, :), derivative(:, :), mass(:, :)
        real(dp), allocatable :: kp(:, :), km(:, :), mp(:, :), mm(:, :)
        real(dp), allocatable :: vectors(:, :), gradient(:)
        real(dp) :: jvp, vjp, scale, error, projection, expected
        integer :: row, column, index, status

        ! One fixed nonuniform FE mesh exercises the forward and both pressure
        ! actions independently of the equilibrium's pressure-sample locations.
        call assemble(equilibrium, 1, 6, .false., stiffness, mass, breaks=breaks)
        call assemble(equilibrium, 1, 6, .false., derivative, mp, direction, breaks)
        allocate (vectors(size(stiffness, 1), 2))
        do row = 1, size(vectors, 1)
            vectors(row, 1) = sin(0.37_dp * real(row, dp))
            vectors(row, 2) = cos(0.29_dp * real(row, dp))
        end do
        jvp = quadratic(vectors(:, 1), derivative) &
            + quadratic(vectors(:, 2), derivative)
        call build_compatible_pressure_trace_gradient(equilibrium, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], &
            1, 2, 64, 8, vectors, gradient, status, radial_breaks=breaks)
        call require(status == compatible_three_component_ok, &
            'nonuniform pressure reverse assembly failed')
        vjp = dot_product(gradient, direction)
        call require(abs(jvp - vjp) < 1.0e-10_dp * max(1.0_dp, abs(jvp)), &
            'nonuniform multivector pressure duality failed')
        scale = maxval(abs(derivative))
        call require(scale > 0.0_dp, 'nonuniform pressure tangent is zero')
        plus = equilibrium
        minus = equilibrium
        do index = 1, size(steps)
            plus%pressure = equilibrium%pressure + steps(index) * direction
            minus%pressure = equilibrium%pressure - steps(index) * direction
            call assemble(plus, 1, 6, .false., kp, mp, breaks=breaks)
            call assemble(minus, 1, 6, .false., km, mm, breaks=breaks)
            error = maxval(abs((kp - km) / (2.0_dp * steps(index)) &
                - derivative)) / scale
            call require(error < 2.0e-6_dp, 'nonuniform pressure FD plateau failed')
            call require(maxval(abs(mp - mass)) == 0.0_dp, &
                'nonuniform positive pressure step changed mass')
            call require(maxval(abs(mm - mass)) == 0.0_dp, &
                'nonuniform negative pressure step changed mass')
        end do

        ! The public spectral actions require an M-orthonormal basis. Use the
        ! independently assembled mass to test that the constructor retains
        ! the explicit mesh for its primal operator and later pressure actions.
        vectors(:, 1) = vectors(:, 1) / sqrt(quadratic(vectors(:, 1), mass))
        projection = 0.0_dp
        do column = 1, size(vectors, 1)
            do row = 1, size(vectors, 1)
                projection = projection + vectors(row, 1) &
                    * mass(row, column) * vectors(column, 2)
            end do
        end do
        vectors(:, 2) = vectors(:, 2) - projection * vectors(:, 1)
        vectors(:, 2) = vectors(:, 2) / sqrt(quadratic(vectors(:, 2), mass))
        call build_fixed_boundary_problem(equilibrium, 1.4_dp, 1.0_dp, &
            1.0e-12_dp, [0, 1, 2], [1, 1, 1], 2, fixed, status, &
            angular_theta=64, angular_zeta=8, radial_breaks=breaks)
        call require(status == fixed_boundary_ok, 'nonuniform fixed problem failed')
        call diagnose_fixed_boundary_energy(fixed, 1, vectors(:, 1), energy, status)
        call require(status == fixed_boundary_ok, 'nonuniform primal diagnosis failed')
        expected = quadratic(vectors(:, 1), stiffness)
        call require(abs(energy%potential_energy - expected) &
            < 1.0e-10_dp * max(1.0_dp, abs(expected)), &
            'fixed constructor did not retain the nonuniform primal mesh')
        expected = quadratic(vectors(:, 1), derivative) &
            + quadratic(vectors(:, 2), derivative)
        call fixed_boundary_pressure_trace_jvp(fixed, equilibrium, 1, vectors, &
            direction, jvp, status)
        call require(status == fixed_boundary_ok, &
            'nonuniform public pressure JVP failed')
        call require(abs(jvp - expected) < 1.0e-10_dp * max(1.0_dp, abs(expected)), &
            'public pressure JVP used a different radial mesh')
        call fixed_boundary_pressure_trace_vjp(fixed, equilibrium, 1, vectors, &
            1.0_dp, gradient, status)
        call require(status == fixed_boundary_ok, &
            'nonuniform public pressure VJP failed')
        vjp = dot_product(gradient, direction)
        call require(abs(vjp - expected) < 1.0e-10_dp * max(1.0_dp, abs(expected)), &
            'public pressure VJP used a different radial mesh')
    end subroutine check_nonuniform_mesh

    subroutine check_period_contraction()
        integer, parameter :: trials = 3, columns = 18, nt = 16, nz = 12
        real(dp) :: fields(nt, nz, 13), fields_dot(nt, nz, 13)
        real(dp) :: drive(nt, nz), drive_dot(nt, nz), gamma_p(nt, nz), jac(nt, nz)
        real(dp) :: h1(2, trials), dh1(2, trials), eta(2, trials), l2(2, trials)
        real(dp) :: vectors(columns, 2), stiffness(columns, columns)
        real(dp) :: terms(columns, columns, 5), slope_bar, gamma_bar, expected
        real(dp) :: angle, scale
        logical :: mask(5)
        integer :: j, k, a, component, selection, status, period_count

        ! Mixed Fourier residue classes exercise both period masks. The
        ! oracle assembles full matrices; the reverse path contracts first.
        fields = 1.0_dp
        fields_dot = 0.0_dp
        jac = 0.13_dp
        gamma_p = 0.7_dp
        do k = 1, nz
            do j = 1, nt
                angle = 2.0_dp * acos(-1.0_dp) &
                    * (real(j - 1, dp) / nt - real(k - 1, dp) / nz)
                fields(j, k, 7) = 1.1_dp + 0.1_dp * cos(angle)
                fields(j, k, 8) = 4.0_dp
                fields(j, k, 9) = 3.0_dp
                fields_dot(j, k, 10) = 0.4_dp * cos(angle)
                fields_dot(j, k, 11) = 0.2_dp
                fields_dot(j, k, 13) = 0.05_dp * sin(angle)
                drive(j, k) = 0.3_dp * cos(angle)
                drive_dot(j, k) = 0.7_dp * sin(angle)
            end do
        end do
        h1 = 0.7_dp
        dh1 = -0.3_dp
        eta = 0.2_dp
        l2 = 0.6_dp
        do a = 1, columns
            vectors(a, 1) = sin(0.37_dp * a)
            vectors(a, 2) = cos(0.29_dp * a)
        end do
        do period_count = 1, 5, 2
            stiffness = 0.0_dp
            terms = 0.0_dp
            call assemble_compatible_pressure_stiffness_tangent(fields, drive, &
                jac, jac, jac, gamma_p, fields_dot, drive_dot, 0.17_dp, &
                [0, 1, 2], [1, 2, 3], [1, 2, 1], period_count, &
                h1, dh1, eta, l2, 0.23_dp, stiffness, terms, status)
            call require(status == 0, 'period tangent oracle failed')
            do selection = 1, 3
                mask = .true.
                if (selection == 2) mask = &
                    [.true., .true., .false., .true., .false.]
                if (selection == 3) mask = &
                    [.false., .false., .true., .false., .true.]
                call contract_compatible_pressure_stiffness(fields, drive, &
                    jac, jac, jac, gamma_p, fields_dot, drive_dot, &
                    [0, 1, 2], [1, 2, 3], [1, 2, 1], period_count, &
                    h1, dh1, eta, l2, 0.23_dp, vectors, mask, &
                    slope_bar, gamma_bar, status)
                call require(status == 0, 'period reverse contraction failed')
                expected = 0.0_dp
                scale = 0.0_dp
                do component = 2, 5
                    if (.not. mask(component)) cycle
                    do a = 1, 2
                        expected = expected + quadratic(vectors(:, a), &
                            terms(:, :, component))
                        scale = scale + abs(quadratic(vectors(:, a), &
                            terms(:, :, component)))
                    end do
                end do
                call require(abs(slope_bar + 0.17_dp * gamma_bar - expected) &
                    < 2.0e-12_dp * max(1.0_dp, scale), &
                    'contracted reverse disagrees with masked period oracle')
            end do
        end do
    end subroutine check_period_contraction

    pure function quadratic(vector, matrix) result(value)
        real(dp), intent(in) :: vector(:), matrix(:, :)
        real(dp) :: value
        integer :: a, b

        value = 0.0_dp
        do b = 1, size(vector)
            do a = 1, size(vector)
                value = value + vector(a) * matrix(a, b) * vector(b)
            end do
        end do
    end function quadratic

    subroutine check_invalid()
        type(gvec_cas3d_equilibrium_t) :: invalid
        type(compatible_three_component_problem_t) :: problem
        real(dp), allocatable :: tangent(:)
        integer :: status

        call build_compatible_pressure_stiffness_tangent(equilibrium, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], 1, &
            2, 32, 8, direction(:size(direction) - 1), problem, status)
        call require(status == compatible_three_component_invalid, &
            'wrong pressure direction shape accepted')
        tangent = direction
        tangent(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        call build_compatible_pressure_stiffness_tangent(equilibrium, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], 1, &
            2, 32, 8, tangent, problem, status)
        call require(status == compatible_three_component_invalid, &
            'nonfinite pressure direction accepted')
        invalid = equilibrium
        invalid%pressure(2) = 0.0_dp
        call build_compatible_pressure_stiffness_tangent(invalid, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], 1, &
            2, 32, 8, direction, problem, status)
        call require(status == compatible_three_component_invalid, &
            'closed pressure parameter domain accepted')
    end subroutine check_invalid

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_pressure_stiffness_tangent
