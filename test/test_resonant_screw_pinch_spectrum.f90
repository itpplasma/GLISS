! Isolated unstable ground eigenvalue of an exact resonant screw pinch.
! References solve the independent perpendicular Euler--Lagrange ODE;
! production matrices use current GLOBAL radial FEEC and surface kernels.
! This does not claim a toroidal primitive import or a gap at zero:
! infinitely many negative modes accumulate at zero, and parallel modes
! are exactly neutral.  Only the most-negative isolated mode is qualified.
program test_resonant_screw_pinch_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dense_spectrum_support, only: refine_dense_eigenpair
    use fixed_boundary_eigen_bracket, only: certify_lowest_bracket
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use resonant_screw_pinch_fixture, only: assemble_resonant_screw_pinch, &
        pinch_blinear, pinch_bz, pinch_lambda_scale, pinch_pressure_edge, &
        pinch_radius
    use physical_constants, only: vacuum_permeability
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_tridiagonal_t
    use variable_generalized_solver, only: variable_eigenvalue_bound, &
        variable_generalized_inertia
    implicit none

    integer, parameter :: meshes(3) = [16, 32, 64]
    real(dp), parameter :: budgets(3) = [1.0e-3_dp, 1.0e-5_dp, 1.0e-7_dp]
    real(dp), parameter :: gammas(2) = [1.0_dp, 5.0_dp / 3.0_dp]
    ! ell=lambda*mu0*rho*a^2/Z^2. Independent 80-digit RK4/Richardson
    ! and DOP853/continuous Pruefer shooting agree below 2e-10 relative.
    ! Regular axis: f~r^(m-1), P~-lambda*rho*r^m/m; wall f(a)=0.
    ! Ground/first overtone phases are pi/2pi, with zero/one interior node.
    real(dp), parameter :: roots(2) = &
        [-0.02435118850778054_dp, -0.02266142767811291_dp]
    real(dp), parameter :: overtones(2) = &
        [-0.0092183751438407_dp, -0.0085296207819971_dp]
    real(dp) :: errors(3, 3, 2, 2), value, error, bound
    integer :: gamma, degree, parity, mesh

    do gamma = 1, size(gammas)
        do degree = 1, size(budgets)
            do parity = 1, 2
                do mesh = 1, size(meshes)
                    call check_case(meshes(mesh), degree, parity, gamma, &
                        value, bound)
                    error = value / roots(gamma) - 1.0_dp
                    errors(mesh, degree, parity, gamma) = error
                    write (*, '(a,4i4,a,3es23.14)') &
                        'gamma/degree/parity/cells ', gamma, degree, parity, &
                        meshes(mesh), ' ell/error/bound ', value, error, bound
                    call require(abs(error) < 0.02_dp, &
                        'coarse unstable ground error exceeds 2%')
                end do
                call require(abs(errors(3, degree, parity, gamma)) &
                    < budgets(degree), 'finest unstable ground accuracy failed')
                do mesh = 2, size(meshes)
                    call require(abs(errors(mesh, degree, parity, gamma)) &
                        <= max(0.8_dp * abs(errors(mesh - 1, degree, parity, &
                        gamma)), 0.1_dp * budgets(degree)), &
                        'unstable ground refinement did not improve')
                end do
            end do
        end do
    end do
    call require(maxval(abs(errors(:, :, 1, :) - errors(:, :, 2, :))) &
        < 1.0e-8_dp, 'unstable ground parity mismatch')
    write (*, '(a)') 'resonant screw pinch global unstable spectrum: PASS'

contains

    subroutine check_case(cells, degree, parity, fixture, value, bound)
        integer, intent(in) :: cells, degree, parity, fixture
        real(dp), intent(out) :: value, bound
        type(variable_block_tridiagonal_t) :: stiffness, mass
        type(fixed_boundary_solver_controls_t) :: controls
        real(dp), allocatable :: dense_k(:, :), dense_m(:, :), seed_k(:, :)
        real(dp), allocatable :: seed_m(:, :), values(:), vectors(:, :), vector(:)
        real(dp) :: residual, resolution, interval, margin, below, above
        real(dp) :: b2, hwall, lower, upper, scale
        integer :: status, normal, tangential, n, widths(1), count

        call assemble_resonant_screw_pinch(cells, degree, parity, &
            gammas(fixture), dense_k, dense_m, normal, tangential)
        ! Only a common units conversion; the assembled physical mass is
        ! retained. Ground ell is O(1e-2), avoiding huge dimensional shifts.
        dense_k = dense_k / pinch_lambda_scale
        n = size(dense_k, 1)
        call require(n == normal + 2 * tangential, 'global dimension mismatch')
        ! With eta=0, nu is purely field parallel. At global resonance
        ! curl(xi cross B)=div(xi)=0: every such column is an exact null.
        scale = maxval(abs(dense_k))
        call require(maxval(abs(dense_k(:, normal + tangential + 1:))) &
            < 1.0e-12_dp * scale, 'parallel resonant neutral space is missing')
        widths(1) = n
        call pack_variable_blocks(dense_k, widths, stiffness, status)
        call require(status == 0, 'resonant stiffness packing failed')
        call pack_variable_blocks(dense_m, widths, mass, status)
        call require(status == 0, 'resonant mass packing failed')
        b2 = (pinch_blinear * pinch_radius / pinch_bz)**2
        hwall = 1.0_dp + b2 + gammas(fixture) * vacuum_permeability &
            * pinch_pressure_edge / pinch_bz**2
        lower = -4.0_dp * b2**2 / hwall * (1.0_dp + 1.0e-8_dp)
        upper = 0.5_dp * (roots(fixture) + overtones(fixture))
        call check_inertia(stiffness, mass, lower, 0)
        call check_inertia(stiffness, mass, upper, 1)
        seed_k = dense_k
        seed_m = dense_m
        call solve_symmetric_generalized_allocated(seed_k, seed_m, values, &
            vectors, status, equilibrate=.true.)
        call require(status == 0, 'resonant equilibrated dense seed failed')
        call require(values(1) < 0.0_dp, 'global instability is missing')
        call refine_dense_eigenpair(stiffness, mass, controls, values, 1, &
            vectors(:, 1), value, vector, residual, resolution, status, interval)
        call require(status == 0, 'resonant indexed eigenpair refinement failed')
        call variable_eigenvalue_bound(stiffness, mass, vector, value, bound, status)
        call require(status == 0, 'resonant mass residual bound failed')
        margin = 0.1_dp * budgets(degree) * abs(roots(fixture))
        call require(bound < margin, 'unstable eigenvalue uncertainty exceeds budget')
        call require(value - bound > lower .and. value + bound < upper, &
            'unstable eigenvalue enclosure crosses independent root window')
        call certify_lowest_bracket(stiffness, mass, 0, value, margin, &
            below, above, status)
        call require(status == 0, 'unstable ground inertia bracket failed')
        call variable_generalized_inertia(stiffness, mass, above, count, status)
        call require(status == 0, 'unstable upper bracket inertia failed')
        call require(count == 1, 'unstable bracket does not isolate ground mode')
    end subroutine check_case

    subroutine check_inertia(stiffness, mass, shift, expected)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        integer, intent(in) :: expected
        integer :: status, count

        call variable_generalized_inertia(stiffness, mass, shift, count, status)
        call require(status == 0, 'independent unstable window inertia failed')
        call require(count == expected, 'independent unstable window count differs')
    end subroutine check_inertia

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, '(a)') trim(message)
        error stop 1
    end subroutine require

end program test_resonant_screw_pinch_spectrum
