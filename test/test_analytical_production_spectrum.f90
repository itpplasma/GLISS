program test_analytical_production_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use exact_cylinder, only: assemble_exact_cylinder, exact_cylinder_parameters_t
    use physical_constants, only: vacuum_permeability
    use symmetric_eigensolver, only: solve_symmetric_generalized_allocated
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use dense_spectrum_support, only: refine_dense_spectrum, refine_dense_eigenpair, &
        certify_dense_spectrum_orthogonality, certify_dense_spectrum_inertia, &
        diagnose_dense_spectrum
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_tridiagonal_t
    use variable_generalized_solver, only: variable_eigenvalue_bound, &
        variable_generalized_inertia, variable_generalized_diagnostics, &
        pencil_roundoff
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Independently tabulated J_3' zeros; checked with intrinsic Bessel functions.
    real(dp), parameter :: roots(3) = &
        [4.201188941210528_dp, 8.015236598375953_dp, 11.345924310743007_dp]
    integer, parameter :: meshes(3) = [8, 16, 32]
    ! Fixed engineering accuracy targets, not variational error bounds.
    real(dp), parameter :: budgets(4) = &
        [1.0e-3_dp, 1.0e-5_dp, 1.0e-7_dp, 1.0e-8_dp]
    type(exact_cylinder_parameters_t) :: physical, changed
    real(dp) :: errors(4, 3, 2, 4), branches(5), baseline(5), factor
    real(dp) :: derivative, worst, case_errors(4), scaling_error
    integer :: mesh, degree, parity, branch, parameter_case

    do branch = 1, size(roots)
        derivative = 0.5_dp*(bessel_jn(2, roots(branch)) &
            - bessel_jn(4, roots(branch)))
        call require(abs(derivative) < 4.0e-15_dp, 'Bessel root is invalid')
    end do
    ! cs = vA and ka = 1 separate the slow branches without a fitted reference.
    physical = exact_cylinder_parameters_t()
    physical%pressure = physical%magnetic_field**2 &
        /(vacuum_permeability*physical%adiabatic_index)
    physical%length = 2.0_dp*pi*physical%radius
    baseline = 0.0_dp
    do degree = 1, 4
        do parity = 1, 2
            do mesh = 1, size(meshes)
                call check_case(meshes(mesh), degree, parity, physical, &
                    branches, errors(:, mesh, parity, degree))
                if (mesh == 2 .and. degree == 3 .and. parity == 1) &
                    baseline = branches
            end do
        end do
        worst = maxval(abs(errors(:, 3, :, degree)))
        write (*, '(a,i0,a,es12.4)') 'degree ', degree, &
            ' finest fast/slow relative error: ', worst
        call require(worst < budgets(degree), &
            'finest branch error exceeds degree accuracy target')
    end do
    call require(maxval(abs(errors)) < 2.0e-2_dp, &
        'coarse branch error exceeds 2%')
    call require(maxval(abs(errors(:, :, 1, :) - errors(:, :, 2, :))) &
        < 1.0e-8_dp, 'parity spectra disagree')
    ! Require useful refinement, without assuming variational upper bounds.
    do degree = 1, 4
        do parity = 1, 2
            do branch = 1, 4
                do mesh = 2, size(meshes)
                    call require(abs(errors(branch, mesh, parity, degree)) &
                        <= max(0.75_dp*abs(errors(branch, mesh - 1, &
                        parity, degree)), 0.1_dp*budgets(degree)), &
                        'radial refinement does not improve branch')
                end do
            end do
        end do
    end do

    call require(all(baseline > 0.0_dp), 'scaling baseline is missing')
    do parameter_case = 1, 6
        changed = physical
        select case (parameter_case)
        case (1)
            changed%density = 2.5_dp*physical%density
            factor = 1.0_dp/2.5_dp
        case (2)
            changed%magnetic_field = 1.5_dp*physical%magnetic_field
            changed%pressure = 1.5_dp**2*physical%pressure
            factor = 1.5_dp**2
        case (3)
            changed%radius = 1.4_dp*physical%radius
            changed%length = 1.4_dp*physical%length
            factor = 1.0_dp/1.4_dp**2
        case (4)
            changed%pressure = 2.0_dp*physical%pressure
            factor = 0.0_dp
        case (5)
            changed%magnetic_field = 1.5_dp*physical%magnetic_field
            factor = 0.0_dp
        case (6)
            changed%length = 1.3_dp*physical%length
            factor = 0.0_dp
        end select
        call check_case(16, 3, 1, changed, branches, case_errors)
        call require(maxval(abs(case_errors)) < budgets(3), &
            'physical parameter case exceeds degree accuracy target')
        write (*, '(a,i0,a,es12.4)') 'physical parameter case ', &
            parameter_case, ' fast/slow relative error: ', &
            maxval(abs(case_errors))
        if (factor > 0.0_dp) then
            scaling_error = maxval(abs(branches/(factor*baseline) - 1.0_dp))
            call require(scaling_error < 1.0e-8_dp, &
                'physical spectral scaling failed')
            write (*, '(a,es12.4)') 'physical scaling relative error: ', &
                scaling_error
        end if
    end do
    write (*, '(a)') 'Analytical production assembly/certifier spectrum passed'

contains

    subroutine references(parameters, slow, fast, alfven)
        type(exact_cylinder_parameters_t), intent(in) :: parameters
        real(dp), intent(out) :: slow(3), fast(3), alfven
        real(dp) :: va2, cs2, k2, q2(3), sum_term(3), product_term(3)
        integer :: i

        va2 = parameters%magnetic_field**2 &
            /(vacuum_permeability*parameters%density)
        cs2 = parameters%adiabatic_index*parameters%pressure/parameters%density
        k2 = (2.0_dp*pi/parameters%length)**2
        do i = 1, size(roots)
            q2(i) = k2 + (roots(i)/parameters%radius)**2
        end do
        sum_term = (va2 + cs2)*q2
        product_term = va2*cs2*k2*q2
        fast = 0.5_dp*(sum_term &
            + sqrt(sum_term**2 - 4.0_dp*product_term))
        slow = product_term/fast
        alfven = va2*k2
    end subroutine references

    subroutine check_case(cells, polynomial_degree, phase, parameters, &
            measured, relative_errors)
        integer, intent(in) :: cells, polynomial_degree, phase
        type(exact_cylinder_parameters_t), intent(in) :: parameters
        real(dp), intent(out) :: measured(5)
        real(dp), optional, intent(out) :: relative_errors(4)
        type(variable_block_tridiagonal_t) :: stiffness, mass
        type(fixed_boundary_solver_controls_t) :: controls
        real(dp), allocatable :: dense_k(:, :), dense_m(:, :), work_k(:, :)
        real(dp), allocatable :: work_m(:, :), values(:), vectors(:, :)
        real(dp), allocatable :: refined_vector(:)
        real(dp) :: slow(3), fast(3), alfven, error(4), lower, upper
        real(dp) :: bound, slow_upper, fast_lower, reference(4)
        real(dp) :: lower_bounds(4), upper_bounds(4)
        real(dp) :: accuracy_target, gap_target, refined_value, residual
        real(dp) :: resolution, interval, quotient
        integer :: info, normal, tangential, nslow, nfast, indices(4), i, n
        integer :: widths(1)
        logical :: replaced

        call assemble_exact_cylinder(cells, polynomial_degree, phase, &
            dense_k, dense_m, parameters=parameters, normal_dofs=normal, &
            tangential_dofs=tangential)
        n = size(dense_k, 1)
        widths(1) = n
        call pack_variable_blocks(dense_k, widths, stiffness, info)
        call require(info == 0, 'stiffness packing failed')
        call pack_variable_blocks(dense_m, widths, mass, info)
        call require(info == 0, 'mass packing failed')
        work_k = dense_k
        work_m = dense_m
        call solve_symmetric_generalized_allocated(work_k, work_m, values, &
            vectors, info, equilibrate=.true.)
        call require(info == 0, 'equilibrated dense seed solve failed')
        call refine_dense_spectrum(stiffness, mass, controls, values, &
            vectors, info)
        call require(info == 0, 'production spectrum refinement failed')
        call certify_spectrum(stiffness, mass, values, vectors)
        call references(parameters, slow, fast, alfven)
        call check_inertia(stiffness, mass, 0.0_dp, 0)
        slow_upper = 0.5_dp*(slow(1) + min(alfven, &
            parameters%adiabatic_index*parameters%pressure &
            /parameters%density*(2.0_dp*pi/parameters%length)**2))
        fast_lower = 0.5_dp*(max(alfven, &
            parameters%adiabatic_index*parameters%pressure &
            /parameters%density*(2.0_dp*pi/parameters%length)**2) + fast(1))
        nslow = count(values < slow_upper)
        nfast = count(values > fast_lower)
        write (*, '(a,3i4,a,3i5)') 'case cells/degree/parity ', cells, &
            polynomial_degree, phase, ' slow/Alfven/fast counts ', &
            nslow, n - nslow - nfast, nfast
        call require(nslow == tangential, 'slow cluster dimension differs')
        ! The m != 0 tangential divergence has full L2 rank. The two
        ! magnetosonic families each have dim(L2); its divergence-free
        ! complement has dim(H1), rather than one Alfven mode per L2 dof.
        call require(nfast == tangential, 'fast cluster dimension differs')
        call require(n - nslow - nfast == normal, &
            'Alfven cluster dimension differs')
        call require(maxval(abs(values(nslow + 1:n - nfast)/alfven &
            - 1.0_dp)) < 1.0e-8_dp, 'Alfven cluster frequency differs')
        indices(1) = nslow
        indices(2) = nslow - 1
        indices(3) = n - nfast + 1
        indices(4) = n - nfast + 2
        reference(:2) = slow(:2)
        reference(3:) = fast(:2)
        lower_bounds(1) = 0.5_dp*(slow(1) + slow(2))
        lower_bounds(2) = 0.5_dp*(slow(2) + slow(3))
        lower_bounds(3) = fast_lower
        lower_bounds(4) = 0.5_dp*(fast(1) + fast(2))
        upper_bounds(1) = slow_upper
        upper_bounds(2) = lower_bounds(1)
        upper_bounds(3) = lower_bounds(4)
        upper_bounds(4) = 0.5_dp*(fast(2) + fast(3))
        replaced = .false.
        do i = 1, 4
            lower = lower_bounds(i)
            upper = upper_bounds(i)
            call check_inertia(stiffness, mass, lower, indices(i) - 1)
            call check_inertia(stiffness, mass, upper, indices(i))
            call variable_eigenvalue_bound(stiffness, mass, &
                vectors(:, indices(i)), values(indices(i)), bound, info)
            call require(info == 0, 'mass residual bound failed')
            accuracy_target = 0.1_dp*budgets(polynomial_degree)*reference(i)
            gap_target = 0.25_dp*min(reference(i) - lower, upper - reference(i))
            if (bound >= min(accuracy_target, gap_target)) then
                ! The dense seed's backward-error stopping rule includes
                ! roundoff/resolution floors. An isolated physical branch
                ! must also meet this stricter mass-residual target.
                write (*, '(a,2i5,a,3es24.15)') 'refine branch/index ', i, &
                    indices(i), ' bound/accuracy/gap ', bound, &
                    accuracy_target, gap_target
                call variable_generalized_diagnostics(stiffness, mass, &
                    vectors(:, indices(i)), values(indices(i)), quotient, &
                    residual, resolution, info)
                call require(info == 0, 'branch seed diagnostics failed')
                write (*, '(a,3es24.15)') 'seed residual/resolution/roundoff ', &
                    residual, resolution, pencil_roundoff(stiffness, mass)
                call refine_dense_eigenpair(stiffness, mass, controls, &
                    values, indices(i), vectors(:, indices(i)), refined_value, &
                    refined_vector, residual, resolution, info, interval)
                call require(info == 0, 'indexed branch refinement failed')
                call variable_eigenvalue_bound(stiffness, mass, refined_vector, &
                    refined_value, bound, info)
                call require(info == 0, 'refined mass residual bound failed')
                write (*, '(a,4es24.15)') &
                    'refined bound/residual/resolution/interval ', &
                    bound, residual, resolution, interval
                values(indices(i)) = refined_value
                vectors(:, indices(i)) = refined_vector
                replaced = .true.
            end if
            call require(bound < gap_target, 'branch uncertainty exceeds gap')
            call require(bound < accuracy_target, &
                'branch uncertainty exceeds accuracy target')
            call require(values(indices(i)) - bound > lower, &
                'branch enclosure crosses lower analytical window')
            call require(values(indices(i)) + bound < upper, &
                'branch enclosure crosses upper analytical window')
            measured(i) = values(indices(i))
        end do
        ! Certificates on the dense seeds do not certify replaced columns.
        if (replaced) call certify_spectrum(stiffness, mass, values, vectors)
        measured(5) = sum(values(nslow + 1:n - nfast))/real(normal, dp)
        error = measured(:4)/reference - 1.0_dp
        call require(maxval(abs(error)) < 2.0e-2_dp, &
            'analytical branch error exceeds 2%')
        if (present(relative_errors)) relative_errors = error
    end subroutine check_case

    subroutine certify_spectrum(stiffness, mass, values, vectors)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: values(:)
        real(dp), contiguous, intent(in) :: vectors(:, :)
        real(dp), allocatable :: quotients(:), residuals(:), resolutions(:)
        integer :: info

        call require(all(values(2:) >= values(:size(values) - 1)), &
            'spectrum order changed after refinement')
        call certify_dense_spectrum_orthogonality(mass, vectors, info)
        call require(info == 0, 'production orthogonality certificate failed')
        call diagnose_dense_spectrum(stiffness, mass, values, vectors, &
            quotients, residuals, resolutions, info)
        call require(info == 0, 'production spectrum diagnostics failed')
        call certify_dense_spectrum_inertia(stiffness, mass, values, &
            residuals, resolutions, info)
        call require(info == 0, 'production gap inertia certificate failed')
    end subroutine certify_spectrum

    subroutine check_inertia(stiffness, mass, shift, expected)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        integer, intent(in) :: expected
        integer :: info, count

        call variable_generalized_inertia(stiffness, mass, shift, count, info)
        call require(info == 0, 'independent branch-window inertia failed')
        call require(count == expected, 'analytical branch identity failed')
    end subroutine check_inertia

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, '(a)') trim(message)
        error stop 1
    end subroutine require
end program test_analytical_production_spectrum
