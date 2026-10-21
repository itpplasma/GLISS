program test_sparse_stability
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_precision, only: precision_policy_t, precision_preset, FP32, &
        STAGE_ESTIMATE, STAGE_CERTIFY
    use gliss_sparse_backend, only: sparse_precision_supported
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t, &
        variable_block_factor_t, pack_variable_blocks, factorize_sparse_blocks, &
        solve_variable_factored
    use variable_generalized_solver, only: variable_generalized_inertia
    implicit none
    integer, parameter :: n = 6
    real(dp), parameter :: d(n) = [-3.0_dp, -1.0_dp, 2.0_dp, &
        4.0_dp, 6.0_dp, 9.0_dp]
    type(variable_block_tridiagonal_t) :: stiffness, mass
    type(variable_block_factor_t) :: factor
    type(precision_policy_t) :: policy
    real(dp) :: l(n, n), ld(n, n), k(n, n), m(n, n), x(n), rhs(n)
    real(dp) :: shifts(0:n)
    integer :: i, info, negative
    logical :: ok

    policy = precision_preset('cpu-fp64', ok)
    call require(ok .and. sparse_precision_supported(policy), 'fp64 policy')
    call policy%set(STAGE_ESTIMATE, FP32, FP32, ok)
    call require(.not. ok, 'estimators must stay fp64')
    call policy%set(STAGE_CERTIFY, FP32, FP32, ok)
    call require(.not. ok, 'certificates must stay fp64')
    policy = precision_preset('cpu-mixed', ok)
    call require(ok .and. .not. sparse_precision_supported(policy), &
        'unsupported candidate precision must fail closed')

    ! Exact dyadic congruence: K=L diag(d) L^T, M=L L^T. Thus
    ! det(K-sigma M)=det(L)^2 product(d-sigma), without a dense oracle.
    ! The same identity is symbolically checked by the pinned LS upstream.
    l = 0.0_dp
    do i = 1, n
        l(i, i) = 1.0_dp
        x(i) = real(i, dp)
    end do
    do i = 2, n
        l(i, i - 1) = 0.25_dp
    end do
    do i = 1, n
        ld(:, i) = l(:, i) * d(i)
    end do
    k = matmul(ld, transpose(l))
    m = matmul(l, transpose(l))
    call pack_variable_blocks(k, [1, 2, 1, 2], stiffness, info)
    call require(info == 0, 'pack K')
    call pack_variable_blocks(m, [1, 2, 1, 2], mass, info)
    call require(info == 0, 'pack M')
    shifts(0) = d(1) - 1.0_dp
    shifts(1:n - 1) = 0.5_dp * (d(:n - 1) + d(2:))
    shifts(n) = d(n) + 1.0_dp
    do i = 0, n
        call variable_generalized_inertia(stiffness, mass, shifts(i), negative, info)
        call require(info == 0 .and. negative == i, 'exact pencil inertia')
    end do
    call factorize_sparse_blocks(stiffness, factor, info)
    call require(info == 0 .and. factor%uses_sparse, 'native sparse factor')
    call require(factor%count_reliable .and. factor%negative_count == 2, &
        'resolved Sylvester count')
    rhs = matmul(k, x)
    call solve_variable_factored(factor, rhs, info)
    call require(info == 0, 'retained solve')
    call require(maxval(abs(rhs - x)) < 1.0e-12_dp, 'manufactured solution')
    stiffness%diagonal(1)%values = 0.0_dp
    stiffness%lower(1)%values = 0.0_dp
    call factorize_sparse_blocks(stiffness, factor, info)
    call require(info /= 0 .or. .not. factor%count_reliable, &
        'singular inertia must not certify a sign')
    print *, 'Sparse inertia, manufactured solve and precision locks PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print *, message
            error stop 1
        end if
    end subroutine require
end program test_sparse_stability
