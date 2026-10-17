program metric_averaging_proof
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, real_valued
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, operator(+), operator(-), &
        operator(*), operator(**)
    implicit none

    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symengine
    type(suite_t) :: checks
    type(expr_t) :: x0, z0, x1, z1, half, quarter
    type(expr_t) :: average_squared_norm, squared_average_norm, jump_squared_norm
    type(expr_t) :: residual
    character(1024) :: report

    report = 'strict-identities.json'
    if (command_argument_count() == 1) call get_command_argument(1, report)
    if (command_argument_count() > 1) error stop 'expected optional report path'
    call arena%init()
    call assumptions%init(arena)
    x0 = sym(arena, 'x0')
    z0 = sym(arena, 'z0')
    x1 = sym(arena, 'x1')
    z1 = sym(arena, 'z1')
    call assume(assumptions, real_valued(x0))
    call assume(assumptions, real_valued(z0))
    call assume(assumptions, real_valued(x1))
    call assume(assumptions, real_valued(z1))
    half = rat(arena, 1_int64, 2_int64)
    quarter = rat(arena, 1_int64, 4_int64)
    native = make_native_engine(arena, assumptions)
    symengine = make_symengine_engine(arena)
    call suite_begin(checks, 'Two-dimensional tangent metric averaging identity')

    average_squared_norm = half*(x0**2 + z0**2 + x1**2 + z1**2)
    squared_average_norm = (half*(x0 + x1))**2 + (half*(z0 + z1))**2
    jump_squared_norm = (x1 - x0)**2 + (z1 - z0)**2
    residual = average_squared_norm - squared_average_norm - quarter*jump_squared_norm
    call check_identity(checks, native, 'quadratic averaging gap [native, real]', residual)
    call check_identity(checks, symengine, 'quadratic averaging gap [SymEngine]', residual)
    if (checks%failed /= 0) error stop 'identity failed or undecided'
    if (checks%by_probe /= 0) error stop 'probe evidence is forbidden'
    if (checks%proved /= checks%total) error stop 'unproved identity'
    call suite_end(checks, trim(report))
end program metric_averaging_proof
