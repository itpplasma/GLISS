program weighted_derivative_proof
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, positive, real_valued
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, exp, log, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    implicit none
    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symengine
    type(suite_t) :: checks
    type(expr_t) :: s, alpha, n, ns, weight, half, b, c, flux, pi, m, logweight, xs
    character(1024) :: report

    report = 'strict-identities.json'
    if (command_argument_count() == 1) call get_command_argument(1, report)
    call arena%init()
    call assumptions%init(arena)
    s = sym(arena, 's')
    alpha = sym(arena, 'alpha')
    b = sym(arena, 'b')
    c = sym(arena, 'c')
    flux = sym(arena, 'T')
    pi = sym(arena, 'pi')
    m = sym(arena, 'm')
    n = sym(arena, 'N')
    ns = sym(arena, 'N_s')
    xs = sym(arena, 'x_s')
    call assume(assumptions, positive(s))
    call assume(assumptions, real_valued(alpha))
    call assume(assumptions, real_valued(n))
    call assume(assumptions, real_valued(ns))
    call assume(assumptions, real_valued(xs))
    call assume(assumptions, real_valued(b))
    call assume(assumptions, real_valued(c))
    call assume(assumptions, positive(flux))
    call assume(assumptions, positive(pi))
    call assume(assumptions, positive(m))
    native = make_native_engine(arena, assumptions)
    symengine = make_symengine_engine(arena)
    call suite_begin(checks, 'Weighted derivative exact identities and kernels')
    half = rat(arena, 1_int64, 2_int64)
    ! On s>0, this is the real definition of s**(-alpha); its
    ! logarithmic form avoids symbolic-exponent branch normalization.
    logweight = -alpha*log(s)
    weight = exp(logweight)
    call check_identity(checks, native, &
        'unfactored weighted product derivative [native]', &
        ds(weight*n) - weight*(ns - alpha*n/s))
    if (checks%failed /= 0) error stop 'unfactored native identity failed'
    ! Divide the derivative identity by its positive real weight,
    ! using d exp(logweight)/ds = weight*d(logweight)/ds.
    call prove('weight-normalized product derivative with real alpha', &
        ds(logweight)*n + ns - (ns - alpha*n/s))
    call prove('alpha one kernel N equals s', diff(s**(-1)*s, s))
    call prove('alpha two kernel N equals s squared', diff(s**(-2)*s**2, s))
    call prove('m1 weight-normalized first-cell derivative and tie factor', &
        diff(-half*log(s), s)*b*s + b - half*b)
    call prove('physical preset axis leading coefficient', &
        (1 - (1 - m*half)) - m*half)
    call prove('zero-curl cylinder eta cancellation for arbitrary x_s', &
        xs + 2*pi*m*(-flux*xs/(2*pi*m))/flux)
    if (checks%by_probe /= 0) error stop 'probe evidence is forbidden'
    if (checks%proved /= checks%total) error stop 'unproved identity'
    call suite_end(checks, trim(report))
contains
    function ds(expression) result(derivative)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: derivative
        derivative = diff(expression, s) + diff(expression, n)*ns
    end function ds

    subroutine prove(label, residual)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: residual
        call check_identity(checks, native, trim(label)//' [native]', residual)
        call check_identity(checks, symengine, trim(label)//' [SymEngine]', residual)
        if (checks%failed /= 0) error stop 'identity failed or undecided'
    end subroutine prove
end program weighted_derivative_proof
