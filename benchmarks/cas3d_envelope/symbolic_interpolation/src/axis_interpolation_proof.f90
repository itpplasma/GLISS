program axis_interpolation_proof
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, positive, real_valued, &
        FACT_POSITIVE, FACT_NONZERO, FACT_REAL
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    implicit none

    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symengine
    type(suite_t) :: checks
    type(expr_t) :: h, a, b, c, s, j, m0, m1, m2, q0, q1, q2
    type(expr_t) :: nak, row, delta, closed_q0, closed_delta, closed_row
    type(expr_t) :: quadratic, exact_q0, exact_q1, exact_q2, axis_estimate
    type(expr_t) :: exact_delta, exact_curvature, imposed_curvature
    type(expr_t) :: first_half, exact_first_half, first_bias, ordinary_half
    type(expr_t) :: left_value, right_value, exact_midpoint, ordinary_bias
    type(expr_t) :: half, quarter
    character(1024) :: report

    report = 'strict-identities.json'
    if (command_argument_count() == 1) call get_command_argument(1, report)
    if (command_argument_count() > 1) error stop 'expected optional report path'
    call arena%init()
    call assumptions%init(arena)
    h = sym(arena, 'h')
    a = sym(arena, 'a')
    b = sym(arena, 'b')
    c = sym(arena, 'c')
    s = sym(arena, 's')
    j = sym(arena, 'j')
    m0 = sym(arena, 'M0')
    m1 = sym(arena, 'M1')
    m2 = sym(arena, 'M2')
    q0 = sym(arena, 'q0')
    q1 = sym(arena, 'q1')
    q2 = sym(arena, 'q2')
    call assume(assumptions, positive(h))
    call assume(assumptions, real_valued(a))
    call assume(assumptions, real_valued(b))
    call assume(assumptions, real_valued(c))
    call assume(assumptions, real_valued(s))
    call assume(assumptions, real_valued(j))
    call assume(assumptions, real_valued(m0))
    call assume(assumptions, real_valued(m1))
    call assume(assumptions, real_valued(m2))
    call assume(assumptions, real_valued(q0))
    call assume(assumptions, real_valued(q1))
    call assume(assumptions, real_valued(q2))
    if (.not. assumptions%has(h, FACT_POSITIVE)) error stop 'missing h>0'
    if (.not. assumptions%has(h, FACT_NONZERO)) error stop 'missing h nonzero'
    if (.not. assumptions%has(h, FACT_REAL)) error stop 'missing h real'
    native = make_native_engine(arena, assumptions)
    symengine = make_symengine_engine(arena)
    call suite_begin(checks, 'Uniform-grid axis and BOOZ interpolation identities')
    half = rat(arena, 1_int64, 2_int64)
    quarter = rat(arena, 1_int64, 4_int64)

    ! K=0 is the uniform not-a-knot row. R=0 is the first cubic spline
    ! row after division by h. These are assumptions of the implication;
    ! the checks below establish the exact elimination identities.
    nak = m0 - 2*m1 + m2
    delta = q2 - 2*q1 + q0
    row = m0 + 4*m1 + m2 - 6*delta/h**2
    call prove('spline row minus not-a-knot row', &
        row - nak - 6*(m1 - delta/h**2))
    call prove('explicit elimination formula for M1', &
        m1 - (row - nak)/6 - delta/h**2)
    closed_q0 = 2*q1 - q2
    closed_delta = q2 - 2*q1 + closed_q0
    closed_row = m0 + 4*m1 + m2 - 6*closed_delta/h**2
    call prove('linear axis closure makes data second difference zero', closed_delta)
    call prove('closed spline equations imply six times M1 is zero', &
        closed_row - nak - 6*m1)

    quadratic = a + b*s + c*s**2
    exact_q0 = a
    exact_q1 = a + b*h + c*h**2
    exact_q2 = a + 2*b*h + 4*c*h**2
    axis_estimate = 2*exact_q1 - exact_q2
    exact_delta = exact_q2 - 2*exact_q1 + exact_q0
    exact_curvature = diff(diff(quadratic, s), s)
    imposed_curvature = (exact_q2 - 2*exact_q1 + axis_estimate)/h**2
    call prove('quadratic true second derivative', exact_curvature - 2*c)
    call prove('quadratic true second difference', exact_delta - 2*c*h**2)
    call prove('true quadratic spline row recovers curvature', &
        exact_delta/h**2 - exact_curvature)
    call prove('quadratic satisfies not-a-knot row', 2*c - 2*(2*c) + 2*c)
    call prove('quadratic satisfies first spline row', &
        2*c + 4*(2*c) + 2*c - 6*exact_delta/h**2)
    call prove('axis quotient extrapolation bias is minus two c h squared', &
        axis_estimate - exact_q0 + 2*c*h**2)
    call prove('imposed first positive-node curvature is zero', imposed_curvature)
    call prove('imposed minus true curvature is minus two c', &
        imposed_curvature - exact_curvature + 2*c)

    ! BOOZ m=1 interpolation operates on q=C/sqrt(s). These are quotient
    ! errors; Cartesian coefficient errors also include sqrt(s_half).
    first_half = (3*exact_q1 - exact_q2)*half
    exact_first_half = a + b*h*half + c*h**2*quarter
    first_bias = first_half - exact_first_half
    call prove('BOOZ first half equals averaging manufactured axis and q1', &
        first_half - (axis_estimate + exact_q1)*half)
    call prove('BOOZ first-half quotient bias is minus three c h squared over four', &
        first_bias + 3*c*h**2*quarter)
    left_value = a + b*j*h + c*(j*h)**2
    right_value = a + b*(j + 1)*h + c*((j + 1)*h)**2
    ordinary_half = (left_value + right_value)*half
    exact_midpoint = a + b*(j + half)*h + c*((j + half)*h)**2
    ordinary_bias = ordinary_half - exact_midpoint
    call prove('ordinary half-grid quotient bias is plus c h squared over four', &
        ordinary_bias - c*h**2*quarter)
    call prove('first-half minus ordinary bias is minus c h squared', &
        first_bias - ordinary_bias + c*h**2)
    call prove('BOOZ first-half rule is exact for a linear quotient', &
        (3*(a + b*h) - (a + 2*b*h))*half - (a + b*h*half))

    if (checks%by_probe /= 0) error stop 'probe evidence is forbidden'
    if (checks%proved /= checks%total) error stop 'unproved identity'
    call suite_end(checks, trim(report))

contains

    subroutine prove(label, residual)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: residual

        call check_identity(checks, native, trim(label)//' [native, h>0]', residual)
        call check_identity(checks, symengine, trim(label)//' [SymEngine]', residual)
        if (checks%failed /= 0) error stop 'identity failed or undecided'
    end subroutine prove

end program axis_interpolation_proof
