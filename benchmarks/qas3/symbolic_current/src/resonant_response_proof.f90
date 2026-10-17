program resonant_response_proof
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, positive, &
        nonnegative, nonzero, real_valued
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_subs, only: subs
    implicit none
    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symengine
    type(suite_t) :: checks
    type(expr_t) :: d, pol, tor, delta, width, scale, kappa, m, n, spread, nfp
    type(expr_t) :: aden, bden, rt, rg, zero, quarter, qasdelta, normalized
    type(expr_t) :: physical_width_squared, rtscaled, rgscaled, path
    character(1024) :: report

    report = 'strict-identities.json'
    if (command_argument_count() == 1) call get_command_argument(1, report)
    call arena%init()
    call assumptions%init(arena)
    d = real_symbol('D')
    pol = real_symbol('P')
    tor = real_symbol('T')
    delta = real_symbol('delta')
    width = real_symbol('W')
    scale = real_symbol('scale')
    kappa = real_symbol('kappa')
    m = real_symbol('m')
    n = real_symbol('n')
    spread = real_symbol('iota_spread')
    nfp = real_symbol('Nfp')
    path = real_symbol('path_slope')
    call assume(assumptions, positive(delta))
    call assume(assumptions, nonnegative(width))
    call assume(assumptions, nonnegative(spread))
    call assume(assumptions, positive(scale))
    call assume(assumptions, nonzero(kappa))
    call assume(assumptions, positive(nfp))
    zero = rat(arena, 0_int64, 1_int64)
    quarter = rat(arena, 1_int64, 4_int64)
    qasdelta = rat(arena, 1_int64, 25_int64)
    aden = d**2 + delta*(d + pol)**2
    bden = d**2 + width**2
    ! This states the common rational-function domain. The separate
    ! source zero-width/gauge branches are recorded in the derivation.
    call assume(assumptions, positive(aden))
    call assume(assumptions, positive(bden))
    rt = d/aden
    rg = d/bden
    native = make_native_engine(arena, assumptions)
    symengine = make_symengine_engine(arena)
    call suite_begin(checks, 'Pinned TERPSICHORE and GLISS resonant responses')

    call prove('TERPSICHORE inverse-equation residual', &
        d*rt - 1 + delta*(d + pol)**2/aden)
    call prove('GLISS inverse-equation residual', d*rg - 1 + width**2/bden)
    call prove('difference of the two rational responses', &
        rt - rg - d*(width**2 - delta*(d + pol)**2)/(aden*bden))
    call prove('equality numerator after common-denominator multiplication', &
        (rt - rg)*aden*bden - d*(width**2 - delta*(d + pol)**2))
    call prove('TERPSICHORE error relative to exact inverse for D nonzero', &
        rt - 1/d + delta*(d + pol)**2/(d*aden))
    call prove('GLISS error relative to exact inverse for D nonzero', &
        rg - 1/d + width**2/(d*bden))
    call prove('TERPSICHORE exact inverse at D plus P zero, D nonzero', &
        subs(rt, pol, -d) - 1/d)
    call prove('GLISS zero-width exact inverse for D nonzero', &
        subs(rg, width, zero) - 1/d)
    call prove('TERPSICHORE zero-delta expression for fixed nonzero D', &
        subs(rt, delta, zero) - 1/d)
    call prove('TERPSICHORE at D zero when P nonzero', subs(rt, d, zero))
    call prove('GLISS at D zero when W positive', subs(rg, d, zero))
    call prove('TERPSICHORE derivative at D zero when P nonzero', &
        subs(diff(rt, d), d, zero) - 1/(delta*pol**2))
    call prove('GLISS derivative at D zero when W positive', &
        subs(diff(rg, d), d, zero) - 1/width**2)
    call prove('TERPSICHORE P zero remains a pole at D zero', &
        subs(rt, pol, zero) - 1/((1 + delta)*d))
    call prove('TERPSICHORE product residual at P zero', &
        d*subs(rt, pol, zero) - 1 + delta/(1 + delta))
    call prove('TERPSICHORE joint approach P equals c D', &
        subs(rt, pol, path*d) - 1/(d*(1 + delta*(1 + path)**2)))
    call prove('public QAS delta point zero four as exact rational', &
        subs(rt, delta, qasdelta) - 25*d/(25*d**2 + (d + pol)**2))

    normalized = ((d/scale)/((d/scale)**2 + (width/scale)**2))/scale
    call prove('GLISS normalized response returns common-chart physical response', &
        normalized - rg)
    physical_width_squared = quarter*m**2*tor**2*spread**2
    call prove('GLISS normalized width squares to physical W squared', &
        scale**2*(quarter*m**2*(tor/scale)**2*spread**2) &
        - physical_width_squared)
    call prove('common-chart harmonic denominator plus P', &
        (m*pol - n*tor) + pol - ((m + 1)*pol - n*tor))
    call prove('width equality condition in common-chart flux variables', &
        subs(width**2 - delta*(d + pol)**2, d, m*pol - n*tor) &
        - (width**2 - delta*((m + 1)*pol - n*tor)**2))

    rtscaled = kappa*d/((kappa*d)**2 + delta*(kappa*d + kappa*pol)**2)
    rgscaled = kappa*d/((kappa*d)**2 + kappa**2*width**2)
    call prove('common flux scaling of TERPSICHORE response', rtscaled - rt/kappa)
    call prove('common flux scaling of GLISS response', rgscaled - rg/kappa)
    call prove('TERPSICHORE dimensionless attenuation survives flux scaling', &
        kappa*d*rtscaled - d*rt)
    call prove('GLISS dimensionless attenuation survives flux scaling', &
        kappa*d*rgscaled - d*rg)
    call prove('flip-poloidal TERPSICHORE denominator mapping', &
        m*(-nfp*pol) - (nfp*n)*(-tor) + nfp*(m*pol - n*tor))
    call prove('flip-poloidal denominator-plus-P mapping', &
        m*(-nfp*pol) - (nfp*n)*(-tor) - nfp*pol &
        + nfp*((m*pol - n*tor) + pol))

    if (checks%by_probe /= 0) error stop 'probe evidence is forbidden'
    if (checks%proved /= checks%total) error stop 'unproved identity'
    call suite_end(checks, trim(report))
contains
    function real_symbol(name) result(expression)
        character(*), intent(in) :: name
        type(expr_t) :: expression
        expression = sym(arena, name)
        call assume(assumptions, real_valued(expression))
    end function real_symbol

    subroutine prove(label, residual)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: residual
        call check_identity(checks, native, trim(label)//' [native]', residual)
        call check_identity(checks, symengine, trim(label)//' [SymEngine]', residual)
        if (checks%failed /= 0) error stop 'identity failed or undecided'
    end subroutine prove
end program resonant_response_proof
