program resonant_pinch_proof
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, positive, negative, &
        nonzero, real_valued, positive_integer, FACT_POSITIVE, FACT_NEGATIVE, &
        FACT_NONZERO
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_subs, only: subs, subs_many
    implicit none

    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symengine
    type(suite_t) :: checks
    type(expr_t) :: r, a, aa, z, mu0, pe, gamma, rho, lambda, m, k
    type(expr_t) :: f, fp, fpp, g, gp, h, hp, hh, hr, pr, prr, ss, sr
    type(expr_t) :: pp, pprime, dtarget, alpha, fc, pc, one
    type(expr_t) :: btheta, bsq, sphys, pressure, prphys, hphys, haxis, hwall
    type(expr_t) :: crossr_imag, crosstheta, crossz, qr_imag, qtheta, qz
    type(expr_t) :: divergence, kres, hperp, dres, qrres, qtres, qzres
    type(expr_t) :: lorentz, raw_energy, reduced_energy, mass, ds, qq, lagrangian
    type(expr_t) :: dfull, full_lagrangian, full_momentum, full_el_g, full_el_h
    type(expr_t) :: momentum, el_f, el_g, gstar, dstar, fpstar, ppstar
    type(expr_t) :: ccoef, ecoef, bcoef, tcoef, scalar_p, scalar_residual
    type(expr_t) :: sturm_residual, indicial_f, indicial_p, bound_ratio
    character(1024) :: report

    report = 'strict-identities.json'
    if (command_argument_count() == 1) call get_command_argument(1, report)
    if (command_argument_count() > 1) error stop 'expected optional report path'
    call arena%init()
    call assumptions%init(arena)
    r = real_symbol('r')
    a = real_symbol('a')
    aa = real_symbol('A')
    z = real_symbol('Z')
    mu0 = real_symbol('mu0')
    pe = real_symbol('p_edge')
    gamma = real_symbol('gamma')
    rho = real_symbol('rho')
    lambda = real_symbol('lambda')
    m = real_symbol('m')
    k = real_symbol('k')
    f = real_symbol('f')
    fp = real_symbol('f_r')
    fpp = real_symbol('f_rr')
    g = real_symbol('g')
    gp = real_symbol('g_r')
    h = real_symbol('h')
    hp = real_symbol('h_r')
    hh = real_symbol('H')
    hr = real_symbol('H_r')
    pr = real_symbol('p_r')
    prr = real_symbol('p_rr')
    ss = real_symbol('S')
    sr = real_symbol('S_r')
    pp = real_symbol('P')
    pprime = real_symbol('P_r')
    dtarget = real_symbol('D_target')
    alpha = real_symbol('alpha')
    fc = real_symbol('F_axis')
    pc = real_symbol('P_axis')
    one = rat(arena, 1_int64, 1_int64)
    call assume(assumptions, positive(r))
    call assume(assumptions, positive(a))
    call assume(assumptions, positive(mu0))
    call assume(assumptions, positive(pe))
    call assume(assumptions, positive(gamma))
    call assume(assumptions, positive(rho))
    call assume(assumptions, positive(hh))
    call assume(assumptions, positive(ss))
    call assume(assumptions, nonzero(z))
    call assume(assumptions, negative(lambda))
    call assume(assumptions, positive_integer(m))
    if (.not. assumptions%has(r, FACT_POSITIVE)) error stop 'missing r>0'
    if (.not. assumptions%has(hh, FACT_POSITIVE)) error stop 'missing H>0'
    if (.not. assumptions%has(rho, FACT_POSITIVE)) error stop 'missing rho>0'
    if (.not. assumptions%has(z, FACT_NONZERO)) error stop 'missing Z nonzero'
    if (.not. assumptions%has(lambda, FACT_NEGATIVE)) error stop 'missing lambda<0'
    native = make_native_engine(arena, assumptions)
    symengine = make_symengine_engine(arena)
    call suite_begin(checks, 'Resonant screw-pinch raw energy and radial ODE')

    btheta = aa*r
    bsq = btheta**2 + z**2
    sphys = bsq/z**2
    pressure = pe + aa**2*(a**2 - r**2)/mu0
    prphys = diff(pressure, r)
    hphys = bsq/mu0 + gamma*pressure
    haxis = z**2/mu0 + gamma*(pe + aa**2*a**2/mu0)
    hwall = (z**2 + aa**2*a**2)/mu0 + gamma*pe
    call prove('equilibrium curl B axial component', dr(r*btheta)/r - 2*aa)
    call prove('equilibrium radial force balance', prphys + 2*aa**2*r/mu0)
    call prove('physical S equals one plus squared transverse field ratio', &
        sphys - 1 - aa**2*r**2/z**2)

    ! Cylindrical curl of xi cross B. Imaginary radial amplitudes are
    ! represented by their real coefficient; all amplitudes are real.
    crossr_imag = z*g - btheta*h
    crosstheta = -z*f
    crossz = btheta*f
    qr_imag = m*crossz/r - k*crosstheta
    qtheta = -k*crossr_imag - dr(crossz)
    qz = (dr(r*crosstheta) + m*crossr_imag)/r
    divergence = dr(r*f)/r - m*g/r - k*h
    call prove('raw radial curl amplitude', qr_imag - (m*aa + k*z)*f)
    call prove('raw theta curl amplitude', &
        qtheta + k*(z*g - aa*r*h) + aa*(f + r*fp))
    call prove('raw axial curl amplitude', &
        qz + z*(fp + f/r) - m*(z*g - aa*r*h)/r)
    call prove('raw cylindrical divergence', divergence - fp - f/r + m*g/r + k*h)
    kres = -m*aa/z
    hperp = -aa*r*g/z
    dfull = subs(divergence, k, kres)
    call prove('resonance alone gives theta curl minus Btheta D', &
        subs(qtheta, k, kres) + btheta*dfull)
    call prove('resonance alone gives axial curl minus Z D', &
        subs(qz, k, kres) + z*dfull)
    call prove('raw energy reduction at resonance before perpendicular restriction', &
        (subs(qtheta, k, kres)**2 + subs(qz, k, kres)**2)/mu0 &
        + 2*aa*f*subs(qtheta, k, kres)/mu0 + prphys*f*dfull &
        + gamma*pressure*dfull**2 - hphys*dfull**2 - 2*prphys*f*dfull)
    full_momentum = hh*dfull + pr*f
    full_lagrangian = r*(hh*dfull**2 + 2*pr*f*dfull &
        - lambda*rho*(f**2 + g**2 + h**2))
    full_el_g = diff(full_lagrangian, g)
    full_el_h = diff(full_lagrangian, h)
    call prove('unrestricted g variation at resonance', &
        full_el_g + 2*(m*full_momentum + lambda*rho*r*g))
    call prove('unrestricted h variation at resonance', &
        full_el_h - 2*r*(m*aa*full_momentum/z - lambda*rho*h))
    call prove('nonzero eigenvalue stationarity forces perpendicular displacement', &
        aa*r*full_el_g + z*full_el_h &
        + 2*lambda*rho*r*(aa*r*g + z*h))
    dres = replace2(divergence, k, h, kres, hperp)
    qrres = replace2(qr_imag, k, h, kres, hperp)
    qtres = replace2(qtheta, k, h, kres, hperp)
    qzres = replace2(qz, k, h, kres, hperp)
    call prove('perpendicular displacement', btheta*g + z*hperp)
    call prove('resonant divergence', dres - fp - f/r + m*sphys*g/r)
    call prove('resonant radial curl is zero', qrres)
    call prove('resonant theta curl equals minus Btheta D', qtres + btheta*dres)
    call prove('resonant axial curl equals minus Z D', qzres + z*dres)

    ! -curl(B).(Q cross conjugate(xi))/mu0 has the radial-curl
    ! product qr_imag*g with positive sign before its overall minus.
    lorentz = -2*aa*(qrres*g - qtres*f)/mu0
    call prove('Lorentz cross term equals pressure-gradient cross term', &
        lorentz - prphys*f*dres)
    raw_energy = (qrres**2 + qtres**2 + qzres**2)/mu0 + lorentz &
        + prphys*f*dres + gamma*pressure*dres**2
    reduced_energy = hphys*dres**2 + 2*prphys*f*dres
    call prove('raw Bernstein energy reduces to H D squared plus two p_r f D', &
        raw_energy - reduced_energy)
    call prove('physical perpendicular mass reduction', &
        rho*(f**2 + g**2 + hperp**2) - rho*(f**2 + sphys*g**2))
    call prove('arbitrary divergence realization by g', &
        subs(dres, g, z**2*(r*fp + f - r*dtarget)/(m*bsq)) - dtarget)

    ds = fp + f/r - m*ss*g/r
    qq = hh*ds**2 + 2*pr*f*ds
    mass = rho*(f**2 + ss*g**2)
    call prove('exact energy square completion', &
        qq - hh*(ds + pr*f/hh)**2 + pr**2*f**2/hh)
    call prove('exact mass-weighted pointwise lower bound decomposition', &
        qq + pr**2*mass/(hh*rho) &
        - hh*(ds + pr*f/hh)**2 - pr**2*ss*g**2/hh)
    call prove('optimal-divergence negative energy value', &
        hh*(-pr*f/hh)**2 + 2*pr*f*(-pr*f/hh) + pr**2*f**2/hh)
    call prove('H equals its regular quadratic radial expression', &
        hphys - haxis - (1 - gamma)*aa**2*r**2/mu0)
    call prove('wall H value', subs(hphys, r, a) - hwall)
    bound_ratio = prphys**2/(rho*hphys)
    call prove('pointwise lower-bound ratio increases toward wall', &
        diff(bound_ratio, r) - 8*aa**4*r*haxis/(mu0**2*rho*hphys**2))
    call prove('wall lower-bound ratio', &
        subs(bound_ratio, r, a) - 4*aa**4*a**2/(mu0**2*rho*hwall))
    call prove('gamma one H is constant', subs(hphys, gamma, one) &
        - (z**2 + mu0*pe + aa**2*a**2)/mu0)

    ! The Euler-Lagrange derivatives are taken from the unreduced
    ! f,g functional. Coefficient derivatives remain independent jets.
    lagrangian = r*(qq - lambda*mass)
    momentum = hh*ds + pr*f
    el_g = diff(lagrangian, g)
    el_f = dr(diff(lagrangian, fp)) - diff(lagrangian, f)
    call prove('canonical radial momentum from independent variation', &
        diff(lagrangian, fp) - 2*r*momentum)
    call prove('algebraic g Euler-Lagrange equation', &
        el_g + 2*ss*(m*momentum + lambda*rho*r*g))
    call prove('differential f Euler-Lagrange equation', &
        el_f - 2*r*(dr(momentum) - pr*ds + lambda*rho*f))
    gstar = -m*pp/(lambda*rho*r)
    dstar = (pp - pr*f)/hh
    fpstar = dstar - f/r + m*ss*gstar/r
    ppstar = pr*dstar - lambda*rho*f
    call prove('first-order f equation derived by eliminating g', &
        fpstar - (1/hh - m**2*ss/(lambda*rho*r**2))*pp &
        + (pr/hh + 1/r)*f)
    call prove('first-order P equation derived by variation', &
        ppstar - (pr/hh)*pp + (pr**2/hh + lambda*rho)*f)
    call prove('first-order reconstruction satisfies momentum definition', &
        replace2(momentum, fp, g, fpstar, gstar) - pp)
    call prove('first-order reconstruction satisfies g variation', &
        replace2(el_g, fp, g, fpstar, gstar))
    call prove('first-order P equation satisfies f variation', &
        ppstar - pr*dstar + lambda*rho*f)

    ! A scalar second-order equation with a positive leading coefficient
    ! on lambda<0. No division by p_r squared/H + lambda*rho is needed.
    ccoef = 1/hh - m**2*ss/(lambda*rho*r**2)
    ecoef = pr/hh + 1/r
    bcoef = pr**2/hh + lambda*rho
    tcoef = 1/ccoef
    scalar_p = tcoef*(fp + ecoef*f)
    scalar_residual = dr(scalar_p) - pr*scalar_p/hh + bcoef*f
    sturm_residual = dr(r*tcoef*fp) &
        + (dr(r*tcoef*ecoef) - r*tcoef*ecoef**2 + r*bcoef)*f
    call prove('scalar ODE momentum reconstruction', ccoef*scalar_p - fp - ecoef*f)
    call prove('scalar radial equation in divergence form', &
        sturm_residual - r*scalar_residual)

    ! Frobenius leading-coefficient algebra; the full existence and
    ! convergence of a Frobenius series is not asserted by zero tests.
    indicial_f = (alpha + 1)*fc + m**2*pc/(lambda*rho)
    indicial_p = (alpha + 1)*pc + lambda*rho*fc
    call prove('axis indicial elimination', &
        (alpha + 1)*indicial_f - m**2*indicial_p/(lambda*rho) &
        - ((alpha + 1)**2 - m**2)*fc)
    call prove('axis indicial factorization', &
        (alpha + 1)**2 - m**2 - (alpha + 1 - m)*(alpha + 1 + m))
    call prove('regular axis f leading coefficient', &
        replace2(indicial_f, alpha, pc, m - 1, -lambda*rho*fc/m))
    call prove('regular axis P leading coefficient', &
        replace2(indicial_p, alpha, pc, m - 1, -lambda*rho*fc/m))
    call prove('regular axis physical f g tie', -m*(-lambda*rho*fc/m)/(lambda*rho) - fc)

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

    function replace2(expression, old1, old2, new1, new2) result(replaced)
        type(expr_t), intent(in) :: expression, old1, old2, new1, new2
        type(expr_t) :: replaced, old(2), new(2)
        old(1) = old1
        old(2) = old2
        new(1) = new1
        new(2) = new2
        replaced = subs_many(expression, old, new)
    end function replace2

    function dr(expression) result(derivative)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: derivative
        derivative = diff(expression, r) + diff(expression, f)*fp &
            + diff(expression, fp)*fpp + diff(expression, g)*gp &
            + diff(expression, h)*hp + diff(expression, hh)*hr &
            + diff(expression, pr)*prr + diff(expression, ss)*sr &
            + diff(expression, pp)*pprime
    end function dr

    subroutine prove(label, residual)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: residual
        call check_identity(checks, native, trim(label)//' [native]', residual)
        call check_identity(checks, symengine, trim(label)//' [SymEngine]', residual)
        if (checks%failed /= 0) error stop 'identity failed or undecided'
    end subroutine prove
end program resonant_pinch_proof
