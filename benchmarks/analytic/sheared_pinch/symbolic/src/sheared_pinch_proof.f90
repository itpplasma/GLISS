program sheared_pinch_proof
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, positive, nonzero, &
        real_valued
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, sqrt, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_subs, only: subs, subs_many
    implicit none

    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symengine
    type(suite_t) :: checks
    type(expr_t) :: r, t, tp, tpp, z, zp, zpp, m, k, mu0, gammap
    type(expr_t) :: f, fp, fpp, g, h, d, tangent, e, mr, k2, b2, ff, gg, jj
    type(expr_t) :: qr, qt, qz, divergence, current_theta, current_z, pressure_r
    type(expr_t) :: lorentz, raw, transformed, gfrom, hfrom, center, completed
    type(expr_t) :: ca, cb, cc, qmin, fn, gn, lagrangian, tangent_trace
    type(expr_t) :: mgg, mgh, mhh, drive, binormal, parallel, magnetic_completed
    type(expr_t) :: aa, shear, rs, zs, wall, pe, pressure, pressure_family_r
    type(expr_t) :: tau, length, flux_t, flux_p, flux_tr, flux_pr, jac, cov_i, cov_j
    type(expr_t) :: jdotb, flux_norm, nu, mu_tangent, coordinate_d, wronskian
    type(expr_t) :: field_eta
    type(expr_t) :: coordinate_binormal, coordinate_parallel, alpha, beta, delta
    type(expr_t) :: fn_family, gn_family, kres, exponent, x, half, critical_shear
    character(1024) :: report

    report = 'strict-identities.json'
    if (command_argument_count() == 1) call get_command_argument(1, report)
    if (command_argument_count() > 1) error stop 'expected optional report path'
    call arena%init()
    call assumptions%init(arena)
    r = real_symbol('r'); t = real_symbol('Btheta'); tp = real_symbol('Btheta_r')
    tpp = real_symbol('Btheta_rr'); z = real_symbol('Bz'); zp = real_symbol('Bz_r')
    zpp = real_symbol('Bz_rr'); m = real_symbol('m'); k = real_symbol('k')
    mu0 = real_symbol('mu0'); gammap = real_symbol('gamma_mu0_pressure')
    f = real_symbol('f'); fp = real_symbol('f_r'); fpp = real_symbol('f_rr')
    g = real_symbol('g'); h = real_symbol('h'); d = real_symbol('D')
    tangent = real_symbol('Ctan'); aa = real_symbol('A'); shear = real_symbol('C')
    rs = real_symbol('r_s'); zs = real_symbol('Z_s'); wall = real_symbol('a')
    pe = real_symbol('p_edge'); tau = real_symbol('two_pi'); length = real_symbol('L')
    exponent = real_symbol('nu_exponent'); x = real_symbol('distance_to_resonance')
    half = rat(arena, 1_int64, 2_int64)
    call assume(assumptions, positive(r))
    call assume(assumptions, positive(z))
    call assume(assumptions, positive(mu0))
    call assume(assumptions, positive(gammap))
    call assume(assumptions, positive(rs))
    call assume(assumptions, positive(zs))
    call assume(assumptions, positive(wall))
    call assume(assumptions, positive(pe))
    call assume(assumptions, positive(tau))
    call assume(assumptions, positive(length))
    call assume(assumptions, nonzero(m))
    call assume(assumptions, nonzero(aa))
    call assume(assumptions, nonzero(shear))
    mr = m/r
    k2 = mr**2 + k**2
    b2 = t**2 + z**2
    ff = mr*t + k*z
    gg = mr*z - k*t
    jj = mr*z + k*t
    call assume(assumptions, nonzero(ff))
    native = make_native_engine(arena, assumptions)
    symengine = make_symengine_engine(arena)
    call suite_begin(checks, &
        'Sheared pinch raw energy, static reduction and local threshold')

    ! Physical cylindrical curl of xi cross B=(i*Ctan,-Bz*f,Btheta*f).
    qr = mr*t*f + k*z*f
    qt = -k*tangent - dr(t*f)
    qz = (m*tangent - dr(r*z*f))/r
    divergence = fp + f/r - mr*g - k*h
    current_theta = -zp
    current_z = tp + t/r
    pressure_r = -t*current_z - z*zp ! This is mu0*p_r.
    lorentz = -current_theta*(qz*f - qr*h) - current_z*(qr*g - qt*f)
    raw = qr**2 + qt**2 + qz**2 + lorentz &
        + pressure_r*f*divergence + gammap*divergence**2
    call prove('physical radial curl amplitude', qr - ff*f)
    call prove('physical theta curl amplitude', qt + k*tangent + tp*f + t*fp)
    call prove('physical axial curl amplitude', qz - mr*tangent + z*fp &
        + (zp + z/r)*f)
    e = fp + f/r - d
    gfrom = (t*e + k*tangent)/ff
    hfrom = (z*e - mr*tangent)/ff
    call prove('invert tangential coordinates Ctan', z*gfrom - t*hfrom - tangent)
    call prove('invert tangential coordinates divergence', &
        replace2(divergence, g, h, gfrom, hfrom) - d)
    call prove('tangent coordinate determinant', -mr*t - k*z + ff)
    call prove('orthogonal wave-field identity', ff**2 + gg**2 - k2*b2)

    center = (gg*fp + jj*f/r)/k2
    ca = ff**2/k2
    cb = (k**2*z**2 - mr**2*t**2)/(r*k2)
    cc = ff**2 + (z**2 - t**2)/r**2 - 2*t*tp/r - jj**2/(r**2*k2)
    qmin = ca*fp**2 + 2*cb*f*fp + cc*f**2
    completed = k2*(tangent - center)**2 + gammap*d**2 + qmin
    transformed = replace2(raw, g, h, gfrom, hfrom)
    call prove('raw force-balanced Bernstein tangent completion', &
        transformed - completed)
    call prove('static Ctan stationarity', &
        diff(transformed, tangent) - 2*k2*(tangent - center))
    call prove('static divergence stationarity', diff(transformed, d) - 2*gammap*d)
    ! Instantiate the already-proved raw completion, retaining its exact DAG.
    call prove('completed static minimum away from F zero', &
        replace2(completed, tangent, d, center, 0*d) - qmin)
    mgg = k2*z**2 + gammap*mr**2
    mgh = -k2*z*t + gammap*mr*k
    mhh = k2*t**2 + gammap*k**2
    tangent_trace = mgg + mhh
    call prove('tangent quadratic trace', tangent_trace - k2*(b2 + gammap))
    call prove('tangent quadratic determinant', mgg*mhh - mgh**2 - gammap*k2*ff**2)
    call prove('static coefficient a', diff(qmin, fp) - 2*(ca*fp + cb*f))
    call prove('static coefficient c', diff(qmin, f) - 2*(cb*fp + cc*f))
    fn = r*ca
    gn = r*cc - dr(r*cb)
    lagrangian = r*qmin
    call prove('cross term integration identity', &
        lagrangian - fn*fp**2 - gn*f**2 - dr(r*cb*f**2))
    call prove('radial Euler Lagrange equation', &
        diff(lagrangian, f) - dr(diff(lagrangian, fp)) &
        - 2*(gn*f - dr(fn*fp)))

    ! This identity ties drive to physical-r force balance, not an imported value.
    drive = dr(pressure_r) + z*zpp + t*(tpp + 2*tp/r) &
        + current_z**2 + zp**2
    call prove('physical radius geometric drive', drive - 2*t*(tp + t/r)/r)
    call prove('linear Btheta drive is four A squared', family(drive) - 4*aa**2)
    binormal = ff*tangent - 2*t*z*f/r
    parallel = gg*tangent - b2*fp + (t**2 - z**2)*f/r
    magnetic_completed = qr**2 + (binormal**2 + parallel**2)/b2 &
        - 2*t*(tp + t/r)*f**2/r + gammap*divergence**2
    call prove('physical completed magnetic vector equals raw energy', &
        subs(raw - magnetic_completed, tangent, z*g - t*h))

    flux_t = tau*r*z; flux_p = length*t
    flux_tr = dr(flux_t); flux_pr = dr(flux_p)
    jac = tau*length*r; cov_i = length*z; cov_j = tau*r*t
    jdotb = -zp*t + current_z*z
    flux_norm = flux_t**2 + flux_p**2
    wronskian = flux_t*flux_pr - flux_tr*flux_p
    field_eta = -(flux_p*tau*m + flux_t*k*length)*tangent/jac
    coordinate_binormal = -(jac*field_eta - wronskian*f + jdotb*jac*f)/jac
    coordinate_parallel = (cov_j*(-k*length*tangent) - cov_i*(-tau*m*tangent) &
        - (flux_t*cov_i + flux_p*cov_j)*fp &
        - (cov_j*flux_pr + cov_i*flux_tr)*f - pressure_r*jac*f)/jac
    call prove('r chart normal curl and angular scale', &
        (flux_p*tau*m + flux_t*k*length)*f/jac - qr)
    call prove('r chart field3 contains axial radial derivative', &
        flux_tr - tau*(z + r*zp))
    call prove('r chart field10 includes axial shear current', &
        family(jdotb) - aa*(2*z - shear*r**2/z))
    call prove('r chart field3 for sheared profile', &
        family(flux_tr) - tau*(z + shear*r**2/z))
    call prove('r chart binormal completed projection', coordinate_binormal - binormal)
    call prove('r chart parallel completed projection', coordinate_parallel - parallel)
    nu = flux_norm*h/z
    ! Generic divergence uses mu; the assembly unknown is nu=mu-J*(P/T)*eta.
    mu_tangent = nu + jac*flux_p*(z*g - t*h)/flux_t
    coordinate_d = fp + f/r &
        + (flux_t*jac*(-tau*m*(z*g - t*h)) &
        - flux_p*jac*(-k*length*(z*g - t*h)) &
        + flux_p*(-tau*m*mu_tangent) &
        + flux_t*(-k*length*mu_tangent))/(jac*flux_norm)
    call prove('r chart eta nu divergence map', coordinate_d - divergence)
    call prove('physical tangent mass decomposition', &
        b2*(g**2 + h**2) - (z*g - t*h)**2 - (t*g + z*h)**2)

    pressure = pe + (2*aa**2 + shear)*(wall**2 - r**2)/(2*mu0)
    pressure_family_r = diff(pressure, r)
    call prove('sheared profile radial force balance', &
        mu0*pressure_family_r - family(pressure_r))
    call prove('squared axial field derivative', family(dr(z**2)) - 2*shear*r)
    fn_family = family(fn)
    gn_family = family(gn)
    kres = -m*aa/zs
    alpha = aa**2*shear**2*rs**5/(zs**2*(zs**2 + aa**2*rs**2))
    beta = -2*aa**2*(2*aa**2 + shear)*rs**3/(zs**2 + aa**2*rs**2)
    delta = 1 - 8*zs**2*(2*aa**2 + shear)/(shear**2*rs**2)
    call prove('simple resonant surface', at_res(family(ff)))
    call prove('simple resonance derivative', &
        at_res(dr_family(family(ff))) + m*aa*shear*rs/zs**2)
    call prove('radial leading coefficient vanishes at resonance', at_res(fn_family))
    call prove('radial leading first derivative vanishes', at_res(dr_family(fn_family)))
    call prove('simple resonance alpha', &
        half*at_res(dr_family(dr_family(fn_family))) - alpha)
    call prove('simple resonance beta', at_res(gn_family) - beta)
    call prove('local discriminant', 1 + 4*beta/alpha - delta)
    call prove('local indicial discriminant form', &
        exponent*(exponent + 1) - beta/alpha - ((2*exponent + 1)**2 - delta)/4)
    call prove('local Hardy square completion with boundary derivative', &
        alpha*x**2*fp**2 + beta*f**2 - alpha*(x*fp + half*f)**2 &
        - (beta + alpha/4)*f**2 + half*alpha*(f**2 + 2*x*f*fp))
    call prove('threshold coefficient is alpha Delta over four', &
        beta + alpha/4 - alpha*delta/4)
    call prove('physical pressure and shear form of threshold', &
        beta + alpha/4 - 2*aa**2*rs**2/(zs**2 + aa**2*rs**2) &
        *(shear**2*rs**3/(8*zs**2) + mu0*at_res(pressure_family_r)))
    call prove('stable fixture Delta', fixture(delta, rat(arena, -9_int64, 5_int64)) &
        - rat(arena, 41_int64, 81_int64))
    critical_shear = 4 - 4*sqrt(rat(arena, 2_int64, 1_int64))
    call prove('critical fixture Delta', fixture(delta, critical_shear))
    call prove('unstable fixture Delta', fixture(delta, -1 + 0*shear) + 7)
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

    function dr(expression) result(derivative)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: derivative
        derivative = diff(expression, r) + diff(expression, t)*tp &
            + diff(expression, tp)*tpp + diff(expression, z)*zp &
            + diff(expression, zp)*zpp + diff(expression, f)*fp &
            + diff(expression, fp)*fpp
    end function dr

    function family(expression) result(replaced)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: replaced, old(5), new(5)
        old(1) = t; old(2) = tp; old(3) = tpp; old(4) = zp; old(5) = zpp
        new(1) = aa*r; new(2) = aa; new(3) = 0*aa
        new(4) = shear*r/z; new(5) = shear/z - shear**2*r**2/z**3
        replaced = subs_many(expression, old, new)
    end function family

    function dr_family(expression) result(derivative)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: derivative
        derivative = diff(expression, r) + diff(expression, z)*shear*r/z
    end function dr_family

    function at_res(expression) result(replaced)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: replaced, old(3), new(3)
        old(1) = r; old(2) = z; old(3) = k
        new(1) = rs; new(2) = zs; new(3) = kres
        replaced = subs_many(expression, old, new)
    end function at_res

    function fixture(expression, shear_value) result(replaced)
        type(expr_t), intent(in) :: expression, shear_value
        type(expr_t) :: replaced, old(4), new(4)
        old(1) = aa; old(2) = zs; old(3) = rs; old(4) = shear
        new(1) = 1 + 0*aa; new(2) = 1 + 0*zs; new(3) = 1 + 0*rs
        new(4) = shear_value
        replaced = subs_many(expression, old, new)
    end function fixture

    function replace2(expression, old1, old2, new1, new2) result(replaced)
        type(expr_t), intent(in) :: expression, old1, old2, new1, new2
        type(expr_t) :: replaced, old(2), new(2)
        old(1) = old1; old(2) = old2
        new(1) = new1; new(2) = new2
        replaced = subs_many(expression, old, new)
    end function replace2

    subroutine prove(label, residual)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: residual
        call check_identity(checks, native, trim(label)//' [native]', residual)
        call check_identity(checks, symengine, trim(label)//' [SymEngine]', residual)
        if (checks%failed /= 0) error stop 'identity failed or undecided'
    end subroutine prove
end program sheared_pinch_proof
