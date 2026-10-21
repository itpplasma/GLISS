program prove_free_boundary
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, positive
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    implicit none
    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: engine
    type(suite_t) :: checks
    type(expr_t) :: r, a, b, t, c, sn, w, e, r0, rw, rt, zw, zt
    type(expr_t) :: phi_edge, q, psi0, dphi_dpsi, dchi_dpsi, ds_dpsi
    type(expr_t) :: potential, response, expected, radial_field, axial, k
    type(expr_t) :: volume, vacuum, total, wall_ratio
    character(1024) :: report
    integer :: m

    report = 'free-boundary-identities.json'
    if (command_argument_count() == 1) call get_command_argument(1, report)
    call arena%init()
    call assumptions%init(arena)
    r = sym(arena, 'R')
    a = sym(arena, 'a')
    b = sym(arena, 'b')
    t = sym(arena, 't')
    w = sym(arena, 'sqrt_psi_n')
    e = sym(arena, 'elongation')
    r0 = sym(arena, 'R0')
    phi_edge = sym(arena, 'Phi_edge')
    q = sym(arena, 'q')
    psi0 = sym(arena, 'psi0')
    radial_field = sym(arena, 'A')
    axial = sym(arena, 'Bz')
    k = sym(arena, 'k')
    wall_ratio = sym(arena, 'wall_radius_squared_over_a_squared')
    call assume(assumptions, positive(r))
    call assume(assumptions, positive(a))
    call assume(assumptions, positive(b))
    call assume(assumptions, positive(w))
    call assume(assumptions, positive(q))
    call assume(assumptions, positive(phi_edge))
    call assume(assumptions, positive(psi0))
    engine = make_native_engine(arena, assumptions)
    call suite_begin(checks, 'Exact free-boundary edge and cylindrical wall identities')

    ! Rational parametrization of the unit circle proves the contour Jacobian
    ! without a numerical trigonometric zero test. R^2=R0^2+2*a*R0*w*c,
    ! Z=e*a*R0*w*sn/R; these partials follow by implicit differentiation.
    c = (1 - t**2)/(1 + t**2)
    sn = 2*t/(1 + t**2)
    call prove('unit circle', c**2 + sn**2 - 1)
    rw = a*r0*c/r
    rt = -a*r0*w*sn/r
    zw = e*a*r0*sn/r - e*a*r0*w*sn*rw/r**2
    zt = e*a*r0*w*c/r - e*a*r0*w*sn*rt/r**2
    call prove('Solovev d(R,Z)/d(psi_n,angle)', &
        (rw*zt - rt*zw)/(2*w) - e*a**2*r0**2/(2*r**2))
    ! dPhi/dpsi_n = F*integral[J/R]dangle = 2*pi*psi0*q.
    ! chi=-2*pi*psi0*psi_n and s=Phi/Phi_edge. Cancel the common 2*pi
    ! scale (psi0 here includes it) before forming the native signed ratio.
    dphi_dpsi = psi0*q
    dchi_dpsi = -psi0
    ds_dpsi = dphi_dpsi/phi_edge
    call prove('true-edge toroidal slope', dphi_dpsi/ds_dpsi - phi_edge)
    call prove('true-edge poloidal slope', dchi_dpsi/ds_dpsi + phi_edge/q)
    call prove('signed chi prime over Phi prime', dchi_dpsi/dphi_dpsi + 1/q)

    ! Current-free annulus, k=0: potential=(r^m+b^(2m)/r^m)*cos(m*theta).
    ! The ideal wall has zero radial derivative. The boundary energy is
    ! Wv=pi*L*a^2*g^2/(2*mu0*m)*(b^(2m)+a^(2m))/(b^(2m)-a^(2m)).
    do m = 1, 3
        potential = r**m + b**(2*m)/r**m
        call prove('annular Laplace radial equation', &
            diff(r*diff(potential, r), r)/r - m**2*potential/r**2)
        call prove('ideal conducting wall Neumann condition', &
            m*b**(m - 1) - m*b**(2*m)/b**(m + 1))
        response = (a**m + b**(2*m)/a**m) &
            /(-m*a**(m - 1) + m*b**(2*m)/a**(m + 1))
        expected = a*(b**(2*m) + a**(2*m))/(m*(b**(2*m) - a**(2*m)))
        call prove('absolute cylindrical Neumann to Dirichlet response', &
            response - expected)
        call prove('wall energy decreases with wall radius', &
            diff(expected, b) + 4*a**(2*m + 1)*b**(2*m - 1) &
            /(b**(2*m) - a**(2*m))**2)
    end do

    ! Rigid helical m=1 translation in a uniform-current screw pinch:
    ! Btheta=A*r, Bz constant, mu0*p=A^2*(a^2-r^2), div(xi)=0.
    ! Q=(A-k*Bz)*(sin(k*z),-cos(k*z),0), curl(B)=2*A*ez.
    volume = (radial_field - k*axial)**2 &
        - 2*radial_field*(radial_field - k*axial)
    call prove('rigid kink raw plasma energy', &
        volume - (k**2*axial**2 - radial_field**2))
    ! Leading long-wavelength vacuum response; q=k*Bz/A, 0<q<1.
    ! Normalize W by pi*L*a^2*xi0^2*A^2/(2*mu0). Corrections of order
    ! (k*a)^2 in the annular vacuum response are NOT included in this limit.
    vacuum = (1 - q)**2*(wall_ratio + 1)/(wall_ratio - 1)
    total = q**2 - 1 + vacuum
    call prove('long wavelength kink wall threshold', &
        total - 2*(1 - q)*(1 - q*wall_ratio)/(wall_ratio - 1))
    call prove('no wall Kruskal Shafranov limit', &
        q**2 - 1 + (1 - q)**2 - 2*q*(q - 1))
    if (checks%by_probe /= 0) error stop 'probe evidence is forbidden'
    if (checks%proved /= checks%total) error stop 'unproved identity'
    call suite_end(checks, trim(report))
contains
    subroutine prove(label, residual)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: residual
        call check_identity(checks, engine, label, residual)
        if (checks%failed /= 0) error stop 'identity failed or undecided'
    end subroutine prove
end program prove_free_boundary
