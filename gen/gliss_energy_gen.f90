program gliss_energy_gen
    ! Source-authored straight-field-coordinate energy principle.
    ! Nuehrenberg, Phys. Plasmas 6,137 (1999), DOI 10.1063/1.873268,
    ! Appendix A1, A5--A12; PROVENANCE.md declares the GLISS sign map.
    ! No Mathematica or Wolfram source/output is an input to this generator.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, rat, sqrt, abs, &
        operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    implicit none

    character(*), parameter :: default_revision = &
        'f5ef0572a524378abf3ae83669ee856684691c0a'
    character(len=32), parameter :: names(27) = [character(len=32) :: &
        'flux_t_slope', 'flux_p_slope', 'flux_t_curve', 'flux_p_curve', &
        'current_i', 'current_j', 'signed_sqrtg', 'bmag', 'grad_s2', &
        'j_dot_b', 'pressure_slope', 'signed_sigma_tilde', 'beta_tilde', &
        'drive_a', 'gamma_pressure_pa', 'xi_s', 'xi_s_radial', 'xi_s_theta', &
        'xi_s_zeta', 'eta_theta', 'eta_zeta', 'sqrtg_xi_s_radial', &
        'sqrtg_eta_theta', 'sqrtg_eta_zeta', 'mu_theta', 'mu_zeta', &
        'vacuum_permeability']
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(suite_t) :: checks
    type(expr_t) :: x(27), perpendicular_inputs(7), components(3), coefficients(18)
    type(expr_t) :: field_xi, field_eta, flux_norm, flux_wronskian
    type(expr_t) :: curl_s, curl_theta, curl_zeta
    type(expr_t) :: two_energy, divergence, three_energy, perpendicular_energy
    type(expr_t) :: normal, binormal, perpendicular_matrix(4), linear
    type(str_t) :: sources(16), symbols(16)
    character(len=1024) :: directory, evidence
    character(len=128) :: revision
    character(len=80) :: label
    integer :: i, row, column, count, unit, status

    directory = 'src/generated'
    evidence = ''
    revision = default_revision
    if (command_argument_count() >= 1) call get_command_argument(1, directory)
    if (command_argument_count() >= 2) call get_command_argument(2, evidence)
    if (command_argument_count() >= 3) call get_command_argument(3, revision)
    if (command_argument_count() > 3) error stop 'expected outputdir [proofJSON] [revision]'
    call arena%init()
    engine = make_symengine_engine(arena)
    call suite_begin(checks, 'GLISS local energy exact identities')
    do i = 1, size(x)
        x(i) = sym(arena, trim(names(i)))
    end do
    count = 0

    ! Let e_i=partial_i r, J=e_s.dot.(e_theta cross e_zeta), T=FT', P=FP',
    ! Q=|grad s|^2 and B=|Bvec|. GLISS uses the signed coordinate convention
    ! Bvec=(P*e_theta+T*e_zeta)/J=beta*grad s+Jc*grad theta+I*grad zeta.
    ! Jc=current_j and I=current_i are covariant B components, not SI current.
    ! At fixed J the paper's flux derivatives are -T,-P; eta is unchanged
    ! and its mu equals -GLISS mu. Define x=xi^s, eta=T*xi^theta-P*xi^zeta.
    ! The covariant components of a=xi cross Bvec are (eta,-T*x,P*x).
    ! Curvilinear curl therefore has contravariant components
    ! u=(P*x_theta+T*x_zeta)/J,
    ! v=(eta_zeta-P*x_s-P'*x)/J, w=(-eta_theta-T*x_s-T'*x)/J.
    ! With jhat=curl Bvec=mu0*j_SI and S=jhat.dot.Bvec, the magnetic vector
    ! is C=curl a+(jhat cross grad s)*x/Q. pressure_slope=mu0*d(p_Pa)/ds.
    ! For b=grad s cross Bvec, b_i=(-sigma*B,-T*Q,P*Q), where
    ! sigma=(I*g_s_theta-Jc*g_s_zeta)/(J*B). Triple products give
    ! (jhat cross grad s).dot.b=-Q*S and
    ! (jhat cross grad s).dot.Bvec=-pressure_slope*Q.
    ! The latter uses force balance jhat cross Bvec=pressure_slope*grad s.
    ! Thus the normal, binormal and field-aligned orthonormal projections are
    ! u/sqrt(Q), (-sigma*B*u-T*Q*v+P*Q*w-S*x)/(sqrt(Q)*B),
    ! and (beta*u+Jc*v+I*w-pressure_slope*x)/B, respectively.
    ! The strict identities below check their expansion independently.
    ! Keep T*P'-T'*P as a Wronskian: T^2*(P/T)' would exclude T=0.
    ! J is signed in differential operators and abs(J) only in integrals.
    ! Both nonzero-J branches are supported; their derivatives cannot cross J=0.
    ! Q>0 and B>0 are required. Divergence also needs T^2+P^2>0.
    ! There is no division by a Fourier resonance m*P-n*T. Axis regularity,
    ! consistent geometry, force balance and the drive_a field remain caller
    ! contracts. Independent input derivatives describe an algebraic extension
    ! until the caller supplies the equilibrium and geometry chain rule.
    field_xi = (x(2)*x(18) + x(1)*x(19))/x(7)
    field_eta = (x(2)*x(20) + x(1)*x(21))/x(7)
    flux_norm = x(1)**2 + x(2)**2
    flux_wronskian = x(1)*x(4) - x(3)*x(2)
    components(1) = field_xi/sqrt(x(9))
    components(2) = -sqrt(x(9))/(x(8)*x(7)) &
        * (x(7)*field_eta - flux_wronskian*x(16) &
        + x(10)*x(7)*x(16)/x(9) + x(12)*x(8)*x(7)*field_xi/x(9))
    components(3) = (x(6)*x(21) - x(5)*x(20) &
        - (x(1)*x(5) + x(2)*x(6))*x(17) &
        - (x(6)*x(4) + x(5)*x(3))*x(16) - x(11)*x(7)*x(16) &
        + x(13)*x(7)*field_xi)/(x(8)*x(7))
    curl_s = (x(2)*x(18) + x(1)*x(19))/x(7)
    curl_theta = (x(21) - x(2)*x(17) - x(4)*x(16))/x(7)
    curl_zeta = (-x(20) - x(1)*x(17) - x(3)*x(16))/x(7)
    call check_identity(checks, engine, 'normal coordinate curl projection', &
        sqrt(x(9))*components(1) - curl_s)
    call check_identity(checks, engine, 'binormal coordinate curl projection', &
        sqrt(x(9))*x(8)*components(2) + x(12)*x(8)*curl_s &
        + x(1)*x(9)*curl_theta - x(2)*x(9)*curl_zeta + x(10)*x(16))
    call check_identity(checks, engine, 'parallel coordinate curl projection', &
        x(8)*components(3) - x(13)*curl_s - x(6)*curl_theta &
        - x(5)*curl_zeta + x(11)*x(16))
    two_energy = (components(1)**2 + components(2)**2 + components(3)**2 &
        - x(14)*x(16)**2)*abs(x(7))

    ! Divergence of a vector density. The three J-weighted derivatives are
    ! independent inputs so radial/angle product rules stay in the caller.
    divergence = x(22)/x(7) + (x(1)*x(23) - x(2)*x(24) &
        + x(2)*x(25) + x(1)*x(26))/(x(7)*flux_norm)
    three_energy = two_energy/x(27) + x(15)*divergence**2*abs(x(7))

    do row = 1, 3
        linear = 0*two_energy
        do column = 1, 6
            i = (column - 1)*3 + row
            coefficients(i) = diff(components(row), x(15 + column))
            linear = linear + coefficients(i)*x(15 + column)
        end do
        write (label, '(a,i0)') 'component linear reconstruction ', row
        call check_identity(checks, engine, trim(label), components(row) - linear)
    end do
    call check_identity(checks, engine, 'gamma pressure derivative', &
        diff(three_energy, x(15)) - divergence**2*abs(x(7)))
    call check_identity(checks, engine, 'zero gamma SI energy', &
        three_energy - x(15)*diff(three_energy, x(15)) - two_energy/x(27))

    call leaf('generated_bending', [components(1)], x([1, 2, 7, 9, 18, 19]))
    call leaf('generated_shear', [components(2)], &
        x([1, 2, 3, 4, 7, 8, 9, 10, 12, 16, 18, 19, 20, 21]))
    call leaf('generated_compression', [components(3)], &
        x([1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 16, 17, 18, 19, 20, 21]))
    call leaf('generated_components', components, [x(:13), x(16:21)])
    call leaf('generated_component_coefficients', coefficients, x(:13), rows=3)
    call leaf('generated_two_energy', [two_energy], [x(:14), x(16:21)])
    call actions('generated_two', two_energy, [x(:14), x(16:21)], 20)
    call leaf('generated_divergence', [divergence], [x([1, 2, 7]), x(22:26)])
    call leaf('generated_three_energy', [three_energy], x)
    call actions('generated_three', three_energy, x, 26)

    ! The kinetic projection retains the normal and binormal displacement
    ! squares of Appendix A8, and removes only the field-parallel square.
    perpendicular_inputs = [x(7), x(8), x(9), x(12), &
        sym(arena, 'density_kg_m3'), x(16), sym(arena, 'eta')]
    normal = perpendicular_inputs(6)/sqrt(perpendicular_inputs(3))
    binormal = perpendicular_inputs(4)*normal &
        + sqrt(perpendicular_inputs(3))*perpendicular_inputs(7)/perpendicular_inputs(2)
    perpendicular_energy = rat(arena, 1_int64, 2_int64) &
        * perpendicular_inputs(5)*abs(perpendicular_inputs(1))*(normal**2 + binormal**2)
    do column = 1, 2
        do row = 1, 2
            i = (column - 1)*2 + row
            perpendicular_matrix(i) = diff(diff(perpendicular_energy, &
                perpendicular_inputs(5 + row)), perpendicular_inputs(5 + column))
        end do
    end do
    call check_identity(checks, engine, 'perpendicular matrix symmetry', &
        perpendicular_matrix(2) - perpendicular_matrix(3))
    call leaf('generated_perpendicular_matrix', perpendicular_matrix, &
        perpendicular_inputs(:5), rows=2)
    call leaf('generated_perpendicular_energy', [perpendicular_energy], perpendicular_inputs)
    call actions('generated_perpendicular', perpendicular_energy, perpendicular_inputs, 7)

    if (len_trim(evidence) > 0) then
        call suite_end(checks, trim(evidence))
    else
        call suite_end(checks)
    end if
    open (newunit=unit, file=trim(directory)//'/gliss_energy_generated.f90', &
        status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'cannot write generated energy module'
    write (unit, '(a)') '! Generated by gen/gliss_energy_gen.f90; do not edit.'
    write (unit, '(a)') '! Domain: J /= 0, B > 0, grad_s2 > 0, FT^2+FP^2 > 0; rho > 0.'
    write (unit, '(a)') 'module gliss_energy_generated'
    write (unit, '(a)') '    implicit none'
    write (unit, '(a)') '    private'
    do i = 1, count
        write (unit, '(a)') '    public :: '//chars(symbols(i))
    end do
    write (unit, '(a)') 'contains'
    do i = 1, count
        write (unit, '(a)') chars(sources(i))
    end do
    write (unit, '(a)') 'end module gliss_energy_generated'
    close (unit)
    write (output_unit, '(a)') 'Wrote '//trim(directory)//'/gliss_energy_generated.f90'

contains

    subroutine actions(prefix, value, variables, active)
        character(*), intent(in) :: prefix
        type(expr_t), intent(in) :: value, variables(:)
        integer, intent(in) :: active
        type(expr_t) :: partials(active), tangent(active), adjoint, forward, reverse(active)
        integer :: entry

        adjoint = sym(arena, 'cotangent')
        do entry = 1, active
            partials(entry) = diff(value, variables(entry))
            tangent(entry) = sym(arena, 'd_'//chars(variables(entry)%name()))
            reverse(entry) = adjoint*partials(entry)
        end do
        forward = partials(1)*tangent(1)
        do entry = 2, active
            forward = forward + partials(entry)*tangent(entry)
        end do
        call leaf(prefix//'_jvp', [value, forward], [variables, tangent])
        call leaf(prefix//'_vjp', [value, reverse], [variables, adjoint])
    end subroutine actions

    subroutine leaf(name, roots, variables, rows)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: roots(:), variables(:)
        integer, optional, intent(in) :: rows
        type(kernel_spec_t) :: spec
        type(expr_t) :: simplified_roots(size(roots))
        type(engine_result_t) :: simplified
        character(:), allocatable :: message
        character(len=80) :: local_label
        integer :: entry, nrows
        logical :: good

        count = count + 1
        symbols(count) = str(name)
        spec%name = str(name)
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%temp_prefix = str('t')
        spec%generator = str('gen/gliss_energy_gen.f90')
        spec%generator_revision = str(trim(revision))
        spec%regenerate_command = str('gliss_energy_gen src/generated')
        allocate (spec%args(size(variables)), spec%outputs(1), spec%output_shapes(1), &
            spec%output_references(size(roots)))
        do entry = 1, size(variables)
            spec%args(entry) = variables(entry)%name()
        end do
        spec%outputs(1) = str('values')
        nrows = 0
        if (present(rows)) nrows = rows
        if (nrows > 0) then
            write (local_label, '(a,i0,a,i0,a)') '(', nrows, ',', size(roots)/nrows, ')'
        else
            write (local_label, '(a,i0,a)') '(', size(roots), ')'
        end if
        spec%output_shapes(1) = str(trim(local_label))
        do entry = 1, size(roots)
            if (nrows > 0) then
                write (local_label, '(a,i0,a,i0,a)') 'values(', &
                    mod(entry - 1, nrows) + 1, ',', (entry - 1)/nrows + 1, ')'
            else
                write (local_label, '(a,i0,a)') 'values(', entry, ')'
            end if
            spec%output_references(entry) = str(trim(local_label))
            simplified = engine%simplify(roots(entry))
            if (.not. simplified%ok) error stop 'native simplification failed'
            simplified_roots(entry) = simplified%value
            write (local_label, '(a,a,i0)') name, ' root ', entry
            call check_identity(checks, engine, trim(local_label), &
                roots(entry) - simplified_roots(entry))
        end do
        sources(count) = emit_kernel(simplified_roots, spec, ok=good, message=message)
        if (.not. good) then
            write (output_unit, '(a)') message
            error stop 'energy emission failed'
        end if
    end subroutine leaf

end program gliss_energy_gen
