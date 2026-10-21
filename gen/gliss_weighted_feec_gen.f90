program gliss_weighted_feec_gen
    ! Local jets of the stored weighted radial coefficient, not a claim of
    ! stability or commuting projection for the assembled displacement space.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use fortsym_arena, only: arena_t
    use fortsym_assume, only: assumption_context_t, assume, positive, real_valued
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, exp, log, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: str_t, str, chars
    implicit none

    character(*), parameter :: default_revision = &
        '37f67f996913252dab8bf7f1b32996ba2fcce636'
    character(len=32), parameter :: names(5) = [character(len=32) :: &
        'coordinate', 'alpha', 'basis_value', &
        'basis_derivative', &
        'basis_second_derivative']
    type(arena_t), target :: arena
    type(assumption_context_t), target :: assumptions
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symengine
    type(suite_t) :: checks
    type(expr_t) :: x(5), direction(5), adjoint(3), weight, log_weight, scale
    type(expr_t) :: jet(3), derived(3), forward(3), reverse(5), objective
    type(expr_t) :: log_jet(3), log_derived(3)
    type(expr_t) :: from_scale(2), residual, lhs, rhs, b, c, q, left, h, t
    type(expr_t) :: affine, coefficient, flux, pi, m, canonical
    type(kernel_spec_t) :: spec
    type(str_t) :: sources(6)
    character(len=1024) :: directory, evidence
    character(len=40) :: revision
    character(len=80) :: label
    integer :: i, j, unit, primal_unit, status

    directory = 'src/generated'
    evidence = ''
    revision = default_revision
    if (command_argument_count() >= 1) call get_command_argument(1, directory)
    if (command_argument_count() >= 2) call get_command_argument(2, evidence)
    if (command_argument_count() >= 3) call get_command_argument(3, revision)
    if (command_argument_count() > 3) then
        error stop 'expected output directory [evidence JSON] [FortSym revision]'
    end if
    call arena%init()
    call assumptions%init(arena)
    do i = 1, size(x)
        x(i) = sym(arena, trim(names(i)))
        direction(i) = sym(arena, 'd_'//trim(names(i)))
        call assume(assumptions, real_valued(x(i)))
        call assume(assumptions, real_valued(direction(i)))
    end do
    call assume(assumptions, positive(x(1)))
    native = make_native_engine(arena, assumptions)
    symengine = make_symengine_engine(arena)
    call suite_begin(checks, 'GLISS weighted FEEC local jets exact identities')

    ! For s>0 and real alpha, real power is exp(-alpha*log(s)). Keeping
    ! power form in generated code preserves the existing custom-power API
    ! and its once-per-mode evaluation in the production caller.
    weight = x(1)**(-x(2))
    derived(1) = weight*x(3)
    derived(2) = ds(derived(1))
    derived(3) = ds(derived(2))
    jet(1) = weight*x(3)
    jet(2) = weight*(x(4) - x(2)*x(3)/x(1))
    jet(3) = weight*(x(5) - 2*x(2)*x(4)/x(1) &
        + x(2)*(x(2) + 1)*x(3)/x(1)**2)
    log_weight = exp(-x(2)*log(x(1)))
    log_derived(1) = log_weight*x(3)
    log_derived(2) = ds(log_derived(1))
    log_derived(3) = ds(log_derived(2))
    log_jet(1) = log_weight*x(3)
    log_jet(2) = log_weight*(x(4) - x(2)*x(3)/x(1))
    log_jet(3) = log_weight*(x(5) - 2*x(2)*x(4)/x(1) &
        + x(2)*(x(2) + 1)*x(3)/x(1)**2)
    do i = 1, size(jet)
        write (label, '(a,i0)') 'weighted total derivative jet ', i - 1
        ! The native engine does not normalize symbolic exponent shifts.
        ! Prove the emitted power form with SymEngine, and independently
        ! prove its defining positive-real logarithmic form with both.
        call prove_power(trim(label), derived(i) - jet(i))
        call prove(trim(label)//' real log definition', log_derived(i) - log_jet(i))
    end do

    adjoint(1) = sym(arena, 'value_cotangent')
    adjoint(2) = sym(arena, 'derivative_cotangent')
    adjoint(3) = sym(arena, 'second_derivative_cotangent')
    do i = 1, size(jet)
        forward(i) = 0*weight
        do j = 1, size(x)
            forward(i) = forward(i) + diff(jet(i), x(j))*direction(j)
        end do
    end do
    objective = adjoint(1)*jet(1) + adjoint(2)*jet(2) + adjoint(3)*jet(3)
    do i = 1, size(x)
        reverse(i) = diff(objective, x(i))
    end do
    lhs = adjoint(1)*forward(1) + adjoint(2)*forward(2) &
        + adjoint(3)*forward(3)
    rhs = reverse(1)*direction(1)
    do i = 2, size(x)
        rhs = rhs + reverse(i)*direction(i)
    end do
    call prove('jet JVP/VJP duality', lhs - rhs)
    call prove('weight alpha derivative retains custom real exponent', &
        diff(jet(1), x(2)) + log(x(1))*jet(1))

    ! A positive scale is an independently supplied value of s**(-alpha).
    ! Emit factored value/first-jet leaves so callers evaluate that power
    ! once per Fourier mode rather than once per active spline basis.
    scale = sym(arena, 'scale')
    from_scale(1) = scale*x(3)
    from_scale(2) = scale*(x(4) - x(2)*x(3)/x(1))
    do i = 1, 2
        write (label, '(a,i0)') 'factored production jet ', i - 1
        call prove(trim(label), weight*from_scale(i) - scale*jet(i))
    end do

    ! Independently specified monomials cover every physical preset power
    ! m=0..6 and the positive-integer kernel that the custom API allows.
    b = sym(arena, 'b')
    c = sym(arena, 'c')
    do i = 0, 6
        coefficient = rat(arena, int(2 - i, int64), 2_int64)
        if (i == 0) coefficient = 0*x(2)
        canonical = x(1)**rat(arena, int(i - 2, int64), 2_int64)
        if (i == 0) canonical = x(1)**0
        write (label, '(a,i0)') 'physical axis preset m=', i
        call prove_power(trim(label), diff(canonical*x(1), x(1)) &
            - canonical*(1 - coefficient))
        canonical = exp(-coefficient*log(x(1)))
        call prove(trim(label)//' real log derivative', diff(canonical*x(1), x(1)) &
            - canonical*(1 - coefficient))
    end do
    call prove('alpha one exact kernel N=s', diff(x(1)**(-1)*x(1), x(1)))
    call prove('alpha two exact kernel N=s squared', &
        diff(x(1)**(-2)*x(1)**2, x(1)))

    ! An interior affine coefficient N=b*s+c has derivative
    ! s**(-alpha)*((1-alpha)*b-alpha*c/s). Its normalized derivative
    ! is not generally constant, unlike the old weighted P0 companion.
    affine = b*x(1) + c
    residual = diff(weight*affine, x(1)) &
        - weight*((1 - x(2))*b - x(2)*c/x(1))
    call prove_power('interior P1 weighted derivative', residual)
    call prove('interior P1 real log derivative', &
        diff(log_weight*affine, x(1)) &
        - log_weight*((1 - x(2))*b - x(2)*c/x(1)))
    call prove('interior P1 normalized derivative obstruction', &
        diff((1 - x(2))*b - x(2)*c/x(1), x(1)) - x(2)*c/x(1)**2)

    ! Reference-cell chain rule: s=left+h*t, h>0. The comparison is
    ! independent of B-spline knot topology and endpoint trace choices.
    left = sym(arena, 'left')
    h = sym(arena, 'h')
    t = sym(arena, 't')
    q = left + h*t
    call prove_power('reference cell weighted affine derivative', &
        diff(q**(-x(2))*(b*q + c), t)/h &
        - q**(-x(2))*((1 - x(2))*b - x(2)*c/q))
    call prove('reference cell real log affine derivative', &
        diff(exp(-x(2)*log(q))*(b*q + c), t)/h &
        - exp(-x(2)*log(q))*((1 - x(2))*b - x(2)*c/q))
    flux = sym(arena, 'flux_t_slope')
    pi = sym(arena, 'pi')
    m = sym(arena, 'mode_m')
    call prove('zero-curl physical eta cancellation', &
        jet(2) + 2*pi*m*(-flux*jet(2)/(2*pi*m))/flux)

    call configure(spec, 'generated_weighted_feec_scale', .false., .false.)
    spec%args = [str('coordinate'), str('alpha')]
    spec%outputs = [str('scale')]
    sources(1) = checked_emission([weight], spec)

    call configure(spec, 'generated_weighted_feec_value_from_scale', .false., .false.)
    spec%args = [str('scale'), str('basis_value')]
    spec%outputs = [str('value')]
    sources(2) = checked_emission([from_scale(1)], spec)

    call configure(spec, 'generated_weighted_feec_first_from_scale', .false., .false.)
    spec%args = [str('coordinate'), str('alpha'), str('scale'), &
        str('basis_value'), str('basis_derivative')]
    spec%outputs = [str('value'), str('derivative')]
    sources(3) = checked_emission(from_scale, spec)

    call configure(spec, 'generated_weighted_feec_jet', .true., .false.)
    spec%outputs = [str('value'), str('derivative'), str('second_derivative')]
    sources(4) = checked_emission(jet, spec)

    call configure(spec, 'generated_weighted_feec_jvp', .true., .true.)
    allocate (spec%outputs(2), spec%output_shapes(2), spec%output_references(6))
    spec%outputs = [str('jet'), str('tangent')]
    spec%output_shapes = [str('(3)'), str('(3)')]
    do i = 1, 3
        write (label, '(a,i0,a)') 'jet(', i, ')'
        spec%output_references(i) = str(trim(label))
        write (label, '(a,i0,a)') 'tangent(', i, ')'
        spec%output_references(3 + i) = str(trim(label))
    end do
    sources(5) = checked_emission([jet, forward], spec)

    call configure(spec, 'generated_weighted_feec_vjp', .true., .false.)
    spec%args = [spec%args, str('value_cotangent'), str('derivative_cotangent'), &
        str('second_derivative_cotangent')]
    allocate (spec%outputs(2), spec%output_shapes(2), spec%output_references(8))
    spec%outputs = [str('jet'), str('gradient')]
    spec%output_shapes = [str('(3)'), str('(5)')]
    do i = 1, 3
        write (label, '(a,i0,a)') 'jet(', i, ')'
        spec%output_references(i) = str(trim(label))
    end do
    do i = 1, 5
        write (label, '(a,i0,a)') 'gradient(', i, ')'
        spec%output_references(3 + i) = str(trim(label))
    end do
    sources(6) = checked_emission([jet, reverse], spec)

    if (checks%by_probe /= 0) error stop 'probe evidence is forbidden'
    if (checks%proved /= checks%total) error stop 'unproved weighted FEEC identity'
    if (len_trim(evidence) > 0) then
        call suite_end(checks, trim(evidence))
    else
        call suite_end(checks)
    end if
    open (newunit=primal_unit, &
        file=trim(directory)//'/gliss_weighted_feec_primal.inc', &
        status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'cannot write generated weighted FEEC primal include'
    write (primal_unit, '(a)') &
        '! Generated by gen/gliss_weighted_feec_gen.f90; do not edit.'
    write (primal_unit, '(a)') '! FortSym revision: '//revision
    do i = 1, 3
        write (primal_unit, '(a)') chars(sources(i))
    end do
    write (primal_unit, '(a)') '! End of generated weighted FEEC primal leaves.'
    close (primal_unit)
    open (newunit=unit, file=trim(directory)//'/gliss_weighted_feec_generated.f90', &
        status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'cannot write generated weighted FEEC module'
    write (unit, '(a)') '! Generated by gen/gliss_weighted_feec_gen.f90; do not edit.'
    write (unit, '(a)') '! FortSym revision: '//revision
    write (unit, '(a)') '! Domain: coordinate > 0 and finite real alpha.'
    write (unit, '(a)') &
        '! Local jets only; topology, traces and assembly remain caller-owned.'
    write (unit, '(a)') 'module gliss_weighted_feec_generated'
    write (unit, '(a)') '    implicit none'
    write (unit, '(a)') '    private'
    write (unit, '(a)') '    public :: generated_weighted_feec_scale'
    write (unit, '(a)') '    public :: generated_weighted_feec_value_from_scale'
    write (unit, '(a)') '    public :: generated_weighted_feec_first_from_scale'
    write (unit, '(a)') '    public :: generated_weighted_feec_jet'
    write (unit, '(a)') &
        '    public :: generated_weighted_feec_jvp, generated_weighted_feec_vjp'
    write (unit, '(a)') 'contains'
    write (unit, '(a)') "    include 'gliss_weighted_feec_primal.inc'"
    do i = 4, size(sources)
        write (unit, '(a)') chars(sources(i))
    end do
    write (unit, '(a)') 'end module gliss_weighted_feec_generated'
    close (unit)
    write (output_unit, '(a)') &
        'Wrote '//trim(directory)//'/gliss_weighted_feec_generated.f90'
    write (output_unit, '(a)') &
        'Wrote '//trim(directory)//'/gliss_weighted_feec_primal.inc'

contains

    function ds(expression) result(derivative)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: derivative

        derivative = diff(expression, x(1)) + diff(expression, x(3))*x(4) &
            + diff(expression, x(4))*x(5)
    end function ds

    subroutine prove(description, expression)
        character(*), intent(in) :: description
        type(expr_t), intent(in) :: expression

        call check_identity(checks, native, &
            trim(description)//' [native]', expression)
        call check_identity(checks, symengine, &
            trim(description)//' [SymEngine]', expression)
        if (checks%failed /= 0) error stop 'weighted FEEC identity failed or undecided'
    end subroutine prove

    subroutine prove_power(description, expression)
        character(*), intent(in) :: description
        type(expr_t), intent(in) :: expression

        call check_identity(checks, symengine, &
            trim(description)//' [power SymEngine]', expression)
        if (checks%failed /= 0) error stop 'weighted power identity failed or undecided'
    end subroutine prove_power

    subroutine configure(configuration, name, inputs, tangent)
        type(kernel_spec_t), intent(out) :: configuration
        character(*), intent(in) :: name
        logical, intent(in) :: inputs, tangent
        integer :: argument, argument_count

        configuration%name = str(name)
        configuration%mode = KERNEL_SUBROUTINE
        configuration%pure_procedure = .true.
        configuration%temp_prefix = str('t')
        configuration%generator = str('gen/gliss_weighted_feec_gen.f90')
        configuration%generator_revision = str(revision)
        configuration%regenerate_command = str('gliss_weighted_feec_gen src/generated')
        if (.not. inputs) return
        argument_count = size(x)
        if (tangent) argument_count = 2*argument_count
        allocate (configuration%args(argument_count))
        do argument = 1, size(x)
            configuration%args(argument) = str(trim(names(argument)))
            if (tangent) configuration%args(size(x) + argument) = &
                str('d_'//trim(names(argument)))
        end do
    end subroutine configure

    function checked_emission(roots, configuration) result(source)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t), intent(in) :: configuration
        type(str_t) :: source
        character(:), allocatable :: message
        logical :: good

        ! Preserve the proved factored power form. In particular, the custom
        ! alpha=1 exact kernel must cancel before multiplication by its weight.
        source = emit_kernel(roots, configuration, ok=good, message=message)
        if (.not. good) then
            write (output_unit, '(a)') message
            error stop 'weighted FEEC emission failed'
        end if
    end function checked_emission

end program gliss_weighted_feec_gen
