program gliss_pressure_adjoint_gen
    ! Reproduce scalar leaves for the fixed-geometry pressure reverse sweep.
    ! The generator derives cotangents from the primal energy and spline maps.
    ! Intervals are positive; all geometry, field and quadrature data are fixed.
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_string, only: str, chars, itoa, words
    implicit none

    type(native_engine_t) :: engine
    type(suite_t) :: suite
    type(expr_t) :: r2, r3, r4, r5, d2, d3, drive_slope, volume, inv_mu0
    type(expr_t) :: pressure_slope, gamma, energy, roots4(4), expected4(4)
    type(expr_t) :: dl, dr, h, yl, yr, ml, mr, value_bar, slope_bar
    type(expr_t) :: value, slope, objective, hl, hr, yc, rhs_bar, rhs
    type(expr_t) :: roots3(3), expected3(3), ho, hi, adjacent, next
    type(expr_t) :: endpoint_bar, endpoint, roots2(2), expected2(2)
    type(kernel_spec_t) :: spec
    character(len=2048) :: output, revision, report
    integer :: output_length, revision_length, report_length, status, unit, ios

    call get_command_argument(1, output, length=output_length, status=status)
    if (status /= 0) error stop 'cannot read output path'
    call get_command_argument(2, revision, length=revision_length, status=status)
    if (status /= 0) error stop 'cannot read FortSym revision'
    report_length = 0
    if (command_argument_count() >= 3) then
        call get_command_argument(3, report, length=report_length, status=status)
        if (status /= 0) error stop 'cannot read identity report path'
    end if
    if (output_length == 0 .or. revision_length == 0) then
        write (output_unit, '(a)') &
            'usage: gliss_pressure_adjoint_gen OUTPUT_PATH FORTSYM_REVISION [IDENTITY_REPORT]'
        error stop 2
    end if

    engine = make_native_engine()
    call suite_begin(suite, 'GLISS pressure reverse scalar identities')

    open (newunit=unit, file=trim(output), status='replace', &
          action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot write generated pressure adjoint module'
    write (unit, '(a)') 'module gliss_pressure_adjoint_generated'
    write (unit, '(a)') '    implicit none'
    write (unit, '(a)') '    private'
    write (unit, '(a)') '    public :: pressure_energy_pullback'
    write (unit, '(a)') '    public :: gliss_spline_evaluation_vjp'
    write (unit, '(a)') '    public :: gliss_spline_rhs_vjp'
    write (unit, '(a)') '    public :: gliss_spline_endpoint_vjp'
    write (unit, '(a)') 'contains'

    r2 = sym('r2')
    r3 = sym('r3')
    r4 = sym('r4')
    r5 = sym('r5')
    d2 = sym('d2')
    d3 = sym('d3')
    drive_slope = sym('drive_slope')
    volume = sym('volume')
    inv_mu0 = sym('inv_mu0')
    pressure_slope = sym('pressure_slope')
    gamma = sym('gamma')
    energy = volume*inv_mu0*(r2**2 + r3**2 &
        - pressure_slope*drive_slope*r4**2) + gamma*volume*r5**2
    roots4 = [diff(energy, r2)*d2, diff(energy, r3)*d3, &
        diff(energy, pressure_slope), diff(energy, gamma)]
    expected4 = [2*volume*inv_mu0*r2*d2, 2*volume*inv_mu0*r3*d3, &
        -drive_slope*volume*inv_mu0*r4**2, volume*r5**2]
    call prove_and_simplify(roots4, expected4, 'pressure energy')
    call initialize_spec('pressure_energy_pullback')
    spec%args = words('r2 r3 r4 r5 d2 d3 drive_slope volume inv_mu0')
    spec%outputs = [str('terms')]
    spec%output_shapes = [str('(4)')]
    spec%output_references = words('terms(1) terms(2) terms(3) terms(4)')
    write (unit, '(a)') chars(emit_kernel(roots4, spec))

    dl = sym('distance_left')
    dr = sym('distance_right')
    h = sym('interval')
    yl = sym('value_left')
    yr = sym('value_right')
    ml = sym('second_left')
    mr = sym('second_right')
    value_bar = sym('value_cotangent')
    slope_bar = sym('slope_cotangent')
    value = ml*dr**3/(6*h) + mr*dl**3/(6*h) &
        + (yl - ml*h**2/6)*dr/h + (yr - mr*h**2/6)*dl/h
    slope = -ml*dr**2/(2*h) + mr*dl**2/(2*h) &
        + (yr - yl)/h - h*(mr - ml)/6
    objective = value_bar*value + slope_bar*slope
    roots4 = [diff(objective, yl), diff(objective, yr), &
        diff(objective, ml), diff(objective, mr)]
    expected4 = [value_bar*dr/h - slope_bar/h, &
        value_bar*dl/h + slope_bar/h, &
        value_bar*(dr**3/(6*h) - h*dr/6) &
        + slope_bar*(-dr**2/(2*h) + h/6), &
        value_bar*(dl**3/(6*h) - h*dl/6) &
        + slope_bar*(dl**2/(2*h) - h/6)]
    call prove_and_simplify(roots4, expected4, 'spline evaluation')
    call initialize_spec('gliss_spline_evaluation_vjp')
    spec%args = words('distance_left distance_right interval value_cotangent slope_cotangent')
    spec%outputs = words('value_left_cotangent value_right_cotangent ' &
        //'second_left_cotangent second_right_cotangent')
    write (unit, '(a)') chars(emit_kernel(roots4, spec))

    hl = sym('interval_left')
    hr = sym('interval_right')
    yc = sym('value_center')
    rhs_bar = sym('rhs_cotangent')
    rhs = 6*((yr - yc)/hr - (yc - yl)/hl)
    objective = rhs_bar*rhs
    roots3 = [diff(objective, yl), diff(objective, yc), diff(objective, yr)]
    expected3 = [6*rhs_bar/hl, -6*rhs_bar*(1/hr + 1/hl), 6*rhs_bar/hr]
    call prove_and_simplify(roots3, expected3, 'spline RHS')
    call initialize_spec('gliss_spline_rhs_vjp')
    spec%args = words('interval_left interval_right rhs_cotangent')
    spec%outputs = words('value_left_cotangent value_center_cotangent value_right_cotangent')
    write (unit, '(a)') chars(emit_kernel(roots3, spec))

    ho = sym('interval_outer')
    hi = sym('interval_inner')
    adjacent = sym('adjacent')
    next = sym('next')
    endpoint_bar = sym('endpoint_cotangent')
    endpoint = ((ho + hi)*adjacent - ho*next)/hi
    objective = endpoint_bar*endpoint
    roots2 = [diff(objective, adjacent), diff(objective, next)]
    expected2 = [endpoint_bar*(ho + hi)/hi, -endpoint_bar*ho/hi]
    call prove_and_simplify(roots2, expected2, 'spline endpoint')
    call initialize_spec('gliss_spline_endpoint_vjp')
    spec%args = words('interval_outer interval_inner endpoint_cotangent')
    spec%outputs = words('adjacent_cotangent next_cotangent')
    write (unit, '(a)') chars(emit_kernel(roots2, spec))

    write (unit, '(a)') 'end module gliss_pressure_adjoint_generated'
    close (unit)
    if (report_length > 0) then
        call suite_end(suite, trim(report))
    else
        call suite_end(suite)
    end if
    write (output_unit, '(a)') 'wrote '//trim(output)

contains

    subroutine prove_and_simplify(roots, expected, label)
        type(expr_t), intent(inout) :: roots(:)
        type(expr_t), intent(in) :: expected(:)
        character(*), intent(in) :: label
        type(engine_result_t) :: simplified
        integer :: k

        do k = 1, size(roots)
            call check_identity(suite, engine, label//' '//itoa(k), &
                roots(k) - expected(k))
            simplified = engine%simplify(roots(k))
            if (.not. simplified%ok) error stop 'FortSym simplification failed'
            call check_identity(suite, engine, label//' simplified '//itoa(k), &
                simplified%value - roots(k))
            roots(k) = simplified%value
        end do
        if (suite%failed /= 0) error stop 'pressure adjoint identity failure'
    end subroutine prove_and_simplify

    subroutine initialize_spec(name)
        character(*), intent(in) :: name

        spec = kernel_spec_t()
        spec%name = str(name)
        spec%pure_procedure = .true.
        spec%temp_prefix = str('t')
        spec%generator = str('gen/gliss_pressure_adjoint_gen.f90')
        spec%generator_revision = str(trim(revision))
        spec%regenerate_command = str( &
            'cmake --build build-codegen --target gliss_regenerate_kernels')
    end subroutine initialize_spec

end program gliss_pressure_adjoint_gen
