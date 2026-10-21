program generate_cartesian_field
    ! Cartesian magnetic field and spatial derivatives from primitive jets.
    ! The defining flux field is differentiated natively by FortSym; staged
    ! emission preserves its quotient/product structure on either J /= 0 chart.
    use, intrinsic :: iso_fortran_env, only: output_unit, int64
    use fortsym_arena, only: arena_t
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, rat, sqrt, &
        operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, CSE_NONE
    use fortsym_string, only: str_t, str, chars
    use fortsym_subs, only: subs_many
    implicit none

    character(*), parameter :: default_revision = &
        'ab86cc4cbf3e3afd6fd3f229ed855f3095e668e3'
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(suite_t) :: checks
    type(expr_t) :: x(32), y(46), roots(46), exact(46), frame(2)
    type(expr_t) :: basis(3, 3), second(3, 3, 3), reciprocal(3, 3)
    type(expr_t) :: cofactors(3, 3)
    type(expr_t) :: field(3), field_q(3, 3), covariant(3), covariant_q(3, 3)
    type(expr_t) :: determinant, contra(2), residual, exterior(3), expanded
    type(kernel_spec_t) :: spec
    type(str_t) :: field_source, frame_source
    character(len=1024) :: directory, report
    character(len=40) :: revision
    character(len=80) :: label
    integer :: i, j, k, n, unit, status

    directory = 'src/generated'
    report = ''
    revision = default_revision
    if (command_argument_count() >= 1) call get_command_argument(1, directory)
    if (command_argument_count() >= 2) call get_command_argument(2, report)
    if (command_argument_count() >= 3) call get_command_argument(3, revision)
    if (command_argument_count() > 3) then
        error stop 'expected output directory [identity JSON] [FortSym revision]'
    end if
    call arena%init()
    engine = make_symengine_engine(arena)
    call suite_begin(checks, 'GLISS Cartesian field exact identities')
    do i = 1, size(x)
        x(i) = indexed_symbol('x', i)
    end do
    do i = 1, size(y)
        y(i) = indexed_symbol('work', i)
    end do
    basis(:, 1) = x(1:3)
    basis(:, 2) = x(4:6)
    basis(:, 3) = x(7:9)
    second(:, 1, 1) = x(10:12)
    second(:, 1, 2) = x(13:15)
    second(:, 2, 1) = x(13:15)
    second(:, 1, 3) = x(16:18)
    second(:, 3, 1) = x(16:18)
    second(:, 2, 2) = x(19:21)
    second(:, 2, 3) = x(22:24)
    second(:, 3, 2) = x(22:24)
    second(:, 3, 3) = x(25:27)
    cofactors(:, 1) = cross3(basis(:, 2), basis(:, 3))
    cofactors(:, 2) = cross3(basis(:, 3), basis(:, 1))
    cofactors(:, 3) = cross3(basis(:, 1), basis(:, 2))
    determinant = dot3(basis(:, 1), cofactors(:, 1))
    contra(1) = -x(29) * x(32) / determinant
    contra(2) = -x(28) / determinant
    do i = 1, 3
        field(i) = contra(1) * basis(i, 2) + contra(2) * basis(i, 3)
        do j = 1, 3
            reciprocal(i, j) = cofactors(i, j) / determinant
        end do
    end do
    do j = 1, 3
        do i = 1, 3
            field_q(i, j) = coordinate_derivative(field(i), j)
        end do
        covariant(j) = dot3(field, basis(:, j))
        do i = 1, 3
            covariant_q(j, i) = coordinate_derivative(covariant(j), i)
        end do
    end do

    ! x(28:32) = Phi_s, chi_s, Phi_ss, chi_ss, inverse_NFP.
    ! Hidden ordered work: J, J_q/J(3), contra(2), contra_q(6),
    ! B_q(9), reciprocal basis(9), B(3), Cartesian gradient(9), curl(3), div.
    exact(1) = determinant
    roots(1) = determinant
    do j = 1, 3
        exact(1 + j) = coordinate_derivative(determinant, j) / determinant
        roots(1 + j) = coordinate_derivative(determinant, j) / y(1)
    end do
    exact(5:6) = contra
    roots(5) = -x(29) * x(32) / y(1)
    roots(6) = -x(28) / y(1)
    do j = 1, 3
        exact(6 + j) = coordinate_derivative(contra(1), j)
        exact(9 + j) = coordinate_derivative(contra(2), j)
        roots(6 + j) = -y(5) * y(1 + j)
        roots(9 + j) = -y(6) * y(1 + j)
        if (j == 1) then
            roots(6 + j) = roots(6 + j) - x(31) * (x(32) / y(1))
            roots(9 + j) = roots(9 + j) - x(30) / y(1)
        end if
        do i = 1, 3
            n = 12 + i + 3 * (j - 1)
            exact(n) = field_q(i, j)
            roots(n) = y(6 + j) * basis(i, 2) + y(5) * second(i, 2, j) &
                + y(9 + j) * basis(i, 3) + y(6) * second(i, 3, j)
            n = 21 + i + 3 * (j - 1)
            exact(n) = reciprocal(i, j)
            roots(n) = cofactors(i, j) / y(1)
        end do
    end do
    exact(31:33) = field
    do i = 1, 3
        roots(30 + i) = y(5) * basis(i, 2) + y(6) * basis(i, 3)
    end do
    do j = 1, 3
        do i = 1, 3
            n = 33 + i + 3 * (j - 1)
            exact(n) = zero()
            roots(n) = zero()
            do k = 1, 3
                exact(n) = exact(n) + field_q(i, k) * reciprocal(j, k)
                roots(n) = roots(n) + y(12 + i + 3 * (k - 1)) &
                    * y(21 + j + 3 * (k - 1))
            end do
        end do
    end do
    exact(43) = exact(39) - exact(41)
    exact(44) = exact(40) - exact(36)
    exact(45) = exact(35) - exact(37)
    exact(46) = exact(34) + exact(38) + exact(42)
    roots(43) = y(39) - y(41)
    roots(44) = y(40) - y(36)
    roots(45) = y(35) - y(37)
    roots(46) = y(34) + y(38) + y(42)
    do i = 1, size(roots)
        expanded = subs_many(roots(i), y, exact)
        write (label, '(a,i0)') 'ordered native derivative stage ', i
        call prove(trim(label), expanded - exact(i))
    end do
    do j = 1, 3
        do i = 1, 3
            write (label, '(a,i0,a,i0)') 'reciprocal basis ', i, ',', j
            residual = dot3(reciprocal(:, i), basis(:, j))
            if (i == j) residual = residual - 1
            call prove(trim(label), residual)
            residual = zero()
            do k = 1, 3
                residual = residual + exact(33 + i + 3 * (k - 1)) * basis(k, j)
            end do
            write (label, '(a,i0,a,i0)') 'Cartesian chain rule ', i, ',', j
            call prove(trim(label), residual - field_q(i, j))
        end do
    end do
    exterior(1) = (covariant_q(3, 2) - covariant_q(2, 3)) / determinant
    exterior(2) = (covariant_q(1, 3) - covariant_q(3, 1)) / determinant
    exterior(3) = (covariant_q(2, 1) - covariant_q(1, 2)) / determinant
    do i = 1, 3
        residual = zero()
        do j = 1, 3
            residual = residual + basis(i, j) * exterior(j)
        end do
        write (label, '(a,i0)') 'covariant exterior derivative curl ', i
        call prove(trim(label), exact(42 + i) - residual)
    end do
    call prove('flux field is divergence free', exact(46))

    ! Frame validation runs before the field kernel divides by signed J.
    frame(1) = determinant
    frame(2) = sqrt(dot3(basis(:, 1), basis(:, 1))) &
        * sqrt(dot3(basis(:, 2), basis(:, 2))) &
        * sqrt(dot3(basis(:, 3), basis(:, 3)))
    call configure(spec, 'generated_cartesian_field_frame', 9, 2, 'values')
    frame_source = verified_emission(frame, spec)
    call configure(spec, 'generated_cartesian_field_work', 32, 46, 'work')
    spec%cse_level = CSE_NONE
    spec%ordered_output_dependencies = .true.
    field_source = verified_emission(roots, spec)
    if (len_trim(report) > 0) then
        call suite_end(checks, trim(report))
    else
        call suite_end(checks)
    end if
    if (checks%failed /= 0) error stop 'Cartesian field strict identity failure'
    open (newunit=unit, &
        file=trim(directory)//'/gliss_cartesian_field_generated.f90', &
        status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'cannot write generated Cartesian field module'
    write (unit, '(a)') '! Generated by gen/generate_cartesian_field.f90; do not edit.'
    write (unit, '(a)') '! FortSym revision: '//trim(revision)
    write (unit, '(a)') '! Domain: finite primitive jets, regular signed J /= 0.'
    write (unit, '(a)') '! x(32): r_s,r_t,r_z,r_ss,r_st,r_sz,r_tt,r_tz,r_zz (3 each),'
    write (unit, '(a)') '! Phi_s,chi_s,Phi_ss,chi_ss,inverse_NFP; positions m, fluxes Wb.'
    write (unit, '(a)') '! Angles are turns; r_z differentiates one-period zeta.'
    write (unit, '(a)') '! values(16): Cartesian B(3), column-major dB_i/dx_j(9),'
    write (unit, '(a)') '! curl B(3), div B; units T and T/m respectively.'
    write (unit, '(a)') '! Parameter JVP/VJP are not emitted by this spatial-jet generator.'
    write (unit, '(a)') 'module gliss_cartesian_field_generated'
    write (unit, '(a)') '    use, intrinsic :: iso_fortran_env, only: dp => real64'
    write (unit, '(a)') '    implicit none'
    write (unit, '(a)') '    private'
    write (unit, '(a)') '    public :: generated_cartesian_field'
    write (unit, '(a)') '    public :: generated_cartesian_field_frame'
    write (unit, '(a)') 'contains'
    write (unit, '(a)') chars(frame_source)
    write (unit, '(a)') chars(field_source)
    write (unit, '(a)') 'pure subroutine generated_cartesian_field(x, values)'
    write (unit, '(a)') '    real(dp), intent(in) :: x(32)'
    write (unit, '(a)') '    real(dp), intent(out) :: values(16)'
    write (unit, '(a)') '    real(dp) :: work(46)'
    write (unit, '(a)') '    call generated_cartesian_field_work(x, work)'
    write (unit, '(a)') '    values = work(31:46)'
    write (unit, '(a)') 'end subroutine generated_cartesian_field'
    write (unit, '(a)') 'end module gliss_cartesian_field_generated'
    close (unit, iostat=status)
    if (status /= 0) error stop 'cannot close generated Cartesian field module'
    write (output_unit, '(a)') &
        'Wrote '//trim(directory)//'/gliss_cartesian_field_generated.f90'

contains

    function indexed_symbol(name, index) result(value)
        character(*), intent(in) :: name
        integer, intent(in) :: index
        type(expr_t) :: value
        character(len=16) :: number

        write (number, '(i0)') index
        value = sym(arena, name//'_'//trim(number))
    end function indexed_symbol

    function zero() result(value)
        type(expr_t) :: value

        value = rat(arena, 0_int64, 1_int64)
    end function zero

    function dot3(left, right) result(value)
        type(expr_t), intent(in) :: left(3), right(3)
        type(expr_t) :: value

        value = left(1) * right(1) + left(2) * right(2) + left(3) * right(3)
    end function dot3

    function cross3(left, right) result(value)
        type(expr_t), intent(in) :: left(3), right(3)
        type(expr_t) :: value(3)

        value(1) = left(2) * right(3) - left(3) * right(2)
        value(2) = left(3) * right(1) - left(1) * right(3)
        value(3) = left(1) * right(2) - left(2) * right(1)
    end function cross3

    function coordinate_derivative(value, coordinate) result(derivative)
        type(expr_t), intent(in) :: value
        integer, intent(in) :: coordinate
        type(expr_t) :: derivative
        integer :: component, axis

        derivative = zero()
        do axis = 1, 3
            do component = 1, 3
                derivative = derivative + diff(value, basis(component, axis)) &
                    * second(component, axis, coordinate)
            end do
        end do
        if (coordinate == 1) derivative = derivative &
            + diff(value, x(28)) * x(30) + diff(value, x(29)) * x(31)
    end function coordinate_derivative

    subroutine prove(name, residual)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: residual

        call check_identity(checks, engine, name, residual)
        if (checks%failed /= 0) then
            error stop 'Cartesian field identity failed or undecided'
        end if
    end subroutine prove

    subroutine configure(configuration, name, inputs, outputs, target)
        type(kernel_spec_t), intent(out) :: configuration
        character(*), intent(in) :: name, target
        integer, intent(in) :: inputs, outputs
        character(len=40) :: reference
        integer :: output

        configuration%name = str(name)
        configuration%pure_procedure = .true.
        configuration%generator = str('gen/generate_cartesian_field.f90')
        configuration%generator_revision = str(revision)
        configuration%regenerate_command = str('generate_cartesian_field src/generated')
        allocate (configuration%args(1), configuration%arg_shapes(1))
        configuration%args(1) = str('x')
        write (reference, '(a,i0,a)') '(', inputs, ')'
        configuration%arg_shapes(1) = str(trim(reference))
        allocate (configuration%outputs(1), configuration%output_shapes(1), &
            configuration%output_references(outputs))
        configuration%outputs(1) = str(target)
        write (reference, '(a,i0,a)') '(', outputs, ')'
        configuration%output_shapes(1) = str(trim(reference))
        do output = 1, outputs
            write (reference, '(a,a,i0,a)') target, '(', output, ')'
            configuration%output_references(output) = str(trim(reference))
        end do
    end subroutine configure

    function verified_emission(values, configuration) result(source)
        type(expr_t), intent(in) :: values(:)
        type(kernel_spec_t), intent(in) :: configuration
        type(str_t) :: source
        type(expr_t) :: symbols(78), references(78), emitted(size(values)), restored
        character(:), allocatable :: message
        character(len=80) :: check_label
        character(len=32) :: reference
        integer :: index
        logical :: good

        symbols(1:32) = x
        symbols(33:78) = y
        do index = 1, 32
            write (reference, '(a,i0,a)') 'x(', index, ')'
            references(index) = sym(arena, trim(reference))
        end do
        do index = 1, 46
            write (reference, '(a,i0,a)') 'work(', index, ')'
            references(32 + index) = sym(arena, trim(reference))
        end do
        do index = 1, size(values)
            emitted(index) = subs_many(values(index), symbols, references)
            restored = subs_many(emitted(index), references, symbols)
            write (check_label, '(a,a,i0)') chars(configuration%name), &
                ' packing ', index
            call prove(trim(check_label), restored - values(index))
        end do
        source = emit_kernel(emitted, configuration, ok=good, message=message)
        if (.not. good) then
            write (output_unit, '(a)') message
            error stop 'Cartesian field emission failed'
        end if
    end function verified_emission

end program generate_cartesian_field
