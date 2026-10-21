program gliss_geometry_gen
    ! Cartesian primitive geometry from source-authored dot/cross identities.
    ! No external CAS expressions or fixtures enter this derivation.
    ! Discrete field-period count and orientation belong to the caller.
    use, intrinsic :: iso_fortran_env, only: output_unit, int64
    use fortsym_arena, only: arena_t
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, sqrt, rat, &
        operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, CSE_NONE
    use fortsym_string, only: str_t, str, chars
    use fortsym_subs, only: subs_many
    implicit none

    character(*), parameter :: default_revision = &
        'f5ef0572a524378abf3ae83669ee856684691c0a'
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(suite_t) :: checks
    type(expr_t) :: x(27), polynomial(13), field(5), curvature(3)
    type(expr_t) :: normal_values(3), contraction(3), composed_partial, partial
    type(expr_t) :: derivative_variables(12), normal_substitution(12)
    type(expr_t) :: angular_cross(3), product_left(3), product_right(3)
    type(expr_t) :: determinant, directional, expected
    type(expr_t) :: orientation, field_squared, normal_squared
    type(str_t) :: sources(12)
    character(len=1024) :: directory, report
    character(len=40) :: revision
    character(len=80) :: label
    integer :: i, j, k, first, second, unit, primal_unit, status

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
    call suite_begin(checks, 'GLISS Cartesian primitive geometry exact identities')
    do i = 1, size(x)
        x(i) = indexed_symbol('x', i)
    end do

    ! x: r_s, r_theta, r_zeta, r_ss, r_s_theta, r_s_zeta,
    ! r_theta_theta, r_theta_zeta, r_zeta_zeta, each a Cartesian 3-vector.
    ! y: g_ss, g_s_theta, g_s_zeta, g_theta_theta, g_theta_zeta,
    ! g_zeta_zeta, J, J_s, J_theta, J_zeta, (r_theta cross r_zeta)(1:3).
    k = 0
    do first = 1, 3
        do second = first, 3
            k = k + 1
            polynomial(k) = dot3(x(3*first - 2:3*first), &
                x(3*second - 2:3*second))
        end do
    end do
    angular_cross = cross3(x(4:6), x(7:9))
    polynomial(7) = dot3(x(1:3), angular_cross)
    polynomial(8:10) = 0*x(1)
    do i = 1, 3
        polynomial(8) = polynomial(8) &
            + diff(polynomial(7), x(i))*x(9 + i) &
            + diff(polynomial(7), x(3 + i))*x(12 + i) &
            + diff(polynomial(7), x(6 + i))*x(15 + i)
        polynomial(9) = polynomial(9) &
            + diff(polynomial(7), x(i))*x(12 + i) &
            + diff(polynomial(7), x(3 + i))*x(18 + i) &
            + diff(polynomial(7), x(6 + i))*x(21 + i)
        polynomial(10) = polynomial(10) &
            + diff(polynomial(7), x(i))*x(15 + i) &
            + diff(polynomial(7), x(3 + i))*x(21 + i) &
            + diff(polynomial(7), x(6 + i))*x(24 + i)
    end do
    polynomial(11:13) = angular_cross
    determinant = polynomial(1)*(polynomial(4)*polynomial(6) &
        - polynomial(5)**2) - polynomial(2)*(polynomial(2)*polynomial(6) &
        - polynomial(3)*polynomial(5)) + polynomial(3) &
        *(polynomial(2)*polynomial(5) - polynomial(3)*polynomial(4))
    call prove('metric determinant equals signed J squared', &
        determinant - polynomial(7)**2)
    call prove('angular area squared equals angular metric determinant', &
        dot3(angular_cross, angular_cross) &
        - polynomial(4)*polynomial(6) + polynomial(5)**2)
    call prove('angular cross perpendicular to theta', dot3(angular_cross, x(4:6)))
    call prove('angular cross perpendicular to zeta', dot3(angular_cross, x(7:9)))
    product_left = cross3(x(13:15), x(7:9))
    product_right = cross3(x(4:6), x(16:18))
    expected = dot3(x(10:12), angular_cross) &
        + dot3(x(1:3), product_left) + dot3(x(1:3), product_right)
    call prove('radial determinant product rule', polynomial(8) - expected)
    product_left = cross3(x(19:21), x(7:9))
    product_right = cross3(x(4:6), x(22:24))
    expected = dot3(x(13:15), angular_cross) &
        + dot3(x(1:3), product_left) + dot3(x(1:3), product_right)
    call prove('theta determinant product rule', polynomial(9) - expected)
    product_left = cross3(x(22:24), x(7:9))
    product_right = cross3(x(4:6), x(25:27))
    expected = dot3(x(16:18), angular_cross) &
        + dot3(x(1:3), product_left) + dot3(x(1:3), product_right)
    call prove('zeta determinant product rule', polynomial(10) - expected)
    do j = 1, size(polynomial)
        directional = 0*x(1)
        do i = 1, size(x)
            directional = directional + x(i)*diff(polynomial(j), x(i))
        end do
        k = 2
        if (j >= 7 .and. j <= 10) k = 3
        write (label, '(a,i0)') 'geometry homogeneous degree ', j
        call prove(trim(label), directional - k*polynomial(j))
    end do
    call emit_family('generated_geometry_polynomial', x, polynomial, sources(1:3))

    ! Field x: g_theta_theta, g_theta_zeta, g_zeta_zeta, signed J,
    ! toroidal_flux_slope, poloidal_flux_slope, inverse_field_periods.
    ! y: B^theta, B^zeta, B_theta, B_zeta, |B|.
    ! J /= 0, positive-definite angular metric, nonzero flux pair.
    field(1) = -x(6)*x(7)/x(4)
    field(2) = -x(5)/x(4)
    field(3) = x(1)*field(1) + x(2)*field(2)
    field(4) = x(2)*field(1) + x(3)*field(2)
    field_squared = field(1)*field(3) + field(2)*field(4)
    field(5) = sqrt(field_squared)
    call prove('poloidal flux sign and field-period factor', &
        -x(4)*field(1) - x(6)*x(7))
    call prove('toroidal flux sign', -x(4)*field(2) - x(5))
    expected = (x(1)*x(6)**2*x(7)**2 &
        + 2*x(2)*x(5)*x(6)*x(7) + x(3)*x(5)**2)/x(4)**2
    call prove('field squared metric contraction', field_squared - expected)
    call prove('positive magnitude squared', field(5)**2 - field_squared)
    do i = 1, 7
        write (label, '(a,i0)') 'magnitude differentiated squared identity ', i
        call prove(trim(label), 2*field(5)*diff(field(5), x(i)) &
            - diff(field_squared, x(i)))
    end do
    call emit_staged_field(x(1:7), field, sources(4:6))

    ! Curvature x: angular_cross(1:3), r_theta_theta(1:3),
    ! r_theta_zeta(1:3), r_zeta_zeta(1:3); y: II_theta_theta,
    ! II_theta_zeta, II_zeta_zeta. Orientation is sign(J) held fixed.
    ! This is the second form with the outward signed-chart normal, not
    ! a derivative across J=0. Cross-product norm must be nonzero.
    orientation = sym(arena, 'orientation')
    normal_squared = dot3(x(1:3), x(1:3))
    do i = 1, 3
        normal_values(i) = orientation*x(i)/sqrt(normal_squared)
    end do
    do i = 1, 3
        curvature(i) = dot3(normal_values, x(3*i + 1:3*i + 3))
        contraction(i) = dot3(x(1:3), x(3*i + 1:3*i + 3))
        directional = 0*x(1)
        do j = 1, 3
            directional = directional + x(j)*diff(curvature(i), x(j))
        end do
        write (label, '(a,i0)') 'normal magnitude scaling invariant ', i
        call prove(trim(label), directional)
        expected = 0*x(1)
        do j = 1, 3
            expected = expected &
                + x(3*i + j)*diff(curvature(i), x(3*i + j))
        end do
        write (label, '(a,i0)') 'second derivative linear contraction ', i
        call prove(trim(label), curvature(i) - expected)
    end do
    derivative_variables = x(1:12)
    normal_substitution = x(1:12)
    normal_substitution(1:3) = normal_values
    do i = 1, 3
        do j = 1, 12
            composed_partial = rat(arena, 0_int64, 1_int64)
            if (j <= 3) then
                do k = 1, 3
                    partial = diff(contraction(i), x(k))
                    partial = subs_many(partial, derivative_variables, &
                        normal_substitution)
                    composed_partial = composed_partial &
                        + partial*diff(normal_values(k), x(j))
                end do
            else
                composed_partial = diff(contraction(i), x(j))
                composed_partial = subs_many(composed_partial, &
                    derivative_variables, normal_substitution)
            end if
            write (label, '(a,i0,a,i0)') 'second form normal composition ', i, ',', j
            call prove(trim(label), composed_partial - diff(curvature(i), x(j)))
        end do
    end do
    call emit_staged_normal(x(1:3), normal_values, sources(7:9))
    call emit_family('generated_geometry_normal_contraction', x(1:12), &
        contraction, sources(10:12))

    ! Strict identities must succeed before replacing production output.
    if (len_trim(report) > 0) then
        call suite_end(checks, trim(report))
    else
        call suite_end(checks)
    end if
    if (checks%failed /= 0) error stop 'geometry strict identity failure'
    open (newunit=primal_unit, &
        file=trim(directory)//'/gliss_geometry_primal.inc', &
        status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'cannot write generated geometry primal include'
    write (primal_unit, '(a)') &
        '! Generated by gen/gliss_geometry_gen.f90; do not edit.'
    write (primal_unit, '(a)') '! FortSym revision: '//trim(revision)
    write (primal_unit, '(a)') &
        '! Shared primal procedures permit inlining within the consuming module.'
    do i = 1, size(sources), 3
        write (primal_unit, '(a)') chars(sources(i))
    end do
    call write_second_form_primal(primal_unit)
    close (primal_unit, iostat=status)
    if (status /= 0) error stop 'cannot close generated geometry primal include'

    open (newunit=unit, file=trim(directory)//'/gliss_geometry_generated.f90', &
        status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'cannot write generated geometry module'
    call write_header(unit)
    write (unit, '(a)') "    include 'gliss_geometry_primal.inc'"
    do i = 1, size(sources)
        if (modulo(i, 3) == 1) cycle
        write (unit, '(a)') chars(sources(i))
    end do
    call write_second_form_products(unit)
    write (unit, '(a)') 'end module gliss_geometry_generated'
    close (unit, iostat=status)
    if (status /= 0) error stop 'cannot close generated geometry module'
    write (output_unit, '(a)') &
        'Wrote '//trim(directory)//'/gliss_geometry_generated.f90'
    write (output_unit, '(a)') &
        'Wrote '//trim(directory)//'/gliss_geometry_primal.inc'

contains

    function indexed_symbol(name, index) result(value)
        character(*), intent(in) :: name
        integer, intent(in) :: index
        type(expr_t) :: value
        character(len=16) :: number

        write (number, '(i0)') index
        value = sym(arena, name//'_'//trim(number))
    end function indexed_symbol

    function dot3(left, right) result(value)
        type(expr_t), intent(in) :: left(3), right(3)
        type(expr_t) :: value

        value = left(1)*right(1) + left(2)*right(2) + left(3)*right(3)
    end function dot3

    function cross3(left, right) result(value)
        type(expr_t), intent(in) :: left(3), right(3)
        type(expr_t) :: value(3)

        value(1) = left(2)*right(3) - left(3)*right(2)
        value(2) = left(3)*right(1) - left(1)*right(3)
        value(3) = left(1)*right(2) - left(2)*right(1)
    end function cross3

    subroutine prove(name, residual)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: residual

        call check_identity(checks, engine, name, residual)
        if (checks%failed /= 0) error stop 'geometry identity failed or undecided'
    end subroutine prove

    subroutine emit_family(name, inputs, values, source, frozen_orientation)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: inputs(:), values(:)
        type(str_t), intent(out) :: source(3)
        logical, intent(in), optional :: frozen_orientation
        type(expr_t) :: tangent(size(inputs)), cotangent(size(values))
        type(expr_t) :: jvp(size(values)), vjp(size(inputs))
        type(expr_t) :: jacobian(size(values), size(inputs)), residual
        type(expr_t) :: forward_roots(2*size(values))
        type(expr_t) :: reverse_roots(size(values) + size(inputs))
        type(kernel_spec_t) :: spec
        character(len=80) :: check_label
        integer :: input, output
        logical :: frozen

        frozen = .false.
        if (present(frozen_orientation)) frozen = frozen_orientation
        do input = 1, size(inputs)
            tangent(input) = indexed_symbol('direction', input)
        end do
        do output = 1, size(values)
            cotangent(output) = indexed_symbol('cotangent', output)
            jvp(output) = 0*inputs(1)
            do input = 1, size(inputs)
                jacobian(output, input) = diff(values(output), inputs(input))
                jvp(output) = jvp(output) &
                    + jacobian(output, input)*tangent(input)
            end do
        end do
        do input = 1, size(inputs)
            vjp(input) = 0*inputs(1)
            do output = 1, size(values)
                vjp(input) = vjp(input) &
                    + cotangent(output)*jacobian(output, input)
            end do
        end do
        residual = 0*inputs(1)
        do output = 1, size(values)
            residual = residual + cotangent(output)*jvp(output)
        end do
        do input = 1, size(inputs)
            residual = residual - tangent(input)*vjp(input)
        end do
        call prove(name//' JVP/VJP duality', residual)
        do output = 1, size(values)
            do input = 1, size(inputs)
                write (check_label, '(a,a,i0,a,i0)') name, &
                    ' mixed derivative ', output, ',', input
                call prove(trim(check_label), &
                    diff(jvp(output), tangent(input)) &
                    - diff(vjp(input), cotangent(output)))
            end do
        end do
        call configure(spec, name, size(inputs), size(values), 0, frozen)
        source(1) = verified_emission(values, spec, .false.)
        forward_roots(1:size(values)) = values
        forward_roots(size(values) + 1:) = jvp
        call configure(spec, name//'_jvp', size(inputs), size(values), 1, frozen)
        source(2) = verified_emission(forward_roots, spec, .false.)
        reverse_roots(1:size(values)) = values
        reverse_roots(size(values) + 1:) = vjp
        call configure(spec, name//'_vjp', size(inputs), size(values), 2, frozen)
        source(3) = verified_emission(reverse_roots, spec, .false.)
    end subroutine emit_family

    subroutine emit_staged_field(inputs, values, source)
        ! Explicit output dependencies preserve the original contra -> cov ->
        ! norm evaluation range. CSE must not hoist a dependent expression above
        ! its output assignment. Strict substitution identities below establish
        ! equivalence to the source map and native FortSym differentiation.
        type(expr_t), intent(in) :: inputs(7), values(5)
        type(str_t), intent(out) :: source(3)
        type(expr_t) :: y(5), d(5), tangent(7), cotangent(5), half
        type(expr_t) :: primal(5), forward(10), reverse(12)
        type(expr_t) :: native_jvp(5), native_vjp(7), expanded, partial
        type(expr_t) :: contra_bar(2), cov_bar(2), total_bar(2), residual
        type(kernel_spec_t) :: spec
        character(len=80) :: check_label
        integer :: input, output

        half = rat(arena, 1_int64, 2_int64)
        do output = 1, 5
            y(output) = indexed_symbol('values', output)
            d(output) = indexed_symbol('derivative', output)
            cotangent(output) = indexed_symbol('cotangent', output)
        end do
        do input = 1, 7
            tangent(input) = indexed_symbol('direction', input)
        end do
        primal(1:2) = values(1:2)
        primal(3) = inputs(1)*y(1) + inputs(2)*y(2)
        primal(4) = inputs(2)*y(1) + inputs(3)*y(2)
        primal(5) = sqrt(y(1)*y(3) + y(2)*y(4))
        forward(1:5) = primal
        forward(6) = -inputs(7)*(tangent(6)/inputs(4)) &
            - inputs(6)*(tangent(7)/inputs(4)) &
            - y(1)*(tangent(4)/inputs(4))
        forward(7) = -tangent(5)/inputs(4) - y(2)*(tangent(4)/inputs(4))
        forward(8) = tangent(1)*y(1) + inputs(1)*d(1) &
            + tangent(2)*y(2) + inputs(2)*d(2)
        forward(9) = tangent(2)*y(1) + inputs(2)*d(1) &
            + tangent(3)*y(2) + inputs(3)*d(2)
        forward(10) = half*(d(1)*(y(3)/y(5)) + d(3)*(y(1)/y(5)) &
            + d(2)*(y(4)/y(5)) + d(4)*(y(2)/y(5)))
        contra_bar(1) = cotangent(1) + half*cotangent(5)*(y(3)/y(5))
        contra_bar(2) = cotangent(2) + half*cotangent(5)*(y(4)/y(5))
        cov_bar(1) = cotangent(3) + half*cotangent(5)*(y(1)/y(5))
        cov_bar(2) = cotangent(4) + half*cotangent(5)*(y(2)/y(5))
        total_bar(1) = contra_bar(1) + inputs(1)*cov_bar(1) &
            + inputs(2)*cov_bar(2)
        total_bar(2) = contra_bar(2) + inputs(2)*cov_bar(1) &
            + inputs(3)*cov_bar(2)
        reverse(1:5) = primal
        reverse(6) = cov_bar(1)*y(1)
        reverse(7) = cov_bar(1)*y(2) + cov_bar(2)*y(1)
        reverse(8) = cov_bar(2)*y(2)
        reverse(9) = -(total_bar(1)*y(1) + total_bar(2)*y(2))/inputs(4)
        reverse(10) = -total_bar(2)/inputs(4)
        reverse(11) = -total_bar(1)*(inputs(7)/inputs(4))
        reverse(12) = -total_bar(1)*(inputs(6)/inputs(4))
        native_jvp = rat(arena, 0_int64, 1_int64)
        native_vjp = rat(arena, 0_int64, 1_int64)
        do output = 1, 5
            do input = 1, 7
                partial = diff(values(output), inputs(input))
                native_jvp(output) = native_jvp(output) + tangent(input)*partial
                native_vjp(input) = native_vjp(input) + cotangent(output)*partial
            end do
        end do
        residual = rat(arena, 0_int64, 1_int64)
        do output = 1, 5
            expanded = subs_many(primal(output), y, values)
            write (check_label, '(a,i0)') 'staged field primal ', output
            call prove(trim(check_label), expanded - values(output))
            expanded = subs_many(forward(5 + output), y, values)
            expanded = subs_many(expanded, d, native_jvp)
            write (check_label, '(a,i0)') 'staged field native JVP ', output
            call prove(trim(check_label), expanded - native_jvp(output))
            residual = residual + cotangent(output)*native_jvp(output)
        end do
        do input = 1, 7
            expanded = subs_many(reverse(5 + input), y, values)
            write (check_label, '(a,i0)') 'staged field native VJP ', input
            call prove(trim(check_label), expanded - native_vjp(input))
            residual = residual - tangent(input)*native_vjp(input)
        end do
        call prove('staged field JVP/VJP duality', residual)
        call configure(spec, 'generated_geometry_field', 7, 5, 0, .false.)
        spec%cse_level = CSE_NONE
        spec%ordered_output_dependencies = .true.
        source(1) = verified_emission(primal, spec, .true.)
        call configure(spec, 'generated_geometry_field_jvp', 7, 5, 1, .false.)
        spec%cse_level = CSE_NONE
        spec%ordered_output_dependencies = .true.
        source(2) = verified_emission(forward, spec, .true.)
        call configure(spec, 'generated_geometry_field_vjp', 7, 5, 2, .false.)
        spec%cse_level = CSE_NONE
        spec%ordered_output_dependencies = .true.
        source(3) = verified_emission(reverse, spec, .true.)
    end subroutine emit_staged_field

    subroutine emit_staged_normal(inputs, values, source)
        ! Project the tangent/cotangent through the unit-normal map. Keeping
        ! the normalized primal as an explicit dependency avoids area**3 and
        ! products of an unnormalized cross with a large second position jet.
        type(expr_t), intent(in) :: inputs(3), values(3)
        type(str_t), intent(out) :: source(3)
        type(expr_t) :: y(3), tangent(3), cotangent(3), area
        type(expr_t) :: forward(6), reverse(6), native_jvp(3), native_vjp(3)
        type(expr_t) :: expanded, partial, residual
        type(kernel_spec_t) :: spec
        character(len=80) :: check_label
        integer :: input, output

        area = sqrt(dot3(inputs, inputs))
        do input = 1, 3
            y(input) = indexed_symbol('values', input)
            tangent(input) = indexed_symbol('direction', input)
            cotangent(input) = indexed_symbol('cotangent', input)
        end do
        forward(1:3) = values
        reverse(1:3) = values
        do input = 1, 3
            forward(3 + input) = orientation*tangent(input)/area &
                - y(input)*(dot3(y, tangent)/(orientation*area))
            reverse(3 + input) = orientation*cotangent(input)/area &
                - y(input)*(dot3(y, cotangent)/(orientation*area))
        end do
        native_jvp = rat(arena, 0_int64, 1_int64)
        native_vjp = rat(arena, 0_int64, 1_int64)
        do output = 1, 3
            do input = 1, 3
                partial = diff(values(output), inputs(input))
                native_jvp(output) = native_jvp(output) + tangent(input)*partial
                native_vjp(input) = native_vjp(input) + cotangent(output)*partial
            end do
        end do
        call prove('signed unit normal magnitude', &
            dot3(values, values) - orientation**2)
        residual = rat(arena, 0_int64, 1_int64)
        do input = 1, 3
            expanded = subs_many(forward(3 + input), y, values)
            write (check_label, '(a,i0)') 'staged unit normal native JVP ', input
            call prove(trim(check_label), expanded - native_jvp(input))
            expanded = subs_many(reverse(3 + input), y, values)
            write (check_label, '(a,i0)') 'staged unit normal native VJP ', input
            call prove(trim(check_label), expanded - native_vjp(input))
            residual = residual + cotangent(input)*native_jvp(input) &
                - tangent(input)*native_vjp(input)
        end do
        call prove('signed unit normal JVP/VJP duality', residual)
        call configure(spec, 'generated_geometry_normal', 3, 3, 0, .true.)
        source(1) = verified_emission(values, spec, .true.)
        call configure(spec, 'generated_geometry_normal_jvp', 3, 3, 1, .true.)
        spec%cse_level = CSE_NONE
        spec%ordered_output_dependencies = .true.
        source(2) = verified_emission(forward, spec, .true.)
        call configure(spec, 'generated_geometry_normal_vjp', 3, 3, 2, .true.)
        spec%cse_level = CSE_NONE
        spec%ordered_output_dependencies = .true.
        source(3) = verified_emission(reverse, spec, .true.)
    end subroutine emit_staged_normal

    subroutine write_second_form_primal(output)
        ! Array packing and the chain of generated maps carry no handwritten
        ! arithmetic. The full packed derivatives are the composition proved
        ! above and checked independently by the native geometry oracle.
        integer, intent(in) :: output

        write (output, '(a)') &
            'pure subroutine generated_geometry_second_form(x, orientation, values)'
        call write_wrapper_declarations(output, 0)
        write (output, '(a)') '    real(dp) :: packed(12)'
        write (output, '(a)') &
            '    call generated_geometry_normal(x(1:3), orientation, packed(1:3))'
        write (output, '(a)') '    packed(4:12) = x(4:12)'
        write (output, '(a)') &
            '    call generated_geometry_normal_contraction(packed, values)'
        write (output, '(a)') 'end subroutine generated_geometry_second_form'
        write (output, '(a)') ''
    end subroutine write_second_form_primal

    subroutine write_second_form_products(output)
        integer, intent(in) :: output

        write (output, '(a)') &
            'pure subroutine generated_geometry_second_form_jvp(x, direction,' &
            //' orientation, values, derivative)'
        call write_wrapper_declarations(output, 1)
        write (output, '(a)') '    real(dp) :: packed(12), packed_direction(12)'
        write (output, '(a)') &
            '    call generated_geometry_normal_jvp(x(1:3), direction(1:3),' &
            //' orientation, &'
        write (output, '(a)') '        packed(1:3), packed_direction(1:3))'
        write (output, '(a)') '    packed(4:12) = x(4:12)'
        write (output, '(a)') '    packed_direction(4:12) = direction(4:12)'
        write (output, '(a)') &
            '    call generated_geometry_normal_contraction_jvp(packed,' &
            //' packed_direction, values, derivative)'
        write (output, '(a)') 'end subroutine generated_geometry_second_form_jvp'
        write (output, '(a)') ''
        write (output, '(a)') &
            'pure subroutine generated_geometry_second_form_vjp(x, cotangent,' &
            //' orientation, values, derivative)'
        call write_wrapper_declarations(output, 2)
        write (output, '(a)') '    real(dp) :: packed(12), packed_cotangent(12)'
        write (output, '(a)') &
            '    call generated_geometry_normal(x(1:3), orientation, packed(1:3))'
        write (output, '(a)') '    packed(4:12) = x(4:12)'
        write (output, '(a)') &
            '    call generated_geometry_normal_contraction_vjp(packed,' &
            //' cotangent, values, packed_cotangent)'
        write (output, '(a)') &
            '    call generated_geometry_normal_vjp(x(1:3), packed_cotangent(1:3),' &
            //' orientation, &'
        write (output, '(a)') '        packed(1:3), derivative(1:3))'
        write (output, '(a)') '    derivative(4:12) = packed_cotangent(4:12)'
        write (output, '(a)') 'end subroutine generated_geometry_second_form_vjp'
    end subroutine write_second_form_products

    subroutine write_wrapper_declarations(output, product)
        integer, intent(in) :: output, product

        write (output, '(a)') &
            '    use, intrinsic :: iso_fortran_env, only: dp => real64'
        write (output, '(a)') '    implicit none'
        write (output, '(a)') '    real(dp), intent(in) :: x(12), orientation'
        write (output, '(a)') '    real(dp), intent(out) :: values(3)'
        if (product == 1) then
            write (output, '(a)') '    real(dp), intent(in) :: direction(12)'
            write (output, '(a)') '    real(dp), intent(out) :: derivative(3)'
        else if (product == 2) then
            write (output, '(a)') '    real(dp), intent(in) :: cotangent(3)'
            write (output, '(a)') '    real(dp), intent(out) :: derivative(12)'
        end if
    end subroutine write_wrapper_declarations

    subroutine configure(spec, name, inputs, outputs, product, frozen)
        type(kernel_spec_t), intent(out) :: spec
        character(*), intent(in) :: name
        integer, intent(in) :: inputs, outputs, product
        logical, intent(in) :: frozen
        character(len=32) :: input_shape, output_shape, index
        integer :: n_arguments, n_outputs, n_roots, i

        write (input_shape, '(a,i0,a)') '(', inputs, ')'
        write (output_shape, '(a,i0,a)') '(', outputs, ')'
        spec = kernel_spec_t()
        spec%name = str(name)
        spec%pure_procedure = .true.
        spec%temp_prefix = str('t')
        spec%generator = str('gen/gliss_geometry_gen.f90')
        spec%generator_revision = str(trim(revision))
        spec%regenerate_command = &
            str('cmake --build build-codegen --target gliss_regenerate_kernels')
        n_arguments = 1
        if (product /= 0) n_arguments = n_arguments + 1
        if (frozen) n_arguments = n_arguments + 1
        allocate (spec%args(n_arguments), spec%arg_shapes(n_arguments))
        spec%args(1) = str('x')
        spec%arg_shapes(1) = str(trim(input_shape))
        if (product == 1) then
            spec%args(2) = str('direction')
            spec%arg_shapes(2) = str(trim(input_shape))
        else if (product == 2) then
            spec%args(2) = str('cotangent')
            spec%arg_shapes(2) = str(trim(output_shape))
        end if
        if (frozen) then
            spec%args(n_arguments) = str('orientation')
            spec%arg_shapes(n_arguments) = str('')
        end if
        n_outputs = 1
        n_roots = outputs
        if (product /= 0) then
            n_outputs = 2
            if (product == 1) n_roots = n_roots + outputs
            if (product == 2) n_roots = n_roots + inputs
        end if
        allocate (spec%outputs(n_outputs), spec%output_shapes(n_outputs))
        allocate (spec%output_references(n_roots))
        spec%outputs(1) = str('values')
        spec%output_shapes(1) = str(trim(output_shape))
        do i = 1, outputs
            write (index, '(a,i0,a)') 'values(', i, ')'
            spec%output_references(i) = str(trim(index))
        end do
        if (product == 0) return
        spec%outputs(2) = str('derivative')
        if (product == 1) spec%output_shapes(2) = str(trim(output_shape))
        if (product == 2) spec%output_shapes(2) = str(trim(input_shape))
        do i = 1, n_roots - outputs
            write (index, '(a,i0,a)') 'derivative(', i, ')'
            spec%output_references(outputs + i) = str(trim(index))
        end do
    end subroutine configure

    function verified_emission(roots, spec, preserve_factors) result(source)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t), intent(in) :: spec
        logical, intent(in) :: preserve_factors
        type(str_t) :: source
        type(expr_t) :: simplified_roots(size(roots))
        type(expr_t) :: emitted_roots(size(roots)), restored
        type(expr_t), allocatable :: symbols(:), references(:)
        type(engine_result_t) :: simplified
        character(:), allocatable :: message
        character(len=80) :: check_label
        character(len=32) :: symbol_name, reference_name, shape
        integer :: root, argument, element, count, total, first, last
        logical :: good

        do root = 1, size(roots)
            simplified = engine%simplify(roots(root))
            if (.not. simplified%ok) error stop 'geometry simplification failed'
            simplified_roots(root) = simplified%value
            write (check_label, '(a,a,i0)') chars(spec%name), ' root ', root
            call prove(trim(check_label), roots(root) - simplified_roots(root))
            if (preserve_factors) simplified_roots(root) = roots(root)
        end do
        ! Keep algebraic symbol identities independent of Fortran subscripting.
        ! Only the verified emission boundary maps x_1 to the array reference
        ! x(1); an inverse substitution proves that no mathematical root changes.
        total = 0
        do argument = 1, size(spec%args)
            if (len(chars(spec%arg_shapes(argument))) == 0) cycle
            shape = chars(spec%arg_shapes(argument))
            read (shape(2:len_trim(shape) - 1), *) count
            total = total + count
        end do
        do argument = 1, size(spec%outputs)
            shape = chars(spec%output_shapes(argument))
            read (shape(2:len_trim(shape) - 1), *) count
            total = total + count
        end do
        allocate (symbols(total), references(total))
        first = 1
        do argument = 1, size(spec%args)
            if (len(chars(spec%arg_shapes(argument))) == 0) cycle
            shape = chars(spec%arg_shapes(argument))
            read (shape(2:len_trim(shape) - 1), *) count
            last = first + count - 1
            do element = 1, count
                write (symbol_name, '(a,a,i0)') chars(spec%args(argument)), &
                    '_', element
                write (reference_name, '(a,a,i0,a)') chars(spec%args(argument)), &
                    '(', element, ')'
                symbols(first + element - 1) = sym(arena, trim(symbol_name))
                references(first + element - 1) = sym(arena, trim(reference_name))
            end do
            first = last + 1
        end do
        do argument = 1, size(spec%outputs)
            shape = chars(spec%output_shapes(argument))
            read (shape(2:len_trim(shape) - 1), *) count
            last = first + count - 1
            do element = 1, count
                write (symbol_name, '(a,a,i0)') chars(spec%outputs(argument)), &
                    '_', element
                write (reference_name, '(a,a,i0,a)') chars(spec%outputs(argument)), &
                    '(', element, ')'
                symbols(first + element - 1) = sym(arena, trim(symbol_name))
                references(first + element - 1) = sym(arena, trim(reference_name))
            end do
            first = last + 1
        end do
        do root = 1, size(roots)
            emitted_roots(root) = subs_many(simplified_roots(root), &
                symbols, references)
            restored = subs_many(emitted_roots(root), references, symbols)
            write (check_label, '(a,a,i0)') chars(spec%name), ' packing ', root
            call prove(trim(check_label), restored - simplified_roots(root))
        end do
        source = emit_kernel(emitted_roots, spec, ok=good, message=message)
        if (.not. good) then
            write (output_unit, '(a)') message
            error stop 'geometry emission failed'
        end if
    end function verified_emission

    subroutine write_header(output)
        integer, intent(in) :: output

        write (output, '(a)') '! Generated by gen/gliss_geometry_gen.f90; do not edit.'
        write (output, '(a)') '! FortSym revision: '//trim(revision)
        write (output, '(a)') &
            '! Inputs are pointwise Cartesian jets in metres; angles/s dimensionless.'
        write (output, '(a)') &
            '! Polynomial x(27): r_s,r_t,r_z,r_ss,r_st,r_sz,r_tt,r_tz,r_zz (3 each).'
        write (output, '(a)') &
            '! Polynomial y(13): g_ss,g_st,g_sz,g_tt,g_tz,g_zz,J,J_s,J_t,J_z,cross(3).'
        write (output, '(a)') &
            '! Field x(7): g_tt,g_tz,g_zz,J,Phi_s,chi_s,inverse_nfp; fluxes in Wb.'
        write (output, '(a)') '! Field y(5): B^theta,B^zeta,B_theta,B_zeta,|B|' &
            //' (SI coordinate components).'
        write (output, '(a)') '! Curvature x(12): angular_cross,r_tt,r_tz,r_zz' &
            //' (3 each); y(3): II_tt,II_tz,II_zz.'
        write (output, '(a)') '! Curvature orientation is sign(J), frozen' &
            //' during JVP/VJP in either nonzero chart.'
        write (output, '(a)') '! Domain: nonzero J and angular area,' &
            //' positive angular metric and nonzero B.'
        write (output, '(a)') '! Field JVP/VJP include the real inverse_nfp input;' &
            //' physical nfp is discrete.'
        write (output, '(a)') '! Callers set its tangent to zero and do not' &
            //' interpret its cotangent as d/dnfp.'
        write (output, '(a)') 'module gliss_geometry_generated'
        write (output, '(a)') '    implicit none'
        write (output, '(a)') '    private'
        write (output, '(a)') '    public :: generated_geometry_polynomial'
        write (output, '(a)') '    public :: generated_geometry_polynomial_jvp,' &
            //' generated_geometry_polynomial_vjp'
        write (output, '(a)') '    public :: generated_geometry_field'
        write (output, '(a)') '    public :: generated_geometry_field_jvp,' &
            //' generated_geometry_field_vjp'
        write (output, '(a)') '    public :: generated_geometry_second_form'
        write (output, '(a)') '    public :: generated_geometry_second_form_jvp,' &
            //' generated_geometry_second_form_vjp'
        write (output, '(a)') 'contains'
    end subroutine write_header

end program gliss_geometry_gen
