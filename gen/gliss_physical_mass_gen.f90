program gliss_physical_mass_gen
    ! Physical kinetic energy in GLISS's (xi^s, eta, mu) variables.
    ! Source: Nuehrenberg, Phys. Plasmas 6 (1999), Appendix A1, A4--A8;
    ! GLISS's flux, eta and signed-Jacobian conventions: PROVENANCE.md.
    ! The three expressions below are orthonormal displacement components.
    ! Their squared norm defines the energy; its Hessian defines the mass.
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
        'bdeb96e0777b1eafcb7c3b4845ca2ba9362a2b00'
    character(len=24), parameter :: names(13) = [character(len=24) :: &
        'flux_t_slope', 'flux_p_slope', 'current_i', 'current_j', &
        'signed_sqrtg', 'bmag', 'grad_s2', 'signed_sigma_tilde', &
        'beta_tilde', 'density_kg_m3', 'xi_s', 'eta', 'mu']
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(suite_t) :: checks
    type(expr_t) :: x(13), direction(13), cotangent, components(3), flux_norm
    type(expr_t) :: energy, mass(3, 3), gradient(13), forward, reverse(13)
    type(expr_t) :: matrix_roots(9), quadratic, linear
    type(kernel_spec_t) :: spec
    type(str_t) :: sources(4)
    character(len=1024) :: directory, evidence
    character(len=40) :: revision
    character(len=40) :: label
    integer :: i, j, k, unit, status

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
    engine = make_symengine_engine(arena)
    call suite_begin(checks, 'GLISS physical mass exact identities')
    do i = 1, size(x)
        x(i) = sym(arena, trim(names(i)))
        direction(i) = sym(arena, 'd_'//trim(names(i)))
    end do
    cotangent = sym(arena, 'cotangent')

    ! Domain: J /= 0, B > 0, grad(s)^2 > 0, FT'^2 + FP'^2 > 0, rho > 0.
    ! Both signs of J are retained. The derivative of |J| is valid on each
    ! nonzero orientation branch; crossing J = 0 is outside this domain.
    flux_norm = x(1)**2 + x(2)**2
    components(1) = x(11)/sqrt(x(7))
    components(2) = x(8)*x(11)/sqrt(x(7)) + sqrt(x(7))*x(12)/x(6)
    components(3) = x(9)*x(11)/x(6) &
        + (x(4)*x(1) - x(3)*x(2))*x(12)/(x(6)*flux_norm) &
        + x(6)*x(13)/flux_norm
    energy = rat(arena, 1_int64, 2_int64)*x(10)*abs(x(5)) &
        * (components(1)**2 + components(2)**2 + components(3)**2)
    do j = 1, 3
        do i = 1, 3
            mass(i, j) = diff(diff(energy, x(10 + i)), x(10 + j))
        end do
    end do
    do i = 1, 13
        gradient(i) = diff(energy, x(i))
        reverse(i) = cotangent*gradient(i)
    end do
    forward = gradient(1)*direction(1)
    do i = 2, 13
        forward = forward + gradient(i)*direction(i)
    end do

    quadratic = 0*energy
    do j = 1, 3
        do i = 1, 3
            quadratic = quadratic &
                + rat(arena, 1_int64, 2_int64)*x(10 + i)*mass(i, j)*x(10 + j)
            write (label, '(a,i0,a,i0)') 'mass symmetry ', i, ',', j
            call check_identity(checks, engine, trim(label), mass(i, j) - mass(j, i))
        end do
    end do
    call check_identity(checks, engine, 'energy equals mass quadratic', energy - quadratic)
    call check_identity(checks, engine, 'linear density scaling', x(10)*gradient(10) - energy)
    do i = 1, 3
        linear = 0*energy
        do j = 1, 3
            linear = linear + mass(i, j)*x(10 + j)
        end do
        write (label, '(a,i0)') 'displacement derivative ', i
        call check_identity(checks, engine, trim(label), gradient(10 + i) - linear)
    end do
    linear = reverse(1)*direction(1)
    do i = 2, 13
        linear = linear + reverse(i)*direction(i)
    end do
    call check_identity(checks, engine, 'JVP/VJP duality', cotangent*forward - linear)

    k = 0
    do j = 1, 3
        do i = 1, 3
            k = k + 1
            matrix_roots(k) = mass(i, j)
        end do
    end do
    call configure(spec, 'generated_physical_mass_matrix', 10, .false., .false.)
    allocate (spec%outputs(1), spec%output_shapes(1), spec%output_references(9))
    spec%outputs(1) = str('mass')
    spec%output_shapes(1) = str('(3,3)')
    k = 0
    do j = 1, 3
        do i = 1, 3
            k = k + 1
            write (label, '(a,i0,a,i0,a)') 'mass(', i, ',', j, ')'
            spec%output_references(k) = str(trim(label))
        end do
    end do
    sources(1) = verified_emission(matrix_roots, spec)

    call configure(spec, 'generated_physical_mass_energy', 13, .false., .false.)
    allocate (spec%outputs(1))
    spec%outputs(1) = str('energy')
    sources(2) = verified_emission([energy], spec)

    call configure(spec, 'generated_physical_mass_jvp', 13, .true., .false.)
    allocate (spec%outputs(2))
    spec%outputs(1) = str('energy')
    spec%outputs(2) = str('derivative')
    sources(3) = verified_emission([energy, forward], spec)

    call configure(spec, 'generated_physical_mass_vjp', 13, .false., .true.)
    allocate (spec%outputs(2), spec%output_shapes(2), spec%output_references(14))
    spec%outputs(1) = str('energy')
    spec%outputs(2) = str('gradient')
    spec%output_shapes(1) = str('')
    spec%output_shapes(2) = str('(13)')
    spec%output_references(1) = str('energy')
    do i = 1, 13
        write (label, '(a,i0,a)') 'gradient(', i, ')'
        spec%output_references(i + 1) = str(trim(label))
    end do
    sources(4) = verified_emission([energy, reverse], spec)

    ! Fail before writing production output if any exact identity was undecided.
    if (len_trim(evidence) > 0) then
        call suite_end(checks, trim(evidence))
    else
        call suite_end(checks)
    end if
    open (newunit=unit, file=trim(directory)//'/gliss_physical_mass_generated.f90', &
        status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'cannot write generated physical mass module'
    write (unit, '(a)') '! Generated by gen/gliss_physical_mass_gen.f90; do not edit.'
    write (unit, '(a)') '! FortSym revision: '//revision
    write (unit, '(a)') '! Domain: J /= 0, B > 0, grad_s2 > 0, FT^2+FP^2 > 0, rho > 0.'
    write (unit, '(a)') 'module gliss_physical_mass_generated'
    write (unit, '(a)') '    implicit none'
    write (unit, '(a)') '    private'
    write (unit, '(a)') '    public :: generated_physical_mass_matrix, generated_physical_mass_energy'
    write (unit, '(a)') '    public :: generated_physical_mass_jvp, generated_physical_mass_vjp'
    write (unit, '(a)') 'contains'
    do i = 1, size(sources)
        write (unit, '(a)') chars(sources(i))
    end do
    write (unit, '(a)') 'end module gliss_physical_mass_generated'
    close (unit)
    write (output_unit, '(a)') 'Wrote '//trim(directory)//'/gliss_physical_mass_generated.f90'

contains

    subroutine configure(configuration, name, inputs, tangent, adjoint)
        type(kernel_spec_t), intent(out) :: configuration
        character(*), intent(in) :: name
        integer, intent(in) :: inputs
        logical, intent(in) :: tangent, adjoint
        integer :: argument_count, argument

        configuration%name = str(name)
        configuration%mode = KERNEL_SUBROUTINE
        configuration%pure_procedure = .true.
        configuration%temp_prefix = str('t')
        configuration%generator = str('gen/gliss_physical_mass_gen.f90')
        configuration%generator_revision = str(revision)
        configuration%regenerate_command = str('gliss_physical_mass_gen src/generated')
        argument_count = inputs
        if (tangent) argument_count = argument_count + inputs
        if (adjoint) argument_count = argument_count + 1
        allocate (configuration%args(argument_count))
        do argument = 1, inputs
            configuration%args(argument) = str(trim(names(argument)))
            if (tangent) configuration%args(inputs + argument) = &
                str('d_'//trim(names(argument)))
        end do
        if (adjoint) configuration%args(argument_count) = str('cotangent')
    end subroutine configure

    function verified_emission(roots, configuration) result(source)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t), intent(in) :: configuration
        type(str_t) :: source
        type(expr_t) :: simplified_roots(size(roots))
        type(engine_result_t) :: simplified
        character(:), allocatable :: message
        character(len=64) :: check_label
        integer :: root
        logical :: good

        do root = 1, size(roots)
            simplified = engine%simplify(roots(root))
            if (.not. simplified%ok) error stop 'native simplification failed'
            simplified_roots(root) = simplified%value
            write (check_label, '(a,a,i0)') chars(configuration%name), ' root ', root
            call check_identity(checks, engine, trim(check_label), &
                roots(root) - simplified_roots(root))
        end do
        source = emit_kernel(simplified_roots, configuration, ok=good, message=message)
        if (.not. good) then
            write (output_unit, '(a)') message
            error stop 'physical mass emission failed'
        end if
    end function verified_emission

end program gliss_physical_mass_gen
