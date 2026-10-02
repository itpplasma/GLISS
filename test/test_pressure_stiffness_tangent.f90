program test_pressure_stiffness_tangent
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_three_component_problem, only: &
        build_compatible_three_component_problem, &
        build_compatible_pressure_stiffness_tangent, &
        compatible_three_component_problem_t, compatible_three_component_ok, &
        compatible_three_component_invalid
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use variable_block_tridiagonal, only: variable_block_to_dense, variable_block_ok
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    real(dp), allocatable :: direction(:)
    character(len=1024) :: directory
    integer :: info, parity, mesh, storage

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) directory = 'test/data'
    call read_gvec_cas3d_file(trim(directory) // '/solovev_q1.045.nc', &
        equilibrium, info)
    call require(info == reader_ok, 'pressure fixture read failed')
    direction = 0.03_dp * equilibrium%pressure * (1.0_dp + equilibrium%s)
    do parity = 0, 2
        do mesh = 1, 2
            do storage = 0, 1
                call check_direction(parity, 6 * mesh, storage == 1)
            end do
        end do
    end do
    call check_duality()
    call check_invalid()
    print *, 'PASS'

contains

    subroutine assemble(data, parity, cells, sparse, matrix, mass, tangent)
        type(gvec_cas3d_equilibrium_t), intent(in) :: data
        integer, intent(in) :: parity, cells
        logical, intent(in) :: sparse
        real(dp), allocatable, intent(out) :: matrix(:, :), mass(:, :)
        real(dp), optional, intent(in) :: tangent(:)
        type(compatible_three_component_problem_t) :: problem
        integer :: status

        if (present(tangent)) then
            call build_compatible_pressure_stiffness_tangent(data, 1.4_dp, 1.0_dp, &
                [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], parity, 2, 64, 8, &
                tangent, problem, status, sparse_storage=sparse, radial_cells=cells)
        else
            call build_compatible_three_component_problem(data, 1.4_dp, 1.0_dp, &
                [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], parity, 2, 64, 8, &
                problem, status, sparse_storage=sparse, radial_cells=cells)
        end if
        call require(status == compatible_three_component_ok, &
            'pressure assembly failed')
        if (sparse) then
            call variable_block_to_dense(problem%sparse_stiffness, matrix, status)
            call require(status == variable_block_ok, &
                'sparse pressure derivative failed')
            call variable_block_to_dense(problem%sparse_mass, mass, status)
            call require(status == variable_block_ok, 'sparse pressure mass failed')
        else
            call move_alloc(problem%stiffness, matrix)
            call move_alloc(problem%mass, mass)
        end if
    end subroutine assemble

    subroutine check_direction(parity, cells, sparse)
        integer, intent(in) :: parity, cells
        logical, intent(in) :: sparse
        real(dp), parameter :: steps(3) = [1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp]
        type(gvec_cas3d_equilibrium_t) :: plus, minus
        real(dp), allocatable :: derivative(:, :), mass(:, :), kp(:, :), km(:, :)
        real(dp), allocatable :: mp(:, :), mm(:, :)
        real(dp) :: error, scale, step
        integer :: index

        call assemble(equilibrium, parity, cells, sparse, derivative, mass, direction)
        scale = maxval(abs(derivative))
        call require(scale > 0.0_dp, 'pressure tangent is identically zero')
        plus = equilibrium
        minus = equilibrium
        do index = 1, size(steps)
            step = steps(index)
            plus%pressure = equilibrium%pressure + step * direction
            minus%pressure = equilibrium%pressure - step * direction
            call assemble(plus, parity, cells, sparse, kp, mp)
            call assemble(minus, parity, cells, sparse, km, mm)
            error = maxval(abs((kp - km) / (2.0_dp * step) - derivative)) / scale
            print '(a,i1,a,i2,a,l1,a,es10.2)', 'parity=', parity, ' cells=', cells, &
                ' sparse=', sparse, ' relative FD=', error
            call require(error < 2.0e-6_dp, 'full pressure K tangent disagrees with FD')
            call require(maxval(abs(mp - mass)) == 0.0_dp, 'pressure changed mass')
            call require(maxval(abs(mm - mass)) == 0.0_dp, 'pressure changed mass')
        end do
    end subroutine check_direction

    subroutine check_duality()
        real(dp), allocatable :: derivative(:, :), mass(:, :), vector(:), basis(:)
        real(dp), allocatable :: gradient(:)
        real(dp) :: jvp, vjp
        integer :: index

        call assemble(equilibrium, 1, 6, .false., derivative, mass, direction)
        allocate (vector(size(derivative, 1)), basis(size(direction)), &
            gradient(size(direction)))
        do index = 1, size(vector)
            vector(index) = sin(0.37_dp * real(index, dp))
        end do
        jvp = quadratic(vector, derivative)
        do index = 1, size(direction)
            basis = 0.0_dp
            basis(index) = 1.0_dp
            call assemble(equilibrium, 1, 6, .false., derivative, mass, basis)
            gradient(index) = quadratic(vector, derivative)
        end do
        vjp = dot_product(gradient, direction)
        call require(abs(jvp - vjp) < 1.0e-10_dp * max(abs(jvp), abs(vjp)), &
            'pressure sample JVP/VJP duality failed')
    end subroutine check_duality

    pure function quadratic(vector, matrix) result(value)
        real(dp), intent(in) :: vector(:), matrix(:, :)
        real(dp) :: value
        integer :: a, b

        value = 0.0_dp
        do b = 1, size(vector)
            do a = 1, size(vector)
                value = value + vector(a) * matrix(a, b) * vector(b)
            end do
        end do
    end function quadratic

    subroutine check_invalid()
        type(gvec_cas3d_equilibrium_t) :: invalid
        type(compatible_three_component_problem_t) :: problem
        real(dp), allocatable :: tangent(:)
        integer :: status

        call build_compatible_pressure_stiffness_tangent(equilibrium, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], 1, &
            2, 32, 8, direction(:size(direction) - 1), problem, status)
        call require(status == compatible_three_component_invalid, &
            'wrong pressure direction shape accepted')
        tangent = direction
        tangent(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        call build_compatible_pressure_stiffness_tangent(equilibrium, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], 1, &
            2, 32, 8, tangent, problem, status)
        call require(status == compatible_three_component_invalid, &
            'nonfinite pressure direction accepted')
        invalid = equilibrium
        invalid%pressure(2) = 0.0_dp
        call build_compatible_pressure_stiffness_tangent(invalid, 1.4_dp, &
            1.0_dp, [0, 1, 2], [1, 1, 1], [0.0_dp, 0.5_dp, 0.0_dp], 1, &
            2, 32, 8, direction, problem, status)
        call require(status == compatible_three_component_invalid, &
            'closed pressure parameter domain accepted')
    end subroutine check_invalid

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_pressure_stiffness_tangent
