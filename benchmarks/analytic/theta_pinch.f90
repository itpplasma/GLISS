program benchmark_theta_pinch
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use exact_cylinder, only: assemble_exact_cylinder
    use cylinder_fixture, only: create_cylinder_fixture
    use compatible_three_component_problem, only: &
        build_compatible_three_component_problem, &
        compatible_three_component_problem_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use symmetric_eigensolver, only: solve_symmetric_generalized
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(compatible_three_component_problem_t) :: problem
    real(dp), allocatable :: values(:), vectors(:, :), gram(:, :), residual(:)
    real(dp), allocatable :: mass_vectors(:, :), stiffness_vectors(:, :)
    real(dp) :: residual_norm, scale, orthogonality
    logical :: exact
    integer :: ns, degree, parity, info, i, unit, mesh, j, dump_unit
    integer, parameter :: meshes(3) = [8, 16, 32]
    character(len=1024) :: filename, option, dump_directory, dump_file

    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) error stop 'Supply output CSV path'
    exact = command_argument_count() == 1
    dump_directory = ''
    call get_command_argument(2, option)
    if (trim(option) == 'dump') then
        exact = .true.
        call get_command_argument(3, dump_directory)
        if (len_trim(dump_directory) == 0) error stop 'Supply dump directory'
    end if
    open (newunit=unit, file=trim(filename), status='replace')
    write (unit, '(a)') 'surfaces,degree,parity,index,omega2_s_minus2,' &
        //'relative_residual,mass_orthogonality'
    do mesh = 1, size(meshes)
        ns = meshes(mesh)
        if (.not. exact) then
            call create_cylinder_fixture('theta_pinch.nc', surfaces=ns, &
                poloidal_scale=0.0_dp)
            call read_gvec_cas3d_file('theta_pinch.nc', equilibrium, info)
            if (info /= 0) error stop 'Fixture read failed'
        end if
        do degree = 1, 4
            do parity = 1, 2
                if (exact) then
                    call assemble_exact_cylinder(ns, degree, parity, &
                        problem%stiffness, problem%mass)
                    info = 0
                else
                    call build_compatible_three_component_problem(equilibrium, &
                        5.0_dp/3.0_dp, 2.0_dp, [3], [1], [-0.5_dp], &
                        parity, degree, 16, 16, problem, info)
                end if
                if (info /= 0) error stop 'Production assembly failed'
                if (len_trim(dump_directory) > 0 .and. parity == 1) then
                    write (dump_file, '(a,a,i0,a,i0,a)') trim(dump_directory), &
                        '/p', degree, '_n', ns, '.dat'
                    open (newunit=dump_unit, file=trim(dump_file), status='replace')
                    write (dump_unit, '(i0)') size(problem%mass, 1)
                    do j = 1, size(problem%mass, 1)
                        do i = 1, size(problem%mass, 1)
                            write (dump_unit, '(2es26.17e3)') &
                                problem%stiffness(i, j), problem%mass(i, j)
                        end do
                    end do
                    close (dump_unit)
                end if
                call solve_symmetric_generalized(problem%stiffness, &
                    problem%mass, values, vectors, info)
                if (info /= 0) error stop 'Spectrum solve failed'
                mass_vectors = matmul(problem%mass, vectors)
                stiffness_vectors = matmul(problem%stiffness, vectors)
                if (allocated(gram)) deallocate (gram)
                if (allocated(residual)) deallocate (residual)
                allocate (gram(size(values), size(values)), residual(size(values)))
                do j = 1, size(values)
                    do i = 1, size(values)
                        gram(i, j) = dot_product(vectors(:, i), mass_vectors(:, j))
                    end do
                    gram(j, j) = gram(j, j) - 1.0_dp
                end do
                orthogonality = maxval(abs(gram))
                do i = 1, size(values)
                    do j = 1, size(values)
                        residual(j) = stiffness_vectors(j, i) &
                            - values(i) * mass_vectors(j, i)
                    end do
                    scale = (norm2(problem%stiffness) &
                        + abs(values(i))*norm2(problem%mass)) &
                        *norm2(vectors(:, i))
                    residual_norm = norm2(residual)/scale
                    write (unit, '(4(i0,a),2(es25.16,a),es25.16)') ns, ',', degree, ',', &
                        parity, ',', i, ',', values(i), ',', residual_norm, &
                        ',', orthogonality
                end do
                write (*, '(a,3i5)') 'Completed surfaces, degree, parity:', &
                    ns, degree, parity
            end do
        end do
    end do
    close (unit)
    if (.not. exact) then
        open (newunit=unit, file='theta_pinch.nc', status='old')
        close (unit, status='delete')
    end if
end program benchmark_theta_pinch
