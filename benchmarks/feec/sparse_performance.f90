program sparse_performance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime, omp_get_num_threads
    use axisymmetric_spectrum, only: build_axisymmetric_mode_table
    use compatible_two_component_problem, only: &
        build_compatible_two_component_problem, compatible_two_component_problem_t
    use gvec_cas3d_reader, only: read_gvec_cas3d_file
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use marginality_spectrum, only: solve_compatible_marginality_problem, &
        marginality_spectrum_result_t
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(compatible_two_component_problem_t) :: problem
    type(marginality_spectrum_result_t) :: result
    integer, allocatable :: mode_m(:), mode_n(:)
    real(dp), allocatable :: powers(:)
    character(len=1024) :: filename, value, message
    real(dp) :: start, assembled, finished
    integer :: info, cells, degree, threads, storage

    if (command_argument_count() /= 4) &
        error stop 'usage: sparse_performance EXPORT CELLS DEGREE SPARSE(0|1)'
    call get_command_argument(1, filename)
    call get_command_argument(2, value)
    read (value, *) cells
    call get_command_argument(3, value)
    read (value, *) degree
    call get_command_argument(4, value)
    read (value, *) storage
    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= 0) error stop 'read failed'
    call build_axisymmetric_mode_table(1, 6, mode_m, mode_n, powers)
    threads = 1
    !$omp parallel
    !$omp single
    threads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    start = omp_get_wtime()
    call build_compatible_two_component_problem(equilibrium, mode_m, mode_n, &
        powers, 1, degree, 64, 8, problem, info, density_kg_m3=1.0_dp, &
        sparse_storage=storage == 1, regular_quotient=.true., radial_cells=cells)
    assembled = omp_get_wtime()
    if (info /= 0) error stop 'assembly failed'
    call solve_compatible_marginality_problem(problem, .true., result, info, message)
    finished = omp_get_wtime()
    if (info /= 0) then
        print *, trim(message)
        error stop 'solve failed'
    end if
    print '(a)', 'threads,cells,degree,sparse,assembly_s,solve_s,lambda,certificate,negative'
    print '(4(i0,","),4(es24.16,","),i0)', threads, cells, degree, storage, &
        assembled - start, finished - assembled, result%lowest_eigenvalue, &
        result%certificate, result%negative_count
end program sparse_performance
