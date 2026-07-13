program test_gliss_spectrum_capi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
        c_f_pointer, c_int, c_loc, c_null_ptr, c_ptr, c_size_t, c_sizeof
    use, intrinsic :: iso_fortran_env, only: error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    implicit none

    integer(c_int), parameter :: status_ok = 0
    integer(c_int), parameter :: status_capacity = 3
    integer(c_int), parameter :: status_invalid_argument = 4
    character(len=*), parameter :: fixture = "spectrum_capi_cylinder.nc"
    real(c_double), parameter :: expected_lowest = &
        -1.2278816610508129e2_c_double
    real(c_double), parameter :: legacy_relative_tolerance = 1.0e-9_c_double

    type, bind(c) :: spectrum_summary_c
        integer(c_size_t) :: struct_size
        integer(c_int) :: has_chart_metric
        integer(c_int) :: has_eigenvector
        integer(c_int) :: field_periods
        integer(c_int) :: parity_class
        integer(c_int) :: radial_quadrature
        integer(c_int) :: angular_theta
        integer(c_int) :: angular_zeta
        integer(c_size_t) :: mode_count
        integer(c_size_t) :: unknowns
        integer(c_size_t) :: normal_unknowns
        integer(c_size_t) :: eta_unknowns
        integer(c_size_t) :: mu_unknowns
        integer(c_size_t) :: negative_count
        integer(c_size_t) :: floor_count
        real(c_double) :: adiabatic_index
        real(c_double) :: density_kg_m3
        real(c_double) :: zero_floor
        real(c_double) :: lowest_eigenvalue
        real(c_double) :: certificate
        real(c_double) :: eigenpair_residual
        real(c_double) :: eigenpair_resolution
        real(c_double) :: inertia_interval
    end type spectrum_summary_c

    type, bind(c) :: energy_terms_c
        integer(c_size_t) :: struct_size
        real(c_double) :: field_line_bending
        real(c_double) :: magnetic_shear
        real(c_double) :: magnetic_compression
        real(c_double) :: pressure_drive
        real(c_double) :: plasma_compressibility
        real(c_double) :: potential_energy
        real(c_double) :: kinetic_energy
        real(c_double) :: rayleigh_quotient
        real(c_double) :: closure_error
        real(c_double) :: closure_tolerance
    end type energy_terms_c

    character(c_char), target :: path(len(fixture)), error_buffer(256)
    integer(c_int), target :: mode_m(2) = [1_c_int, 2_c_int]
    integer(c_int), target :: mode_n(2) = [1_c_int, 1_c_int]
    type(c_ptr), target :: equilibrium, problem, rejected
    integer(c_size_t), target :: unknowns, written
    integer(c_size_t), target :: eigenvalues_written, eigenvectors_written
    type(spectrum_summary_c), target :: summary
    type(energy_terms_c), target :: energy
    real(c_double), allocatable, target :: eigenvector(:), sentinel(:)
    real(c_double), allocatable, target :: eigenvalues(:), eigenvectors(:)
    real(c_double), allocatable, target :: residuals(:), resolutions(:)
    real(c_double), allocatable, target :: rayleigh_quotients(:)
    real(c_double), allocatable, target :: value_sentinel(:), vector_sentinel(:)
    integer(c_int) :: status

    interface
        function equilibrium_create(path_pointer, path_length, handle, &
                error_pointer, error_capacity) &
                bind(c, name="gliss_equilibrium_create") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: path_pointer, handle, error_pointer
            integer(c_size_t), value :: path_length, error_capacity
            integer(c_int) :: result
        end function equilibrium_create

        function equilibrium_destroy(handle, error_pointer, error_capacity) &
                bind(c, name="gliss_equilibrium_destroy") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, error_pointer
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function equilibrium_destroy

        function problem_create(equilibrium_handle, gamma, density, floor, &
                mode_count, poloidal, toroidal, radial_quadrature, handle, &
                error_pointer, error_capacity) &
                bind(c, name="gliss_stability_problem_create") result(result)
            import c_double, c_int, c_ptr, c_size_t
            type(c_ptr), value :: equilibrium_handle, poloidal, toroidal
            type(c_ptr), value :: handle, error_pointer
            real(c_double), value :: gamma, density, floor
            integer(c_size_t), value :: mode_count, error_capacity
            integer(c_int), value :: radial_quadrature
            integer(c_int) :: result
        end function problem_create

        function problem_destroy(handle, error_pointer, error_capacity) &
                bind(c, name="gliss_stability_problem_destroy") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, error_pointer
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function problem_destroy

        function problem_unknown_count(handle, parity_class, unknown_count, &
                error_pointer, error_capacity) &
                bind(c, name="gliss_stability_problem_unknown_count") &
                result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, unknown_count, error_pointer
            integer(c_int), value :: parity_class
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function problem_unknown_count

        function problem_solve_class(handle, parity_class, capacity, vector, &
                written_pointer, summary_pointer, error_pointer, &
                error_capacity) bind(c, name="gliss_stability_problem_solve_class") &
                result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, vector, written_pointer
            type(c_ptr), value :: summary_pointer, error_pointer
            integer(c_int), value :: parity_class
            integer(c_size_t), value :: capacity, error_capacity
            integer(c_int) :: result
        end function problem_solve_class

        function problem_full_spectrum(handle, parity_class, value_capacity, &
                values, residuals, resolutions, rayleigh_quotients, &
                vector_capacity, vectors, values_written, &
                vectors_written, error_pointer, error_capacity) bind(c, &
                name="gliss_stability_problem_full_spectrum") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, values, residuals, resolutions
            type(c_ptr), value :: rayleigh_quotients
            type(c_ptr), value :: vectors, values_written
            type(c_ptr), value :: vectors_written, error_pointer
            integer(c_int), value :: parity_class
            integer(c_size_t), value :: value_capacity, vector_capacity
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function problem_full_spectrum

        function problem_energy(handle, parity_class, vector_count, vector, &
                terms, error_pointer, error_capacity) &
                bind(c, name="gliss_stability_problem_energy") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, vector, terms, error_pointer
            integer(c_int), value :: parity_class
            integer(c_size_t), value :: vector_count, error_capacity
            integer(c_int) :: result
        end function problem_energy
    end interface

    call create_cylinder_fixture(fixture)
    call copy_chars(fixture, path)
    equilibrium = c_null_ptr
    problem = c_null_ptr
    status = equilibrium_create(c_loc(path), int(size(path), c_size_t), &
        c_loc(equilibrium), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "equilibrium creation failed")
    rejected = c_null_ptr
    status = problem_create(c_null_ptr, 5.0_c_double / 3.0_c_double, &
        2.0_c_double, 1.0_c_double, 2_c_size_t, c_loc(mode_m), &
        c_loc(mode_n), 1_c_int, c_loc(rejected), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null equilibrium handle was accepted")
    call require(.not. c_associated(rejected), &
        "rejected create returned a problem handle")
    status = problem_create(equilibrium, 5.0_c_double / 3.0_c_double, &
        2.0_c_double, 1.0_c_double, 0_c_size_t, c_null_ptr, c_null_ptr, &
        1_c_int, c_loc(rejected), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "empty mode arrays were accepted")
    status = problem_create(equilibrium, 5.0_c_double / 3.0_c_double, &
        2.0_c_double, 1.0_c_double, 2_c_size_t, c_loc(mode_m), &
        c_loc(mode_n), 1_c_int, c_null_ptr, c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null problem output pointer was accepted")
    status = problem_create(equilibrium, 5.0_c_double / 3.0_c_double, &
        2.0_c_double, 1.0_c_double, 2_c_size_t, c_loc(mode_m), &
        c_loc(mode_n), 1_c_int, c_loc(problem), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "problem creation failed")
    call require(c_associated(problem), "problem handle is null")
    status = equilibrium_destroy(c_loc(equilibrium), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "early equilibrium destroy failed")

    status = problem_unknown_count(problem, 1_c_int, c_null_ptr, &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null unknown-count output was accepted")
    status = problem_unknown_count(problem, 1_c_int, c_loc(unknowns), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok, "unknown-count query failed")
    call require(unknowns == 196_c_size_t, "unknown count is wrong")
    allocate (eigenvector(unknowns), sentinel(unknowns))
    sentinel = 7.0_c_double
    eigenvector = sentinel
    summary%struct_size = c_sizeof(summary)
    status = problem_solve_class(problem, 1_c_int, unknowns - 1_c_size_t, &
        c_loc(eigenvector), c_loc(written), c_loc(summary), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_capacity, "small capacity was not rejected")
    call require(written == unknowns, "required vector capacity is wrong")
    call require(all(eigenvector == sentinel), &
        "capacity failure modified the eigenvector")

    summary%struct_size = c_sizeof(summary)
    status = problem_solve_class(problem, 1_c_int, unknowns, &
        c_loc(eigenvector), c_loc(written), c_loc(summary), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok, "class solve failed")
    call require(written == unknowns, "written vector size is wrong")
    call require(summary%parity_class == 1, "summary parity is wrong")
    call require(summary%unknowns == unknowns, "summary size is wrong")
    call require(summary%normal_unknowns == 64_c_size_t, &
        "summary normal size is wrong")
    call require(summary%eta_unknowns == 66_c_size_t, &
        "summary eta size is wrong")
    call require(summary%mu_unknowns == 66_c_size_t, &
        "summary mu size is wrong")
    call require(summary%negative_count == 1_c_size_t, &
        "negative inertia count is wrong")
    call require(summary%floor_count == 5_c_size_t, &
        "floor count is wrong")
    call require(abs(summary%lowest_eigenvalue - expected_lowest) &
        <= legacy_relative_tolerance * abs(expected_lowest), &
        "C API changed the legacy eigenvalue")
    call require(all(ieee_is_finite(eigenvector)), &
        "C API eigenvector contains nonfinite values")
    call require(summary%certificate == summary%inertia_interval &
        + summary%eigenpair_residual + summary%eigenpair_resolution, &
        "C API certificate components do not close")

    energy%struct_size = c_sizeof(energy)
    status = problem_energy(problem, 1_c_int, unknowns, c_loc(eigenvector), &
        c_loc(energy), c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok, "C API energy decomposition failed")
    call require(abs(energy%kinetic_energy - 1.0_c_double) < 1.0e-11_c_double, &
        "C API energy vector is not mass normalized")
    call require(abs(energy%rayleigh_quotient - summary%lowest_eigenvalue) &
        <= 1.0e-10_c_double * abs(summary%lowest_eigenvalue), &
        "C API energy quotient differs from the eigenvalue")
    call require(energy%closure_error <= energy%closure_tolerance, &
        "C API energy terms do not close")
    status = problem_energy(problem, 1_c_int, unknowns - 1_c_size_t, &
        c_loc(eigenvector), c_loc(energy), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "C API energy accepted the wrong vector size")
    status = problem_energy(problem, 1_c_int, unknowns, c_null_ptr, &
        c_loc(energy), c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "C API energy accepted a null vector")
    sentinel = eigenvector
    eigenvector(1) = ieee_value(0.0_c_double, ieee_quiet_nan)
    status = problem_energy(problem, 1_c_int, unknowns, c_loc(eigenvector), &
        c_loc(energy), c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "C API energy accepted a nonfinite vector")
    eigenvector = 0.0_c_double
    status = problem_energy(problem, 1_c_int, unknowns, c_loc(eigenvector), &
        c_loc(energy), c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "C API energy accepted a zero-norm vector")
    eigenvector = sentinel
    status = problem_energy(problem, 0_c_int, unknowns, c_loc(eigenvector), &
        c_loc(energy), c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "C API energy accepted an invalid parity class")
    energy%struct_size = 0_c_size_t
    status = problem_energy(problem, 1_c_int, unknowns, c_loc(eigenvector), &
        c_loc(energy), c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "C API energy accepted an incompatible struct size")
    energy%struct_size = c_sizeof(energy)
    status = problem_energy(problem, 1_c_int, unknowns, c_loc(eigenvector), &
        c_null_ptr, c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "C API energy accepted a null output")

    allocate (eigenvalues(unknowns), residuals(unknowns), &
        resolutions(unknowns), rayleigh_quotients(unknowns), &
        value_sentinel(unknowns))
    allocate (eigenvectors(unknowns * unknowns), &
        vector_sentinel(unknowns * unknowns))
    value_sentinel = 8.0_c_double
    vector_sentinel = 9.0_c_double
    eigenvalues = value_sentinel
    residuals = value_sentinel
    resolutions = value_sentinel
    rayleigh_quotients = value_sentinel
    eigenvectors = vector_sentinel
    status = problem_full_spectrum(problem, 1_c_int, unknowns - 1_c_size_t, &
        c_loc(eigenvalues), c_loc(residuals), c_loc(resolutions), &
        c_loc(rayleigh_quotients), unknowns * unknowns, c_loc(eigenvectors), &
        c_loc(eigenvalues_written), &
        c_loc(eigenvectors_written), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_capacity, &
        "small full-spectrum eigenvalue capacity was not rejected")
    call require(eigenvalues_written == unknowns, &
        "required full-spectrum eigenvalue count is wrong")
    call require(eigenvectors_written == unknowns * unknowns, &
        "required full-spectrum eigenvector count is wrong")
    call require(all(eigenvalues == value_sentinel), &
        "full-spectrum capacity failure modified eigenvalues")
    call require(all(residuals == value_sentinel), &
        "full-spectrum capacity failure modified residuals")
    call require(all(resolutions == value_sentinel), &
        "full-spectrum capacity failure modified resolutions")
    call require(all(rayleigh_quotients == value_sentinel), &
        "full-spectrum capacity failure modified Rayleigh quotients")
    call require(all(eigenvectors == vector_sentinel), &
        "full-spectrum capacity failure modified eigenvectors")
    status = problem_full_spectrum(problem, 1_c_int, unknowns, &
        c_loc(eigenvalues), c_loc(residuals), c_loc(resolutions), &
        c_loc(rayleigh_quotients), unknowns * unknowns - 1_c_size_t, &
        c_loc(eigenvectors), &
        c_loc(eigenvalues_written), &
        c_loc(eigenvectors_written), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_capacity, &
        "small full-spectrum eigenvector capacity was not rejected")
    status = problem_full_spectrum(problem, 1_c_int, unknowns, &
        c_loc(eigenvalues), c_loc(residuals), c_loc(resolutions), &
        c_loc(rayleigh_quotients), unknowns * unknowns, c_loc(eigenvectors), &
        c_loc(eigenvalues_written), &
        c_loc(eigenvectors_written), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "full-spectrum C API solve failed")
    call require(eigenvalues_written == unknowns, &
        "full-spectrum C API wrote the wrong eigenvalue count")
    call require(eigenvectors_written == unknowns * unknowns, &
        "full-spectrum C API wrote the wrong eigenvector count")
    call require(all(ieee_is_finite(eigenvalues)), &
        "full-spectrum C API returned nonfinite eigenvalues")
    call require(all(ieee_is_finite(residuals)) .and. all(residuals >= 0.0), &
        "full-spectrum C API returned invalid residuals")
    call require(all(ieee_is_finite(resolutions)) &
        .and. all(resolutions >= 0.0), &
        "full-spectrum C API returned invalid resolutions")
    call require(all(ieee_is_finite(rayleigh_quotients)), &
        "full-spectrum C API returned nonfinite Rayleigh quotients")
    call require(all(ieee_is_finite(eigenvectors)), &
        "full-spectrum C API returned nonfinite eigenvectors")
    call require(all(eigenvalues(2:) >= &
        eigenvalues(:size(eigenvalues) - 1)), &
        "full-spectrum C API eigenvalues are not sorted")
    call require(abs(eigenvalues(merge(1_c_size_t, summary%floor_count &
        + 1_c_size_t, summary%negative_count > 0_c_size_t)) &
        - summary%lowest_eigenvalue) &
        <= summary%certificate + 1.0e-12_c_double &
        * abs(summary%lowest_eigenvalue), &
        "full-spectrum C API disagrees with certified active eigenvalue")
    call require(abs(dot_product(eigenvectors(:unknowns), eigenvector)) &
        / sqrt(dot_product(eigenvectors(:unknowns), eigenvectors(:unknowns)) &
        * dot_product(eigenvector, eigenvector)) &
        > 1.0_c_double - 1.0e-10_c_double, &
        "full-spectrum C API component order is wrong")

    status = problem_full_spectrum(problem, 1_c_int, unknowns, c_null_ptr, &
        c_loc(residuals), c_loc(resolutions), c_loc(rayleigh_quotients), &
        unknowns * unknowns, c_loc(eigenvectors), &
        c_loc(eigenvalues_written), c_loc(eigenvectors_written), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null full-spectrum eigenvalue output was accepted")
    status = problem_full_spectrum(problem, 1_c_int, unknowns, &
        c_loc(eigenvalues), c_null_ptr, c_loc(resolutions), &
        c_loc(rayleigh_quotients), unknowns * unknowns, c_loc(eigenvectors), &
        c_loc(eigenvalues_written), c_loc(eigenvectors_written), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null full-spectrum residual output was accepted")
    status = problem_full_spectrum(problem, 1_c_int, unknowns, &
        c_loc(eigenvalues), c_loc(residuals), c_loc(resolutions), &
        c_loc(rayleigh_quotients), unknowns * unknowns, c_loc(eigenvectors), &
        c_null_ptr, &
        c_loc(eigenvectors_written), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null full-spectrum written output was accepted")
    status = problem_full_spectrum(problem, 0_c_int, 0_c_size_t, c_null_ptr, &
        c_null_ptr, c_null_ptr, c_null_ptr, 0_c_size_t, c_null_ptr, &
        c_loc(eigenvalues_written), c_loc(eigenvectors_written), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "full-spectrum C API accepted an invalid parity class")

    status = problem_solve_class(problem, 1_c_int, unknowns, c_null_ptr, &
        c_loc(written), c_loc(summary), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null eigenvector output was accepted")
    status = problem_solve_class(problem, 1_c_int, unknowns, &
        c_loc(eigenvector), c_null_ptr, c_loc(summary), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null written output was accepted")
    status = problem_solve_class(problem, 1_c_int, unknowns, &
        c_loc(eigenvector), c_loc(written), c_null_ptr, c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null summary output was accepted")

    summary%struct_size = 0_c_size_t
    status = problem_solve_class(problem, 1_c_int, unknowns, &
        c_loc(eigenvector), c_loc(written), c_loc(summary), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "wrong summary size was accepted")
    status = problem_unknown_count(problem, 0_c_int, c_loc(unknowns), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "invalid parity was accepted")

    status = problem_destroy(c_loc(problem), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "problem destroy failed")
    status = problem_destroy(c_loc(problem), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "repeated problem destroy failed")
    status = problem_destroy(c_null_ptr, c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null problem handle pointer was accepted")
    status = equilibrium_destroy(c_loc(equilibrium), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "equilibrium destroy failed")
    call delete_fixture()
    write (*, "(a)") "PASS"

contains

    subroutine copy_chars(source, destination)
        character(len=*), intent(in) :: source
        character(c_char), intent(out) :: destination(:)
        integer :: i

        do i = 1, size(destination)
            destination(i) = source(i:i)
        end do
    end subroutine copy_chars

    subroutine delete_fixture()
        integer :: unit, info

        open (newunit=unit, file=fixture, status="old", iostat=info)
        call require(info == 0, "failed to open fixture for deletion")
        close (unit, status="delete", iostat=info)
        call require(info == 0, "failed to delete fixture")
    end subroutine delete_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_gliss_spectrum_capi
