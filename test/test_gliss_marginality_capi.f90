program test_gliss_marginality_capi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_loc, &
        c_null_ptr, c_ptr, c_size_t, c_sizeof
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cylinder_fixture, only: create_cylinder_fixture
    implicit none

    integer(c_int), parameter :: status_ok = 0
    integer(c_int), parameter :: status_invalid_argument = 4
    character(len=*), parameter :: base_file = "marginality_capi_base.nc"
    character(len=*), parameter :: perturbed_file = &
        "marginality_capi_perturbed.nc"

    type, bind(c) :: marginality_result_c
        integer(c_size_t) :: struct_size
        integer(c_int) :: has_eigenpair
        integer(c_int) :: field_periods
        integer(c_size_t) :: mode_count
        integer(c_size_t) :: radial_surfaces
        integer(c_int) :: parity_class
        integer(c_int) :: radial_quadrature
        integer(c_int) :: angular_theta
        integer(c_int) :: angular_zeta
        integer(c_size_t) :: negative_count
        real(c_double) :: lowest_eigenvalue
        real(c_double) :: certificate
        real(c_double) :: eigenpair_residual
        real(c_double) :: force_balance_residual
    end type marginality_result_c

    type, bind(c) :: axisymmetric_result_c
        integer(c_size_t) :: struct_size
        integer(c_int) :: has_eigenpair
        integer(c_int) :: field_periods
        integer(c_int) :: toroidal_mode
        integer(c_int) :: poloidal_max
        integer(c_size_t) :: mode_count
        integer(c_size_t) :: radial_surfaces
        integer(c_int) :: parity_class
        integer(c_int) :: radial_quadrature
        integer(c_size_t) :: negative_count
        real(c_double) :: lowest_eigenvalue
        real(c_double) :: certificate
        real(c_double) :: eigenpair_residual
        real(c_double) :: force_balance_residual
    end type axisymmetric_result_c

    integer(c_int), target :: mode_m(7) = &
        [0, 1, 1, 2, 2, 3, 3]
    integer(c_int), target :: mode_n(7) = &
        [1, -1, 1, -1, 1, -1, 1]
    integer(c_int), target :: coupled_m(2) = [1, 1]
    integer(c_int), target :: coupled_n(2) = [0, 1]
    integer(c_int), target :: duplicate_m(2) = [1, 1]
    integer(c_int), target :: duplicate_n(2) = [1, 1]
    integer(c_int), target :: phase_direct_m(3) = [1, 2, 0]
    integer(c_int), target :: phase_direct_n(3) = [1, 1, 1]
    integer(c_int), target :: envelope_m(2) = [0, 1]
    integer(c_int), target :: envelope_n(2) = [0, 0]
    integer(c_int), target :: collision_m(3) = [0, 0, 0]
    integer(c_int), target :: collision_n(3) = [0, -1, 1]
    integer(c_int), target :: bad_envelope_m(1) = [1]
    integer(c_int), target :: bad_envelope_n(1) = [0]
    character(c_char), target :: error_buffer(256)
    type(c_ptr), target :: base, perturbed
    type(marginality_result_c), target :: general, inertia, shifted
    type(marginality_result_c), target :: phase_direct, phase_prefix
    type(marginality_result_c), target :: phase_collision
    type(axisymmetric_result_c), target :: convenience
    integer(c_int) :: status

    interface
        function equilibrium_create(path, path_length, handle, error, &
                error_capacity) bind(c, name="gliss_equilibrium_create") &
                result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: path, handle, error
            integer(c_size_t), value :: path_length, error_capacity
            integer(c_int) :: result
        end function equilibrium_create

        function equilibrium_destroy(handle, error, error_capacity) bind(c, &
                name="gliss_equilibrium_destroy") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, error
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function equilibrium_destroy

        function cas3d_marginality(equilibrium, mode_count, poloidal, &
                toroidal, parity_class, radial_quadrature, angular_theta, &
                angular_zeta, solve_eigenpair, result_pointer, error, &
                error_capacity) bind(c, name="gliss_cas3d_marginality") &
                result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: equilibrium, poloidal, toroidal
            type(c_ptr), value :: result_pointer, error
            integer(c_size_t), value :: mode_count, error_capacity
            integer(c_int), value :: parity_class, radial_quadrature
            integer(c_int), value :: angular_theta, angular_zeta
            integer(c_int), value :: solve_eigenpair
            integer(c_int) :: result
        end function cas3d_marginality

        function cas3d_phase_envelope(equilibrium, base_m, base_n, &
                envelope_count, poloidal, toroidal, parity_class, &
                radial_quadrature, angular_theta, angular_zeta, &
                solve_eigenpair, result_pointer, error, error_capacity) &
                bind(c, name="gliss_cas3d_phase_envelope") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: equilibrium, poloidal, toroidal
            type(c_ptr), value :: result_pointer, error
            integer(c_int), value :: base_m, base_n
            integer(c_size_t), value :: envelope_count, error_capacity
            integer(c_int), value :: parity_class, radial_quadrature
            integer(c_int), value :: angular_theta, angular_zeta
            integer(c_int), value :: solve_eigenpair
            integer(c_int) :: result
        end function cas3d_phase_envelope

        function axisymmetric_spectrum(equilibrium, toroidal_mode, &
                poloidal_max, radial_quadrature, solve_eigenpair, &
                result_pointer, error, error_capacity) bind(c, &
                name="gliss_axisymmetric_spectrum") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: equilibrium, result_pointer, error
            integer(c_int), value :: toroidal_mode, poloidal_max
            integer(c_int), value :: radial_quadrature, solve_eigenpair
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function axisymmetric_spectrum
    end interface

    call create_cylinder_fixture(base_file, chart_shift=0.0_dp, surfaces=32)
    call create_cylinder_fixture(perturbed_file, chart_shift=0.0_dp, &
        surfaces=32, toroidal_perturbation=0.02_dp)
    call load_equilibrium(base_file, base)
    call load_equilibrium(perturbed_file, perturbed)

    general%struct_size = c_sizeof(general)
    status = cas3d_marginality(base, size(mode_m, kind=c_size_t), &
        c_loc(mode_m), c_loc(mode_n), 1_c_int, 1_c_int, 64_c_int, 8_c_int, &
        1_c_int, c_loc(general), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "forced-general solve failed")
    convenience%struct_size = c_sizeof(convenience)
    status = axisymmetric_spectrum(base, 1_c_int, 3_c_int, 1_c_int, &
        1_c_int, c_loc(convenience), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "axisymmetric convenience solve failed")
    call require_same_result(general, convenience)

    phase_direct%struct_size = c_sizeof(phase_direct)
    status = cas3d_marginality(base, &
        size(phase_direct_m, kind=c_size_t), c_loc(phase_direct_m), &
        c_loc(phase_direct_n), 1_c_int, 1_c_int, 64_c_int, 16_c_int, &
        1_c_int, c_loc(phase_direct), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "direct phase-prefix solve failed")
    phase_prefix%struct_size = c_sizeof(phase_prefix)
    status = cas3d_phase_envelope(base, 1_c_int, 1_c_int, &
        size(envelope_m, kind=c_size_t), c_loc(envelope_m), &
        c_loc(envelope_n), 1_c_int, 1_c_int, 64_c_int, 16_c_int, 1_c_int, &
        c_loc(phase_prefix), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "phase-prefix solve failed")
    call require_same_marginality(phase_direct, phase_prefix)

    phase_collision%struct_size = c_sizeof(phase_collision)
    status = cas3d_phase_envelope(base, 1_c_int, 1_c_int, &
        size(collision_m, kind=c_size_t), c_loc(collision_m), &
        c_loc(collision_n), 1_c_int, 1_c_int, 64_c_int, 16_c_int, 0_c_int, &
        c_loc(phase_collision), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "colliding phase envelope failed")
    call require(phase_collision%mode_count == 5_c_size_t, &
        "colliding phase envelope lost labeled sidebands")

    inertia%struct_size = c_sizeof(inertia)
    status = cas3d_marginality(base, size(mode_m, kind=c_size_t), &
        c_loc(mode_m), c_loc(mode_n), 2_c_int, 1_c_int, 64_c_int, 8_c_int, &
        0_c_int, c_loc(inertia), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "parity-two inertia failed")
    call require(inertia%has_eigenpair == 0, &
        "inertia call reported an eigenpair")
    call require(inertia%parity_class == 2, "parity metadata is wrong")
    call require(ieee_is_nan(inertia%lowest_eigenvalue), &
        "inertia eigenvalue is not NaN")

    call solve_coupled(base, general)
    call solve_coupled(perturbed, shifted)
    call require(general%lowest_eigenvalue /= shifted%lowest_eigenvalue, &
        "nonaxisymmetric coupling did not change the spectrum")

    general%mode_count = 999_c_size_t
    status = cas3d_marginality(base, size(duplicate_m, kind=c_size_t), &
        c_loc(duplicate_m), c_loc(duplicate_n), 1_c_int, 1_c_int, 64_c_int, &
        8_c_int, 1_c_int, c_loc(general), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_invalid_argument, &
        "duplicate modes were accepted")
    call require(general%mode_count == 999_c_size_t, &
        "failed mode validation modified the result")
    status = cas3d_marginality(base, size(mode_m, kind=c_size_t), &
        c_loc(mode_m), c_loc(mode_n), 1_c_int, 1_c_int, 64_c_int, 8_c_int, &
        2_c_int, c_loc(general), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_invalid_argument, &
        "invalid solve_eigenpair was accepted")
    call require(general%mode_count == 999_c_size_t, &
        "failed solve validation modified the result")

    phase_collision%mode_count = 777_c_size_t
    status = cas3d_phase_envelope(base, 1_c_int, 1_c_int, &
        size(bad_envelope_m, kind=c_size_t), c_loc(bad_envelope_m), &
        c_loc(bad_envelope_n), 1_c_int, 1_c_int, 64_c_int, 16_c_int, &
        0_c_int, c_loc(phase_collision), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_invalid_argument, &
        "phase envelope without origin was accepted")
    call require(phase_collision%mode_count == 777_c_size_t, &
        "failed phase-envelope validation modified the result")

    status = equilibrium_destroy(c_loc(base), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "base equilibrium destroy failed")
    status = equilibrium_destroy(c_loc(perturbed), c_loc(error_buffer), &
        size(error_buffer, kind=c_size_t))
    call require(status == status_ok, "perturbed equilibrium destroy failed")
    call delete_file(base_file)
    call delete_file(perturbed_file)
    write (*, "(a)") "PASS"

contains

    subroutine load_equilibrium(filename, equilibrium)
        character(len=*), intent(in) :: filename
        type(c_ptr), intent(out), target :: equilibrium
        character(c_char), target :: path(len(filename))
        integer :: index

        do index = 1, len(filename)
            path(index) = filename(index:index)
        end do
        equilibrium = c_null_ptr
        status = equilibrium_create(c_loc(path), &
            int(size(path), c_size_t), c_loc(equilibrium), &
            c_loc(error_buffer), size(error_buffer, kind=c_size_t))
        call require(status == status_ok, "equilibrium creation failed")
    end subroutine load_equilibrium

    subroutine solve_coupled(equilibrium, result)
        type(c_ptr), intent(in) :: equilibrium
        type(marginality_result_c), intent(out), target :: result

        result%struct_size = c_sizeof(result)
        status = cas3d_marginality(equilibrium, &
            size(coupled_m, kind=c_size_t), c_loc(coupled_m), &
            c_loc(coupled_n), 1_c_int, 1_c_int, 64_c_int, 16_c_int, &
            1_c_int, c_loc(result), c_loc(error_buffer), &
            size(error_buffer, kind=c_size_t))
        call require(status == status_ok, "coupled general solve failed")
        call require(ieee_is_finite(result%lowest_eigenvalue), &
            "coupled result is not finite")
    end subroutine solve_coupled

    subroutine require_same_result(general_result, axisymmetric_result)
        type(marginality_result_c), intent(in) :: general_result
        type(axisymmetric_result_c), intent(in) :: axisymmetric_result

        call require(general_result%has_eigenpair == &
            axisymmetric_result%has_eigenpair, "eigenpair flags differ")
        call require(general_result%field_periods == &
            axisymmetric_result%field_periods, "field periods differ")
        call require(general_result%mode_count == &
            axisymmetric_result%mode_count, "mode counts differ")
        call require(general_result%radial_surfaces == &
            axisymmetric_result%radial_surfaces, "radial counts differ")
        call require(general_result%parity_class == &
            axisymmetric_result%parity_class, "parity classes differ")
        call require(general_result%radial_quadrature == &
            axisymmetric_result%radial_quadrature, "quadrature differs")
        call require(general_result%negative_count == &
            axisymmetric_result%negative_count, "inertia counts differ")
        call require(general_result%lowest_eigenvalue == &
            axisymmetric_result%lowest_eigenvalue, "eigenvalues differ")
        call require(general_result%certificate == &
            axisymmetric_result%certificate, "certificates differ")
        call require(general_result%eigenpair_residual == &
            axisymmetric_result%eigenpair_residual, "residuals differ")
        call require(general_result%force_balance_residual == &
            axisymmetric_result%force_balance_residual, &
            "force residuals differ")
    end subroutine require_same_result

    subroutine require_same_marginality(first, second)
        type(marginality_result_c), intent(in) :: first, second

        call require(first%has_eigenpair == second%has_eigenpair, &
            "phase-prefix eigenpair flags differ")
        call require(first%field_periods == second%field_periods, &
            "phase-prefix field periods differ")
        call require(first%mode_count == second%mode_count, &
            "phase-prefix mode counts differ")
        call require(first%radial_surfaces == second%radial_surfaces, &
            "phase-prefix radial counts differ")
        call require(first%parity_class == second%parity_class, &
            "phase-prefix parity classes differ")
        call require(first%radial_quadrature == second%radial_quadrature, &
            "phase-prefix quadratures differ")
        call require(first%negative_count == second%negative_count, &
            "phase-prefix inertia counts differ")
        call require(first%lowest_eigenvalue == second%lowest_eigenvalue, &
            "phase-prefix eigenvalues differ")
        call require(first%certificate == second%certificate, &
            "phase-prefix certificates differ")
        call require(first%eigenpair_residual == second%eigenpair_residual, &
            "phase-prefix residuals differ")
        call require(first%force_balance_residual == &
            second%force_balance_residual, &
            "phase-prefix force-balance residuals differ")
    end subroutine require_same_marginality

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: io_status, unit

        open (newunit=unit, file=filename, status="old", iostat=io_status)
        if (io_status == 0) close (unit, status="delete")
    end subroutine delete_file

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_gliss_marginality_capi
