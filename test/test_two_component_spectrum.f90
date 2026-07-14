program test_two_component_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use axisymmetric_spectrum, only: axisymmetric_spectrum_ok, &
        axisymmetric_spectrum_result_t, build_axisymmetric_mode_table, &
        compute_axisymmetric_spectrum
    use cylinder_fixture, only: create_cylinder_fixture
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: equilibrium_is_axisymmetric, &
        gvec_cas3d_equilibrium_t
    use two_component_spectrum, only: compute_two_component_spectrum, &
        two_component_spectrum_invalid, two_component_spectrum_ok, &
        two_component_spectrum_result_t
    implicit none

    character(len=*), parameter :: fixture = "two_component_spectrum.nc"
    type(axisymmetric_spectrum_result_t) :: axisymmetric
    type(two_component_spectrum_result_t) :: direct
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    integer, allocatable :: mode_m(:), mode_n(:)
    real(dp), allocatable :: stored_power(:)
    character(len=256) :: message
    integer :: info

    call create_cylinder_fixture(fixture, chart_shift=0.0_dp, surfaces=32)
    call read_gvec_cas3d_file(fixture, equilibrium, info)
    call delete_file(fixture)
    call require(info == reader_ok, "cylinder equilibrium was rejected")
    call require(equilibrium_is_axisymmetric(equilibrium), &
        "rotating Cartesian position harmonics broke axisymmetry detection")
    call build_axisymmetric_mode_table(1, 3, mode_m, mode_n, stored_power)
    call compute_axisymmetric_spectrum(equilibrium, 1, 3, 1, .true., &
        axisymmetric, info, message)
    call require(info == axisymmetric_spectrum_ok, &
        "axisymmetric convenience solve failed: " // trim(message))
    call compute_two_component_spectrum(equilibrium, mode_m, mode_n, &
        stored_power, 1, 1, 64, 8, .true., direct, info, message)
    call require(info == two_component_spectrum_ok, &
        "forced-general solve failed: " // trim(message))
    call require_same_result(axisymmetric, direct)
    call check_count_only(equilibrium, mode_m, mode_n, stored_power)
    call check_nonaxisymmetric(equilibrium)
    call compute_two_component_spectrum(equilibrium, mode_m, mode_n, &
        stored_power, 1, 2, 64, 8, .false., direct, info, message)
    call require(info == two_component_spectrum_invalid, &
        "unsupported quadrature was accepted")
    write (*, "(a)") "PASS"

contains

    subroutine check_count_only(state, poloidal, toroidal, powers)
        type(gvec_cas3d_equilibrium_t), intent(in) :: state
        integer, intent(in) :: poloidal(:), toroidal(:)
        real(dp), intent(in) :: powers(:)
        type(two_component_spectrum_result_t) :: result
        character(len=256) :: local_message
        integer :: status

        call compute_two_component_spectrum(state, poloidal, toroidal, &
            powers, 1, 1, 64, 8, .false., result, status, local_message)
        call require(status == two_component_spectrum_ok, &
            "forced-general inertia solve failed")
        call require(.not. result%has_eigenpair, &
            "count-only solve returned an eigenpair")
        call require(ieee_is_nan(result%lowest_eigenvalue), &
            "count-only eigenvalue is not NaN")
    end subroutine check_count_only

    subroutine check_nonaxisymmetric(axisymmetric_state)
        type(gvec_cas3d_equilibrium_t), intent(in) :: axisymmetric_state
        type(gvec_cas3d_equilibrium_t) :: state
        type(two_component_spectrum_result_t) :: reference, result
        real(dp) :: scale
        character(len=256) :: local_message
        integer :: harmonic, status

        call compute_two_component_spectrum(axisymmetric_state, [1, 1], &
            [0, 1], [0.0_dp, 0.0_dp], 1, 1, 64, 16, .true., reference, &
            status, local_message)
        call require(status == two_component_spectrum_ok, &
            "nonaxisymmetric reference solve failed")
        state = axisymmetric_state
        harmonic = findloc(state%toroidal_modes, 1, dim=1)
        call require(harmonic > 0, "fixture lacks a toroidal harmonic")
        scale = max(1.0_dp, maxval(abs(state%mod_b%cosine)))
        state%mod_b%cosine(:, 1, harmonic) = 0.02_dp * scale
        call require(.not. equilibrium_is_axisymmetric(state), &
            "toroidal perturbation remained axisymmetric")
        call compute_two_component_spectrum(state, [1, 1], [0, 1], &
            [0.0_dp, 0.0_dp], 1, 1, 64, 16, .true., result, status, &
            local_message)
        call require(status == two_component_spectrum_ok, &
            "forced-3D perturbed solve failed")
        call require(ieee_is_finite(result%lowest_eigenvalue), &
            "forced-3D result is not finite")
        call require(result%lowest_eigenvalue /= reference%lowest_eigenvalue, &
            "toroidal coupling did not change the spectrum")
    end subroutine check_nonaxisymmetric

    subroutine require_same_result(wrapper, general)
        type(axisymmetric_spectrum_result_t), intent(in) :: wrapper
        type(two_component_spectrum_result_t), intent(in) :: general

        call require(wrapper%has_eigenpair .eqv. general%has_eigenpair, &
            "eigenpair flags differ")
        call require(wrapper%field_periods == general%field_periods, &
            "field periods differ")
        call require(wrapper%mode_count == general%mode_count, &
            "mode counts differ")
        call require(wrapper%radial_surfaces == general%radial_surfaces, &
            "radial counts differ")
        call require(wrapper%parity_class == general%parity_class, &
            "parity classes differ")
        call require(wrapper%radial_quadrature == general%radial_quadrature, &
            "quadrature metadata differ")
        call require(wrapper%negative_count == general%negative_count, &
            "inertia counts differ")
        call require(wrapper%lowest_eigenvalue == general%lowest_eigenvalue, &
            "eigenvalues differ")
        call require(wrapper%certificate == general%certificate, &
            "certificates differ")
        call require(wrapper%eigenpair_residual == &
            general%eigenpair_residual, "residuals differ")
        call require(wrapper%force_balance_residual == &
            general%force_balance_residual, "force residuals differ")
    end subroutine require_same_result

    subroutine delete_file(path)
        character(len=*), intent(in) :: path
        integer :: status, unit

        open (newunit=unit, file=path, status="old", iostat=status)
        if (status == 0) close (unit, status="delete")
    end subroutine delete_file

    subroutine require(condition, text)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: text

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // text
        error stop 1
    end subroutine require

end program test_two_component_spectrum
