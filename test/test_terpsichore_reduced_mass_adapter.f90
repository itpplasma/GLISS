program test_terpsichore_reduced_mass_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: dynamic_family_layout_t
    use fourier_phase_kind, only: phase_sine
    use terpsichore_matrix_fixture, only: &
        read_terpsichore_fixed_boundary_fixture, &
        read_terpsichore_fixed_boundary_potential_fixture, &
        read_terpsichore_potential_fixture, &
        terpsichore_matrix_fixture_invalid, terpsichore_matrix_fixture_ok, &
        terpsichore_matrix_fixture_t, &
        terpsichore_potential_metadata_is_valid
    use terpsichore_reduced_mass_adapter, only: &
        assemble_terpsichore_fixture_reduced_mass, &
        assemble_terpsichore_fixture_reduced_mass_free_boundary, &
        terpsichore_reduced_adapter_ok
    use terpsichore_reduced_mass_family_assembly, only: &
        assemble_terpsichore_reduced_fixed_boundary_mass
    implicit none

    call test_fixture_reader()
    call test_potential_fixture_reader()
    call test_vacuum_potential_fixture_reader()
    call test_fixture_adapter()
    call test_free_boundary_adapter()
    call test_fixture_rejections()
    write (*, "(a)") "PASS"

contains

    subroutine test_fixture_reader()
        type(terpsichore_matrix_fixture_t) :: fixture
        integer :: info, unit

        open (newunit=unit, status="scratch", form="unformatted", &
            access="sequential")
        call write_fixture(unit)
        rewind (unit)
        call read_terpsichore_fixed_boundary_fixture(unit, 0, fixture, info)
        close (unit)
        call require(info == terpsichore_matrix_fixture_ok, &
            "valid TERPSICHORE matrix fixture was rejected")
        call require(fixture%intervals == 3 .and. fixture%poloidal_points == 2 &
            .and. fixture%toroidal_points == 2 .and. fixture%modes == 2, &
            "TERPSICHORE fixture dimensions are wrong")
        call require(all(fixture%mode_m == [1, 2]) .and. &
            all(fixture%mode_n == [0, 1]), &
            "TERPSICHORE fixture modes are wrong")
        call require(maxval(abs(fixture%s - [0.0_dp, 0.2_dp, 0.6_dp, &
            1.0_dp])) < 1.0e-7_dp, "TERPSICHORE radial grid is wrong")
        call require(abs(fixture%flux_t_slope(2) - 1.1_dp) < 1.0e-7_dp, &
            "TERPSICHORE toroidal flux slope is wrong")
        call require(abs(fixture%signed_bjac(4, 3) + 1.6_dp) < 1.0e-7_dp, &
            "TERPSICHORE BJAC ordering is wrong")
    end subroutine test_fixture_reader

    subroutine write_fixture(unit, modelk, vacuum_intervals)
        integer, intent(in) :: unit
        integer, intent(in), optional :: modelk, vacuum_intervals
        real(dp), allocatable :: bjac(:, :), full_surface(:, :), s(:)
        real(dp) :: pth(3), fpp(4), ftp(4)
        real(dp) :: current_i(3), current_j(3), pressure_slope(3)
        real(dp) :: parity, ql(2), surface_field(4, 0:3)
        real(dp) :: cell(4, 3), radial_surface(4), radial_cell(3)
        real(dp) :: mode_surface(2, 0:3)
        integer :: mode_m(2), mode_n(2), point, surface_index
        integer :: legacy_modelk, vacuum

        vacuum = 0
        if (present(vacuum_intervals)) vacuum = vacuum_intervals
        allocate (s(0:3 + vacuum), bjac(4, 0:3 + vacuum), &
            full_surface(4, 0:3 + vacuum))
        s(:3) = [0.0_dp, 0.2_dp, 0.6_dp, 1.0_dp]
        if (vacuum > 0) s(4:) = [(1.0_dp + 0.3_dp * &
            real(surface_index, dp), surface_index = 1, vacuum)]
        pth = 0.0_dp
        fpp = [0.2_dp, 0.25_dp, 0.3_dp, 0.35_dp]
        ftp = [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp]
        current_i = 0.1_dp
        current_j = 0.6_dp
        pressure_slope = -0.01_dp
        legacy_modelk = 0
        if (present(modelk)) legacy_modelk = modelk
        parity = 0.0_dp
        mode_m = [1, 2]
        mode_n = [0, 1]
        ql = [0.5_dp, 0.0_dp]
        do surface_index = 0, 3 + vacuum
            do point = 1, 4
                bjac(point, surface_index) = -1.0_dp &
                    - 0.1_dp * real(point - 1, dp) &
                    - 0.1_dp * real(surface_index, dp)
            end do
        end do
        bjac(:, 0) = 0.0_dp
        write (unit) 3, 2, 2, 1, 1, 2
        write (unit) s, pth, fpp, ftp, current_j, current_i, pressure_slope, &
            parity
        write (unit) mode_m, mode_n, ql
        write (unit) 0.0_dp
        surface_field = 1.25_dp
        full_surface = 1.25_dp
        write (unit) bjac, surface_field, full_surface + 1.0_dp, &
            full_surface + 2.0_dp, full_surface + 3.0_dp
        radial_surface = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
        radial_cell = [0.5_dp, 0.6_dp, 0.7_dp]
        mode_surface = 0.0_dp
        write (unit) 2, mode_m, mode_n, mode_surface, surface_field + 4.0_dp, &
            radial_surface, radial_surface + 0.5_dp, radial_cell, radial_cell, &
            1.0_dp, 2.0_dp
        cell = 0.5_dp
        write (unit) cell, cell + 0.1_dp, cell + 1.0_dp, cell + 0.2_dp, &
            cell + 0.3_dp, cell + 2.0_dp, 0.75_dp, legacy_modelk
    end subroutine write_fixture

    subroutine test_potential_fixture_reader()
        type(terpsichore_matrix_fixture_t) :: fixture
        integer :: info, modelk, unit

        do modelk = 0, 3
            open (newunit=unit, status="scratch", form="unformatted", &
                access="sequential")
            call write_fixture(unit, modelk)
            rewind (unit)
            call read_terpsichore_fixed_boundary_potential_fixture(unit, 0, &
                fixture, info)
            close (unit)
            call require(info == terpsichore_matrix_fixture_ok, &
                "valid TERPSICHORE potential fixture was rejected")
            call require(fixture%legacy_modelk == modelk, &
                "raw TERPSICHORE MODELK is wrong")
        end do
        call require(abs(fixture%flux_p_slope(2) - 0.25_dp) < 1.0e-14_dp, &
            "TERPSICHORE poloidal flux slope is wrong")
        call require(abs(fixture%current_i(2) - 0.1_dp) < 1.0e-14_dp, &
            "TERPSICHORE I current profile is wrong")
        call require(abs(fixture%current_j(2) - 0.6_dp) < 1.0e-14_dp, &
            "TERPSICHORE J current profile is wrong")
        call require(abs(fixture%pressure_slope(2) + 0.01_dp) < 1.0e-14_dp, &
            "TERPSICHORE pressure slope is wrong")
        call require(abs(fixture%flux_t_curve(3) - 0.3_dp) < 1.0e-14_dp, &
            "TERPSICHORE toroidal flux derivative is wrong")
        call require(abs(fixture%flux_p_curve(3) - 0.8_dp) < 1.0e-14_dp, &
            "TERPSICHORE poloidal flux derivative is wrong")
        call require(abs(fixture%signed_bjac_radial(2, 3) - 5.25_dp) &
            < 1.0e-14_dp, "TERPSICHORE Jacobian derivative is wrong")
        call require(abs(fixture%sigma_b_s(1, 0) - 1.25_dp) < 1.0e-14_dp, &
            "TERPSICHORE sigma-Bs field is wrong")
        call require(abs(fixture%metric_ss_over_jacobian(1, 1) - 2.25_dp) &
            < 1.0e-14_dp, "TERPSICHORE GSSL field is wrong")
        call require(abs(fixture%metric_st_over_jacobian(1, 1) - 3.25_dp) &
            < 1.0e-14_dp, "TERPSICHORE GSTL field is wrong")
        call require(abs(fixture%metric_tt_over_jacobian(1, 3) - 4.25_dp) &
            < 1.0e-14_dp, "TERPSICHORE GTTL field is wrong")
        call require(abs(fixture%sigma_b(4, 2) - 1.5_dp) < 1.0e-14_dp, &
            "TERPSICHORE sigma-B field is wrong")
        call require(abs(fixture%parallel_current(1, 1) - 2.5_dp) &
            < 1.0e-14_dp, "TERPSICHORE parallel current is wrong")
        call require(abs(fixture%current_factor - 0.75_dp) < 1.0e-14_dp, &
            "TERPSICHORE current factor is wrong")
    end subroutine test_potential_fixture_reader

    subroutine test_vacuum_potential_fixture_reader()
        type(terpsichore_matrix_fixture_t) :: fixture
        integer :: info, unit

        open (newunit=unit, status="scratch", form="unformatted", &
            access="sequential")
        call write_fixture(unit, vacuum_intervals=1)
        rewind (unit)
        call read_terpsichore_potential_fixture(unit, 1, fixture, info)
        close (unit)
        call require(info == terpsichore_matrix_fixture_ok, &
            "valid vacuum TERPSICHORE fixture was rejected")
        call require(all(shape(fixture%signed_bjac) == [4, 4]) &
            .and. all(shape(fixture%metric_ss_over_jacobian) == [4, 4]), &
            "vacuum fixture retained non-plasma geometry surfaces")
        call require(maxval(abs(fixture%s - [0.0_dp, 0.2_dp, 0.6_dp, &
            1.0_dp])) < 1.0e-14_dp, &
            "vacuum radial extension shifted the plasma grid")
        call require(abs(fixture%flux_p_slope(1) - 0.2_dp) < 1.0e-14_dp &
            .and. abs(fixture%sigma_b_s(1, 0) - 1.25_dp) < 1.0e-14_dp &
            .and. abs(fixture%metric_ss_over_jacobian(1, 0) - 2.25_dp) &
            < 1.0e-14_dp, &
            "vacuum record lengths shifted plasma fields")
    end subroutine test_vacuum_potential_fixture_reader

    subroutine test_fixture_adapter()
        type(terpsichore_matrix_fixture_t) :: fixture
        type(dynamic_family_layout_t) :: actual_layout, expected_layout
        real(dp), allocatable :: actual(:, :), expected(:, :)
        real(dp), allocatable :: normal_phase(:, :, :), tangent_phase(:, :, :)
        real(dp), allocatable :: radial_factor(:, :), radial_weight(:)
        integer :: info

        call build_adapter_fixture(fixture)
        call assemble_terpsichore_fixture_reduced_mass(fixture, actual, &
            actual_layout, info)
        call require(info == terpsichore_reduced_adapter_ok, &
            "TERPSICHORE fixture adapter failed")
        call build_expected_inputs(fixture, normal_phase, tangent_phase, &
            radial_factor, radial_weight)
        call assemble_terpsichore_reduced_fixed_boundary_mass( &
            fixture%signed_bjac(:, 1:), &
            fixture%flux_t_slope(:fixture%intervals), &
            normal_phase, tangent_phase, radial_factor, radial_weight, &
            fixture%mode_m, fixture%mode_n, [phase_sine, phase_sine], &
            expected, expected_layout, info)
        call require(info == 0, "direct reduced fixture assembly failed")
        call require(actual_layout%total_unknowns &
            == expected_layout%total_unknowns, "fixture layouts differ")
        call require(maxval(abs(actual - expected)) < 2.0e-14_dp, &
            "TERPSICHORE fixture adapter matrix is wrong")
    end subroutine test_fixture_adapter

    subroutine test_free_boundary_adapter()
        type(terpsichore_matrix_fixture_t) :: fixture
        type(dynamic_family_layout_t) :: fixed_layout, free_layout
        real(dp), allocatable :: fixed(:, :), free(:, :)
        integer, allocatable :: retained(:)
        integer :: info

        call build_adapter_fixture(fixture)
        call assemble_terpsichore_fixture_reduced_mass(fixture, fixed, &
            fixed_layout, info)
        call require(info == terpsichore_reduced_adapter_ok, &
            "fixed mass for free-boundary comparison failed")
        call assemble_terpsichore_fixture_reduced_mass_free_boundary( &
            fixture, free, free_layout, info)
        call require(info == terpsichore_reduced_adapter_ok, &
            "free-boundary reduced mass failed")
        call require(free_layout%outer_normal_retained &
            .and. free_layout%normal_unknowns == 6, &
            "free-boundary reduced mass layout is wrong")
        retained = [1, 2, 3, 4, 7, 8, 9, 10, 11, 12]
        call require(maxval(abs(free(retained, retained) - fixed)) &
            < 2.0e-14_dp, &
            "fixed reduced mass is not the constrained free mass")
        call require(maxval(abs(free - transpose(free))) < 2.0e-14_dp, &
            "free-boundary reduced mass is not symmetric")
    end subroutine test_free_boundary_adapter

    subroutine test_fixture_rejections()
        type(terpsichore_matrix_fixture_t) :: fixture
        type(dynamic_family_layout_t) :: layout
        real(dp), allocatable :: mass(:, :)
        integer :: info, unit

        open (newunit=unit, status="scratch", form="unformatted", &
            access="sequential")
        call write_fixture(unit)
        rewind (unit)
        call read_terpsichore_fixed_boundary_fixture(unit, 1, fixture, info)
        close (unit)
        call require(info == terpsichore_matrix_fixture_invalid, &
            "nonzero IVAC was accepted by the fixed-boundary reader")
        call require_header_rejected([997, 2, 2, 1, 1, 2], &
            "oversized external dimensions were accepted")
        call require_header_rejected([2, 999, 999, 1, 1, 3], &
            "oversized phase storage was accepted")
        call require_header_rejected([996, 1, 1, 1, 1, 3], &
            "oversized dense matrix was accepted")
        fixture = terpsichore_matrix_fixture_t(intervals=64, &
            poloidal_points=300, toroidal_points=100, stability_periods=1, &
            field_periods=1, modes=2)
        call require(.not. terpsichore_potential_metadata_is_valid(fixture), &
            "oversized aggregate potential storage was accepted")
        call build_adapter_fixture(fixture)
        fixture%stability_periods = 2
        call assemble_terpsichore_fixture_reduced_mass(fixture, mass, &
            layout, info)
        call require(info /= terpsichore_reduced_adapter_ok, &
            "incommensurate stability periods were accepted")
        fixture%stability_periods = 1
        fixture%s(fixture%intervals) = 0.9_dp
        call assemble_terpsichore_fixture_reduced_mass(fixture, mass, &
            layout, info)
        call require(info /= terpsichore_reduced_adapter_ok, &
            "incomplete normalized radial domain was accepted")
    end subroutine test_fixture_rejections

    subroutine require_header_rejected(header, message)
        integer, intent(in) :: header(6)
        character(len=*), intent(in) :: message
        type(terpsichore_matrix_fixture_t) :: fixture
        integer :: info, unit

        open (newunit=unit, status="scratch", form="unformatted", &
            access="sequential")
        write (unit) header
        rewind (unit)
        call read_terpsichore_fixed_boundary_fixture(unit, 0, fixture, info)
        close (unit)
        call require(info == terpsichore_matrix_fixture_invalid, message)
    end subroutine require_header_rejected

    subroutine build_adapter_fixture(fixture)
        type(terpsichore_matrix_fixture_t), intent(out) :: fixture
        integer :: point, surface

        fixture%intervals = 3
        fixture%poloidal_points = 4
        fixture%toroidal_points = 3
        fixture%stability_periods = 1
        fixture%field_periods = 1
        fixture%modes = 2
        fixture%parity = 0.0_dp
        allocate (fixture%s(0:3))
        fixture%s = [0.0_dp, 0.2_dp, 0.6_dp, 1.0_dp]
        fixture%flux_t_slope = [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp]
        fixture%mode_m = [1, 2]
        fixture%mode_n = [0, 1]
        fixture%radial_power = [0.5_dp, 0.0_dp]
        allocate (fixture%signed_bjac(12, 0:3))
        do surface = 0, 3
            do point = 1, 12
                fixture%signed_bjac(point, surface) = -1.0_dp &
                    - 0.01_dp * real(point + surface, dp)
            end do
        end do
        fixture%signed_bjac(:, 0) = 0.0_dp
    end subroutine build_adapter_fixture

    subroutine build_expected_inputs(fixture, normal_phase, tangent_phase, &
            radial_factor, radial_weight)
        type(terpsichore_matrix_fixture_t), intent(in) :: fixture
        real(dp), allocatable, intent(out) :: normal_phase(:, :, :)
        real(dp), allocatable, intent(out) :: tangent_phase(:, :, :)
        real(dp), allocatable, intent(out) :: radial_factor(:, :)
        real(dp), allocatable, intent(out) :: radial_weight(:)
        real(dp) :: angle, midpoint, theta, zeta
        integer :: interval, j, k, mode, point

        allocate (normal_phase(2, 12, 3), tangent_phase(2, 12, 3))
        allocate (radial_factor(2, 3), radial_weight(3))
        do interval = 1, 3
            midpoint = 0.5_dp * (fixture%s(interval - 1) + fixture%s(interval))
            radial_factor(:, interval) = midpoint**(-fixture%radial_power)
            radial_weight(interval) = 3.0_dp &
                * (fixture%s(interval) - fixture%s(interval - 1))
            do k = 1, 3
                zeta = 2.0_dp * acos(-1.0_dp) * real(k - 1, dp) / 3.0_dp
                do j = 1, 4
                    theta = 2.0_dp * acos(-1.0_dp) * real(j - 1, dp) / 4.0_dp
                    point = j + 4 * (k - 1)
                    do mode = 1, 2
                        angle = real(fixture%mode_m(mode), dp) * theta &
                            - real(fixture%mode_n(mode), dp) * zeta
                        normal_phase(mode, point, interval) = sin(angle)
                        tangent_phase(mode, point, interval) = cos(angle)
                    end do
                end do
            end do
        end do
    end subroutine build_expected_inputs

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_reduced_mass_adapter
