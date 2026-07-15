program test_terpsichore_reduced_mass_family_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: build_dynamic_element_map, &
        dynamic_family_layout_t
    use terpsichore_reduced_mass, only: &
        assemble_terpsichore_reduced_mass_element, &
        terpsichore_reduced_element_energy
    use terpsichore_reduced_mass_family_assembly, only: &
        assemble_terpsichore_reduced_fixed_boundary_mass, &
        assemble_terpsichore_reduced_family_mass_fixed_layout, &
        terpsichore_reduced_family_invalid, &
        terpsichore_reduced_family_ok
    implicit none

    integer, parameter :: intervals = 4, modes = 2, points = 32
    integer, parameter :: trial_m(modes) = [0, 0]
    integer, parameter :: trial_n(modes) = [0, 0]
    integer, parameter :: trial_parity(modes) = [1, 2]
    real(dp) :: signed_bjac(points, intervals), flux_t_slope(intervals)
    real(dp) :: normal_phase(modes, points, intervals)
    real(dp) :: tangential_phase(modes, points, intervals)
    real(dp) :: radial_factor(modes, intervals), radial_weight(intervals)
    real(dp), allocatable :: mass(:, :)
    type(dynamic_family_layout_t) :: layout
    integer :: info

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*)
            real(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    call build_fixture()
    call assemble_terpsichore_reduced_fixed_boundary_mass(signed_bjac, &
        flux_t_slope, normal_phase, tangential_phase, radial_factor, &
        radial_weight, trial_m, trial_n, trial_parity, mass, layout, info)
    call require(info == terpsichore_reduced_family_ok, &
        "valid reduced family assembly failed")
    call require(layout%normal_unknowns == intervals - 1, &
        "reduced family normal boundary elimination is wrong")
    call require(layout%eta_unknowns == intervals, &
        "reduced family tangential ownership is wrong")
    call require(layout%mu_unknowns == 0, &
        "reduced family retained the absent compressional component")
    call require(maxval(abs(mass - transpose(mass))) < 2.0e-14_dp, &
        "reduced family mass is not symmetric")
    call check_independent_gather(mass, layout)
    call check_positive_definite(mass)
    call check_fixed_rejection(layout)
    call check_invalid_partition(mass)

    write (*, "(a)") "PASS"

contains

    subroutine build_fixture()
        real(dp) :: angle, s
        integer :: interval, point

        do interval = 1, intervals
            s = (real(interval, dp) - 0.5_dp) / real(intervals, dp)
            flux_t_slope(interval) = 1.1_dp + 0.1_dp * s
            radial_weight(interval) = 0.4_dp + 0.2_dp * s
            radial_factor(1, interval) = s**0.25_dp
            radial_factor(2, interval) = s**0.5_dp
            do point = 1, points
                angle = 2.0_dp * acos(-1.0_dp) * real(point - 1, dp) &
                    / real(points, dp)
                signed_bjac(point, interval) = -(1.2_dp + 0.1_dp * s &
                    + 0.07_dp * cos(angle) - 0.03_dp * sin(2.0_dp * angle))
                normal_phase(:, point, interval) = [1.0_dp, 0.0_dp]
                tangential_phase(:, point, interval) = [0.0_dp, 1.0_dp]
            end do
        end do
    end subroutine build_fixture

    subroutine check_independent_gather(actual, local_layout)
        real(dp), intent(in) :: actual(:, :)
        type(dynamic_family_layout_t), intent(in) :: local_layout
        real(dp) :: expected(size(actual, 1), size(actual, 2))
        real(dp) :: displacement(size(actual, 1)), local_displacement(3 * modes)
        real(dp) :: direct_energy, global_energy, image
        real(dp), allocatable :: element(:, :)
        integer, allocatable :: element_map(:, :)
        integer, parameter :: expected_map(3 * modes, intervals) = reshape([ &
            0, 0, 1, 0, 0, 4, &
            1, 0, 2, 0, 0, 5, &
            2, 0, 3, 0, 0, 6, &
            3, 0, 0, 0, 0, 7], [3 * modes, intervals])
        integer :: a, b, ga, gb, interval, local_info, row

        expected = 0.0_dp
        call build_dynamic_element_map(local_layout, element_map, local_info)
        call require(local_info == 0, "reduced family element map failed")
        call require(all(element_map(:3 * modes, :) == expected_map), &
            "reduced family activity or boundary map is wrong")
        do interval = 1, intervals
            call assemble_terpsichore_reduced_mass_element( &
                signed_bjac(:, interval), flux_t_slope(interval), &
                normal_phase(:, :, interval), &
                tangential_phase(:, :, interval), &
                radial_factor(:, interval), radial_weight(interval), &
                element, local_info)
            call require(local_info == 0, "reduced family oracle element failed")
            do b = 1, 3 * modes
                gb = element_map(b, interval)
                if (gb == 0) cycle
                do a = 1, 3 * modes
                    ga = element_map(a, interval)
                    if (ga == 0) cycle
                    expected(ga, gb) = expected(ga, gb) + element(a, b)
                end do
            end do
        end do
        call require(maxval(abs(actual - expected)) < 2.0e-14_dp, &
            "reduced family gather disagrees with the element oracle")
        do a = 1, size(displacement)
            displacement(a) = 0.1_dp * real(a, dp)
        end do
        direct_energy = 0.0_dp
        do interval = 1, intervals
            local_displacement = 0.0_dp
            do a = 1, 3 * modes
                ga = expected_map(a, interval)
                if (ga > 0) local_displacement(a) = displacement(ga)
            end do
            direct_energy = direct_energy + terpsichore_reduced_element_energy( &
                signed_bjac(:, interval), flux_t_slope(interval), &
                normal_phase(:, :, interval), &
                tangential_phase(:, :, interval), &
                radial_factor(:, interval), radial_weight(interval), &
                local_displacement)
        end do
        global_energy = 0.0_dp
        do row = 1, size(actual, 1)
            image = 0.0_dp
            do a = 1, size(actual, 2)
                image = image + actual(row, a) * displacement(a)
            end do
            global_energy = global_energy + displacement(row) * image
        end do
        global_energy = 0.5_dp * global_energy
        call require(abs(global_energy - direct_energy) < 2.0e-13_dp, &
            "reduced family energy disagrees with the cell sum")
    end subroutine check_independent_gather

    subroutine check_positive_definite(matrix)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: copy(size(matrix, 1), size(matrix, 2))
        real(dp) :: eigenvalues(size(matrix, 1)), work(8 * size(matrix, 1))
        integer :: local_info

        copy = matrix
        call dsyev("N", "U", size(copy, 1), copy, size(copy, 1), &
            eigenvalues, work, size(work), local_info)
        call require(local_info == 0, "reduced family eigensolve failed")
        call require(eigenvalues(1) > 1.0e-10_dp, &
            "boundary-constrained reduced family mass is not positive definite")
    end subroutine check_positive_definite

    subroutine check_fixed_rejection(local_layout)
        type(dynamic_family_layout_t), intent(in) :: local_layout
        real(dp) :: rejected(local_layout%total_unknowns, &
            local_layout%total_unknowns)
        real(dp) :: bad_slope(intervals), bad_weight(intervals)
        integer, allocatable :: full_map(:, :)
        integer :: bad_map(3 * modes, intervals), local_info

        call build_dynamic_element_map(local_layout, full_map, local_info)
        call require(local_info == 0, "reduced family rejection map failed")
        bad_map = full_map(:3 * modes, :)
        bad_slope = flux_t_slope
        bad_slope(2) = 0.0_dp
        call assemble_terpsichore_reduced_family_mass_fixed_layout( &
            signed_bjac, bad_slope, normal_phase, tangential_phase, &
            radial_factor, radial_weight, bad_map, rejected, local_info)
        call require(local_info == terpsichore_reduced_family_invalid, &
            "zero flux slope was accepted by the family gather")
        bad_slope = flux_t_slope
        bad_weight = radial_weight
        bad_weight(1) = -1.0_dp
        call assemble_terpsichore_reduced_family_mass_fixed_layout( &
            signed_bjac, bad_slope, normal_phase, tangential_phase, &
            radial_factor, bad_weight, bad_map, rejected, local_info)
        call require(local_info == terpsichore_reduced_family_invalid, &
            "negative radial weight was accepted by the family gather")
        bad_weight = radial_weight
        bad_map(1, 1) = -1
        call assemble_terpsichore_reduced_family_mass_fixed_layout( &
            signed_bjac, bad_slope, normal_phase, tangential_phase, &
            radial_factor, bad_weight, bad_map, rejected, local_info)
        call require(local_info == terpsichore_reduced_family_invalid, &
            "malformed map was accepted by the family gather")
    end subroutine check_fixed_rejection

    subroutine check_invalid_partition(original)
        real(dp), intent(in) :: original(:, :)
        real(dp), allocatable :: rejected(:, :)
        type(dynamic_family_layout_t) :: rejected_layout
        integer :: local_info

        call assemble_terpsichore_reduced_fixed_boundary_mass( &
            signed_bjac(:, :1), &
            flux_t_slope(:1), normal_phase(:, :, :1), &
            tangential_phase(:, :, :1), radial_factor(:, :1), &
            radial_weight(:1), trial_m, trial_n, trial_parity, rejected, &
            rejected_layout, local_info)
        call require(local_info /= terpsichore_reduced_family_ok, &
            "single-cell reduced family was accepted")
        call require(size(original, 1) == layout%total_unknowns, &
            "valid reduced family matrix shape changed")
    end subroutine check_invalid_partition

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_reduced_mass_family_assembly
