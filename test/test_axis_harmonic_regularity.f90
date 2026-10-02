program test_axis_harmonic_regularity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use axis_regular_harmonic_spline, only: axis_regular_harmonic_field_t, &
        axis_regular_harmonic_ok, evaluate_axis_regular_harmonics, &
        fit_axis_regular_harmonics
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        radial_cubic_spline_grid_t, radial_cubic_spline_ok
    implicit none

    call check_axis_jets()
    call check_convergence()
    call check_roundoff()
    print *, 'PASS'

contains

    subroutine check_axis_jets()
        integer, parameter :: modes(6) = [4, 5, 9, 10, 24, 36]
        real(dp), parameter :: nodes(6) = [0.04_dp, 0.17_dp, 0.36_dp, &
            0.58_dp, 0.79_dp, 0.96_dp]
        type(radial_cubic_spline_grid_t) :: grid
        type(axis_regular_harmonic_field_t) :: field
        real(dp) :: samples(6, 6), values(6), slopes(6), seconds(6), exponent
        real(dp) :: left(6, 3), right(6, 3), offset
        integer :: column, info

        call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
        call require(info == radial_cubic_spline_ok, 'axis grid failed')
        do column = 1, size(modes)
            exponent = 0.5_dp * real(modes(column), dp)
            samples(:, column) = nodes**exponent * (1.0_dp + nodes**3)
        end do
        call fit_axis_regular_harmonics(grid, modes, samples, field, info)
        call require(info == axis_regular_harmonic_ok, 'axis fit failed')
        call evaluate_axis_regular_harmonics(grid, field, 0.0_dp, &
            values, slopes, seconds, info)
        call require(info == axis_regular_harmonic_ok, 'finite axis jets rejected')
        call require(all(values == 0.0_dp), 'nonzero high-mode axis value')
        call require(all(slopes == 0.0_dp), 'nonzero smooth harmonic axis slope')
        call require(all(seconds(2:) == 0.0_dp), 'nonzero high-mode axis curvature')
        ! A C2 join makes the primitive metric and its first radial derivative
        ! continuous.  Check both sides with an offset far below the data scale.
        offset = 1.0e-10_dp * nodes(1)
        call evaluate_axis_regular_harmonics(grid, field, nodes(1) - offset, &
            left(:, 1), left(:, 2), left(:, 3), info)
        call evaluate_axis_regular_harmonics(grid, field, nodes(1) + offset, &
            right(:, 1), right(:, 2), right(:, 3), info)
        call require(maxval(abs(left - right)) < 1.0e-8_dp, 'axis patch is not C2')
    end subroutine check_axis_jets

    subroutine check_convergence()
        integer, parameter :: modes(4) = [4, 10, 24, 36]
        integer, parameter :: meshes(3) = [64, 128, 256]
        type(radial_cubic_spline_grid_t) :: grid
        type(axis_regular_harmonic_field_t) :: field
        real(dp), allocatable :: nodes(:), samples(:, :)
        real(dp) :: values(4), slopes(4), seconds(4), exact(3), errors(4, 3, 3)
        real(dp) :: exponent, coordinate, rates(3), minimum_rates(3)
        integer :: mesh, column, query, info, n

        errors = 0.0_dp
        minimum_rates = [3.3_dp, 2.3_dp, 1.3_dp]
        do mesh = 1, size(meshes)
            n = meshes(mesh)
            allocate (nodes(n), samples(n, size(modes)))
            do query = 1, n
                nodes(query) = (real(query, dp) - 0.5_dp) / real(n, dp)
            end do
            do column = 1, size(modes)
                exponent = 0.5_dp * real(modes(column), dp)
                samples(:, column) = nodes**exponent * (1.0_dp + nodes**3)
            end do
            call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
            call fit_axis_regular_harmonics(grid, modes, samples, field, info)
            call require(info == axis_regular_harmonic_ok, 'convergence fit failed')
            do query = 1, 4 * n
                coordinate = real(query, dp) / real(4 * n, dp)
                call evaluate_axis_regular_harmonics(grid, field, coordinate, &
                    values, slopes, seconds, info)
                call require(info == axis_regular_harmonic_ok, 'jet evaluation failed')
                do column = 1, size(modes)
                    exponent = 0.5_dp * real(modes(column), dp)
                    ! Independent derivatives of r^m(1+r^6), s=r^2.
                    exact(1) = coordinate**exponent * (1.0_dp + coordinate**3)
                    exact(2) = exponent * coordinate**(exponent - 1.0_dp) &
                        + (exponent + 3.0_dp) * coordinate**(exponent + 2.0_dp)
                    exact(3) = exponent * (exponent - 1.0_dp) &
                        * coordinate**(exponent - 2.0_dp) &
                        + (exponent + 3.0_dp) * (exponent + 2.0_dp) &
                        * coordinate**(exponent + 1.0_dp)
                    errors(column, :, mesh) = max(errors(column, :, mesh), &
                        abs([values(column), slopes(column), seconds(column)] - exact))
                end do
            end do
            deallocate (nodes, samples)
        end do
        do column = 1, size(modes)
            do mesh = 1, size(meshes) - 1
                rates = log(errors(column, :, mesh) / errors(column, :, mesh + 1)) &
                    / log(2.0_dp)
                print '(a,i2,a,3f8.3)', 'm=', modes(column), ' jet orders ', rates
                call require(all(rates > minimum_rates), 'harmonic jet loses convergence')
            end do
        end do
    end subroutine check_convergence

    subroutine check_roundoff()
        integer, parameter :: n = 64, modes(4) = [16, 24, 28, 36]
        type(radial_cubic_spline_grid_t) :: grid
        type(axis_regular_harmonic_field_t) :: field
        real(dp) :: nodes(n), samples(n, 4), values(4), slopes(4), seconds(4)
        real(dp) :: maxima(3), coordinate
        integer :: node, column, query, info

        do node = 1, n
            nodes(node) = (real(node, dp) - 0.5_dp) / real(n, dp)
            do column = 1, size(modes)
                samples(node, column) = 1.0e-16_dp &
                    * sin(real(7 * node + column, dp))
            end do
        end do
        call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
        call fit_axis_regular_harmonics(grid, modes, samples, field, info)
        maxima = 0.0_dp
        do query = 0, 2048
            coordinate = real(query, dp) / 2048.0_dp
            call evaluate_axis_regular_harmonics(grid, field, coordinate, &
                values, slopes, seconds, info)
            call require(info == axis_regular_harmonic_ok, 'roundoff jet rejected')
            maxima = max(maxima, [maxval(abs(values)), &
                maxval(abs(slopes)) / real(n, dp), &
                maxval(abs(seconds)) / real(n, dp)**2])
        end do
        print '(a,3es12.3)', 'scaled roundoff jets ', maxima
        call require(maxima(1) < 1.0e-15_dp, 'roundoff value amplified')
        call require(maxval(maxima(2:)) < 1.0e-13_dp, 'roundoff derivatives amplified')
    end subroutine check_roundoff

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_axis_harmonic_regularity
