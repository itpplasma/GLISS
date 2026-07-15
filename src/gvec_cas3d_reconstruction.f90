module gvec_cas3d_reconstruction
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_types, only: harmonic_pair_t
    implicit none
    private

    integer, parameter, public :: reconstruction_ok = 0
    integer, parameter, public :: reconstruction_invalid_surface = 1
    integer, parameter, public :: reconstruction_shape_mismatch = 2
    integer, parameter, public :: reconstruction_nonfinite_input = 3
    integer, parameter, public :: reconstruction_invalid_spacing = 4

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

    public :: reconstruct_harmonic_grid
    public :: periodic_sixth_order_derivatives
    public :: project_harmonic_grid

    interface
        pure subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(dp), intent(in) :: alpha, beta
            complex(dp), intent(in) :: a(lda, *), b(ldb, *)
            complex(dp), intent(inout) :: c(ldc, *)
        end subroutine zgemm
    end interface

contains

    pure subroutine reconstruct_harmonic_grid(pair, radial_surface, &
            poloidal_modes, toroidal_modes, theta, zeta_period, values, &
            derivative_theta, derivative_zeta_period, info, &
            derivative_theta_theta, derivative_theta_zeta_period, &
            derivative_zeta_period_zeta_period)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: radial_surface
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        real(dp), intent(in) :: theta(:), zeta_period(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        real(dp), allocatable, intent(out) :: derivative_theta(:, :)
        real(dp), allocatable, intent(out) :: derivative_zeta_period(:, :)
        integer, intent(out) :: info
        real(dp), allocatable, intent(out), optional :: &
            derivative_theta_theta(:, :)
        real(dp), allocatable, intent(out), optional :: &
            derivative_theta_zeta_period(:, :)
        real(dp), allocatable, intent(out), optional :: &
            derivative_zeta_period_zeta_period(:, :)
        complex(dp), allocatable :: coefficients(:, :), theta_phase(:, :)
        complex(dp), allocatable :: theta_derivative_phase(:, :)
        complex(dp), allocatable :: zeta_phase(:, :)
        complex(dp), allocatable :: zeta_derivative_phase(:, :)
        complex(dp), allocatable :: theta_partial(:, :), work_partial(:, :)
        complex(dp), allocatable :: complex_grid(:, :)

        call validate_inputs(pair, radial_surface, poloidal_modes, &
            toroidal_modes, theta, zeta_period, info)
        if (info /= reconstruction_ok) return
        coefficients = cmplx(pair%cosine(radial_surface, :, :), &
            -pair%sine(radial_surface, :, :), kind=dp)
        call build_phase_matrices(poloidal_modes, toroidal_modes, theta, &
            zeta_period, theta_phase, theta_derivative_phase, zeta_phase, &
            zeta_derivative_phase)
        allocate (theta_partial(size(theta), size(toroidal_modes)))
        allocate (work_partial(size(theta), size(toroidal_modes)))
        allocate (complex_grid(size(theta), size(zeta_period)))
        call multiply_complex(theta_phase, coefficients, theta_partial)
        call multiply_complex(theta_partial, zeta_phase, complex_grid)
        call copy_real_grid(complex_grid, values)
        call multiply_complex(theta_derivative_phase, coefficients, &
            work_partial)
        call multiply_complex(work_partial, zeta_phase, complex_grid)
        call copy_real_grid(complex_grid, derivative_theta)
        call multiply_complex(theta_partial, zeta_derivative_phase, &
            complex_grid)
        call copy_real_grid(complex_grid, derivative_zeta_period)
        call reconstruct_second_angular_derivatives(coefficients, &
            poloidal_modes, toroidal_modes, theta_phase, &
            theta_derivative_phase, zeta_phase, zeta_derivative_phase, &
            theta_partial, work_partial, complex_grid, &
            derivative_theta_theta, derivative_theta_zeta_period, &
            derivative_zeta_period_zeta_period)
        info = reconstruction_ok
    end subroutine reconstruct_harmonic_grid

    pure subroutine reconstruct_second_angular_derivatives(coefficients, &
            poloidal_modes, toroidal_modes, theta_phase, &
            theta_derivative_phase, zeta_phase, zeta_derivative_phase, &
            theta_partial, work_partial, complex_grid, &
            derivative_theta_theta, derivative_theta_zeta_period, &
            derivative_zeta_period_zeta_period)
        complex(dp), contiguous, intent(in) :: coefficients(:, :)
        complex(dp), contiguous, intent(in) :: theta_phase(:, :)
        complex(dp), contiguous, intent(in) :: theta_derivative_phase(:, :)
        complex(dp), contiguous, intent(in) :: zeta_phase(:, :)
        complex(dp), contiguous, intent(in) :: zeta_derivative_phase(:, :)
        complex(dp), contiguous, intent(in) :: theta_partial(:, :)
        complex(dp), contiguous, intent(inout) :: work_partial(:, :)
        complex(dp), contiguous, intent(inout) :: complex_grid(:, :)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        real(dp), allocatable, intent(out), optional :: &
            derivative_theta_theta(:, :)
        real(dp), allocatable, intent(out), optional :: &
            derivative_theta_zeta_period(:, :)
        real(dp), allocatable, intent(out), optional :: &
            derivative_zeta_period_zeta_period(:, :)
        complex(dp), allocatable :: second_phase(:, :)
        integer :: mode

        if (present(derivative_theta_theta)) then
            allocate (second_phase, mold=theta_phase)
            do mode = 1, size(poloidal_modes)
                second_phase(:, mode) = -(two_pi &
                    * real(poloidal_modes(mode), dp))**2 &
                    * theta_phase(:, mode)
            end do
            call multiply_complex(second_phase, coefficients, work_partial)
            call multiply_complex(work_partial, zeta_phase, complex_grid)
            call copy_real_grid(complex_grid, derivative_theta_theta)
            deallocate (second_phase)
        end if
        if (present(derivative_theta_zeta_period)) then
            call multiply_complex(theta_derivative_phase, coefficients, &
                work_partial)
            call multiply_complex(work_partial, zeta_derivative_phase, &
                complex_grid)
            call copy_real_grid(complex_grid, derivative_theta_zeta_period)
        end if
        if (present(derivative_zeta_period_zeta_period)) then
            allocate (second_phase, mold=zeta_phase)
            do mode = 1, size(toroidal_modes)
                second_phase(mode, :) = -(two_pi &
                    * real(toroidal_modes(mode), dp))**2 &
                    * zeta_phase(mode, :)
            end do
            call multiply_complex(theta_partial, second_phase, complex_grid)
            call copy_real_grid(complex_grid, &
                derivative_zeta_period_zeta_period)
        end if
    end subroutine reconstruct_second_angular_derivatives

    pure subroutine multiply_complex(left, right, product)
        complex(dp), contiguous, intent(in) :: left(:, :), right(:, :)
        complex(dp), contiguous, intent(out) :: product(:, :)
        complex(dp), parameter :: one = cmplx(1.0_dp, 0.0_dp, dp)
        complex(dp), parameter :: zero = cmplx(0.0_dp, 0.0_dp, dp)

        call zgemm("N", "N", size(left, 1), size(right, 2), &
            size(left, 2), one, left, size(left, 1), right, size(right, 1), &
            zero, product, size(product, 1))
    end subroutine multiply_complex

    pure subroutine copy_real_grid(source, target)
        complex(dp), intent(in) :: source(:, :)
        real(dp), allocatable, intent(out) :: target(:, :)
        integer :: column, row

        allocate (target(size(source, 1), size(source, 2)))
        do column = 1, size(source, 2)
            do row = 1, size(source, 1)
                target(row, column) = real(source(row, column), dp)
            end do
        end do
    end subroutine copy_real_grid

    pure subroutine project_harmonic_grid(values, poloidal_modes, &
            toroidal_modes, theta, zeta_period, cosine, sine)
        real(dp), intent(in) :: values(:, :)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        real(dp), intent(in) :: theta(:), zeta_period(:)
        real(dp), intent(out) :: cosine(:, :), sine(:, :)
        real(dp) :: phase, cosine_sum, sine_sum, weight
        integer :: mode_m, mode_n, m, n, j, k

        do mode_n = 1, size(toroidal_modes)
            n = toroidal_modes(mode_n)
            do mode_m = 1, size(poloidal_modes)
                m = poloidal_modes(mode_m)
                cosine(mode_m, mode_n) = 0.0_dp
                sine(mode_m, mode_n) = 0.0_dp
                if (m == 0 .and. n < 0) cycle
                cosine_sum = 0.0_dp
                sine_sum = 0.0_dp
                do k = 1, size(zeta_period)
                    do j = 1, size(theta)
                        phase = two_pi * (real(m, dp) * theta(j) &
                            - real(n, dp) * zeta_period(k))
                        cosine_sum = cosine_sum + values(j, k) * cos(phase)
                        sine_sum = sine_sum + values(j, k) * sin(phase)
                    end do
                end do
                weight = 2.0_dp
                if (m == 0 .and. n == 0) weight = 1.0_dp
                weight = weight / real(size(theta) * size(zeta_period), dp)
                cosine(mode_m, mode_n) = weight * cosine_sum
                sine(mode_m, mode_n) = weight * sine_sum
            end do
        end do
    end subroutine project_harmonic_grid

    pure subroutine build_phase_matrices(poloidal_modes, toroidal_modes, &
            theta, zeta_period, theta_phase, theta_derivative_phase, &
            zeta_phase, zeta_derivative_phase)
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        real(dp), intent(in) :: theta(:), zeta_period(:)
        complex(dp), allocatable, intent(out) :: theta_phase(:, :)
        complex(dp), allocatable, intent(out) :: theta_derivative_phase(:, :)
        complex(dp), allocatable, intent(out) :: zeta_phase(:, :)
        complex(dp), allocatable, intent(out) :: zeta_derivative_phase(:, :)
        complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)
        integer :: mode, point

        allocate (theta_phase(size(theta), size(poloidal_modes)))
        allocate (theta_derivative_phase(size(theta), size(poloidal_modes)))
        allocate (zeta_phase(size(toroidal_modes), size(zeta_period)))
        allocate (zeta_derivative_phase(size(toroidal_modes), &
            size(zeta_period)))
        do mode = 1, size(poloidal_modes)
            do point = 1, size(theta)
                theta_phase(point, mode) = exp(imaginary_unit * two_pi * &
                    real(poloidal_modes(mode), dp) * theta(point))
                theta_derivative_phase(point, mode) = imaginary_unit * &
                    two_pi * real(poloidal_modes(mode), dp) * &
                    theta_phase(point, mode)
            end do
        end do
        do point = 1, size(zeta_period)
            do mode = 1, size(toroidal_modes)
                zeta_phase(mode, point) = exp(-imaginary_unit * two_pi * &
                    real(toroidal_modes(mode), dp) * zeta_period(point))
                zeta_derivative_phase(mode, point) = -imaginary_unit * &
                    two_pi * real(toroidal_modes(mode), dp) * &
                    zeta_phase(mode, point)
            end do
        end do
    end subroutine build_phase_matrices

    pure subroutine validate_inputs(pair, radial_surface, poloidal_modes, &
            toroidal_modes, theta, zeta_period, info)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: radial_surface
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        real(dp), intent(in) :: theta(:), zeta_period(:)
        integer, intent(out) :: info

        info = reconstruction_shape_mismatch
        if (.not. allocated(pair%cosine)) return
        if (.not. allocated(pair%sine)) return
        if (any(shape(pair%cosine) /= shape(pair%sine))) return
        if (size(pair%cosine, 2) /= size(poloidal_modes)) return
        if (size(pair%cosine, 3) /= size(toroidal_modes)) return
        if (size(poloidal_modes) < 1 .or. size(toroidal_modes) < 1) return
        if (size(theta) < 1) return
        if (size(zeta_period) < 1) return
        info = reconstruction_invalid_surface
        if (radial_surface < 1) return
        if (radial_surface > size(pair%cosine, 1)) return
        info = reconstruction_nonfinite_input
        if (.not. all(ieee_is_finite(pair%cosine))) return
        if (.not. all(ieee_is_finite(pair%sine))) return
        if (.not. all(ieee_is_finite(theta))) return
        if (.not. all(ieee_is_finite(zeta_period))) return
        info = reconstruction_ok
    end subroutine validate_inputs

    pure subroutine periodic_sixth_order_derivatives(values, theta_spacing, &
            zeta_spacing, derivative_theta, derivative_zeta, info)
        real(dp), intent(in) :: values(:, :), theta_spacing, zeta_spacing
        real(dp), allocatable, intent(out) :: derivative_theta(:, :)
        real(dp), allocatable, intent(out) :: derivative_zeta(:, :)
        integer, intent(out) :: info
        integer :: theta, zeta
        integer :: tm3, tm2, tm1, tp1, tp2, tp3
        integer :: zm3, zm2, zm1, zp1, zp2, zp3

        call validate_periodic_derivative_inputs(values, theta_spacing, &
            zeta_spacing, info)
        if (info /= reconstruction_ok) return
        allocate (derivative_theta(size(values, 1), size(values, 2)))
        allocate (derivative_zeta(size(values, 1), size(values, 2)))
        do theta = 1, size(values, 1)
            tm3 = modulo(theta - 4, size(values, 1)) + 1
            tm2 = modulo(theta - 3, size(values, 1)) + 1
            tm1 = modulo(theta - 2, size(values, 1)) + 1
            tp1 = modulo(theta, size(values, 1)) + 1
            tp2 = modulo(theta + 1, size(values, 1)) + 1
            tp3 = modulo(theta + 2, size(values, 1)) + 1
            do zeta = 1, size(values, 2)
                derivative_theta(theta, zeta) = &
                    (-values(tm3, zeta) + 9.0_dp * values(tm2, zeta) &
                    - 45.0_dp * values(tm1, zeta) &
                    + 45.0_dp * values(tp1, zeta) &
                    - 9.0_dp * values(tp2, zeta) &
                    + values(tp3, zeta)) / (60.0_dp * theta_spacing)
            end do
        end do
        do zeta = 1, size(values, 2)
            zm3 = modulo(zeta - 4, size(values, 2)) + 1
            zm2 = modulo(zeta - 3, size(values, 2)) + 1
            zm1 = modulo(zeta - 2, size(values, 2)) + 1
            zp1 = modulo(zeta, size(values, 2)) + 1
            zp2 = modulo(zeta + 1, size(values, 2)) + 1
            zp3 = modulo(zeta + 2, size(values, 2)) + 1
            do theta = 1, size(values, 1)
                derivative_zeta(theta, zeta) = &
                    (-values(theta, zm3) + 9.0_dp * values(theta, zm2) &
                    - 45.0_dp * values(theta, zm1) &
                    + 45.0_dp * values(theta, zp1) &
                    - 9.0_dp * values(theta, zp2) &
                    + values(theta, zp3)) / (60.0_dp * zeta_spacing)
            end do
        end do
        info = reconstruction_ok
    end subroutine periodic_sixth_order_derivatives

    pure subroutine validate_periodic_derivative_inputs(values, theta_spacing, &
            zeta_spacing, info)
        real(dp), intent(in) :: values(:, :), theta_spacing, zeta_spacing
        integer, intent(out) :: info

        info = reconstruction_shape_mismatch
        if (size(values, 1) < 7) return
        if (size(values, 2) < 7) return
        info = reconstruction_nonfinite_input
        if (.not. all(ieee_is_finite(values))) return
        if (.not. ieee_is_finite(theta_spacing)) return
        if (.not. ieee_is_finite(zeta_spacing)) return
        info = reconstruction_invalid_spacing
        if (theta_spacing <= 0.0_dp) return
        if (zeta_spacing <= 0.0_dp) return
        info = reconstruction_ok
    end subroutine validate_periodic_derivative_inputs

end module gvec_cas3d_reconstruction
