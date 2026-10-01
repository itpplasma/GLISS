! Sine-parity forms from the cosine-parity angular products.
!
! A sine-parity trial is its cosine-parity partner turned a quarter period:
! the normal component forward and the tangential ones back. The stiffness,
! its five terms and the physical mass of the sine-parity trials, returned
! alongside the cosine-parity ones, must equal an independent assembly with
! sine-parity trials. Checked on real kernel fields of the Solov'ev tokamak
! and of the two-period quasi-axisymmetric stellarator fixture, with modes
! whose toroidal numbers fall in the same and the opposite period class,
! including (0, 0), whose sine-parity normal component vanishes.
program test_sine_parity_pairing
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface
    use compatible_physical_mass_assembly, only: &
        assemble_compatible_physical_mass_surface
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use phase_assembly_policy, only: phase_assembly_transformed
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
    implicit none

    integer, parameter :: trials = 6, h1_count = 3, l2_count = 2
    integer, parameter :: mode_m(trials) = [0, 0, 1, 1, 2, 3]
    integer, parameter :: mode_n(trials) = [0, 1, -1, 1, 2, -3]
    character(len=1024) :: directory

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop "supply the fixture directory"
    call check_fixture(trim(directory)//"/solovev_q1.035.nc", 0.37_dp)
    call check_fixture(trim(directory)//"/qa_lowres.nc", 0.61_dp)
    write (*, "(a)") "PASS"

contains

    subroutine check_fixture(path, coordinate)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: coordinate
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        type(primitive_equilibrium_spline_t) :: spline
        real(dp), allocatable :: fields(:, :, :), drive(:, :)
        real(dp), allocatable :: jacobian_s(:, :), jacobian_t(:, :)
        real(dp), allocatable :: jacobian_z(:, :), gamma_p(:, :)
        real(dp) :: h1(h1_count, trials), dh1(h1_count, trials)
        real(dp) :: eta(l2_count, trials), l2(l2_count, trials)
        real(dp) :: theta(24), zeta(10), pressure
        real(dp), allocatable :: k_cos(:, :), t_cos(:, :, :), m_cos(:, :)
        real(dp), allocatable :: k_sin(:, :), t_sin(:, :, :), m_sin(:, :)
        real(dp), allocatable :: k_ref(:, :), t_ref(:, :, :), m_ref(:, :)
        integer :: columns, i, info, sine(trials), cosine(trials)

        call read_gvec_cas3d_file(path, equilibrium, info)
        call require(info == reader_ok, "fixture read failed")
        call fit_primitive_equilibrium(equilibrium, spline, info)
        call require(info == primitive_equilibrium_ok, "spline fit failed")
        do i = 1, size(theta)
            theta(i) = real(i - 1, dp) / real(size(theta), dp)
        end do
        do i = 1, size(zeta)
            zeta(i) = real(i - 1, dp) / real(size(zeta), dp)
        end do
        call evaluate_primitive_kernel_surface(spline, coordinate, theta, &
            zeta, fields, drive, info, jacobian_s, jacobian_t, jacobian_z, &
            pressure)
        call require(info == primitive_kernel_ok, "kernel surface failed")
        allocate (gamma_p(size(theta), size(zeta)), &
            source=5.0_dp / 3.0_dp * pressure)
        call fill_basis(h1, dh1, eta, l2)
        cosine = 1
        sine = 2
        columns = trials * (h1_count + 2 * l2_count)
        allocate (k_cos(columns, columns), t_cos(columns, columns, 5), &
            m_cos(columns, columns), k_sin(columns, columns), &
            t_sin(columns, columns, 5), m_sin(columns, columns), &
            k_ref(columns, columns), t_ref(columns, columns, 5), &
            m_ref(columns, columns))
        k_cos = 0.0_dp
        t_cos = 0.0_dp
        m_cos = 0.0_dp
        k_sin = 0.0_dp
        t_sin = 0.0_dp
        m_sin = 0.0_dp
        k_ref = 0.0_dp
        t_ref = 0.0_dp
        m_ref = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, mode_n, &
            cosine, spline%field_periods, h1, dh1, eta, l2, 0.3_dp, &
            phase_assembly_transformed, k_cos, info, t_cos, k_sin, t_sin)
        call require(info == 0, "paired stiffness failed")
        call assemble_compatible_compressible_stiffness_surface(fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, mode_n, &
            sine, spline%field_periods, h1, dh1, eta, l2, 0.3_dp, &
            phase_assembly_transformed, k_ref, info, t_ref)
        call require(info == 0, "sine stiffness failed")
        call require(close(k_sin, k_ref), "sine stiffness differs")
        do i = 1, 5
            call require(close(t_sin(:, :, i), t_ref(:, :, i)), &
                "sine stiffness term differs")
        end do
        call assemble_compatible_physical_mass_surface(fields, 1.5_dp, &
            mode_m, mode_n, cosine, spline%field_periods, h1, eta, l2, &
            0.3_dp, phase_assembly_transformed, m_cos, info, m_sin)
        call require(info == 0, "paired mass failed")
        call assemble_compatible_physical_mass_surface(fields, 1.5_dp, &
            mode_m, mode_n, sine, spline%field_periods, h1, eta, l2, &
            0.3_dp, phase_assembly_transformed, m_ref, info)
        call require(info == 0, "sine mass failed")
        call require(close(m_sin, m_ref), "sine mass differs")
        ! Pairing is defined for cosine-parity trials only.
        call assemble_compatible_physical_mass_surface(fields, 1.5_dp, &
            mode_m, mode_n, sine, spline%field_periods, h1, eta, l2, &
            0.3_dp, phase_assembly_transformed, m_ref, info, m_sin)
        call require(info /= 0, "pairing of sine-parity trials was accepted")
    end subroutine check_fixture

    ! Distinct radial basis values per function and trial.
    subroutine fill_basis(h1, dh1, eta, l2)
        real(dp), intent(out) :: h1(:, :), dh1(:, :), eta(:, :), l2(:, :)
        integer :: basis, trial

        do trial = 1, size(h1, 2)
            do basis = 1, size(h1, 1)
                h1(basis, trial) = 0.3_dp + 0.1_dp * basis + 0.07_dp * trial
                dh1(basis, trial) = -0.8_dp + 0.5_dp * basis - 0.03_dp * trial
            end do
            do basis = 1, size(l2, 1)
                eta(basis, trial) = 0.9_dp - 0.2_dp * basis + 0.05_dp * trial
                l2(basis, trial) = 0.4_dp + 0.3_dp * basis - 0.02_dp * trial
            end do
        end do
    end subroutine fill_basis

    logical function close(actual, expected)
        real(dp), intent(in) :: actual(:, :), expected(:, :)

        close = all(abs(actual - expected) <= 1.0e-12_dp &
            * max(maxval(abs(expected)), tiny(1.0_dp)))
    end function close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") "FAIL: "//message
            error stop 1
        end if
    end subroutine require

end program test_sine_parity_pairing
