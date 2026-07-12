program gliss_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int
    use compressible_geometry, only: build_compressible_geometry, &
        compressible_geometry_ok
    use compressible_stiffness_family_assembly, only: &
        assemble_compressible_family_stiffness
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        dynamic_family_layout_t, dynamic_layout_ok
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mass_density_policy, only: mass_density_ok, &
        mass_density_profile_t, validate_mass_density_profile
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    use mode_topology, only: build_mode_family, mode_family_t
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_mass_family_assembly, only: assemble_physical_family_mass
    use radial_space_policy, only: radial_space_config_t
    use variable_block_tridiagonal, only: pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_generalized_ok, &
        variable_generalized_inertia
    use variable_spectrum_analysis, only: analyze_variable_spectrum, &
        variable_spectrum_ok, variable_spectrum_summary_t
    implicit none

    integer, parameter :: n_theta = 64, n_zeta = 64
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(mode_family_t) :: family
    character(len=1024) :: filename, token
    real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
    real(dp), allocatable :: jacobian_s(:, :, :), jacobian_t(:, :, :)
    real(dp), allocatable :: jacobian_z(:, :, :), gamma_p(:, :, :)
    integer, allocatable :: mode_m(:), mode_n(:)
    real(dp) :: adiabatic_index, density_kg_m3, zero_floor
    integer :: info, i, j, arguments, comma
    integer :: family_index, poloidal_max, toroidal_max
    logical :: generated_family

    interface
        subroutine terminate_process(status) bind(C, name="exit")
            import c_int
            integer(c_int), value :: status
        end subroutine terminate_process
    end interface

    arguments = command_argument_count()
    if (arguments < 5) then
        call fail_usage("missing required arguments")
    end if
    call read_argument(1, "EXPORT_FILE", filename)
    call read_real_argument(2, "GAMMA", adiabatic_index)
    call read_real_argument(3, "DENSITY", density_kg_m3)
    call read_real_argument(4, "FLOOR", zero_floor)
    if (adiabatic_index <= 0.0_dp) call fail_usage("GAMMA must be positive")
    if (density_kg_m3 <= 0.0_dp) call fail_usage("DENSITY must be positive")
    if (zero_floor <= 0.0_dp) call fail_usage("FLOOR must be positive")
    call read_argument(5, "mode or --family", token)
    generated_family = trim(token) == "--family"
    if (generated_family) then
        if (arguments /= 8) then
            call fail_usage("--family requires exactly INDEX MMAX NMAX")
        end if
        call read_integer_argument(6, "INDEX", family_index)
        call read_integer_argument(7, "MMAX", poloidal_max)
        call read_integer_argument(8, "NMAX", toroidal_max)
        if (family_index < 0) call fail_usage("INDEX must be nonnegative")
        if (poloidal_max < 0) call fail_usage("MMAX must be nonnegative")
        if (toroidal_max < 0) call fail_usage("NMAX must be nonnegative")
    else
        allocate (mode_m(arguments - 4), mode_n(arguments - 4))
        do i = 5, arguments
            call read_argument(i, "mode", token)
            comma = index(token, ",")
            if (comma <= 1 .or. comma == len_trim(token)) &
                call fail_usage("modes must be given as m,n")
            if (index(token(comma + 1:), ",") > 0) &
                call fail_usage("modes must contain exactly one comma")
            call parse_integer(token(:comma - 1), "poloidal mode", &
                mode_m(i - 4))
            call parse_integer(token(comma + 1:), "toroidal mode", &
                mode_n(i - 4))
            if (mode_m(i - 4) < 0) &
                call fail_usage("poloidal mode m must be nonnegative")
            if (mode_m(i - 4) == 0 .and. mode_n(i - 4) < 0) &
                call fail_usage("axis modes require nonnegative n")
            do j = 1, i - 5
                if (mode_m(j) == mode_m(i - 4) .and. &
                    mode_n(j) == mode_n(i - 4)) &
                    call fail_usage("duplicate mode")
            end do
        end do
    end if

    call read_gvec_cas3d_file(trim(filename), equilibrium, info)
    if (info /= reader_ok) then
        write (error_unit, "(a, i0)") "reader error ", info
        error stop 1
    end if
    if (generated_family) then
        call build_mode_family(equilibrium%field_periods, family_index, &
            poloidal_max, toroidal_max, family, info)
        if (info /= 0) then
            write (error_unit, "(a)") "invalid mode-family configuration"
            error stop 1
        end if
        mode_m = family%poloidal
        mode_n = family%toroidal
    end if
    call build_kernel_geometry(equilibrium, n_theta, n_zeta, fields, &
        drive, info)
    if (info /= mercier_ok) then
        write (error_unit, "(a, i0)") "geometry error ", info
        error stop 1
    end if
    call build_compressible_geometry(equilibrium, n_theta, n_zeta, &
        adiabatic_index, jacobian_s, jacobian_t, jacobian_z, gamma_p, info)
    if (info /= compressible_geometry_ok) then
        write (error_unit, "(a, i0)") "compressible geometry error ", info
        error stop 1
    end if

    write (*, "(a)") "chart_metric,field_periods,modes,parity_class," // &
        "adiabatic_index,density_kg_m3,unknowns,negative_count," // &
        "floor_count,lowest_eigenvalue,certificate,eigenpair_residual," // &
        "eigenpair_resolution,inertia_interval"
    do i = 1, 2
        call report_class(i)
    end do

contains

    subroutine fail_usage(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") "gliss_spectrum: " // trim(message)
        write (error_unit, "(a)") &
            "usage: gliss_spectrum EXPORT_FILE GAMMA DENSITY FLOOR " // &
            "m,n [m,n ...] | --family INDEX MMAX NMAX"
        call terminate_process(2_c_int)
    end subroutine fail_usage

    subroutine read_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        integer :: status

        call get_command_argument(position, value, status=status)
        if (status /= 0) call fail_usage(trim(name) // " is too long")
        if (len_trim(value) == 0) call fail_usage(trim(name) // " is empty")
    end subroutine read_argument

    subroutine read_real_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: value

        call read_argument(position, name, token)
        call parse_real(token, name, value)
    end subroutine read_real_argument

    subroutine read_integer_argument(position, name, value)
        integer, intent(in) :: position
        character(len=*), intent(in) :: name
        integer, intent(out) :: value

        call read_argument(position, name, token)
        call parse_integer(token, name, value)
    end subroutine read_integer_argument

    subroutine parse_real(text, name, value)
        character(len=*), intent(in) :: text, name
        real(dp), intent(out) :: value
        integer :: status

        read (text, *, iostat=status) value
        if (status /= 0) call fail_usage(trim(name) // " must be a number")
        if (.not. ieee_is_finite(value)) &
            call fail_usage(trim(name) // " must be finite")
    end subroutine parse_real

    subroutine parse_integer(text, name, value)
        character(len=*), intent(in) :: text, name
        integer, intent(out) :: value
        integer :: status

        read (text, *, iostat=status) value
        if (status /= 0) call fail_usage(trim(name) // " must be an integer")
    end subroutine parse_integer

    subroutine report_class(parity_class)
        integer, intent(in) :: parity_class
        type(dynamic_family_layout_t) :: layout, mass_layout
        type(mass_density_profile_t) :: density
        type(radial_space_config_t) :: radial_space
        type(variable_block_tridiagonal_t) :: k_blocks, m_blocks
        type(variable_spectrum_summary_t) :: summary
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        integer, allocatable :: trial_parity(:), permutation(:), widths(:)
        real(dp), allocatable :: stored_power(:)
        real(dp) :: step, lowest, certificate, residual, resolution
        real(dp) :: inertia_interval
        integer :: ns, status

        ns = size(equilibrium%s)
        step = 1.0_dp / real(ns, dp)
        allocate (trial_parity(size(mode_m)), source=parity_class)
        allocate (stored_power(size(mode_m)), source=0.0_dp)
        density%s = [0.0_dp, 1.0_dp]
        density%kilograms_per_cubic_metre = [density_kg_m3, density_kg_m3]
        call validate_mass_density_profile(density, status)
        if (status /= mass_density_ok) then
            write (error_unit, "(a)") "invalid density profile"
            error stop 1
        end if
        call assemble_compressible_family_stiffness(fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p, mode_m, mode_n, &
            trial_parity, stored_power, equilibrium%field_periods, &
            radial_space, step, phase_assembly_transformed, stiffness, &
            layout, status)
        if (status /= 0) then
            write (error_unit, "(a, i0)") "stiffness assembly error ", status
            error stop 1
        end if
        call assemble_physical_family_mass(fields, density, mode_m, mode_n, &
            trial_parity, stored_power, equilibrium%field_periods, &
            radial_space, step, phase_assembly_transformed, mass, &
            mass_layout, status)
        if (status /= 0) then
            write (error_unit, "(a, i0)") "mass assembly error ", status
            error stop 1
        end if
        if (mass_layout%total_unknowns /= layout%total_unknowns) then
            write (error_unit, "(a)") "stiffness and mass layouts differ"
            error stop 1
        end if
        call build_dynamic_block_permutation(layout, widths, permutation, &
            status)
        if (status /= dynamic_layout_ok) then
            write (error_unit, "(a, i0)") "permutation error ", status
            error stop 1
        end if
        call pack_permuted_variable_blocks(stiffness, permutation, widths, &
            k_blocks, status)
        if (status /= variable_block_ok) then
            write (error_unit, "(a, i0)") "stiffness packing error ", status
            error stop 1
        end if
        call pack_permuted_variable_blocks(mass, permutation, widths, &
            m_blocks, status)
        if (status /= variable_block_ok) then
            write (error_unit, "(a, i0)") "mass packing error ", status
            error stop 1
        end if
        call analyze_variable_spectrum(k_blocks, m_blocks, zero_floor, &
            summary, status)
        if (status /= variable_spectrum_ok) then
            write (error_unit, "(a, i0)") "spectrum analysis error ", status
            error stop 1
        end if
        call resolve_lowest(k_blocks, m_blocks, summary, lowest, &
            inertia_interval, residual, resolution)
        certificate = inertia_interval + residual + resolution
        write (*, "(l1, a, i0, a, i0, a, i0, a, es9.2, a, es9.2, a, i0, " // &
            "a, i0, a, i0, a, es24.16, 4(a, es9.2))") &
            equilibrium%has_chart_metric, ",", equilibrium%field_periods, &
            ",", size(mode_m), ",", parity_class, ",", adiabatic_index, &
            ",", density_kg_m3, ",", layout%total_unknowns, ",", &
            summary%negative_count, ",", summary%zero_count, ",", lowest, &
            ",", certificate, ",", residual, ",", resolution, ",", &
            inertia_interval
    end subroutine report_class

    subroutine resolve_lowest(k_blocks, m_blocks, summary, lowest, &
            inertia_interval, residual, resolution)
        type(variable_block_tridiagonal_t), intent(in) :: k_blocks, m_blocks
        type(variable_spectrum_summary_t), intent(in) :: summary
        real(dp), intent(out) :: lowest, inertia_interval, residual, resolution
        real(dp), allocatable :: vector(:)
        real(dp) :: lower, upper, middle, shift
        integer :: status, count, iteration

        if (summary%negative_count == 0) then
            if (.not. summary%has_positive) then
                lowest = 0.0_dp
                inertia_interval = summary%zero_floor
                residual = 0.0_dp
                resolution = 0.0_dp
                return
            end if
            shift = 0.5_dp * (summary%first_positive_lower &
                + summary%first_positive_upper)
            inertia_interval = summary%first_positive_upper &
                - summary%first_positive_lower
        else
            lower = -2.0_dp * summary%zero_floor
            do iteration = 1, 200
                call variable_generalized_inertia(k_blocks, m_blocks, &
                    lower, count, status)
                if (status /= variable_generalized_ok) then
                    lower = lower * (1.0_dp + 1.0e-8_dp)
                    cycle
                end if
                if (count == 0) exit
                lower = 2.0_dp * lower
            end do
            upper = -summary%zero_floor
            do iteration = 1, 200
                middle = 0.5_dp * (lower + upper)
                if (upper - lower <= 1.0e-9_dp * abs(middle) &
                    + 1.0e-3_dp * summary%zero_floor) exit
                call variable_generalized_inertia(k_blocks, m_blocks, &
                    middle, count, status)
                if (status /= variable_generalized_ok) then
                    middle = middle * (1.0_dp + 1.0e-8_dp)
                    cycle
                end if
                if (count == 0) then
                    lower = middle
                else
                    upper = middle
                end if
            end do
            shift = lower
            inertia_interval = upper - lower
        end if
        call iterate_variable_generalized_eigenvalue(k_blocks, m_blocks, &
            shift, lowest, vector, residual, resolution, status)
        if (status /= variable_generalized_ok) then
            write (error_unit, "(a, i0)") "inverse iteration error ", status
            error stop 1
        end if
    end subroutine resolve_lowest

end program gliss_spectrum
