program enzyme_compressible_stiffness_family_gradient
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use, intrinsic :: iso_fortran_env, only: error_unit
    use enzyme_compressible_stiffness_family_fixture, only: &
        global_compressible_stiffness_energy
    implicit none

    interface
        function enzyme_autodiff(function_pointer, x) result(gradient) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: function_pointer
            real(c_double), value :: x
            real(c_double) :: gradient
        end function enzyme_autodiff
    end interface

    real(c_double), parameter :: point = 0.1_c_double
    real(c_double), parameter :: step = 1.0e-5_c_double
    real(c_double) :: centered, gradient

    gradient = enzyme_autodiff(c_funloc(global_compressible_stiffness_energy), &
        point)
    centered = (global_compressible_stiffness_energy(point + step) &
        - global_compressible_stiffness_energy(point - step)) &
        / (2.0_c_double * step)
    if (abs(centered) <= 1.0e-8_c_double) then
        write (error_unit, "(a)") &
            "FAIL inactive global compressible stiffness gradient"
        error stop 1
    end if
    if (abs(gradient - centered) > 1.0e-8_c_double &
        * max(1.0_c_double, abs(centered))) then
        write (error_unit, "(a,2es24.16)") &
            "FAIL global compressible stiffness gradient ", gradient, centered
        error stop 1
    end if
    write (*, "(a)") "PASS"
end program enzyme_compressible_stiffness_family_gradient
