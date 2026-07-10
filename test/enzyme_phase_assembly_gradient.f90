program enzyme_phase_assembly_gradient
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use, intrinsic :: iso_fortran_env, only: error_unit
    use enzyme_phase_fixture, only: transformed_phase_energy
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

    real(c_double), parameter :: step = 1.0e-6_c_double
    real(c_double) :: centered, gradient

    gradient = enzyme_autodiff(c_funloc(transformed_phase_energy), 0.0_c_double)
    centered = (transformed_phase_energy(step) &
        - transformed_phase_energy(-step)) / (2.0_c_double * step)
    if (abs(gradient - centered) > 1.0e-7_c_double &
        * max(1.0_c_double, abs(centered))) then
        write (error_unit, "(a,2es24.16)") &
            "FAIL transformed phase gradient ", gradient, centered
        error stop 1
    end if
    write (*, "(a)") "PASS"
end program enzyme_phase_assembly_gradient
