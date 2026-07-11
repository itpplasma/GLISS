program enzyme_helical_cylinder_gradient
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use, intrinsic :: iso_fortran_env, only: error_unit
    use helical_cylinder_limit, only: benchmark_helical_vertical_margin
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

    gradient = enzyme_autodiff(c_funloc(benchmark_helical_vertical_margin), &
        0.5_c_double)
    centered = (benchmark_helical_vertical_margin(0.5_c_double + step) &
        - benchmark_helical_vertical_margin(0.5_c_double - step)) &
        / (2.0_c_double * step)
    if (abs(gradient - 5.0_c_double) > 1.0e-12_c_double) then
        write (error_unit, "(a,es24.16)") &
            "FAIL helical-cylinder gradient ", gradient
        error stop 1
    end if
    if (abs(gradient - centered) > 1.0e-9_c_double) then
        write (error_unit, "(a,2es24.16)") &
            "FAIL helical-cylinder finite-difference gradient ", &
            gradient, centered
        error stop 1
    end if
    write (*, "(a)") "PASS"
end program enzyme_helical_cylinder_gradient
