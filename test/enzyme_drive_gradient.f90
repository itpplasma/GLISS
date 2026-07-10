program enzyme_drive_gradient
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use, intrinsic :: iso_fortran_env, only: error_unit
    use local_mode_model, only: benchmark_mode_energy
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

    real(c_double) :: gradient

    gradient = enzyme_autodiff(c_funloc(benchmark_mode_energy), 0.1_c_double)
    if (abs(gradient + 0.5_c_double) > 1.0e-12_c_double) then
        write (error_unit, "(a,es24.16)") "FAIL drive gradient ", gradient
        error stop 1
    end if
    write (*, "(a)") "PASS"
end program enzyme_drive_gradient
