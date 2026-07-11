program enzyme_compressible_stiffness_gradient
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use, intrinsic :: iso_fortran_env, only: error_unit
    use enzyme_compressible_stiffness_fixture, only: &
        transformed_compressible_stiffness_energy, &
        transformed_radial_coordinate_energy, transformed_radial_step_energy, &
        transformed_stored_power_energy
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

    gradient = enzyme_autodiff( &
        c_funloc(transformed_compressible_stiffness_energy), point)
    centered = centered_geometry()
    call require_gradient(gradient, centered, "geometry")
    gradient = enzyme_autodiff(c_funloc(transformed_radial_coordinate_energy), &
        point)
    centered = centered_radial_coordinate()
    call require_gradient(gradient, centered, "radial coordinate")
    gradient = enzyme_autodiff(c_funloc(transformed_radial_step_energy), point)
    centered = centered_radial_step()
    call require_gradient(gradient, centered, "radial step")
    gradient = enzyme_autodiff(c_funloc(transformed_stored_power_energy), point)
    centered = centered_stored_power()
    call require_gradient(gradient, centered, "stored power")
    write (*, "(a)") "PASS"

contains

    subroutine require_gradient(gradient, centered, label)
        real(c_double), intent(in) :: gradient, centered
        character(len=*), intent(in) :: label

        if (abs(centered) <= 1.0e-8_c_double) then
            write (error_unit, "(a)") &
                "FAIL inactive compressible stiffness " // label // " gradient"
            error stop 1
        end if
        if (abs(gradient - centered) > 1.0e-8_c_double &
            * max(1.0_c_double, abs(centered))) then
            write (error_unit, "(a,2es24.16)") &
                "FAIL compressible stiffness " // label // " gradient ", &
                gradient, centered
            error stop 1
        end if
    end subroutine require_gradient

    function centered_geometry() result(centered)
        real(c_double) :: centered

        centered = (transformed_compressible_stiffness_energy(point + step) &
            - transformed_compressible_stiffness_energy(point - step)) &
            / (2.0_c_double * step)
    end function centered_geometry

    function centered_radial_coordinate() result(centered)
        real(c_double) :: centered

        centered = (transformed_radial_coordinate_energy(point + step) &
            - transformed_radial_coordinate_energy(point - step)) &
            / (2.0_c_double * step)
    end function centered_radial_coordinate

    function centered_radial_step() result(centered)
        real(c_double) :: centered

        centered = (transformed_radial_step_energy(point + step) &
            - transformed_radial_step_energy(point - step)) &
            / (2.0_c_double * step)
    end function centered_radial_step

    function centered_stored_power() result(centered)
        real(c_double) :: centered

        centered = (transformed_stored_power_energy(point + step) &
            - transformed_stored_power_energy(point - step)) &
            / (2.0_c_double * step)
    end function centered_stored_power

end program enzyme_compressible_stiffness_gradient
