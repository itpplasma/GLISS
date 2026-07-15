module fixed_boundary_solver_controls
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    type, public :: fixed_boundary_solver_controls_t
        real(dp) :: eigenvalue_relative = 1.0e-13_dp
        real(dp) :: residual_relative = 1.0e-12_dp
        real(dp) :: negative_bracket_relative = 1.0e-9_dp
        real(dp) :: negative_bracket_floor = 1.0e-3_dp
        integer :: inverse_iteration_limit = 500
        integer :: bracket_iteration_limit = 200
    end type fixed_boundary_solver_controls_t

    public :: valid_fixed_boundary_solver_controls

contains

    pure function valid_fixed_boundary_solver_controls(controls) result(valid)
        type(fixed_boundary_solver_controls_t), intent(in) :: controls
        logical :: valid
        real(dp) :: tolerances(4)

        tolerances(1) = controls%eigenvalue_relative
        tolerances(2) = controls%residual_relative
        tolerances(3) = controls%negative_bracket_relative
        tolerances(4) = controls%negative_bracket_floor
        valid = all(ieee_is_finite(tolerances))
        if (.not. valid) return
        valid = all(tolerances > 0.0_dp)
        if (.not. valid) return
        valid = controls%inverse_iteration_limit >= 1 &
            .and. controls%bracket_iteration_limit >= 1
    end function valid_fixed_boundary_solver_controls

end module fixed_boundary_solver_controls
