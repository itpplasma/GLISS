module physical_constants
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter, public :: vacuum_permeability = &
        1.2566370614359172954e-6_dp

end module physical_constants
