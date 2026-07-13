module gliss_c_contexts
    use fixed_boundary_spectrum, only: fixed_boundary_problem_t
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none
    private

    type, public :: equilibrium_context_t
        type(gvec_cas3d_equilibrium_t) :: equilibrium
    end type equilibrium_context_t

    type, public :: stability_problem_context_t
        type(fixed_boundary_problem_t) :: problem
    end type stability_problem_context_t

end module gliss_c_contexts
