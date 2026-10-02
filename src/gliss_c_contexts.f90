module gliss_c_contexts
    use fixed_boundary_spectrum, only: fixed_boundary_problem_t
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    implicit none
    private

    type, public :: equilibrium_context_t
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        integer :: references = 1
    end type equilibrium_context_t

    type, public :: stability_problem_context_t
        type(fixed_boundary_problem_t) :: problem
        type(equilibrium_context_t), pointer :: equilibrium => null()
    end type stability_problem_context_t

    public :: retain_equilibrium_context, release_equilibrium_context

contains

    subroutine retain_equilibrium_context(context, info)
        type(equilibrium_context_t), pointer, intent(inout) :: context
        integer, intent(out) :: info

        info = 1
        if (.not. associated(context)) return
        if (context%references < 1) return
        if (context%references == huge(context%references)) return
        context%references = context%references + 1
        info = 0
    end subroutine retain_equilibrium_context

    subroutine release_equilibrium_context(context, info)
        type(equilibrium_context_t), pointer, intent(inout) :: context
        integer, intent(out) :: info

        info = 0
        if (.not. associated(context)) return
        if (context%references < 1) then
            info = 1
            return
        end if
        if (context%references == 1) then
            deallocate (context, stat=info)
            if (info /= 0) return
        else
            context%references = context%references - 1
        end if
        nullify (context)
    end subroutine release_equilibrium_context

end module gliss_c_contexts
