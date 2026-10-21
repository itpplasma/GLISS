module gliss_sparse_backend
    use, intrinsic :: iso_fortran_env, only: dp => real64
#ifdef GLISS_HAVE_FORTSPARSE
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_inertia, fortsparse_status_t, i8, FORTSPARSE_OK, &
        FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_INDEFINITE, &
        FORTSPARSE_INERTIA_RELIABLE
    use fortnum_precision, only: precision_policy_t, precision_pair_t, FP64, &
        STAGE_EVALUATE, STAGE_CERTIFY
#endif
    implicit none
    private
#ifdef GLISS_HAVE_FORTSPARSE
    logical, parameter, public :: gliss_sparse_enabled = .true.
#else
    logical, parameter, public :: gliss_sparse_enabled = .false.
#endif

    type, public :: gliss_sparse_factor_t
        private
        logical :: ready = .false.
#ifdef GLISS_HAVE_FORTSPARSE
        type(sparse_solver_t) :: solver
        ! The initial consumer supports cpu-fp64 only. In particular, inertia,
        ! estimators and accepted residuals never use narrowed storage/math.
        type(precision_policy_t) :: precision
#endif
    end type gliss_sparse_factor_t

    public :: factor_sparse_csc, solve_sparse_factored
#ifdef GLISS_HAVE_FORTSPARSE
    public :: sparse_precision_supported
#endif
contains

#ifdef GLISS_HAVE_FORTSPARSE
    pure logical function sparse_precision_supported(policy) result(supported)
        type(precision_policy_t), intent(in) :: policy
        type(precision_pair_t) :: pair
        integer :: stage

        supported = .false.
        do stage = STAGE_EVALUATE, STAGE_CERTIFY
            pair = policy%get(stage)
            if (pair%storage /= FP64 .or. pair%compute /= FP64) return
        end do
        supported = .true.
    end function sparse_precision_supported
#endif

    subroutine factor_sparse_csc(n, pointers, rows, values, factor, &
            negative, reliable, info)
        integer, intent(in) :: n, pointers(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(gliss_sparse_factor_t), intent(out) :: factor
        integer, intent(out) :: negative, info
        logical, intent(out) :: reliable
#ifdef GLISS_HAVE_FORTSPARSE
        type(fortsparse_status_t) :: status
        integer(i8) :: nnegative, nzero, npositive
        integer :: state
#endif

        info = -1
        negative = -1
        reliable = .false.
        if (n < 1 .or. size(pointers) /= n + 1) return
        if (size(rows) /= size(values)) return
#ifdef GLISS_HAVE_FORTSPARSE
        if (.not. sparse_precision_supported(factor%precision)) return
        factor%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        factor%solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        ! Raw CSC avoids a triplet sort and any dense n-by-n allocation.
        call sparse_factor(factor%solver, n, n, size(values), pointers, &
            rows, values, status)
        if (status%code /= FORTSPARSE_OK) return
        call sparse_inertia(factor%solver, nnegative, nzero, npositive, &
            state, status)
        if (status%code /= FORTSPARSE_OK) return
        negative = int(nnegative)
        reliable = state == FORTSPARSE_INERTIA_RELIABLE .and. nzero == 0_i8
        if (nnegative + nzero + npositive /= int(n, i8)) return
        factor%ready = .true.
        info = 0
#endif
    end subroutine factor_sparse_csc

    subroutine solve_sparse_factored(factor, rhs, info)
        type(gliss_sparse_factor_t), intent(inout) :: factor
        real(dp), intent(inout) :: rhs(:)
        integer, intent(out) :: info
#ifdef GLISS_HAVE_FORTSPARSE
        type(fortsparse_status_t) :: status
#endif

        info = -1
        if (.not. factor%ready .or. size(rhs) < 1) return
#ifdef GLISS_HAVE_FORTSPARSE
        call sparse_solve(factor%solver, rhs, status)
        if (status%code == FORTSPARSE_OK) info = 0
#endif
    end subroutine solve_sparse_factored
end module gliss_sparse_backend
