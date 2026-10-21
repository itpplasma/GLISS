subroutine mhd_cpu_reference(family, n, x, y) bind(C)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use ideal_mhd_batch
    implicit none
    integer(c_int), value, intent(in) :: family, n
    real(c_double), intent(in) :: x(n, *)
    real(c_double), intent(out) :: y(n, *)
    select case (family)
    case (1)
        call mhd_two_fp64(n, x(:, :41), y(:, :22))
    case (2)
        call mhd_three_fp64(n, x(:, :55), y(:, :29))
    case (3)
        call mhd_mass_fp64(n, x(:, :15), y(:, :9))
    case (4)
        call mhd_xi_fp64(n, x(:, :47), y(:, :25))
    end select
end subroutine mhd_cpu_reference
