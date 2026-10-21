program test_ideal_mhd_batches
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    use physical_constants, only: vacuum_permeability
    use ideal_mhd_batch, only: mhd_xi_fp64, mhd_xi_fp32, mhd_xi_mixed
    use two_component_kernel, only: two_component_energy_batch, &
        two_component_energy_jvp, two_component_energy_vjp
    use three_component_kernel, only: three_component_energy_batch, &
        three_component_energy_jvp, three_component_energy_vjp
    use perpendicular_kinetic_kernel, only: perpendicular_kinetic_energy_batch, &
        perpendicular_kinetic_energy_jvp, perpendicular_kinetic_energy_vjp
    implicit none
    integer, parameter :: dp = real64, sp = real32
    integer :: n = 257, repeats = 0, family = 0
    character(32) :: arg
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg); read (arg, *) family
        call get_command_argument(2, arg); read (arg, *) n
        call get_command_argument(3, arg); read (arg, *) repeats
    end if
    if (n < 1 .or. repeats < 0) error stop 'invalid batch dimensions'
    call check_dense_scalar()
    if (family == 0 .or. family == 1) call check_two()
    if (family == 0 .or. family == 2) call check_three()
    if (family == 0 .or. family == 3) call check_mass()
    if (family == 0 .or. family == 4) call check_xi()
    print *, "PASS shared MHD batch analytical/derivative/kind oracles"
contains
    subroutine check_dense_scalar()
        ! Nondyadic geometry/background inputs exercise every derivative entry.
        integer, parameter :: count=17
        real(dp) :: a(count,55), y(count,29), gradient(27), value, tangent
        real(dp) :: inputs27(27), direction27(27)
        real(sp) :: single(count,55), low(count,29)
        integer :: j,k
        a=0
        do j=1,count
            do k=1,27
                a(j,k)=0.1_dp*sin(real(7*j+3*k,dp))
                a(j,27+k)=cos(real(2*j+k,dp))/100
            end do
            a(j,1)=1.1_dp; a(j,2)=0.23_dp; a(j,7)=1.2_dp
            if(mod(j,2)==0) a(j,7)=-1.2_dp
            a(j,8)=1.7_dp; a(j,9)=1.3_dp; a(j,15)=0.3_dp
            a(j,27)=vacuum_permeability; a(j,54)=0; a(j,55)=1
        end do
        call three_component_energy_batch(count,a,y)
        do j=1,count
            inputs27=a(j,1:27); direction27=a(j,28:54)
            call three_component_energy_jvp(inputs27(:26),direction27(:26), &
                value,tangent)
            call three_component_energy_vjp(inputs27(:26),1.0_dp,value, &
                gradient(:26))
            if(abs(y(j,1)-value)>2e-12_dp*(1+abs(value))) error stop 'dense three value'
            if(abs(y(j,2)-tangent)>2e-12_dp*(1+abs(tangent))) error stop 'dense three JVP'
            if(maxval(abs(y(j,3:28)-gradient(:26))/(1+abs(gradient(:26))))>2e-12_dp) &
                error stop 'dense three VJP'
        end do
        single=real(a,sp)
        call three_component_energy_batch(count,single,low)
        if(maxval(abs(real(low,dp)-y)/(1+abs(y)))>512*epsilon(1.0_sp)) error stop 'dense three fp32'
        ! Compress the two-component physical input order (gamma is omitted).
        do j=1,count
            inputs27=a(j,1:27); direction27=a(j,28:54)
            a(j,1:14)=inputs27(:14)
            a(j,15:20)=inputs27(16:21)
            a(j,21:34)=direction27(:14)
            a(j,35:40)=direction27(16:21)
            a(j,41)=1
        end do
        call two_component_energy_batch(count,a(:,:41),y(:,:22))
        do j=1,count
            inputs27(:20)=a(j,:20); direction27(:20)=a(j,21:40)
            call two_component_energy_jvp(inputs27(:20),direction27(:20), &
                value,tangent)
            call two_component_energy_vjp(inputs27(:20),1.0_dp,value,gradient(:20))
            if(abs(y(j,1)-value)>2e-12_dp*(1+abs(value))) error stop 'dense two value'
            if(abs(y(j,2)-tangent)>2e-12_dp*(1+abs(tangent))) error stop 'dense two JVP'
            if(maxval(abs(y(j,3:22)-gradient(:20))/(1+abs(gradient(:20))))>2e-12_dp) &
                error stop 'dense two VJP'
        end do
        do j=1,count
            a(j,:7)=[1.2_dp,1.7_dp,1.3_dp,0.1_dp,2.3_dp,0.4_dp,-0.3_dp]
            if(mod(j,2)==0) a(j,1)=-1.2_dp
            do k=1,7
                a(j,7+k)=cos(real(2*j+k,dp))/100
            end do
            a(j,15)=1
        end do
        call perpendicular_kinetic_energy_batch(count,a(:,:15),y(:,:9))
        do j=1,count
            inputs27(:7)=a(j,:7); direction27(:7)=a(j,8:14)
            call perpendicular_kinetic_energy_jvp(inputs27(:7),direction27(:7), &
                value,tangent)
            call perpendicular_kinetic_energy_vjp(inputs27(:7),1.0_dp,value, &
                gradient(:7))
            if(abs(y(j,1)-value)>2e-12_dp*(1+abs(value))) error stop 'dense mass value'
            if(abs(y(j,2)-tangent)>2e-12_dp*(1+abs(tangent))) error stop 'dense mass JVP'
            if(maxval(abs(y(j,3:9)-gradient(:7))/(1+abs(gradient(:7))))>2e-12_dp) &
                error stop 'dense mass VJP'
        end do
    end subroutine check_dense_scalar

    subroutine check_two()
        real(dp), allocatable :: x(:, :), y(:, :), z(:, :), yp(:, :), ym(:, :), xm(:, :)
        real(sp), allocatable :: xs(:, :), ys(:, :)
        real(dp) :: u, v, du, dv, expected, forward, h, defect, pairing
        real(dp) :: elapsed, started, checksum, max_error
        integer :: j, k, sample, mode
        allocate (x(n, 41), y(n, 22), z(n, 22), &
                  yp(n, 22), ym(n, 22), xm(n, 41), &
                  xs(n, 41), ys(n, 22))
        x = 0
        x(:, 41) = 1
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            x(j, 1) = 1; x(j, 7) = 1; x(j, 8) = 1; x(j, 9) = 1
            if (mod(j, 2) == 0) x(j, 7) = -1
            x(j, 18) = u; x(j, 20) = v
            x(j, 38) = du; x(j, 40) = dv
        end do
        call two_component_energy_batch(n, x, y)
        if (.not. all(ieee_is_finite(y))) error stop 'two nonfinite fp64'
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            expected = u*u + v*v; forward = 2*(u*du + v*dv)
            if (abs(y(j, 1) - expected) > 2e-12_dp*(1 + abs(expected))) &
                error stop 'two analytical value'
            if (abs(y(j, 2) - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'two analytical JVP'
            pairing = sum(y(j, 3:)*x(j, 21:40))
            if (abs(pairing - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'two VJP pairing'
        end do
        ! Finite differences activate every background and displacement input.
        do k = 1, 20
            x(:, 20 + k) = real((-1)**k, dp)/64
        end do
        call two_component_energy_batch(n, x, y)
        h = 2.0_dp**(-18)
        xm = x
        xm(:, :20) = x(:, :20) + h*x(:, 21:40)
        call two_component_energy_batch(n, xm, yp)
        xm(:, :20) = x(:, :20) - h*x(:, 21:40)
        call two_component_energy_batch(n, xm, ym)
        defect = maxval(abs((yp(:, 1) - ym(:, 1))/(2*h) - y(:, 2)))
        if (defect > 2e-8_dp*(1 + maxval(abs(y(:, 2))))) &
            error stop 'two centered directional difference'
        do j = 1, n
            pairing = sum(y(j, 3:)*x(j, 21:40))
            if (abs(pairing - y(j, 2)) > 2e-12_dp*(1 + abs(y(j, 2)))) &
                error stop 'two full VJP pairing'
        end do
        xs = real(x, sp)
        call two_component_energy_batch(n, xs, ys)
        call two_component_energy_batch(n, xs, z)
        if (.not. all(ieee_is_finite(ys)) .or. .not. all(ieee_is_finite(z))) &
            error stop 'two nonfinite low precision'
        max_error = max(maxval(abs(real(ys, dp) - y)), maxval(abs(z - y)))
        ! Well-conditioned dyadic samples, bounded away from all denominators.
        ! 512*u32 bounds the operation chains and input sensitivity here;
        ! this regression budget is empirical admission, not an enclosure.
        if (max_error > 512*epsilon(1.0_sp)*(1 + maxval(abs(y)))) &
            error stop 'two low precision rounding budget'
        if (repeats == 0) return
        do mode = 1, 3
            do sample = -1, repeats
                started = omp_get_wtime()
                select case (mode)
                case (1)
                    call two_component_energy_batch(n, x, y)
                case (2)
                    call two_component_energy_batch(n, xs, ys)
                case (3)
                    call two_component_energy_batch(n, xs, z)
                end select
                elapsed = omp_get_wtime() - started
                select case (mode)
                case (1)
                    checksum = sum(y(:, 1))
                case (2)
                    checksum = sum(real(ys(:, 1), dp))
                case (3)
                    checksum = sum(z(:, 1))
                end select
                if (sample > 0) print '(a,4(",",i0),4(",",es24.16))', &
                    'two', omp_get_max_threads(), n, mode, sample, &
                    elapsed, max_error, checksum, real(size(x)*8 + size(y)*8, dp)
            end do
        end do
    end subroutine check_two
    subroutine check_three()
        real(dp), allocatable :: x(:, :), y(:, :), z(:, :), yp(:, :), ym(:, :), xm(:, :)
        real(sp), allocatable :: xs(:, :), ys(:, :)
        real(dp) :: u, v, du, dv, expected, forward, h, defect, pairing
        real(dp) :: elapsed, started, checksum, max_error
        integer :: j, k, sample, mode
        allocate (x(n, 55), y(n, 29), z(n, 29), &
                  yp(n, 29), ym(n, 29), xm(n, 55), &
                  xs(n, 55), ys(n, 29))
        x = 0
        x(:, 55) = 1
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            x(j, 1) = 1; x(j, 7) = 1; x(j, 8) = 1; x(j, 9) = 1
            if (mod(j, 2) == 0) x(j, 7) = -1
            x(j, 19) = u; x(j, 21) = v; x(j, 27) = 2
            x(j, 15) = 3; x(j, 22) = 0.5_dp
            x(j, 46) = du; x(j, 48) = dv
        end do
        call three_component_energy_batch(n, x, y)
        if (.not. all(ieee_is_finite(y))) error stop 'three nonfinite fp64'
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            expected = (u*u + v*v)/2 + 0.75_dp; forward = u*du + v*dv
            if (abs(y(j, 1) - expected) > 2e-12_dp*(1 + abs(expected))) &
                error stop 'three analytical value'
            if (abs(y(j, 2) - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'three analytical JVP'
            pairing = sum(y(j, 3:)*x(j, 28:54))
            if (abs(pairing - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'three VJP pairing'
        end do
        ! Finite differences activate every background and displacement input.
        do k = 1, 27
            x(:, 27 + k) = real((-1)**k, dp)/64
        end do
        call three_component_energy_batch(n, x, y)
        h = 2.0_dp**(-18)
        xm = x
        xm(:, :27) = x(:, :27) + h*x(:, 28:54)
        call three_component_energy_batch(n, xm, yp)
        xm(:, :27) = x(:, :27) - h*x(:, 28:54)
        call three_component_energy_batch(n, xm, ym)
        defect = maxval(abs((yp(:, 1) - ym(:, 1))/(2*h) - y(:, 2)))
        if (defect > 2e-8_dp*(1 + maxval(abs(y(:, 2))))) &
            error stop 'three centered directional difference'
        do j = 1, n
            pairing = sum(y(j, 3:)*x(j, 28:54))
            if (abs(pairing - y(j, 2)) > 2e-12_dp*(1 + abs(y(j, 2)))) &
                error stop 'three full VJP pairing'
        end do
        xs = real(x, sp)
        call three_component_energy_batch(n, xs, ys)
        call three_component_energy_batch(n, xs, z)
        if (.not. all(ieee_is_finite(ys)) .or. .not. all(ieee_is_finite(z))) &
            error stop 'three nonfinite low precision'
        max_error = max(maxval(abs(real(ys, dp) - y)), maxval(abs(z - y)))
        ! Well-conditioned dyadic samples, bounded away from all denominators.
        ! 512*u32 bounds the operation chains and input sensitivity here;
        ! this regression budget is empirical admission, not an enclosure.
        if (max_error > 512*epsilon(1.0_sp)*(1 + maxval(abs(y)))) &
            error stop 'three low precision rounding budget'
        if (repeats == 0) return
        do mode = 1, 3
            do sample = -1, repeats
                started = omp_get_wtime()
                select case (mode)
                case (1)
                    call three_component_energy_batch(n, x, y)
                case (2)
                    call three_component_energy_batch(n, xs, ys)
                case (3)
                    call three_component_energy_batch(n, xs, z)
                end select
                elapsed = omp_get_wtime() - started
                select case (mode)
                case (1)
                    checksum = sum(y(:, 1))
                case (2)
                    checksum = sum(real(ys(:, 1), dp))
                case (3)
                    checksum = sum(z(:, 1))
                end select
                if (sample > 0) print '(a,4(",",i0),4(",",es24.16))', &
                    'three', omp_get_max_threads(), n, mode, sample, &
                    elapsed, max_error, checksum, real(size(x)*8 + size(y)*8, dp)
            end do
        end do
    end subroutine check_three
    subroutine check_mass()
        real(dp), allocatable :: x(:, :), y(:, :), z(:, :), yp(:, :), ym(:, :), xm(:, :)
        real(sp), allocatable :: xs(:, :), ys(:, :)
        real(dp) :: u, v, du, dv, expected, forward, h, defect, pairing
        real(dp) :: elapsed, started, checksum, max_error
        integer :: j, k, sample, mode
        allocate (x(n, 15), y(n, 9), z(n, 9), &
                  yp(n, 9), ym(n, 9), xm(n, 15), &
                  xs(n, 15), ys(n, 9))
        x = 0
        x(:, 15) = 1
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            x(j, 1) = 1; x(j, 2) = 2; x(j, 3) = 4; x(j, 5) = 3
            if (mod(j, 2) == 0) x(j, 1) = -1
            x(j, 6) = u; x(j, 7) = v
            x(j, 13) = du; x(j, 14) = dv
        end do
        call perpendicular_kinetic_energy_batch(n, x, y)
        if (.not. all(ieee_is_finite(y))) error stop 'mass nonfinite fp64'
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            expected = 1.5_dp*(u*u/4 + v*v); forward = 3*(u*du/4 + v*dv)
            if (abs(y(j, 1) - expected) > 2e-12_dp*(1 + abs(expected))) &
                error stop 'mass analytical value'
            if (abs(y(j, 2) - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'mass analytical JVP'
            pairing = sum(y(j, 3:)*x(j, 8:14))
            if (abs(pairing - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'mass VJP pairing'
        end do
        ! Finite differences activate every background and displacement input.
        do k = 1, 7
            x(:, 7 + k) = real((-1)**k, dp)/64
        end do
        call perpendicular_kinetic_energy_batch(n, x, y)
        h = 2.0_dp**(-18)
        xm = x
        xm(:, :7) = x(:, :7) + h*x(:, 8:14)
        call perpendicular_kinetic_energy_batch(n, xm, yp)
        xm(:, :7) = x(:, :7) - h*x(:, 8:14)
        call perpendicular_kinetic_energy_batch(n, xm, ym)
        defect = maxval(abs((yp(:, 1) - ym(:, 1))/(2*h) - y(:, 2)))
        if (defect > 2e-8_dp*(1 + maxval(abs(y(:, 2))))) &
            error stop 'mass centered directional difference'
        do j = 1, n
            pairing = sum(y(j, 3:)*x(j, 8:14))
            if (abs(pairing - y(j, 2)) > 2e-12_dp*(1 + abs(y(j, 2)))) &
                error stop 'mass full VJP pairing'
        end do
        xs = real(x, sp)
        call perpendicular_kinetic_energy_batch(n, xs, ys)
        call perpendicular_kinetic_energy_batch(n, xs, z)
        if (.not. all(ieee_is_finite(ys)) .or. .not. all(ieee_is_finite(z))) &
            error stop 'mass nonfinite low precision'
        max_error = max(maxval(abs(real(ys, dp) - y)), maxval(abs(z - y)))
        ! Well-conditioned dyadic samples, bounded away from all denominators.
        ! 512*u32 bounds the operation chains and input sensitivity here;
        ! this regression budget is empirical admission, not an enclosure.
        if (max_error > 512*epsilon(1.0_sp)*(1 + maxval(abs(y)))) &
            error stop 'mass low precision rounding budget'
        if (repeats == 0) return
        do mode = 1, 3
            do sample = -1, repeats
                started = omp_get_wtime()
                select case (mode)
                case (1)
                    call perpendicular_kinetic_energy_batch(n, x, y)
                case (2)
                    call perpendicular_kinetic_energy_batch(n, xs, ys)
                case (3)
                    call perpendicular_kinetic_energy_batch(n, xs, z)
                end select
                elapsed = omp_get_wtime() - started
                select case (mode)
                case (1)
                    checksum = sum(y(:, 1))
                case (2)
                    checksum = sum(real(ys(:, 1), dp))
                case (3)
                    checksum = sum(z(:, 1))
                end select
                if (sample > 0) print '(a,4(",",i0),4(",",es24.16))', &
                    'mass', omp_get_max_threads(), n, mode, sample, &
                    elapsed, max_error, checksum, real(size(x)*8 + size(y)*8, dp)
            end do
        end do
    end subroutine check_mass
    subroutine check_xi()
        real(dp), allocatable :: x(:, :), y(:, :), z(:, :), yp(:, :), ym(:, :), xm(:, :)
        real(sp), allocatable :: xs(:, :), ys(:, :)
        real(dp) :: u, v, du, dv, expected, forward, h, defect, pairing
        real(dp) :: elapsed, started, checksum, max_error
        integer :: j, k, sample, mode
        allocate (x(n, 47), y(n, 25), z(n, 25), &
                  yp(n, 25), ym(n, 25), xm(n, 47), &
                  xs(n, 47), ys(n, 25))
        x = 0
        x(:, 47) = 1
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            x(j, 3) = 1; x(j, 5) = 1
            x(j, 12) = u; x(j, 21) = v
            x(j, 35) = du; x(j, 44) = dv
        end do
        call mhd_xi_fp64(n, x, y)
        if (.not. all(ieee_is_finite(y))) error stop 'xi nonfinite fp64'
        do j = 1, n
            u = 0.25_dp + real(mod(j, 97), dp)/128
            v = -0.125_dp + real(mod(j, 53), dp)/64
            du = 0.125_dp; dv = -0.25_dp
            expected = u*v; forward = du*v + u*dv
            if (abs(y(j, 1) - expected) > 2e-12_dp*(1 + abs(expected))) &
                error stop 'xi analytical value'
            if (abs(y(j, 2) - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'xi analytical JVP'
            pairing = sum(y(j, 3:)*x(j, 24:46))
            if (abs(pairing - forward) > 2e-12_dp*(1 + abs(forward))) &
                error stop 'xi VJP pairing'
        end do
        ! Finite differences activate every background and displacement input.
        do k = 1, 23
            x(:, 23 + k) = real((-1)**k, dp)/64
        end do
        call mhd_xi_fp64(n, x, y)
        h = 2.0_dp**(-18)
        xm = x
        xm(:, :23) = x(:, :23) + h*x(:, 24:46)
        call mhd_xi_fp64(n, xm, yp)
        xm(:, :23) = x(:, :23) - h*x(:, 24:46)
        call mhd_xi_fp64(n, xm, ym)
        defect = maxval(abs((yp(:, 1) - ym(:, 1))/(2*h) - y(:, 2)))
        if (defect > 2e-8_dp*(1 + maxval(abs(y(:, 2))))) &
            error stop 'xi centered directional difference'
        do j = 1, n
            pairing = sum(y(j, 3:)*x(j, 24:46))
            if (abs(pairing - y(j, 2)) > 2e-12_dp*(1 + abs(y(j, 2)))) &
                error stop 'xi full VJP pairing'
        end do
        xs = real(x, sp)
        call mhd_xi_fp32(n, xs, ys)
        call mhd_xi_mixed(n, xs, z)
        if (.not. all(ieee_is_finite(ys)) .or. .not. all(ieee_is_finite(z))) &
            error stop 'xi nonfinite low precision'
        max_error = max(maxval(abs(real(ys, dp) - y)), maxval(abs(z - y)))
        ! Well-conditioned dyadic samples, bounded away from all denominators.
        ! 512*u32 bounds the operation chains and input sensitivity here;
        ! this regression budget is empirical admission, not an enclosure.
        if (max_error > 512*epsilon(1.0_sp)*(1 + maxval(abs(y)))) &
            error stop 'xi low precision rounding budget'
        if (repeats == 0) return
        do mode = 1, 3
            do sample = -1, repeats
                started = omp_get_wtime()
                select case (mode)
                case (1)
                    call mhd_xi_fp64(n, x, y)
                case (2)
                    call mhd_xi_fp32(n, xs, ys)
                case (3)
                    call mhd_xi_mixed(n, xs, z)
                end select
                elapsed = omp_get_wtime() - started
                select case (mode)
                case (1)
                    checksum = sum(y(:, 1))
                case (2)
                    checksum = sum(real(ys(:, 1), dp))
                case (3)
                    checksum = sum(z(:, 1))
                end select
                if (sample > 0) print '(a,4(",",i0),4(",",es24.16))', &
                    'xi', omp_get_max_threads(), n, mode, sample, &
                    elapsed, max_error, checksum, real(size(x)*8 + size(y)*8, dp)
            end do
        end do
    end subroutine check_xi
end program test_ideal_mhd_batches
