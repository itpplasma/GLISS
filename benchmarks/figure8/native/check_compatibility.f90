program check_compatibility
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        evaluate_primitive_equilibrium, primitive_equilibrium_spline_t
    use primitive_geometry_grid, only: primitive_geometry_grid_t
    use plasma_vacuum_boundary, only: build_plasma_boundary_surface
    implicit none
    type(gvec_cas3d_equilibrium_t) :: eq
    type(primitive_equilibrium_spline_t) :: spline
    type(primitive_geometry_grid_t) :: geometry
    real(dp), allocatable :: xyz(:,:,:), normal(:,:,:)
    real(dp) :: theta(37), zeta(73), pressure, slope, fp, ft
    real(dp), parameter :: radial(9)=[.0001_dp,.001_dp,.01_dp,.05_dp, &
        .15_dp,.37_dp,.63_dp,.89_dp,.99_dp]
    integer, parameter :: nu(3)=[25,33,49], nv(3)=[64,96,128]
    character(len=2048) :: path
    integer :: info,i,k
    call get_command_argument(1,path)
    call read_gvec_cas3d_file(trim(path),eq,info)
    write(*,'(a,i0)') 'reader_info=',info
    if(info/=reader_ok) stop 1
    write(*,'(a,i0)') 'field_periods=',eq%field_periods
    write(*,'(a,i0)') 'boundary_schema=',eq%boundary%schema_version
    if(eq%field_periods/=2 .or. eq%boundary%schema_version/=1) stop 8
    call fit_primitive_equilibrium(eq,spline,info)
    write(*,'(a,i0)') 'spline_info=',info
    if(info/=0) stop 2
    do i=1,size(theta)
        theta(i)=(real(i-1,dp)+.371_dp)/real(size(theta),dp)
    enddo
    do k=1,size(zeta)
        zeta(k)=(real(k-1,dp)+.293_dp)/real(size(zeta),dp)
    enddo
    do i=1,size(radial)
        call evaluate_primitive_equilibrium(spline,radial(i),theta,zeta, &
            geometry,pressure,slope,info)
        if(info/=0) then
            write(*,'(a,es20.12,a,i0)') 'radial_s=',radial(i),' geometry_info=',info
            stop 3
        else
            write(*,'(a,es20.12,a,2es20.12,a,es20.12)') 'radial_s=',radial(i), &
                ' signed_J_minmax=',minval(geometry%signed_jacobian), &
                maxval(geometry%signed_jacobian), &
                ' min_B=',minval(geometry%mod_b)
            if(.not.all(ieee_is_finite(geometry%signed_jacobian))) stop 9
            if(.not.all(ieee_is_finite(geometry%mod_b))) stop 10
            if(.not.ieee_is_finite(pressure) .or. .not.ieee_is_finite(slope)) stop 11
            if(maxval(geometry%signed_jacobian)>=0.0_dp) stop 4
            if(minval(geometry%mod_b)<=0.0_dp) stop 5
            if(abs(pressure)>tiny(1.0_dp)) stop 6
        endif
    enddo
    do i=1,3
        call build_plasma_boundary_surface(spline,nu(i),nv(i),xyz,normal,fp,ft,info)
        write(*,'(a,2i5,a,i0,a,2es20.12)') 'boundary_nu_nv=',nu(i),nv(i), &
            ' info=',info,' fp_ft=',fp,ft
        if(info/=0) stop 7
        if(.not.all(ieee_is_finite(xyz)) .or. .not.all(ieee_is_finite(normal))) stop 12
        if(.not.ieee_is_finite(fp) .or. .not.ieee_is_finite(ft)) stop 13
        if(abs(ft-0.03141592654_dp)>1.e-12_dp .or. abs(fp+0.7_dp*ft)>1.e-12_dp) stop 14
    enddo
end program check_compatibility
