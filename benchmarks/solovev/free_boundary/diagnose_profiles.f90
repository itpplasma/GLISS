! Inspect the edge quantities actually used by the native vacuum operator.
program benchmark_boundary_profiles
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use plasma_vacuum_boundary, only: build_plasma_boundary_surface, &
        plasma_vacuum_ok
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use radial_cubic_spline, only: evaluate_radial_cubic_spline, &
        evaluate_radial_cubic_spline_field, fit_radial_cubic_spline, &
        radial_cubic_spline_ok, radial_cubic_spline_t
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(primitive_equilibrium_spline_t) :: spline
    type(radial_cubic_spline_t) :: iota
    character(len=1024) :: path
    real(dp), allocatable :: surface(:, :, :), normal(:, :, :)
    real(dp) :: fp, ft, iota_edge, discard
    real(dp) :: profiles(3), slopes(3), seconds(3)
    integer :: info

    if (command_argument_count() /= 1) &
        error stop 'usage: benchmark_boundary_profiles EXPORT'
    call get_command_argument(1, path)
    call read_gvec_cas3d_file(trim(path), equilibrium, info)
    if (info /= reader_ok) error stop 'read failed'
    call fit_primitive_equilibrium(equilibrium, spline, info)
    if (info /= primitive_equilibrium_ok) error stop 'fit failed'
    call build_plasma_boundary_surface(spline, 128, 128, surface, normal, &
        fp, ft, info)
    if (info /= plasma_vacuum_ok) error stop 'boundary failed'
    call fit_radial_cubic_spline(spline%radial_grid, &
        equilibrium%rotational_transform, iota, info)
    if (info /= radial_cubic_spline_ok) error stop 'iota fit failed'
    call evaluate_radial_cubic_spline(spline%radial_grid, iota, 1.0_dp, &
        iota_edge, discard, info)
    if (info /= radial_cubic_spline_ok) error stop 'iota evaluation failed'
    call evaluate_radial_cubic_spline_field(spline%radial_grid, &
        spline%profiles, 1.0_dp, profiles, slopes, seconds, info)
    if (info /= radial_cubic_spline_ok) error stop 'profile evaluation failed'
    write (*, '(a,es26.18)') 'native PhiPrime=', -ft
    write (*, '(a,es26.18)') 'native chiPrime/PhiPrime=', fp / ft
    write (*, '(a,es26.18)') 'stored iota edge spline=', iota_edge
    write (*, '(a,es26.18)') 'ratio minus stored iota=', fp / ft - iota_edge
    write (*, '(a,es26.18)') 'native edge pressure Pa=', profiles(3)
end program benchmark_boundary_profiles
