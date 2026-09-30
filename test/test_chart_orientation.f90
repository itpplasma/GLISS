! Chart orientation is a property of the reconstructed map.  Reflecting the
! positions (z -> -z) reverses the handedness of (s, theta, zeta) but leaves
! the physics unchanged, so the n=1 stability count of the Solov'ev fixture
! must be the same for both orientations.
program test_chart_orientation
    use axisymmetric_spectrum, only: axisymmetric_spectrum_ok, &
        axisymmetric_spectrum_result_t, compute_axisymmetric_spectrum
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use primitive_kernel_geometry, only: primitive_chart_orientation, &
        primitive_kernel_ok
    implicit none
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(axisymmetric_spectrum_result_t) :: original, mirrored
    character(len=1024) :: directory
    character(len=256) :: message
    integer :: orientation, info

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop 'supply the fixture directory'
    call read_gvec_cas3d_file(trim(directory) // '/solovev_q1.035.nc', &
        equilibrium, info)
    if (info /= reader_ok) error stop 'fixture read failed'
    call primitive_chart_orientation(equilibrium, orientation, info)
    if (info /= primitive_kernel_ok .or. orientation /= -1) &
        error stop 'GVEC export is not reported left-handed'
    call compute_axisymmetric_spectrum(equilibrium, 1, 6, 2, .false., &
        original, info, message)
    if (info /= axisymmetric_spectrum_ok) error stop 'original spectrum failed'
    equilibrium%zhat%cosine = -equilibrium%zhat%cosine
    equilibrium%zhat%sine = -equilibrium%zhat%sine
    call primitive_chart_orientation(equilibrium, orientation, info)
    if (info /= primitive_kernel_ok .or. orientation /= 1) &
        error stop 'mirrored export is not reported right-handed'
    call compute_axisymmetric_spectrum(equilibrium, 1, 6, 2, .false., &
        mirrored, info, message)
    if (info /= axisymmetric_spectrum_ok) error stop 'mirrored spectrum failed'
    if (mirrored%negative_count /= original%negative_count) &
        error stop 'stability count depends on chart orientation'
    write (*, '(a,i0)') 'chart orientation checks passed; count ', &
        original%negative_count
end program test_chart_orientation
