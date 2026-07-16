program test_cas3d_phase_envelope_transform
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use cas3d_phase_envelope_transform, only: &
        apply_cas3d_phase_envelope_congruence, &
        build_cas3d_phase_envelope_map, cas3d_phase_envelope_map_t, &
        cas3d_phase_transform_invalid, cas3d_phase_transform_ok
    use fourier_phase_kind, only: phase_cosine
    implicit none

    type(cas3d_phase_envelope_map_t) :: collision_map, map
    real(dp), allocatable :: mass(:, :), stiffness(:, :)
    real(dp), allocatable :: stiffness_terms(:, :, :)
    real(dp), parameter :: diagonal_values(3) = [2.0_dp, 4.0_dp, 6.0_dp]
    real(dp) :: expected_stiffness(3, 3), first_null(5), image(5)
    real(dp) :: second_null(5)
    integer :: info

    call build_cas3d_phase_envelope_map([3, 4, 2], [2, 2, 2], &
        [1, 1, 1], [3, 4, 2], [2, 2, 2], phase_cosine, map, info)
    call require(info == cas3d_phase_transform_ok, &
        "noncolliding carrier-plus-pair map failed")
    call require(all(map%physical_row(:, 1) == [1, 0]), &
        "carrier physical row differs")
    call require(all(map%physical_row(:, 2) == [2, 3]), &
        "even-envelope physical rows differ")
    call require(all(map%normal_weight(:, 2) == [0.5_dp, 0.5_dp]), &
        "normal even-envelope weights differ")
    call require(all(map%normal_weight(:, 3) == [-0.5_dp, 0.5_dp]), &
        "normal odd-envelope weights differ")
    call require(all(map%eta_weight(:, 2) == [0.5_dp, 0.5_dp]), &
        "eta even-envelope weights differ")
    call require(all(map%eta_weight(:, 3) == [0.5_dp, -0.5_dp]), &
        "eta odd-envelope weights differ")
    call build_cas3d_phase_envelope_map([3, 4, 2], [2, 2, 2], &
        [1, -1, 1], [3, 4, 2], [2, 2, 2], phase_cosine, map, info)
    call require(info == cas3d_phase_transform_ok, &
        "canonical-orientation map failed")
    call require(all(map%normal_weight(:, 2) == [0.5_dp, 0.5_dp]), &
        "cosine canonicalization changed a normal coefficient")
    call require(all(map%eta_weight(:, 2) == [-0.5_dp, 0.5_dp]), &
        "sine canonicalization did not reverse the eta coefficient")
    call build_cas3d_phase_envelope_map([0, 1, 1], [0, 0, 0], &
        [1, 1, -1], [0, 1], [0, 0], phase_cosine, map, info)
    call require(info == cas3d_phase_transform_ok, &
        "self-colliding sideband pair map failed")
    call require(all(map%physical_row(:, 2) == [2, 0]) &
        .and. all(map%physical_row(:, 3) == [2, 0]), &
        "self-colliding support was not reduced")
    call require(all(map%normal_weight(:, 2) == [1.0_dp, 0.0_dp]) &
        .and. all(map%normal_weight(:, 3) == [0.0_dp, 0.0_dp]), &
        "self-colliding normal weights differ")
    call require(all(map%eta_weight(:, 2) == [0.0_dp, 0.0_dp]) &
        .and. all(map%eta_weight(:, 3) == [1.0_dp, 0.0_dp]), &
        "self-colliding eta weights differ")
    call build_cas3d_phase_envelope_map([3, 4, 2], [2, 2, 2], &
        [1, 1, 1], [3, 4, 2, 5], [2, 2, 2, 2], phase_cosine, map, info)
    call require(info == cas3d_phase_transform_invalid, &
        "map with an unspanned physical mode was accepted")

    call build_cas3d_phase_envelope_map([3, 4, 2], [2, 2, 2], &
        [1, -1, 1], [3, 4, 2], [2, 2, 2], phase_cosine, map, info)
    call require(info == cas3d_phase_transform_ok, &
        "valid map rebuild failed after rejected input")

    call diagonal_fixture(diagonal_values, stiffness, &
        stiffness_terms, mass)
    call apply_cas3d_phase_envelope_congruence(map, 1, 0, 0.25_dp, &
        stiffness, stiffness_terms, mass, info)
    call require(info == cas3d_phase_transform_ok, &
        "noncolliding coefficient congruence failed")
    expected_stiffness = reshape([2.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 2.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 2.5_dp], [3, 3])
    call require(all(stiffness == expected_stiffness), &
        "carrier-plus-pair stiffness congruence differs")
    call require(all(stiffness_terms(:, :, 1) == stiffness), &
        "termwise congruence differs from total stiffness")
    call require_diagonal_mass(mass, 0.25_dp, &
        "printed envelope coefficient mass differs")

    call build_cas3d_phase_envelope_map([3, 4, 2, 2, 4], &
        [2, 2, 2, 2, 2], [1, 1, 1, 1, 1], [3, 4, 2], [2, 2, 2], &
        phase_cosine, collision_map, info)
    call require(info == cas3d_phase_transform_ok, &
        "colliding labeled map failed")
    call diagonal_fixture(diagonal_values, stiffness, &
        stiffness_terms, mass)
    call apply_cas3d_phase_envelope_congruence(collision_map, 1, 0, &
        0.25_dp, stiffness, stiffness_terms, mass, info)
    call require(info == cas3d_phase_transform_ok, &
        "colliding labeled congruence failed")
    first_null = [0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp]
    second_null = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
    image = matmul(stiffness, first_null)
    call require(maxval(abs(image)) == 0.0_dp, &
        "first repeated-label null direction is not exact")
    image = matmul(stiffness, second_null)
    call require(maxval(abs(image)) == 0.0_dp, &
        "second repeated-label null direction is not exact")
    call require_diagonal_mass(mass, 0.25_dp, &
        "colliding labels lost positive coefficient mass")
    write (*, "(a)") "PASS"

contains

    subroutine diagonal_fixture(diagonal, matrix, terms, local_mass)
        real(dp), intent(in) :: diagonal(:)
        real(dp), allocatable, intent(out) :: matrix(:, :), terms(:, :, :)
        real(dp), allocatable, intent(out) :: local_mass(:, :)
        integer :: index

        allocate (matrix(size(diagonal), size(diagonal)), source=0.0_dp)
        allocate (terms(size(diagonal), size(diagonal), 1), source=0.0_dp)
        allocate (local_mass(size(diagonal), size(diagonal)), source=0.0_dp)
        do index = 1, size(diagonal)
            matrix(index, index) = diagonal(index)
            terms(index, index, 1) = diagonal(index)
            local_mass(index, index) = 1.0_dp
        end do
    end subroutine diagonal_fixture

    subroutine require_diagonal_mass(matrix, diagonal, message)
        real(dp), intent(in) :: matrix(:, :), diagonal
        character(len=*), intent(in) :: message
        integer :: column, row

        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                if (row == column) then
                    call require(matrix(row, column) == diagonal, message)
                else
                    call require(matrix(row, column) == 0.0_dp, message)
                end if
            end do
        end do
    end subroutine require_diagonal_mass

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_cas3d_phase_envelope_transform
