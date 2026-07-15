program test_generalized_block_solver
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use block_tridiagonal, only: apply_block_tridiagonal_into, &
        block_tridiagonal_t
    use family_assembly, only: iterate_block_eigenvalue
    use generalized_block_solver, only: generalized_inertia, &
        generalized_eigenpair_diagnostics, generalized_invalid, &
        generalized_mass_not_spd, generalized_ok, &
        iterate_generalized_eigenvalue
    use symmetric_eigensolver, only: solve_symmetric_generalized
    implicit none

    type(block_tridiagonal_t) :: stiffness, mass, identity_mass, corrupt
    real(dp), allocatable :: dense_k(:, :), dense_m(:, :), eigenvalues(:)
    real(dp), allocatable :: eigenvectors(:, :), vector(:, :), bad_vector(:, :)
    real(dp) :: eigenvalue, quotient, residual, shift, standard
    integer :: info, count, i

    call build_fixture(stiffness, mass)
    call block_to_dense(stiffness, dense_k)
    call block_to_dense(mass, dense_m)
    call solve_symmetric_generalized(dense_k, dense_m, eigenvalues, &
        eigenvectors, info)
    call require(info == 0, "dense generalized oracle failed")
    call require(all(eigenvalues(2:) > eigenvalues(:size(eigenvalues) - 1)), &
        "generalized fixture spectrum is not simple")
    call check_dense_eigenpairs(dense_k, dense_m, eigenvalues, eigenvectors)

    call generalized_inertia(stiffness, mass, eigenvalues(1) - 1.0_dp, &
        count, info)
    call require(info == generalized_ok .and. count == 0, &
        "generalized inertia below the spectrum is wrong")
    do i = 1, size(eigenvalues) - 1
        shift = 0.5_dp * (eigenvalues(i) + eigenvalues(i + 1))
        call generalized_inertia(stiffness, mass, shift, count, info)
        call require(info == generalized_ok .and. count == i, &
            "generalized inertia between eigenvalues is wrong")
    end do
    call generalized_inertia(stiffness, mass, &
        eigenvalues(size(eigenvalues)) + 1.0_dp, count, info)
    call require(info == generalized_ok .and. count == size(eigenvalues), &
        "generalized inertia above the spectrum is wrong")

    shift = eigenvalues(1) - 0.25_dp &
        * (eigenvalues(2) - eigenvalues(1))
    call iterate_generalized_eigenvalue(stiffness, mass, shift, eigenvalue, &
        vector, residual, info)
    call require(info == generalized_ok, "generalized iteration failed")
    call require(abs(eigenvalue - eigenvalues(1)) < 1.0e-11_dp, &
        "generalized iteration disagrees with the dense oracle")
    call require(abs(mass_norm(vector, mass) - 1.0_dp) < 2.0e-14_dp, &
        "generalized eigenvector is not mass-normalized")
    call require(residual < 1.0e-11_dp, &
        "generalized eigenvector residual is too large")
    call generalized_eigenpair_diagnostics(stiffness, mass, vector, &
        eigenvalue, quotient, standard, info)
    call require(info == generalized_ok, &
        "generalized eigenpair diagnostics failed")
    call require(abs(quotient - eigenvalue) < 1.0e-14_dp, &
        "generalized Rayleigh quotient is inconsistent")
    call require(abs(standard - residual) < 1.0e-14_dp, &
        "generalized residual is inconsistent")

    allocate (bad_vector(1, 1), source=1.0_dp)
    call generalized_eigenpair_diagnostics(stiffness, mass, bad_vector, &
        eigenvalue, quotient, residual, info)
    call require(info == generalized_invalid, &
        "wrong-shaped diagnostic vector was accepted")
    deallocate (bad_vector)
    allocate (bad_vector(size(vector, 1), size(vector, 2)))
    bad_vector = ieee_value(0.0_dp, ieee_quiet_nan)
    call generalized_eigenpair_diagnostics(stiffness, mass, bad_vector, &
        eigenvalue, quotient, residual, info)
    call require(info == generalized_invalid, &
        "nonfinite diagnostic vector was accepted")
    vector = 0.0_dp
    call generalized_eigenpair_diagnostics(stiffness, mass, vector, &
        eigenvalue, quotient, residual, info)
    call require(info == generalized_invalid, &
        "zero diagnostic vector was accepted")

    call build_identity_mass(stiffness, identity_mass)
    call iterate_block_eigenvalue(stiffness, 1.0_dp, shift, standard, info)
    call require(info == 0, "standard block iteration failed")
    call iterate_generalized_eigenvalue(stiffness, identity_mass, shift, &
        eigenvalue, vector, residual, info)
    call require(info == generalized_ok, "identity-mass iteration failed")
    call require(abs(eigenvalue - standard) < 1.0e-11_dp, &
        "identity mass changes the standard block eigenvalue")

    corrupt = stiffness
    corrupt%diag(1, 2, 1) = corrupt%diag(1, 2, 1) + 0.1_dp
    call generalized_inertia(corrupt, mass, 0.0_dp, count, info)
    call require(info == generalized_invalid, &
        "nonsymmetric stiffness was accepted")
    corrupt = mass
    corrupt%diag(1, 1, 1) = -1.0_dp
    call generalized_inertia(stiffness, corrupt, 0.0_dp, count, info)
    call require(info == generalized_mass_not_spd, &
        "indefinite mass was accepted")
    write (*, "(a)") "PASS"

contains

    subroutine build_fixture(stiffness, mass)
        type(block_tridiagonal_t), intent(out) :: stiffness, mass
        integer :: block

        allocate (stiffness%diag(2, 2, 3), stiffness%off(2, 2, 2))
        allocate (mass%diag(2, 2, 3), mass%off(2, 2, 2))
        do block = 1, 3
            stiffness%diag(1, 1, block) = 2.0_dp + 0.3_dp * block
            stiffness%diag(2, 1, block) = 0.2_dp
            stiffness%diag(1, 2, block) = 0.2_dp
            stiffness%diag(2, 2, block) = 3.0_dp + 0.4_dp * block
            mass%diag(1, 1, block) = 1.5_dp + 0.1_dp * block
            mass%diag(2, 1, block) = 0.1_dp
            mass%diag(1, 2, block) = 0.1_dp
            mass%diag(2, 2, block) = 1.2_dp + 0.05_dp * block
        end do
        stiffness%off(:, :, 1) = reshape([&
            -0.4_dp, 0.02_dp, 0.05_dp, -0.3_dp], [2, 2])
        stiffness%off(:, :, 2) = reshape([&
            -0.35_dp, 0.03_dp, 0.04_dp, -0.25_dp], [2, 2])
        mass%off(:, :, 1) = reshape([&
            0.08_dp, 0.0_dp, 0.01_dp, 0.05_dp], [2, 2])
        mass%off(:, :, 2) = reshape([&
            0.06_dp, 0.01_dp, 0.0_dp, 0.04_dp], [2, 2])
    end subroutine build_fixture

    subroutine check_dense_eigenpairs(stiffness, mass, eigenvalues, vectors)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenvalues(:)
        real(dp), intent(in) :: vectors(:, :)
        real(dp) :: mass_image(size(eigenvalues))
        real(dp) :: residual_vector(size(eigenvalues))
        real(dp) :: stiffness_image(size(eigenvalues)), product, residual
        integer :: i, j

        do i = 1, size(eigenvalues)
            call matrix_vector_product(stiffness, vectors(:, i), &
                stiffness_image)
            call matrix_vector_product(mass, vectors(:, i), mass_image)
            residual_vector = stiffness_image - eigenvalues(i) * mass_image
            residual = norm2(residual_vector)
            call require(residual < 1.0e-12_dp, &
                "dense generalized residual is too large")
            do j = 1, size(eigenvalues)
                call matrix_vector_product(mass, vectors(:, j), mass_image)
                product = dot_product(vectors(:, i), mass_image)
                if (i == j) product = product - 1.0_dp
                call require(abs(product) < 1.0e-12_dp, &
                    "dense eigenvectors are not mass-orthonormal")
            end do
        end do
    end subroutine check_dense_eigenpairs

    pure subroutine matrix_vector_product(matrix, vector, image)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp), intent(out) :: image(:)
        integer :: column, row

        do row = 1, size(matrix, 1)
            image(row) = 0.0_dp
            do column = 1, size(matrix, 2)
                image(row) = image(row) + matrix(row, column) * vector(column)
            end do
        end do
    end subroutine matrix_vector_product

    subroutine block_to_dense(blocks, dense)
        type(block_tridiagonal_t), intent(in) :: blocks
        real(dp), allocatable, intent(out) :: dense(:, :)
        integer :: block, width, first, next

        width = size(blocks%diag, 1)
        allocate (dense(width * size(blocks%diag, 3), &
            width * size(blocks%diag, 3)), source=0.0_dp)
        do block = 1, size(blocks%diag, 3)
            first = (block - 1) * width + 1
            dense(first:first + width - 1, first:first + width - 1) = &
                blocks%diag(:, :, block)
            if (block == size(blocks%diag, 3)) cycle
            next = first + width
            dense(next:next + width - 1, first:first + width - 1) = &
                blocks%off(:, :, block)
            dense(first:first + width - 1, next:next + width - 1) = &
                transpose(blocks%off(:, :, block))
        end do
    end subroutine block_to_dense

    subroutine build_identity_mass(stiffness, mass)
        type(block_tridiagonal_t), intent(in) :: stiffness
        type(block_tridiagonal_t), intent(out) :: mass
        integer :: block, entry

        allocate (mass%diag, mold=stiffness%diag)
        allocate (mass%off, mold=stiffness%off)
        mass%diag = 0.0_dp
        mass%off = 0.0_dp
        do block = 1, size(mass%diag, 3)
            do entry = 1, size(mass%diag, 1)
                mass%diag(entry, entry, block) = 1.0_dp
            end do
        end do
    end subroutine build_identity_mass

    function mass_norm(vector, mass) result(squared_norm)
        real(dp), contiguous, intent(in) :: vector(:, :)
        type(block_tridiagonal_t), intent(in) :: mass
        real(dp) :: squared_norm
        real(dp) :: image(size(vector, 1), size(vector, 2))

        call apply_block_tridiagonal_into(mass, vector, image)
        squared_norm = sum(vector * image)
    end function mass_norm

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_generalized_block_solver
