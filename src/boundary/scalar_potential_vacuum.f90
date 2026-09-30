module scalar_potential_vacuum
    ! Energy of the current-free field outside a toroidal plasma edge, inside
    ! an optional ideal wall.
    !
    ! The field B = grad(Phi) of the vacuum region V solves the Neumann
    ! problem with prescribed normal flux on the plasma edge P and none on
    ! the wall W. Green's identity on the boundary,
    !     c(x) Phi(x) + int Phi dG/dn_V = int G dPhi/dn_V,   G = 1/(4 pi r),
    ! with n_V the outward normal of V, is collocated at the centroids of flat
    ! triangles carrying a constant Phi. The single layer is the exact
    ! potential of a uniform triangle and the double layer its exact signed
    ! solid angle. The jump c(x) is each surface's discrete solid-angle sum,
    ! 1/2 in the continuum. The energy is
    !     integral_V |B|^2 dV = int_{P} Phi dPhi/dn_V dA,
    ! so only the vacuum region contributes. A current-sheet representation
    ! of the same Neumann data also stores field energy inside the plasma.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boundary_mesh_geometry, only: meshes_intersect, point_inside_mesh
    implicit none
    private

    integer, parameter, public :: vacuum_bie_ok = 0
    integer, parameter, public :: vacuum_bie_invalid = 1
    integer, parameter, public :: vacuum_bie_degenerate = 2
    integer, parameter, public :: vacuum_bie_not_nested = 3
    integer, parameter, public :: vacuum_bie_singular = 4

    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Relative node mismatch below which a rotation is a mesh symmetry.
    real(dp), parameter :: symmetry_tolerance = 1.0e-10_dp

    public :: assemble_exterior_vacuum, edge_triangle_centroids

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv

        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine zgesv

        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta
            real(dp), intent(in) :: a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm

        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                beta, c, ldc)
            import :: dp
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(dp), intent(in) :: alpha, beta
            complex(dp), intent(in) :: a(lda, *), b(ldb, *)
            complex(dp), intent(inout) :: c(ldc, *)
        end subroutine zgemm
    end interface

contains

    ! Centroids (u, v) of the 2 nu nv edge triangles in assembly order; each
    ! triangle covers 1 / (2 nu nv) of the unit (u, v) square.
    pure subroutine edge_triangle_centroids(nu, nv, centroid)
        integer, intent(in) :: nu, nv
        real(dp), intent(out) :: centroid(:, :)
        real(dp) :: du, dv, u0, v0
        integer :: i, k, triangle

        du = 1.0_dp / real(nu, dp)
        dv = 1.0_dp / real(nv, dp)
        triangle = 0
        do k = 1, nv
            v0 = real(k - 1, dp) * dv
            do i = 1, nu
                u0 = real(i - 1, dp) * du
                triangle = triangle + 1
                centroid(1, triangle) = u0 + du / 3.0_dp
                centroid(2, triangle) = v0 + 2.0_dp * dv / 3.0_dp
                triangle = triangle + 1
                centroid(1, triangle) = u0 + 2.0_dp * du / 3.0_dp
                centroid(2, triangle) = v0 + dv / 3.0_dp
            end do
        end do
    end subroutine edge_triangle_centroids

    ! flux(t, k) is B.n dA / (du dv) of Neumann datum k at the centroid of
    ! edge triangle t, for either consistent orientation of n. energy(k, l)
    ! is the bilinear form of integral_V B_k . B_l dV (units mu0 = 1).
    !
    ! When the edge and wall meshes are invariant under a rotation by 2 pi / P
    ! about the z axis that shifts the toroidal node index by nv / P, the
    ! collocation matrices are block circulant over the P sectors: only the
    ! rows of the first sector are assembled, and a discrete Fourier
    ! transform over sectors leaves P independent systems of one sector's
    ! size. symmetric = .false. solves the full system instead; both give the
    ! same energy to rounding.
    subroutine assemble_exterior_vacuum(plasma, flux, energy, info, wall, &
            symmetric)
        real(dp), contiguous, intent(in) :: plasma(:, :, :)
        real(dp), intent(in) :: flux(:, :)
        real(dp), allocatable, intent(out) :: energy(:, :)
        integer, intent(out) :: info
        real(dp), contiguous, intent(in), optional :: wall(:, :, :)
        logical, intent(in), optional :: symmetric
        real(dp), allocatable :: area(:), centre(:, :), corners(:, :, :)
        real(dp), allocatable :: normal(:, :), wall_corners(:, :, :)
        real(dp), allocatable :: datum(:, :, :), double_layer(:, :, :)
        real(dp), allocatable :: single_layer(:, :, :), harmonic(:, :, :)
        integer, allocatable :: order(:, :), owner(:), status(:)
        real(dp) :: average, weight
        integer :: block, column, data, harmonics, j, plasma_block
        integer :: plasma_count, row, sectors, total, wall_block

        info = vacuum_bie_invalid
        if (.not. valid_surface(plasma)) return
        plasma_count = 2 * size(plasma, 2) * size(plasma, 3)
        if (size(flux, 1) /= plasma_count .or. size(flux, 2) < 1) return
        if (.not. all(ieee_is_finite(flux))) return
        if (present(wall)) then
            if (.not. valid_surface(wall)) return
        end if
        call build_surface(plasma, corners, info)
        if (info /= vacuum_bie_ok) return
        total = plasma_count
        if (present(wall)) then
            call build_surface(wall, wall_corners, info)
            if (info /= vacuum_bie_ok) return
            total = total + size(wall_corners, 3)
        end if
        sectors = 1
        if (present(symmetric)) then
            if (symmetric) sectors = rotational_sectors(plasma, wall)
        else
            sectors = rotational_sectors(plasma, wall)
        end if
        plasma_block = plasma_count / sectors
        wall_block = (total - plasma_count) / sectors
        block = plasma_block + wall_block
        if (present(wall)) then
            ! Rotation carries the first sector onto every other one.
            if (.not. nested(plasma, wall, corners, wall_corners, sectors)) &
                then
                info = vacuum_bie_not_nested
                return
            end if
        end if
        allocate (centre(3, total), normal(3, total), area(total), &
            owner(total), order(block, sectors))
        call triangle_frames(corners, 1, centre, normal, area, owner, 1)
        if (present(wall)) call triangle_frames(wall_corners, 2, centre, &
            normal, area, owner, plasma_count + 1)
        ! n_V points into the plasma on its edge and out through the wall.
        call orient(centre, normal, area, owner, 1, -1.0_dp)
        if (present(wall)) call orient(centre, normal, area, owner, 2, 1.0_dp)
        ! Unknown row of sector s: its plasma triangles, then its wall ones.
        do j = 1, sectors
            do row = 1, plasma_block
                order(row, j) = (j - 1) * plasma_block + row
            end do
            do row = 1, wall_block
                order(plasma_block + row, j) = plasma_count &
                    + (j - 1) * wall_block + row
            end do
        end do

        data = size(flux, 2)
        allocate (datum(plasma_block, sectors, data))
        average = 1.0_dp / real(plasma_count, dp)
        do column = 1, data
            do j = 1, sectors
                do row = 1, plasma_block
                    ! Mean dPhi/dn_V over the triangle; the sign convention
                    ! cancels in the quadratic energy.
                    datum(row, j, column) = flux(order(row, j), column) &
                        * average / area(order(row, j))
                end do
            end do
        end do
        allocate (double_layer(block, block, sectors), &
            single_layer(block, plasma_block, sectors))
        call assemble_sector_rows(corners, wall_corners, centre, normal, &
            order, plasma_count, plasma_block, double_layer, single_layer)

        harmonics = sectors / 2 + 1
        allocate (harmonic(data, data, harmonics), status(harmonics))
        !$omp parallel do schedule(dynamic) if (harmonics > 1)
        do j = 1, harmonics
            call solve_harmonic(j - 1, sectors, double_layer, single_layer, &
                datum, area(1:plasma_block), present(wall), total, &
                harmonic(:, :, j), status(j))
        end do
        !$omp end parallel do
        info = vacuum_bie_singular
        if (any(status /= 0)) return
        ! Parseval over sectors; harmonics j and P - j are conjugate.
        allocate (energy(data, data))
        energy = 0.0_dp
        do j = 1, harmonics
            weight = 2.0_dp
            if (j == 1 .or. 2 * (j - 1) == sectors) weight = 1.0_dp
            energy = energy + weight / real(sectors, dp) * harmonic(:, :, j)
        end do
        do column = 1, data
            do row = 1, column - 1
                average = 0.5_dp * (energy(row, column) + energy(column, row))
                energy(row, column) = average
                energy(column, row) = average
            end do
        end do
        if (.not. all(ieee_is_finite(energy))) return
        info = vacuum_bie_ok
    end subroutine assemble_exterior_vacuum

    ! Rows of the first sector against the columns of every sector:
    ! double_layer(:, :, d) and single_layer(:, :, d) couple the collocation
    ! points of sector 1 to the triangles of sector d. The self term of the
    ! double layer is the jump coefficient, each surface's discrete
    ! solid-angle sum: +1/2 on the plasma edge, -1/2 on the wall in the
    ! continuum, whose n_V is the outward normal of the enclosed body.
    subroutine assemble_sector_rows(corners, wall_corners, centre, normal, &
            order, plasma_count, plasma_block, double_layer, single_layer)
        real(dp), intent(in) :: corners(:, :, :)
        real(dp), allocatable, intent(in) :: wall_corners(:, :, :)
        real(dp), intent(in) :: centre(:, :), normal(:, :)
        integer, intent(in) :: order(:, :), plasma_count, plasma_block
        real(dp), intent(out) :: double_layer(:, :, :), single_layer(:, :, :)
        real(dp) :: direction(3), jump, point(3), sense, vertices(3, 3)
        integer :: block, column, index, row, sector, sectors

        block = size(order, 1)
        sectors = size(order, 2)
        !$omp parallel do collapse(2) private(direction, index, point, &
        !$omp sense, vertices, row)
        do sector = 1, sectors
            do column = 1, block
                index = order(column, sector)
                if (index <= plasma_count) then
                    vertices = corners(:, :, index)
                else
                    vertices = wall_corners(:, :, index - plasma_count)
                end if
                direction = normal(:, index)
                sense = orientation(vertices, direction) / (4.0_dp * pi)
                do row = 1, block
                    point = centre(:, order(row, 1))
                    if (index == order(row, 1)) then
                        double_layer(row, column, sector) = 0.0_dp
                    else
                        double_layer(row, column, sector) = &
                            -solid_angle(vertices, point) * sense
                    end if
                    if (column <= plasma_block) &
                        single_layer(row, column, sector) = &
                        triangle_potential(vertices, point) / (4.0_dp * pi)
                end do
            end do
        end do
        !$omp end parallel do
        do row = 1, block
            jump = 0.0_dp
            do sector = 1, sectors
                if (row <= plasma_block) then
                    do column = 1, plasma_block
                        jump = jump + double_layer(row, column, sector)
                    end do
                else
                    do column = plasma_block + 1, block
                        jump = jump + double_layer(row, column, sector)
                    end do
                end if
            end do
            double_layer(row, row, 1) = abs(jump)
        end do
    end subroutine assemble_sector_rows

    ! Energy of Fourier harmonic j over the sectors: with
    ! x^_j = sum_s x_s exp(-2 pi i j s / P), the block-circulant system
    ! sum_b A_(b-a) x_b = sum_b S_(b-a) d_b becomes A^_j x^_j = S^_j d^_j with
    ! A^_j = sum_d A_d exp(2 pi i j d / P), and the energy sum_s d_s . W x_s
    ! is (1/P) sum_j Re(conj(d^_j) . W x^_j). Harmonics 0 and P/2 are real.
    ! The walled region is bounded, so harmonic 0 carries the free constant
    ! of the Neumann potential, fixed by a zero mean.
    subroutine solve_harmonic(j, sectors, double_layer, single_layer, datum, &
            area, walled, total, energy, info)
        integer, intent(in) :: j, sectors, total
        real(dp), intent(in) :: double_layer(:, :, :), single_layer(:, :, :)
        real(dp), intent(in) :: datum(:, :, :), area(:)
        logical, intent(in) :: walled
        real(dp), intent(out) :: energy(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: system(:, :), rhs(:, :), transformed(:, :)
        real(dp), allocatable :: layer(:, :)
        complex(dp), allocatable :: system_c(:, :), rhs_c(:, :)
        complex(dp), allocatable :: transformed_c(:, :), layer_c(:, :)
        integer, allocatable :: pivots(:)
        real(dp) :: angle, phase(sectors)
        complex(dp) :: rotation(sectors)
        integer :: block, column, count, data, plasma_block, row, sector

        block = size(double_layer, 1)
        plasma_block = size(single_layer, 2)
        data = size(datum, 3)
        do sector = 1, sectors
            angle = 2.0_dp * pi * real(modulo(j * (sector - 1), sectors), dp) &
                / real(sectors, dp)
            phase(sector) = cos(angle)
            rotation(sector) = cmplx(cos(angle), sin(angle), dp)
        end do
        if (j == 0 .or. 2 * j == sectors) then
            if (2 * j == sectors .and. j > 0) then
                do sector = 1, sectors
                    phase(sector) = real(1 - 2 * modulo(sector - 1, 2), dp)
                end do
            else
                phase = 1.0_dp
            end if
            count = block
            if (walled .and. j == 0) count = block + 1
            allocate (system(count, count), rhs(count, data), &
                transformed(plasma_block, data), layer(block, plasma_block), &
                pivots(count))
            system = 0.0_dp
            layer = 0.0_dp
            transformed = 0.0_dp
            do sector = 1, sectors
                system(1:block, 1:block) = system(1:block, 1:block) &
                    + phase(sector) * double_layer(:, :, sector)
                layer = layer + phase(sector) * single_layer(:, :, sector)
                transformed = transformed + phase(sector) * datum(:, sector, :)
            end do
            rhs = 0.0_dp
            call dgemm("N", "N", block, data, plasma_block, 1.0_dp, layer, &
                block, transformed, plasma_block, 0.0_dp, rhs, count)
            if (count > block) then
                system(1:block, count) = 1.0_dp
                system(count, 1:block) = real(sectors, dp) / real(total, dp)
            end if
            call dgesv(count, data, system, count, pivots, rhs, count, info)
            if (info /= 0) return
            do column = 1, data
                do row = 1, data
                    energy(row, column) = sum(transformed(:, row) &
                        * area * rhs(1:plasma_block, column))
                end do
            end do
        else
            allocate (system_c(block, block), rhs_c(block, data), &
                transformed_c(plasma_block, data), &
                layer_c(block, plasma_block), pivots(block))
            system_c = (0.0_dp, 0.0_dp)
            layer_c = (0.0_dp, 0.0_dp)
            transformed_c = (0.0_dp, 0.0_dp)
            do sector = 1, sectors
                system_c = system_c + rotation(sector) &
                    * double_layer(:, :, sector)
                layer_c = layer_c + rotation(sector) &
                    * single_layer(:, :, sector)
                transformed_c = transformed_c + conjg(rotation(sector)) &
                    * datum(:, sector, :)
            end do
            call zgemm("N", "N", block, data, plasma_block, &
                (1.0_dp, 0.0_dp), layer_c, block, transformed_c, &
                plasma_block, (0.0_dp, 0.0_dp), rhs_c, block)
            call zgesv(block, data, system_c, block, pivots, rhs_c, block, &
                info)
            if (info /= 0) return
            do column = 1, data
                do row = 1, data
                    energy(row, column) = real(sum(conjg( &
                        transformed_c(:, row)) * area &
                        * rhs_c(1:plasma_block, column)), dp)
                end do
            end do
        end if
        info = 0
    end subroutine solve_harmonic

    ! Largest P dividing the toroidal node counts such that rotating every
    ! node by 2 pi / P about the z axis (either sense) gives the node nv / P
    ! toroidal indices on, on the edge and on the wall.
    integer function rotational_sectors(plasma, wall) result(sectors)
        real(dp), intent(in) :: plasma(:, :, :)
        real(dp), intent(in), optional :: wall(:, :, :)
        real(dp) :: sense
        integer :: candidate, turn

        sectors = 1
        do candidate = size(plasma, 3), 2, -1
            if (modulo(size(plasma, 3), candidate) /= 0) cycle
            if (present(wall)) then
                if (modulo(size(wall, 3), candidate) /= 0) cycle
            end if
            do turn = 1, 2
                sense = real(3 - 2 * turn, dp)
                if (.not. rotation_invariant(plasma, candidate, sense)) cycle
                if (present(wall)) then
                    if (.not. rotation_invariant(wall, candidate, sense)) cycle
                end if
                sectors = candidate
                return
            end do
        end do
    end function rotational_sectors

    pure logical function rotation_invariant(surface, sectors, sense) &
            result(invariant)
        real(dp), intent(in) :: surface(:, :, :), sense
        integer, intent(in) :: sectors
        real(dp) :: angle, cosine, rotated(3), scale, sine
        integer :: i, k, nv, shifted

        nv = size(surface, 3)
        angle = sense * 2.0_dp * pi / real(sectors, dp)
        cosine = cos(angle)
        sine = sin(angle)
        scale = maxval(abs(surface))
        invariant = .false.
        do k = 1, nv
            shifted = modulo(k - 1 + nv / sectors, nv) + 1
            do i = 1, size(surface, 2)
                rotated(1) = cosine * surface(1, i, k) - sine * surface(2, i, k)
                rotated(2) = sine * surface(1, i, k) + cosine * surface(2, i, k)
                rotated(3) = surface(3, i, k)
                if (norm2(rotated - surface(:, i, shifted)) &
                    > symmetry_tolerance * scale) return
            end do
        end do
        invariant = .true.
    end function rotation_invariant

    ! +1 when n_V agrees with the vertex-order normal of the triangle.
    pure real(dp) function orientation(vertices, normal_v) result(sign_value)
        real(dp), intent(in) :: vertices(3, 3), normal_v(3)
        real(dp) :: edge12(3), edge13(3), product(3)

        edge12 = vertices(:, 2) - vertices(:, 1)
        edge13 = vertices(:, 3) - vertices(:, 1)
        call cross_product(edge12, edge13, product)
        sign_value = sign(1.0_dp, dot_product(product, normal_v))
    end function orientation

    subroutine triangle_frames(corners, surface, centre, normal, area, &
            owner, first)
        real(dp), intent(in) :: corners(:, :, :)
        integer, intent(in) :: surface, first
        real(dp), intent(inout) :: centre(:, :), normal(:, :), area(:)
        integer, intent(inout) :: owner(:)
        real(dp) :: edge12(3), edge13(3), vector(3)
        integer :: t, index

        do t = 1, size(corners, 3)
            index = first + t - 1
            centre(:, index) = (corners(:, 1, t) + corners(:, 2, t) &
                + corners(:, 3, t)) / 3.0_dp
            edge12 = corners(:, 2, t) - corners(:, 1, t)
            edge13 = corners(:, 3, t) - corners(:, 1, t)
            call cross_product(edge12, edge13, vector)
            area(index) = 0.5_dp * norm2(vector)
            normal(:, index) = vector / norm2(vector)
            owner(index) = surface
        end do
    end subroutine triangle_frames

    ! Flip the normals of one surface so that they point along direction
    ! times its outward normal (divergence theorem: sum centre . n A > 0).
    subroutine orient(centre, normal, area, owner, surface, direction)
        real(dp), intent(in) :: centre(:, :), area(:)
        real(dp), intent(inout) :: normal(:, :)
        integer, intent(in) :: owner(:), surface
        real(dp), intent(in) :: direction
        real(dp) :: volume
        integer :: t

        volume = 0.0_dp
        do t = 1, size(area)
            if (owner(t) == surface) volume = volume + area(t) &
                * (centre(1, t) * normal(1, t) + centre(2, t) * normal(2, t) &
                + centre(3, t) * normal(3, t))
        end do
        do t = 1, size(area)
            if (owner(t) == surface) normal(:, t) = sign(1.0_dp, volume) &
                * direction * normal(:, t)
        end do
    end subroutine orient

    pure logical function valid_surface(surface) result(valid)
        real(dp), intent(in) :: surface(:, :, :)

        valid = size(surface, 1) == 3 .and. size(surface, 2) >= 3 &
            .and. size(surface, 3) >= 3
        if (valid) valid = all(ieee_is_finite(surface))
    end function valid_surface

    subroutine build_surface(surface, corners, info)
        real(dp), contiguous, intent(in) :: surface(:, :, :)
        real(dp), allocatable, intent(out) :: corners(:, :, :)
        integer, intent(out) :: info
        real(dp) :: vertices(3, 3)
        integer :: i, i1, k, k1, nu, nv, t

        nu = size(surface, 2)
        nv = size(surface, 3)
        allocate (corners(3, 3, 2 * nu * nv))
        t = 0
        do k = 1, nv
            k1 = modulo(k, nv) + 1
            do i = 1, nu
                i1 = modulo(i, nu) + 1
                t = t + 1
                corners(:, 1, t) = surface(:, i, k1)
                corners(:, 2, t) = surface(:, i1, k1)
                corners(:, 3, t) = surface(:, i, k)
                t = t + 1
                corners(:, 1, t) = surface(:, i1, k)
                corners(:, 2, t) = surface(:, i, k)
                corners(:, 3, t) = surface(:, i1, k1)
            end do
        end do
        info = vacuum_bie_degenerate
        do t = 1, size(corners, 3)
            vertices = corners(:, :, t)
            if (degenerate(vertices)) return
        end do
        info = vacuum_bie_ok
    end subroutine build_surface

    pure logical function degenerate(vertices)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: edge12(3), edge13(3), edge23(3), product(3), scale

        edge12 = vertices(:, 2) - vertices(:, 1)
        edge13 = vertices(:, 3) - vertices(:, 1)
        edge23 = vertices(:, 3) - vertices(:, 2)
        scale = max(norm2(edge12), norm2(edge13), norm2(edge23))
        call cross_product(edge12, edge13, product)
        degenerate = norm2(product) <= 256.0_dp * epsilon(1.0_dp) * scale**2
    end function degenerate

    ! The meshes do not intersect, every edge node lies inside the wall and
    ! no wall node inside the edge. With rotational symmetry the first of the
    ! sectors decides for all of them.
    function nested(plasma, wall, corners, wall_corners, sectors)
        real(dp), contiguous, intent(in) :: plasma(:, :, :), wall(:, :, :)
        real(dp), contiguous, intent(in) :: corners(:, :, :)
        real(dp), contiguous, intent(in) :: wall_corners(:, :, :)
        integer, intent(in) :: sectors
        logical :: nested
        real(dp) :: point(3)
        integer :: i, k

        nested = .false.
        if (meshes_intersect(corners(:, :, 1:size(corners, 3) / sectors), &
            wall_corners)) return
        do k = 1, size(plasma, 3) / sectors
            do i = 1, size(plasma, 2)
                point = plasma(:, i, k)
                if (.not. point_inside_mesh(point, wall_corners)) return
            end do
        end do
        do k = 1, size(wall, 3) / sectors
            do i = 1, size(wall, 2)
                point = wall(:, i, k)
                if (point_inside_mesh(point, corners)) return
            end do
        end do
        nested = .true.
    end function nested

    ! Signed solid angle of a flat triangle seen from point, positive when
    ! the vertex-order normal points away from it (Van Oosterom-Strackee).
    pure real(dp) function solid_angle(vertices, point) result(angle)
        real(dp), intent(in) :: vertices(3, 3), point(3)
        real(dp) :: a(3), b(3), c(3), la, lb, lc, numerator, denominator
        real(dp) :: product(3)

        a = vertices(:, 1) - point
        b = vertices(:, 2) - point
        c = vertices(:, 3) - point
        la = norm2(a)
        lb = norm2(b)
        lc = norm2(c)
        call cross_product(b, c, product)
        numerator = dot_product(a, product)
        denominator = la * lb * lc + dot_product(a, b) * lc &
            + dot_product(a, c) * lb + dot_product(b, c) * la
        angle = 2.0_dp * atan2(numerator, denominator)
    end function solid_angle

    ! Integral of 1/|point - y| over a flat triangle (exact).
    pure real(dp) function triangle_potential(vertices, point) result(value)
        real(dp), intent(in) :: vertices(3, 3), point(3)
        real(dp) :: a, d_minus, d_plus, denominator, edge(3), edge_length
        real(dp) :: h, l_minus, l_plus, next_offset(3), normal(3)
        real(dp) :: edge12(3), edge13(3), moment(3), offset(3), ratio
        integer :: i, next

        edge12 = vertices(:, 2) - vertices(:, 1)
        edge13 = vertices(:, 3) - vertices(:, 1)
        call cross_product(edge12, edge13, normal)
        normal = normal / norm2(normal)
        value = 0.0_dp
        do i = 1, 3
            next = modulo(i, 3) + 1
            edge = vertices(:, next) - vertices(:, i)
            edge_length = norm2(edge)
            next_offset = vertices(:, next) - point
            offset = vertices(:, i) - point
            l_plus = norm2(next_offset) / edge_length
            l_minus = norm2(offset) / edge_length
            ratio = (l_plus + l_minus + 1.0_dp) / (l_plus + l_minus - 1.0_dp)
            call cross_product(offset, edge, moment)
            a = dot_product(moment, normal) / edge_length**2
            h = abs(dot_product(offset, normal)) / edge_length
            value = value + edge_length * a * log(ratio)
            if (h > 64.0_dp * epsilon(1.0_dp)) then
                d_plus = dot_product(next_offset, edge) / edge_length**2
                d_minus = dot_product(offset, edge) / edge_length**2
                denominator = a * a + h * h
                value = value - edge_length * h &
                    * (atan(a * d_plus / (denominator + l_plus * h)) &
                    - atan(a * d_minus / (denominator + l_minus * h)))
            end if
        end do
    end function triangle_potential

    pure subroutine cross_product(a, b, c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(out) :: c(3)

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end subroutine cross_product

end module scalar_potential_vacuum
