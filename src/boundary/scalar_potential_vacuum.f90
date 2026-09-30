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

    public :: assemble_exterior_vacuum, edge_triangle_centroids

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
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
    subroutine assemble_exterior_vacuum(plasma, flux, energy, info, wall)
        real(dp), contiguous, intent(in) :: plasma(:, :, :)
        real(dp), intent(in) :: flux(:, :)
        real(dp), allocatable, intent(out) :: energy(:, :)
        integer, intent(out) :: info
        real(dp), contiguous, intent(in), optional :: wall(:, :, :)
        real(dp), allocatable :: area(:), centre(:, :), corners(:, :, :)
        real(dp), allocatable :: normal(:, :), system(:, :), rhs(:, :)
        real(dp), allocatable :: wall_corners(:, :, :), datum(:, :)
        integer, allocatable :: owner(:), pivots(:)
        real(dp) :: symmetric, uv_area
        integer :: column, count, data, lapack_info, plasma_count, row, total

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
            if (.not. nested(plasma, wall, corners, wall_corners)) then
                info = vacuum_bie_not_nested
                return
            end if
            total = total + size(wall_corners, 3)
        end if
        allocate (centre(3, total), normal(3, total), area(total), &
            owner(total))
        call triangle_frames(corners, 1, centre, normal, area, owner, 1)
        if (present(wall)) call triangle_frames(wall_corners, 2, centre, &
            normal, area, owner, plasma_count + 1)
        ! n_V points into the plasma on its edge and out through the wall.
        call orient(centre, normal, area, owner, 1, -1.0_dp)
        if (present(wall)) call orient(centre, normal, area, owner, 2, 1.0_dp)

        count = total
        if (present(wall)) count = total + 1
        data = size(flux, 2)
        allocate (system(count, count), rhs(count, data), &
            datum(plasma_count, data), pivots(count))
        uv_area = 1.0_dp / real(plasma_count, dp)
        do row = 1, plasma_count
            ! Mean dPhi/dn_V over the triangle; the sign convention cancels
            ! in the quadratic energy.
            datum(row, :) = flux(row, :) * uv_area / area(row)
        end do
        call assemble_system(corners, wall_corners, centre, normal, owner, &
            plasma_count, total, datum, system, rhs, present(wall))
        call dgesv(count, data, system, count, pivots, rhs, count, &
            lapack_info)
        if (lapack_info /= 0) then
            info = vacuum_bie_singular
            return
        end if
        allocate (energy(data, data))
        do column = 1, data
            do row = 1, data
                energy(row, column) = sum(datum(:, row) &
                    * rhs(1:plasma_count, column) * area(1:plasma_count))
            end do
        end do
        do column = 1, data
            do row = 1, column - 1
                symmetric = 0.5_dp * (energy(row, column) + energy(column, row))
                energy(row, column) = symmetric
                energy(column, row) = symmetric
            end do
        end do
        info = vacuum_bie_singular
        if (.not. all(ieee_is_finite(energy))) return
        info = vacuum_bie_ok
    end subroutine assemble_exterior_vacuum

    subroutine assemble_system(corners, wall_corners, centre, normal, owner, &
            plasma_count, total, datum, system, rhs, walled)
        real(dp), intent(in) :: corners(:, :, :)
        real(dp), allocatable, intent(in) :: wall_corners(:, :, :)
        real(dp), intent(in) :: centre(:, :), normal(:, :)
        integer, intent(in) :: owner(:), plasma_count, total
        real(dp), intent(in) :: datum(:, :)
        real(dp), intent(out) :: system(:, :), rhs(:, :)
        logical, intent(in) :: walled
        real(dp) :: double_layer, jump(2), point(3), single, vertices(3, 3)
        real(dp) :: direction(3)
        integer :: column, row

        system = 0.0_dp
        rhs = 0.0_dp
        !$omp parallel do private(column, direction, double_layer, jump, &
        !$omp point, single, vertices)
        do row = 1, total
            jump = 0.0_dp
            point = centre(:, row)
            do column = 1, total
                if (column <= plasma_count) then
                    vertices = corners(:, :, column)
                else
                    vertices = wall_corners(:, :, column - plasma_count)
                end if
                double_layer = 0.0_dp
                direction = normal(:, column)
                if (column /= row) double_layer = -solid_angle(vertices, &
                    point) * orientation(vertices, direction) / (4.0_dp * pi)
                system(row, column) = double_layer
                if (owner(column) == owner(row)) jump(owner(row)) = &
                    jump(owner(row)) + double_layer
                if (column <= plasma_count) then
                    single = triangle_potential(vertices, point) &
                        / (4.0_dp * pi)
                    rhs(row, :) = rhs(row, :) + single * datum(column, :)
                end if
            end do
            ! Continuum: +1/2 on the plasma edge, -1/2 on the wall, whose
            ! n_V is the outward normal of the enclosed body.
            system(row, row) = abs(jump(owner(row)))
        end do
        !$omp end parallel do
        if (walled) then
            ! V is bounded: fix the free constant of the Neumann potential.
            system(:, total + 1) = 0.0_dp
            system(1:total, total + 1) = 1.0_dp
            system(total + 1, :) = 0.0_dp
            system(total + 1, 1:total) = 1.0_dp / real(total, dp)
            rhs(total + 1, :) = 0.0_dp
        end if
    end subroutine assemble_system

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

    function nested(plasma, wall, corners, wall_corners)
        real(dp), contiguous, intent(in) :: plasma(:, :, :), wall(:, :, :)
        real(dp), contiguous, intent(in) :: corners(:, :, :)
        real(dp), contiguous, intent(in) :: wall_corners(:, :, :)
        logical :: nested
        real(dp) :: point(3)
        integer :: i, k

        nested = .false.
        if (meshes_intersect(corners, wall_corners)) return
        do k = 1, size(plasma, 3)
            do i = 1, size(plasma, 2)
                point = plasma(:, i, k)
                if (.not. point_inside_mesh(point, wall_corners)) return
            end do
        end do
        do k = 1, size(wall, 3)
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
