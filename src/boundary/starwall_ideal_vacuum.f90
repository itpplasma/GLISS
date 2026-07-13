module starwall_ideal_vacuum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use starwall_mesh_geometry, only: meshes_intersect, point_inside_mesh
    implicit none
    private

    integer, parameter, public :: starwall_ok = 0
    integer, parameter, public :: starwall_invalid_input = 1
    integer, parameter, public :: starwall_degenerate_surface = 2
    integer, parameter, public :: starwall_surfaces_not_nested = 3
    integer, parameter, public :: starwall_not_spd = 4

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: triangle_t
        real(dp) :: r(3, 3)
        real(dp) :: uv(2, 3)
        real(dp) :: edge_current(3, 3)
        real(dp) :: cross_norm
        integer :: node(3)
        integer :: surface
    end type triangle_t

    public :: assemble_starwall_ideal_vacuum

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf

        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
    end interface

contains

    subroutine assemble_starwall_ideal_vacuum(plasma, fp, ft, stiffness, &
            response, info, wall)
        real(dp), intent(in) :: plasma(:, :, :), fp, ft
        real(dp), allocatable, intent(out) :: stiffness(:, :), response(:, :)
        integer, intent(out) :: info
        real(dp), intent(in), optional :: wall(:, :, :)
        type(triangle_t), allocatable :: triangles(:), plasma_triangles(:)
        type(triangle_t), allocatable :: wall_triangles(:)
        real(dp), allocatable :: coupling(:, :), inductance(:, :)
        integer :: current_count, lapack_info, plasma_nodes, wall_nodes

        info = starwall_invalid_input
        if (.not. valid_surface_array(plasma)) return
        if (.not. ieee_is_finite(fp) .or. .not. ieee_is_finite(ft)) return
        plasma_nodes = size(plasma, 2) * size(plasma, 3)
        wall_nodes = 0
        if (present(wall)) then
            if (.not. valid_surface_array(wall)) return
            wall_nodes = size(wall, 2) * size(wall, 3)
        end if

        call build_triangles(plasma, 1, plasma_triangles, info)
        if (info /= starwall_ok) return
        if (present(wall)) then
            call build_triangles(wall, 2, wall_triangles, info)
            if (info /= starwall_ok) return
            if (.not. nested_surfaces(plasma, wall, plasma_triangles, &
                wall_triangles)) then
                info = starwall_surfaces_not_nested
                return
            end if
            allocate (triangles(size(plasma_triangles) + size(wall_triangles)))
            triangles(:size(plasma_triangles)) = plasma_triangles
            triangles(size(plasma_triangles) + 1:) = wall_triangles
        else
            allocate (triangles, source=plasma_triangles)
        end if

        current_count = plasma_nodes + wall_nodes
        allocate (inductance(current_count, current_count), &
            coupling(current_count, plasma_nodes))
        call assemble_inductance(triangles, plasma_nodes, inductance)
        call assemble_coupling(plasma_triangles, size(plasma, 2), &
            size(plasma, 3), fp, ft, coupling)

        response = -coupling
        call dpotrf("L", current_count, inductance, current_count, lapack_info)
        if (lapack_info /= 0) then
            info = starwall_not_spd
            return
        end if
        call dpotrs("L", current_count, plasma_nodes, inductance, &
            current_count, response, current_count, lapack_info)
        if (lapack_info /= 0) then
            info = starwall_invalid_input
            return
        end if
        stiffness = -matmul(transpose(coupling), response)
        stiffness = 0.5_dp * (stiffness + transpose(stiffness))
        info = starwall_ok
    end subroutine assemble_starwall_ideal_vacuum

    pure function valid_surface_array(surface) result(valid)
        real(dp), intent(in) :: surface(:, :, :)
        logical :: valid

        valid = size(surface, 1) == 3 .and. size(surface, 2) >= 3 .and. &
            size(surface, 3) >= 3 .and. all(ieee_is_finite(surface))
    end function valid_surface_array

    subroutine build_triangles(surface, surface_index, triangles, info)
        real(dp), intent(in) :: surface(:, :, :)
        integer, intent(in) :: surface_index
        type(triangle_t), allocatable, intent(out) :: triangles(:)
        integer, intent(out) :: info
        integer :: i, i1, k, k1, nu, nv, triangle_index
        real(dp) :: u0, u1, v0, v1

        nu = size(surface, 2)
        nv = size(surface, 3)
        allocate (triangles(2 * nu * nv))
        triangle_index = 0
        do k = 1, nv
            k1 = modulo(k, nv) + 1
            v0 = real(k - 1, dp) / real(nv, dp)
            v1 = real(k, dp) / real(nv, dp)
            do i = 1, nu
                i1 = modulo(i, nu) + 1
                u0 = real(i - 1, dp) / real(nu, dp)
                u1 = real(i, dp) / real(nu, dp)
                triangle_index = triangle_index + 1
                call set_triangle(triangles(triangle_index), &
                    reshape([surface(:, i, k1), surface(:, i1, k1), &
                    surface(:, i, k)], [3, 3]), &
                    reshape([u0, v1, u1, v1, u0, v0], [2, 3]), &
                    [node_index(i, k1, nu), node_index(i1, k1, nu), &
                    node_index(i, k, nu)], surface_index, info)
                if (info /= starwall_ok) return
                triangle_index = triangle_index + 1
                call set_triangle(triangles(triangle_index), &
                    reshape([surface(:, i1, k), surface(:, i, k), &
                    surface(:, i1, k1)], [3, 3]), &
                    reshape([u1, v0, u0, v0, u1, v1], [2, 3]), &
                    [node_index(i1, k, nu), node_index(i, k, nu), &
                    node_index(i1, k1, nu)], surface_index, info)
                if (info /= starwall_ok) return
            end do
        end do
        info = starwall_ok
    end subroutine build_triangles

    subroutine set_triangle(triangle, points, uv, nodes, surface, info)
        type(triangle_t), intent(out) :: triangle
        real(dp), intent(in) :: points(3, 3), uv(2, 3)
        integer, intent(in) :: nodes(3), surface
        integer, intent(out) :: info
        real(dp) :: cross_value(3), edge_scale

        triangle%r = points
        triangle%uv = uv
        triangle%node = nodes
        triangle%surface = surface
        cross_value = cross_product(points(:, 2) - points(:, 1), &
            points(:, 3) - points(:, 1))
        triangle%cross_norm = norm2(cross_value)
        edge_scale = maxval([norm2(points(:, 2) - points(:, 1)), &
            norm2(points(:, 3) - points(:, 2)), &
            norm2(points(:, 1) - points(:, 3))])
        if (triangle%cross_norm <= 256.0_dp * epsilon(1.0_dp) * edge_scale**2) then
            info = starwall_degenerate_surface
            return
        end if
        triangle%edge_current(:, 1) = (points(:, 2) - points(:, 3)) &
            / triangle%cross_norm
        triangle%edge_current(:, 2) = (points(:, 3) - points(:, 1)) &
            / triangle%cross_norm
        triangle%edge_current(:, 3) = (points(:, 1) - points(:, 2)) &
            / triangle%cross_norm
        info = starwall_ok
    end subroutine set_triangle

    pure integer function node_index(i, k, nu) result(index)
        integer, intent(in) :: i, k, nu

        index = i + nu * (k - 1)
    end function node_index

    subroutine assemble_inductance(triangles, plasma_nodes, matrix)
        type(triangle_t), intent(in) :: triangles(:)
        integer, intent(in) :: plasma_nodes
        real(dp), intent(out) :: matrix(:, :)
        real(dp) :: local(3, 3), scalar
        integer :: first, second

        matrix = 0.0_dp
        do second = 1, size(triangles)
            do first = 1, second
                if (first == second) then
                    call self_inductance(triangles(first), local)
                else
                    scalar = triangle_pair_integral(triangles(first), &
                        triangles(second))
                    call mutual_inductance(triangles(first), triangles(second), &
                        scalar, local)
                end if
                call scatter_inductance(triangles(first), triangles(second), &
                    local, plasma_nodes, first /= second, matrix)
            end do
        end do
    end subroutine assemble_inductance

    subroutine self_inductance(triangle, local)
        type(triangle_t), intent(in) :: triangle
        real(dp), intent(out) :: local(3, 3)
        real(dp) :: edge(3), factor, perimeter
        integer :: i, j

        edge = [norm2(triangle%r(:, 2) - triangle%r(:, 1)), &
            norm2(triangle%r(:, 3) - triangle%r(:, 2)), &
            norm2(triangle%r(:, 1) - triangle%r(:, 3))]
        perimeter = sum(edge)
        factor = triangle%cross_norm**2 / (12.0_dp * pi) * &
            sum(log(perimeter / (perimeter - 2.0_dp * edge)) / edge)
        do j = 1, 3
            do i = 1, 3
                local(i, j) = dot_product(triangle%edge_current(:, i), &
                    triangle%edge_current(:, j)) * factor
            end do
        end do
    end subroutine self_inductance

    subroutine mutual_inductance(first, second, scalar, local)
        type(triangle_t), intent(in) :: first, second
        real(dp), intent(in) :: scalar
        real(dp), intent(out) :: local(3, 3)
        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                local(i, j) = dot_product(first%edge_current(:, i), &
                    second%edge_current(:, j)) * scalar / (4.0_dp * pi)
            end do
        end do
    end subroutine mutual_inductance

    subroutine scatter_inductance(first, second, local, plasma_nodes, mirror, &
            matrix)
        type(triangle_t), intent(in) :: first, second
        real(dp), intent(in) :: local(3, 3)
        integer, intent(in) :: plasma_nodes
        logical, intent(in) :: mirror
        real(dp), intent(inout) :: matrix(:, :)
        integer :: a, b, ia(2), ib(2), na, nb, p, q
        real(dp) :: ca(2), cb(2)

        do b = 1, 3
            call local_current_map(second, b, plasma_nodes, ib, cb, nb)
            do a = 1, 3
                call local_current_map(first, a, plasma_nodes, ia, ca, na)
                do q = 1, nb
                    do p = 1, na
                        matrix(ia(p), ib(q)) = matrix(ia(p), ib(q)) &
                            + ca(p) * cb(q) * local(a, b)
                        if (mirror) matrix(ib(q), ia(p)) = &
                            matrix(ib(q), ia(p)) + ca(p) * cb(q) * local(a, b)
                    end do
                end do
            end do
        end do
    end subroutine scatter_inductance

    subroutine local_current_map(triangle, vertex, plasma_nodes, index, &
            coefficient, count)
        type(triangle_t), intent(in) :: triangle
        integer, intent(in) :: vertex, plasma_nodes
        integer, intent(out) :: index(2), count
        real(dp), intent(out) :: coefficient(2)
        integer :: offset

        offset = merge(0, plasma_nodes, triangle%surface == 1)
        count = 1
        index(1) = offset + 1
        coefficient(1) = triangle%uv(1, vertex)
        if (triangle%node(vertex) > 1) then
            count = 2
            index(2) = offset + triangle%node(vertex)
            coefficient(2) = 1.0_dp
        end if
    end subroutine local_current_map

    subroutine assemble_coupling(triangles, nu, nv, fp, ft, coupling)
        type(triangle_t), intent(in) :: triangles(:)
        integer, intent(in) :: nu, nv
        real(dp), intent(in) :: fp, ft
        real(dp), intent(out) :: coupling(:, :)
        real(dp) :: area, derivative(3), determinant, d_du(3), d_dv(3)
        integer :: a, b, current_index, t

        coupling = 0.0_dp
        do t = 1, size(triangles)
            determinant = (triangles(t)%uv(1, 2) - triangles(t)%uv(1, 1)) &
                * (triangles(t)%uv(2, 3) - triangles(t)%uv(2, 1)) &
                - (triangles(t)%uv(1, 3) - triangles(t)%uv(1, 1)) &
                * (triangles(t)%uv(2, 2) - triangles(t)%uv(2, 1))
            area = 0.5_dp * abs(determinant)
            d_du = [triangles(t)%uv(2, 2) - triangles(t)%uv(2, 3), &
                triangles(t)%uv(2, 3) - triangles(t)%uv(2, 1), &
                triangles(t)%uv(2, 1) - triangles(t)%uv(2, 2)] / determinant
            d_dv = [triangles(t)%uv(1, 3) - triangles(t)%uv(1, 2), &
                triangles(t)%uv(1, 1) - triangles(t)%uv(1, 3), &
                triangles(t)%uv(1, 2) - triangles(t)%uv(1, 1)] / determinant
            derivative = fp * d_du + ft * d_dv
            do a = 1, 3
                if (triangles(t)%node(a) == 1) cycle
                current_index = triangles(t)%node(a)
                do b = 1, 3
                    coupling(current_index, triangles(t)%node(b)) = &
                        coupling(current_index, triangles(t)%node(b)) &
                        + area * derivative(b) / 3.0_dp
                end do
            end do
        end do
        coupling(1, :) = coupling(1, :) - fp / real(nu * nv, dp)
    end subroutine assemble_coupling

    function triangle_pair_integral(source, target) result(value)
        type(triangle_t), intent(in) :: source, target
        real(dp) :: value
        integer :: refinement

        refinement = merge(3, 0, triangles_touch(source, target))
        value = integrate_triangle(source, target%r, refinement)
    end function triangle_pair_integral

    recursive function integrate_triangle(source, target, level) result(value)
        type(triangle_t), intent(in) :: source
        real(dp), intent(in) :: target(3, 3)
        integer, intent(in) :: level
        real(dp) :: value
        real(dp), parameter :: root15 = sqrt(15.0_dp)
        real(dp), parameter :: zeta(7) = [1.0_dp / 3.0_dp, &
            (6.0_dp + root15) / 21.0_dp, (9.0_dp - 2.0_dp * root15) / 21.0_dp, &
            (6.0_dp + root15) / 21.0_dp, (6.0_dp - root15) / 21.0_dp, &
            (9.0_dp + 2.0_dp * root15) / 21.0_dp, (6.0_dp - root15) / 21.0_dp]
        real(dp), parameter :: eta(7) = [1.0_dp / 3.0_dp, &
            (6.0_dp + root15) / 21.0_dp, (6.0_dp + root15) / 21.0_dp, &
            (9.0_dp - 2.0_dp * root15) / 21.0_dp, (6.0_dp - root15) / 21.0_dp, &
            (6.0_dp - root15) / 21.0_dp, (9.0_dp + 2.0_dp * root15) / 21.0_dp]
        real(dp), parameter :: weight(7) = [270.0_dp / 2400.0_dp, &
            (155.0_dp + root15) / 2400.0_dp, &
            (155.0_dp + root15) / 2400.0_dp, &
            (155.0_dp + root15) / 2400.0_dp, &
            (155.0_dp - root15) / 2400.0_dp, &
            (155.0_dp - root15) / 2400.0_dp, &
            (155.0_dp - root15) / 2400.0_dp]
        real(dp) :: point(3)
        real(dp) :: midpoint12(3), midpoint23(3), midpoint31(3)
        integer :: q

        if (level > 0) then
            midpoint12 = 0.5_dp * (target(:, 1) + target(:, 2))
            midpoint23 = 0.5_dp * (target(:, 2) + target(:, 3))
            midpoint31 = 0.5_dp * (target(:, 3) + target(:, 1))
            value = integrate_triangle(source, reshape([target(:, 1), &
                midpoint12, midpoint31], [3, 3]), level - 1) &
                + integrate_triangle(source, reshape([midpoint12, target(:, 2), &
                midpoint23], [3, 3]), level - 1) &
                + integrate_triangle(source, reshape([midpoint31, midpoint23, &
                target(:, 3)], [3, 3]), level - 1) &
                + integrate_triangle(source, reshape([midpoint12, midpoint23, &
                midpoint31], [3, 3]), level - 1)
            return
        end if
        value = 0.0_dp
        do q = 1, 7
            point = target(:, 1) + zeta(q) * &
                (target(:, 2) - target(:, 1)) + eta(q) * &
                (target(:, 3) - target(:, 1))
            value = value + weight(q) * triangle_potential(source, point)
        end do
        value = norm2(cross_product(target(:, 2) - target(:, 1), &
            target(:, 3) - target(:, 1))) * value
    end function integrate_triangle

    pure function triangles_touch(first, second) result(touch)
        type(triangle_t), intent(in) :: first, second
        logical :: touch
        integer :: i, j

        touch = .false.
        if (first%surface /= second%surface) return
        do j = 1, 3
            do i = 1, 3
                touch = touch .or. first%node(i) == second%node(j)
            end do
        end do
    end function triangles_touch

    function triangle_potential(triangle, point) result(value)
        type(triangle_t), intent(in) :: triangle
        real(dp), intent(in) :: point(3)
        real(dp) :: value
        real(dp) :: a, d_minus, d_plus, edge(3), edge_length, h
        real(dp) :: l_minus, l_plus, normal(3), ratio, denominator
        integer :: i, next

        normal = cross_product(triangle%r(:, 2) - triangle%r(:, 1), &
            triangle%r(:, 3) - triangle%r(:, 1)) / triangle%cross_norm
        value = 0.0_dp
        do i = 1, 3
            next = modulo(i, 3) + 1
            edge = triangle%r(:, next) - triangle%r(:, i)
            edge_length = norm2(edge)
            l_plus = norm2(triangle%r(:, next) - point) / edge_length
            l_minus = norm2(triangle%r(:, i) - point) / edge_length
            ratio = (l_plus + l_minus + 1.0_dp) &
                / (l_plus + l_minus - 1.0_dp)
            a = dot_product(cross_product(triangle%r(:, i) - point, edge), &
                normal) / edge_length**2
            h = abs(dot_product(triangle%r(:, i) - point, normal)) / edge_length
            value = value + edge_length * a * log(ratio)
            if (h > 64.0_dp * epsilon(1.0_dp)) then
                d_plus = dot_product(triangle%r(:, next) - point, edge) &
                    / edge_length**2
                d_minus = dot_product(triangle%r(:, i) - point, edge) &
                    / edge_length**2
                denominator = a * a + h * h
                value = value - edge_length * h * &
                    (atan(a * d_plus / (denominator + l_plus * h)) &
                    - atan(a * d_minus / (denominator + l_minus * h)))
            end if
        end do
    end function triangle_potential

    function nested_surfaces(plasma, wall, plasma_triangles, wall_triangles) &
            result(nested)
        real(dp), intent(in) :: plasma(:, :, :), wall(:, :, :)
        type(triangle_t), intent(in) :: plasma_triangles(:), wall_triangles(:)
        logical :: nested
        real(dp), allocatable :: plasma_mesh(:, :, :), wall_mesh(:, :, :)
        integer :: i, k, triangle

        nested = .true.
        allocate (plasma_mesh(3, 3, size(plasma_triangles)), &
            wall_mesh(3, 3, size(wall_triangles)))
        do triangle = 1, size(plasma_triangles)
            plasma_mesh(:, :, triangle) = plasma_triangles(triangle)%r
        end do
        do triangle = 1, size(wall_triangles)
            wall_mesh(:, :, triangle) = wall_triangles(triangle)%r
        end do
        if (meshes_intersect(plasma_mesh, wall_mesh)) then
            nested = .false.
            return
        end if
        do k = 1, size(plasma, 3)
            do i = 1, size(plasma, 2)
                if (.not. point_inside_mesh(plasma(:, i, k), wall_mesh)) then
                    nested = .false.
                    return
                end if
            end do
        end do
        do k = 1, size(wall, 3)
            do i = 1, size(wall, 2)
                if (point_inside_mesh(wall(:, i, k), plasma_mesh)) then
                    nested = .false.
                    return
                end if
            end do
        end do
    end function nested_surfaces

    pure function cross_product(first, second) result(value)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: value(3)

        value = [first(2) * second(3) - first(3) * second(2), &
            first(3) * second(1) - first(1) * second(3), &
            first(1) * second(2) - first(2) * second(1)]
    end function cross_product

end module starwall_ideal_vacuum
