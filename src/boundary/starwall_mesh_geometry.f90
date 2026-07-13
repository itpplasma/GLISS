module starwall_mesh_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: meshes_intersect, point_inside_mesh

contains

    function meshes_intersect(first, second) result(intersect)
        real(dp), intent(in) :: first(:, :, :), second(:, :, :)
        logical :: intersect
        integer :: i, j

        intersect = .false.
        do j = 1, size(second, 3)
            do i = 1, size(first, 3)
                if (triangles_intersect(first(:, :, i), second(:, :, j))) then
                    intersect = .true.
                    return
                end if
            end do
        end do
    end function meshes_intersect

    function point_inside_mesh(point, triangles) result(inside)
        real(dp), intent(in) :: point(3), triangles(:, :, :)
        logical :: inside
        real(dp) :: a(3), b(3), c(3), denominator, solid_angle
        integer :: triangle

        solid_angle = 0.0_dp
        do triangle = 1, size(triangles, 3)
            a = triangles(:, 1, triangle) - point
            b = triangles(:, 2, triangle) - point
            c = triangles(:, 3, triangle) - point
            denominator = norm2(a) * norm2(b) * norm2(c) &
                + dot_product(a, b) * norm2(c) &
                + dot_product(b, c) * norm2(a) &
                + dot_product(c, a) * norm2(b)
            solid_angle = solid_angle + 2.0_dp * atan2( &
                dot_product(a, cross_product(b, c)), denominator)
        end do
        inside = abs(solid_angle) > 2.0_dp * pi
    end function point_inside_mesh

    function triangles_intersect(first, second) result(intersect)
        real(dp), intent(in) :: first(3, 3), second(3, 3)
        logical :: intersect
        real(dp) :: distance_first(3), distance_second(3), normal_first(3)
        real(dp) :: normal_second(3), normal_scale, scale, tolerance
        integer :: edge, next

        normal_first = cross_product(first(:, 2) - first(:, 1), &
            first(:, 3) - first(:, 1))
        normal_second = cross_product(second(:, 2) - second(:, 1), &
            second(:, 3) - second(:, 1))
        scale = max(max_edge(first), max_edge(second))
        tolerance = 4096.0_dp * epsilon(1.0_dp) * scale
        distance_second = plane_distances(second, first(:, 1), normal_first)
        distance_first = plane_distances(first, second(:, 1), normal_second)
        if (same_strict_sign(distance_second, tolerance * norm2(normal_first)) &
            .or. same_strict_sign(distance_first, &
            tolerance * norm2(normal_second))) then
            intersect = .false.
            return
        end if

        normal_scale = norm2(normal_first) * norm2(normal_second)
        if (norm2(cross_product(normal_first, normal_second)) &
            <= 4096.0_dp * epsilon(1.0_dp) * normal_scale) then
            if (maxval(abs(distance_second)) > tolerance * norm2(normal_first)) then
                intersect = .false.
            else
                intersect = coplanar_triangles_intersect(first, second, normal_first)
            end if
            return
        end if

        intersect = .false.
        do edge = 1, 3
            next = modulo(edge, 3) + 1
            intersect = intersect .or. segment_intersects_triangle( &
                first(:, edge), first(:, next), second, tolerance)
            intersect = intersect .or. segment_intersects_triangle( &
                second(:, edge), second(:, next), first, tolerance)
        end do
    end function triangles_intersect

    pure function plane_distances(points, origin, normal) result(distance)
        real(dp), intent(in) :: points(3, 3), origin(3), normal(3)
        real(dp) :: distance(3)
        integer :: i

        do i = 1, 3
            distance(i) = dot_product(points(:, i) - origin, normal)
        end do
    end function plane_distances

    pure function same_strict_sign(values, tolerance) result(same)
        real(dp), intent(in) :: values(3), tolerance
        logical :: same

        same = all(values > tolerance) .or. all(values < -tolerance)
    end function same_strict_sign

    function segment_intersects_triangle(first, second, triangle, tolerance) &
            result(intersect)
        real(dp), intent(in) :: first(3), second(3), triangle(3, 3), tolerance
        logical :: intersect
        real(dp) :: determinant, direction(3), edge1(3), edge2(3), h(3)
        real(dp) :: inverse, q(3), s(3), parameter, u, v

        direction = second - first
        edge1 = triangle(:, 2) - triangle(:, 1)
        edge2 = triangle(:, 3) - triangle(:, 1)
        h = cross_product(direction, edge2)
        determinant = dot_product(edge1, h)
        if (abs(determinant) <= tolerance * norm2(edge1) * norm2(edge2)) then
            intersect = .false.
            return
        end if
        inverse = 1.0_dp / determinant
        s = first - triangle(:, 1)
        u = inverse * dot_product(s, h)
        q = cross_product(s, edge1)
        v = inverse * dot_product(direction, q)
        parameter = inverse * dot_product(edge2, q)
        intersect = u >= -tolerance .and. v >= -tolerance &
            .and. u + v <= 1.0_dp + tolerance &
            .and. parameter >= -tolerance .and. parameter <= 1.0_dp + tolerance
    end function segment_intersects_triangle

    function coplanar_triangles_intersect(first, second, normal) result(intersect)
        real(dp), intent(in) :: first(3, 3), second(3, 3), normal(3)
        logical :: intersect
        real(dp) :: first_2d(2, 3), second_2d(2, 3), tolerance
        integer :: drop, edge_first, edge_second, next_first, next_second

        drop = maxloc(abs(normal), dim=1)
        call project_triangle(first, drop, first_2d)
        call project_triangle(second, drop, second_2d)
        tolerance = 4096.0_dp * epsilon(1.0_dp) &
            * max(maxval(abs(first_2d)), maxval(abs(second_2d)), 1.0_dp)
        intersect = point_in_triangle_2d(first_2d(:, 1), second_2d, tolerance) &
            .or. point_in_triangle_2d(second_2d(:, 1), first_2d, tolerance)
        do edge_second = 1, 3
            next_second = modulo(edge_second, 3) + 1
            do edge_first = 1, 3
                next_first = modulo(edge_first, 3) + 1
                intersect = intersect .or. segments_intersect_2d( &
                    first_2d(:, edge_first), first_2d(:, next_first), &
                    second_2d(:, edge_second), second_2d(:, next_second), &
                    tolerance)
            end do
        end do
    end function coplanar_triangles_intersect

    pure subroutine project_triangle(points, drop, projected)
        real(dp), intent(in) :: points(3, 3)
        integer, intent(in) :: drop
        real(dp), intent(out) :: projected(2, 3)
        integer :: coordinate, output

        output = 0
        do coordinate = 1, 3
            if (coordinate == drop) cycle
            output = output + 1
            projected(output, :) = points(coordinate, :)
        end do
    end subroutine project_triangle

    pure function point_in_triangle_2d(point, triangle, tolerance) result(inside)
        real(dp), intent(in) :: point(2), triangle(2, 3), tolerance
        logical :: inside
        real(dp) :: orientation(3)

        orientation = [orient2d(triangle(:, 1), triangle(:, 2), point), &
            orient2d(triangle(:, 2), triangle(:, 3), point), &
            orient2d(triangle(:, 3), triangle(:, 1), point)]
        inside = all(orientation >= -tolerance) &
            .or. all(orientation <= tolerance)
    end function point_in_triangle_2d

    pure function segments_intersect_2d(a, b, c, d, tolerance) result(intersect)
        real(dp), intent(in) :: a(2), b(2), c(2), d(2), tolerance
        logical :: intersect
        real(dp) :: ab_c, ab_d, cd_a, cd_b

        ab_c = orient2d(a, b, c)
        ab_d = orient2d(a, b, d)
        cd_a = orient2d(c, d, a)
        cd_b = orient2d(c, d, b)
        intersect = ab_c * ab_d <= tolerance**2 &
            .and. cd_a * cd_b <= tolerance**2 &
            .and. max(min(a(1), b(1)), min(c(1), d(1))) &
            <= min(max(a(1), b(1)), max(c(1), d(1))) + tolerance &
            .and. max(min(a(2), b(2)), min(c(2), d(2))) &
            <= min(max(a(2), b(2)), max(c(2), d(2))) + tolerance
    end function segments_intersect_2d

    pure real(dp) function orient2d(a, b, c) result(value)
        real(dp), intent(in) :: a(2), b(2), c(2)

        value = (b(1) - a(1)) * (c(2) - a(2)) &
            - (b(2) - a(2)) * (c(1) - a(1))
    end function orient2d

    pure real(dp) function max_edge(triangle) result(value)
        real(dp), intent(in) :: triangle(3, 3)

        value = max(norm2(triangle(:, 2) - triangle(:, 1)), &
            norm2(triangle(:, 3) - triangle(:, 2)), &
            norm2(triangle(:, 1) - triangle(:, 3)))
    end function max_edge

    pure function cross_product(first, second) result(value)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: value(3)

        value = [first(2) * second(3) - first(3) * second(2), &
            first(3) * second(1) - first(1) * second(3), &
            first(1) * second(2) - first(2) * second(1)]
    end function cross_product

end module starwall_mesh_geometry
