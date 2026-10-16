module boundary_mesh_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: meshes_intersect, point_inside_mesh
    public :: periodic_surface_is_embedded, triangle_pair_is_embedded

contains

    ! Admit the piecewise-linear realization of the periodic parameter mesh.
    ! Every face pair must intersect only in its common abstract simplex.
    ! This says nothing about a smooth interpolant between the sampled nodes.
    logical function periodic_surface_is_embedded(surface, orientation) &
            result(embedded)
        real(dp), intent(in) :: surface(:, :, :)
        integer, intent(out) :: orientation
        real(dp), allocatable :: corners(:, :, :), boxes(:, :, :)
        integer, allocatable :: nodes(:, :)
        real(dp) :: origin(3), scale, volume, magnitude, contribution, edges
        real(dp) :: normal(3), product(3)
        integer :: nu, nv, i, i1, k, k1, t, first, triangles

        embedded = .false.
        orientation = 0
        if (size(surface, 1) /= 3) return
        if (.not. all(ieee_is_finite(surface))) return
        nu = size(surface, 2)
        nv = size(surface, 3)
        if (nu < 3 .or. nv < 3) return
        if (nu > huge(1) / nv / 2) return
        origin = surface(:, 1, 1)
        scale = 0.0_dp
        do k = 1, nv
            do i = 1, nu
                scale = max(scale, norm2(surface(:, i, k) - origin))
            end do
        end do
        if (.not. ieee_is_finite(scale) .or. scale <= tiny(1.0_dp)) return
        triangles = 2 * nu * nv
        allocate (corners(3, 3, triangles), nodes(3, triangles), &
            boxes(3, 2, triangles))
        t = 0
        do k = 1, nv
            k1 = modulo(k, nv) + 1
            do i = 1, nu
                i1 = modulo(i, nu) + 1
                ! Match scalar_potential_vacuum's diagonal and orientation.
                t = t + 1
                nodes(1, t) = i + nu * (k1 - 1)
                nodes(2, t) = i1 + nu * (k1 - 1)
                nodes(3, t) = i + nu * (k - 1)
                corners(:, 1, t) = (surface(:, i, k1) - origin) / scale
                corners(:, 2, t) = (surface(:, i1, k1) - origin) / scale
                corners(:, 3, t) = (surface(:, i, k) - origin) / scale
                t = t + 1
                nodes(1, t) = i1 + nu * (k - 1)
                nodes(2, t) = i + nu * (k - 1)
                nodes(3, t) = i1 + nu * (k1 - 1)
                corners(:, 1, t) = (surface(:, i1, k) - origin) / scale
                corners(:, 2, t) = (surface(:, i, k) - origin) / scale
                corners(:, 3, t) = (surface(:, i1, k1) - origin) / scale
            end do
        end do
        volume = 0.0_dp
        magnitude = 0.0_dp
        do t = 1, triangles
            edges = max_edge(corners(:, :, t))
            normal = face_normal(corners(:, :, t))
            if (norm2(normal) &
                <= 8192.0_dp * epsilon(1.0_dp) * edges**2) return
            product = vector_cross(corners(:, 2, t), corners(:, 3, t))
            contribution = dot_product(corners(:, 1, t), product)
            volume = volume + contribution
            magnitude = magnitude + abs(contribution)
        end do
        if (abs(volume) <= 8192.0_dp * epsilon(1.0_dp) * magnitude) return
        call bounding_boxes(corners, boxes)
        embedded = .true.
        !$omp parallel do default(none) schedule(dynamic) &
        !$omp shared(corners, nodes, boxes, triangles) private(first) &
        !$omp reduction(.and.:embedded)
        do t = 2, triangles
            do first = 1, t - 1
                if (.not. embedded) cycle
                if (any(boxes(:, 1, first) > boxes(:, 2, t)) &
                    .or. any(boxes(:, 1, t) > boxes(:, 2, first))) cycle
                embedded = triangle_pair_is_embedded(corners(:, :, first), &
                    corners(:, :, t), nodes(:, first), nodes(:, t))
            end do
        end do
        !$omp end parallel do
        if (embedded) orientation = int(sign(1.0_dp, volume))
    end function periodic_surface_is_embedded

    ! Closed cones include boundary rays: an extra edge contact is invalid.
    ! Near-coplanar overlapping cones are rejected conservatively.
    pure logical function triangle_pair_is_embedded(first_input, second_input, &
            first_nodes, second_nodes) result(embedded)
        real(dp), intent(in) :: first_input(3, 3), second_input(3, 3)
        integer, intent(in) :: first_nodes(3), second_nodes(3)
        integer :: shared_first(3), shared_second(3), count, i, j, other
        real(dp) :: first(3, 3), second(3, 3)
        real(dp) :: a(3), b(3), c(3), d(3), n1(3), n2(3), direction(3), reverse(3)
        real(dp) :: tolerance, normal_scale, scale, cone_tolerance, condition

        embedded = .false.
        if (.not. all(ieee_is_finite(first_input)) &
            .or. .not. all(ieee_is_finite(second_input))) return
        scale = max(max_edge(first_input), max_edge(second_input))
        tolerance = 8192.0_dp * epsilon(1.0_dp)
        if (.not. ieee_is_finite(scale) .or. scale <= tiny(1.0_dp)) return
        do i = 1, 3
            first(:, i) = (first_input(:, i) - first_input(:, 1)) / scale
            second(:, i) = (second_input(:, i) - first_input(:, 1)) / scale
        end do
        if (.not. all(ieee_is_finite(first)) &
            .or. .not. all(ieee_is_finite(second))) return
        n1 = face_normal(first)
        n2 = face_normal(second)
        if (norm2(n1) <= tolerance * max_edge(first)**2) return
        if (norm2(n2) <= tolerance * max_edge(second)**2) return
        do i = 1, 3
            do j = i + 1, 3
                if (first_nodes(i) == first_nodes(j) &
                    .or. second_nodes(i) == second_nodes(j)) return
            end do
        end do
        count = 0
        do i = 1, 3
            do j = 1, 3
                if (first_nodes(i) /= second_nodes(j)) cycle
                if (norm2(first(:, i) - second(:, j)) > tolerance) return
                count = count + 1
                if (count > 3) return
                shared_first(count) = i
                shared_second(count) = j
            end do
        end do
        select case (count)
        case (0)
            embedded = .not. triangles_intersect(first, second)
        case (1)
            i = shared_first(1)
            j = shared_second(1)
            a = first(:, modulo(i, 3) + 1) - first(:, i)
            b = first(:, modulo(i + 1, 3) + 1) - first(:, i)
            c = second(:, modulo(j, 3) + 1) - second(:, j)
            d = second(:, modulo(j + 1, 3) + 1) - second(:, j)
            a = a / norm2(a)
            b = b / norm2(b)
            c = c / norm2(c)
            d = d / norm2(d)
            n1 = vector_cross(a, b)
            n2 = vector_cross(c, d)
            direction = vector_cross(n1, n2)
            reverse = -direction
            normal_scale = norm2(n1) * norm2(n2)
            condition = norm2(a) * norm2(b) / norm2(n1) &
                + norm2(c) * norm2(d) / norm2(n2)
            cone_tolerance = tolerance * condition
            if (norm2(direction) > cone_tolerance * normal_scale) then
                cone_tolerance = cone_tolerance * normal_scale / norm2(direction)
                embedded = .not. ((cone_contains(direction, a, b, cone_tolerance) &
                    .and. cone_contains(direction, c, d, cone_tolerance)) &
                    .or. (cone_contains(reverse, a, b, cone_tolerance) &
                    .and. cone_contains(reverse, c, d, cone_tolerance)))
            else
                embedded = .not. (cone_contains(a, c, d, cone_tolerance) &
                    .or. cone_contains(b, c, d, cone_tolerance) &
                    .or. cone_contains(c, a, b, cone_tolerance) &
                    .or. cone_contains(d, a, b, cone_tolerance))
            end if
        case (2)
            a = first(:, shared_first(2)) - first(:, shared_first(1))
            other = 6 - sum(shared_first(:2))
            b = first(:, other) - first(:, shared_first(1))
            other = 6 - sum(shared_second(:2))
            c = second(:, other) - second(:, shared_second(1))
            a = a / norm2(a)
            b = b / norm2(b)
            c = c / norm2(c)
            n1 = vector_cross(a, b)
            n2 = vector_cross(a, c)
            normal_scale = norm2(n1) * norm2(n2)
            condition = 1.0_dp / norm2(n1) + 1.0_dp / norm2(n2)
            cone_tolerance = tolerance * condition
            direction = vector_cross(n1, n2)
            embedded = norm2(direction) > cone_tolerance * normal_scale
            if (.not. embedded) embedded = dot_product(n1, n2) &
                < -cone_tolerance * normal_scale
        end select
    end function triangle_pair_is_embedded

    pure logical function cone_contains(direction, first, second, angular_tolerance) &
            result(contains)
        real(dp), intent(in) :: direction(3), first(3), second(3), angular_tolerance
        real(dp) :: a(3), b(3), d(3), normal(3), tolerance
        real(dp) :: first_product(3), second_product(3)

        a = first / norm2(first)
        b = second / norm2(second)
        d = direction / norm2(direction)
        normal = vector_cross(a, b)
        tolerance = 8192.0_dp * epsilon(1.0_dp) &
            + angular_tolerance * norm2(normal)
        ! The caller supplies an intersection-line direction, or conservatively
        ! projects a ray from a numerically coplanar face. Testing its plane
        ! residual here could hide a crossing after cancellation in n1 x n2.
        first_product = vector_cross(d, b)
        second_product = vector_cross(a, d)
        contains = dot_product(first_product, normal) >= -tolerance &
            .and. dot_product(second_product, normal) >= -tolerance
    end function cone_contains

    pure function vector_cross(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        call cross_product(first, second, product)
    end function vector_cross

    pure function face_normal(triangle) result(normal)
        real(dp), intent(in) :: triangle(3, 3)
        real(dp) :: normal(3), first_edge(3), second_edge(3)

        first_edge = triangle(:, 2) - triangle(:, 1)
        second_edge = triangle(:, 3) - triangle(:, 1)
        normal = vector_cross(first_edge, second_edge)
    end function face_normal

    ! Triangles whose axis-aligned bounding boxes, widened by the contact
    ! tolerance scale, are disjoint cannot intersect and skip the exact test.
    function meshes_intersect(first, second) result(intersect)
        real(dp), contiguous, intent(in) :: first(:, :, :), second(:, :, :)
        logical :: intersect
        real(dp) :: first_box(3, 2, size(first, 3))
        real(dp) :: second_box(3, 2, size(second, 3))
        integer :: i, j

        intersect = .false.
        call bounding_boxes(first, first_box)
        call bounding_boxes(second, second_box)
        do j = 1, size(second, 3)
            do i = 1, size(first, 3)
                if (any(first_box(:, 1, i) > second_box(:, 2, j)) &
                    .or. any(second_box(:, 1, j) > first_box(:, 2, i))) cycle
                if (triangles_intersect(first(:, :, i), second(:, :, j))) then
                    intersect = .true.
                    return
                end if
            end do
        end do
    end function meshes_intersect

    pure subroutine bounding_boxes(triangles, boxes)
        real(dp), intent(in) :: triangles(:, :, :)
        real(dp), intent(out) :: boxes(:, :, :)
        real(dp) :: edges, margin, vertices(3, 3)
        integer :: component, t

        do t = 1, size(triangles, 3)
            vertices = triangles(:, :, t)
            edges = max_edge(vertices)
            ! Wider than the tolerance of the exact test of any pair.
            margin = 8192.0_dp * epsilon(1.0_dp) &
                * max(maxval(abs(vertices)), edges)
            do component = 1, 3
                boxes(component, 1, t) = min(vertices(component, 1), &
                    vertices(component, 2), vertices(component, 3)) - margin
                boxes(component, 2, t) = max(vertices(component, 1), &
                    vertices(component, 2), vertices(component, 3)) + margin
            end do
        end do
    end subroutine bounding_boxes

    function point_inside_mesh(point, triangles) result(inside)
        real(dp), intent(in) :: point(3), triangles(:, :, :)
        logical :: inside
        real(dp) :: a(3), b(3), c(3), cross(3), denominator, solid_angle
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
            call cross_product(b, c, cross)
            solid_angle = solid_angle + 2.0_dp &
                * atan2(dot_product(a, cross), denominator)
        end do
        inside = abs(solid_angle) > 2.0_dp * pi
    end function point_inside_mesh

    ! A separating projection certifies disjointness. Face normals and edge
    ! cross products cover noncoplanar triangles; in-plane face-edge normals
    ! also cover the coplanar case. Uncertain contacts remain intersections.
    pure function triangles_intersect(first, second) result(intersect)
        real(dp), intent(in) :: first(3, 3), second(3, 3)
        logical :: intersect
        real(dp) :: a(3, 3), b(3, 3), first_edges(3, 3), second_edges(3, 3)
        real(dp) :: n1(3), n2(3), first_axis(3), second_axis(3), scale
        integer :: i, j, next

        intersect = .true.
        scale = max(max_edge(first), max_edge(second))
        if (.not. ieee_is_finite(scale) .or. scale <= tiny(1.0_dp)) return
        do i = 1, 3
            a(:, i) = (first(:, i) - first(:, 1)) / scale
            b(:, i) = (second(:, i) - first(:, 1)) / scale
        end do
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        do i = 1, 3
            next = modulo(i, 3) + 1
            first_edges(:, i) = a(:, next) - a(:, i)
            second_edges(:, i) = b(:, next) - b(:, i)
            if (norm2(first_edges(:, i)) <= tiny(1.0_dp) &
                .or. norm2(second_edges(:, i)) <= tiny(1.0_dp)) return
            first_edges(:, i) = first_edges(:, i) / norm2(first_edges(:, i))
            second_edges(:, i) = second_edges(:, i) / norm2(second_edges(:, i))
        end do
        n1 = vector_cross(first_edges(:, 1), first_edges(:, 2))
        n2 = vector_cross(second_edges(:, 1), second_edges(:, 2))
        if (axis_separates(n1, a, b) .or. axis_separates(n2, a, b)) then
            intersect = .false.
            return
        end if
        do i = 1, 3
            first_axis = vector_cross(n1, first_edges(:, i))
            second_axis = vector_cross(n2, second_edges(:, i))
            if (axis_separates(first_axis, a, b) &
                .or. axis_separates(second_axis, a, b)) then
                intersect = .false.
                return
            end if
            do j = 1, 3
                first_axis = vector_cross(first_edges(:, i), second_edges(:, j))
                if (axis_separates(first_axis, a, b)) then
                    intersect = .false.
                    return
                end if
            end do
        end do
    end function triangles_intersect

    pure logical function axis_separates(axis, first, second) result(separates)
        real(dp), intent(in) :: axis(3), first(3, 3), second(3, 3)
        real(dp) :: direction(3), a(3), b(3), length, tolerance
        integer :: i

        separates = .false.
        length = norm2(axis)
        if (length <= 8192.0_dp * epsilon(1.0_dp)) return
        direction = axis / length
        do i = 1, 3
            a(i) = dot_product(direction, first(:, i))
            b(i) = dot_product(direction, second(:, i))
        end do
        tolerance = 8192.0_dp * epsilon(1.0_dp) &
            * max(maxval(abs(a)), maxval(abs(b)), 1.0_dp)
        separates = minval(a) > maxval(b) + tolerance &
            .or. minval(b) > maxval(a) + tolerance
    end function axis_separates

    pure real(dp) function max_edge(triangle) result(value)
        real(dp), intent(in) :: triangle(3, 3)

        value = max(norm2(triangle(:, 2) - triangle(:, 1)), &
            norm2(triangle(:, 3) - triangle(:, 2)), &
            norm2(triangle(:, 1) - triangle(:, 3)))
    end function max_edge

    pure subroutine cross_product(first, second, value)
        real(dp), intent(in) :: first(3), second(3)
        real(dp), intent(out) :: value(3)

        value(1) = first(2) * second(3) - first(3) * second(2)
        value(2) = first(3) * second(1) - first(1) * second(3)
        value(3) = first(1) * second(2) - first(2) * second(1)
    end subroutine cross_product

end module boundary_mesh_geometry
