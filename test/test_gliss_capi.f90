program test_gliss_capi
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
        c_int, c_loc, c_null_char, c_null_ptr, c_ptr, c_size_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: error_unit
    use cylinder_fixture, only: create_cylinder_fixture, fixture_ns
    implicit none

    integer(c_int), parameter :: status_ok = 0
    integer(c_int), parameter :: status_read_error = 1
    integer(c_int), parameter :: status_capacity = 3
    integer(c_int), parameter :: status_invalid_argument = 4
    character(len=*), parameter :: fixture = "capi_cylinder.nc"
    character(len=*), parameter :: output_file = "capi_roundtrip.nc"

    character(c_char), target :: path(len(fixture)), error_buffer(256)
    character(c_char), target :: output_path(len(output_file))
    character(c_char), target :: embedded_nul(3), missing_path(10)
    real(c_double), allocatable, target :: s_values(:), d_mercier(:)
    real(c_double), allocatable :: legacy_s(:), legacy_d(:)
    type(c_ptr), target :: context, roundtrip_context
    integer(c_size_t), target :: surfaces, written
    integer(c_int), target :: schema_version
    integer(c_int) :: status, legacy_status, legacy_surfaces

    interface
        function equilibrium_create(path_pointer, path_length, handle, &
                error_pointer, error_capacity) &
                bind(c, name="gliss_equilibrium_create") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: path_pointer
            integer(c_size_t), value :: path_length
            type(c_ptr), value :: handle
            type(c_ptr), value :: error_pointer
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function equilibrium_create

        function equilibrium_destroy(handle, error_pointer, error_capacity) &
                bind(c, name="gliss_equilibrium_destroy") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle
            type(c_ptr), value :: error_pointer
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function equilibrium_destroy

        function equilibrium_surface_count(handle, surfaces, error_pointer, &
                error_capacity) bind(c, name="gliss_equilibrium_surface_count") &
                result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle
            type(c_ptr), value :: surfaces
            type(c_ptr), value :: error_pointer
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function equilibrium_surface_count

        function equilibrium_schema_version(handle, version, error_pointer, &
                error_capacity) bind(c, &
                name="gliss_equilibrium_schema_version") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, version, error_pointer
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function equilibrium_schema_version

        function equilibrium_write(handle, path_pointer, path_length, &
                error_pointer, error_capacity) bind(c, &
                name="gliss_equilibrium_write") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle, path_pointer, error_pointer
            integer(c_size_t), value :: path_length, error_capacity
            integer(c_int) :: result
        end function equilibrium_write

        function mercier_profile_context(handle, n_theta, n_zeta, capacity, &
                s_pointer, d_pointer, written, error_pointer, error_capacity) &
                bind(c, name="gliss_mercier_profile_context") result(result)
            import c_int, c_ptr, c_size_t
            type(c_ptr), value :: handle
            integer(c_int), value :: n_theta, n_zeta
            integer(c_size_t), value :: capacity
            type(c_ptr), value :: s_pointer, d_pointer
            type(c_ptr), value :: written
            type(c_ptr), value :: error_pointer
            integer(c_size_t), value :: error_capacity
            integer(c_int) :: result
        end function mercier_profile_context

        subroutine mercier_profile_legacy(path_value, path_length, n_theta, &
                n_zeta, capacity, surfaces, s_values, d_mercier, status) &
                bind(c, name="gliss_mercier_profile")
            import c_char, c_double, c_int
            character(c_char), intent(in) :: path_value(*)
            integer(c_int), value, intent(in) :: path_length, n_theta, n_zeta
            integer(c_int), value, intent(in) :: capacity
            integer(c_int), intent(out) :: surfaces, status
            real(c_double), intent(out) :: s_values(*), d_mercier(*)
        end subroutine mercier_profile_legacy
    end interface

    call create_cylinder_fixture(fixture)
    call copy_chars(fixture, path)
    context = c_null_ptr
    status = equilibrium_create(c_loc(path), int(size(path), c_size_t), &
        c_loc(context), c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok, "valid equilibrium create failed")
    call require(c_associated(context), "create returned a null handle")
    call require(len_trim(error_text(error_buffer)) == 0, &
        "success left stale error text")
    schema_version = -1_c_int
    status = equilibrium_schema_version(context, c_loc(schema_version), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok .and. schema_version == 0_c_int, &
        "legacy schema version query failed")
    status = equilibrium_schema_version(context, c_null_ptr, &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null schema-version output was accepted")
    call copy_chars(output_file, output_path)
    schema_version = -1_c_int
    status = equilibrium_schema_version(c_null_ptr, c_loc(schema_version), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument .and. &
        schema_version == 0_c_int, "null schema-version handle was accepted")
    status = equilibrium_write(c_null_ptr, c_loc(output_path), &
        int(size(output_path), c_size_t), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null equilibrium write handle was accepted")
    status = equilibrium_write(context, c_null_ptr, 1_c_size_t, &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null equilibrium write path was accepted")
    status = equilibrium_write(context, c_loc(output_path), 0_c_size_t, &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "empty equilibrium write path was accepted")
    status = equilibrium_write(context, c_loc(output_path), &
        int(size(output_path), c_size_t), c_null_ptr, 1_c_size_t)
    call require(status == status_invalid_argument, &
        "null nonempty write error buffer was accepted")
    status = equilibrium_write(context, c_loc(output_path), &
        int(size(output_path), c_size_t), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "equilibrium write failed")
    status = equilibrium_write(context, c_loc(output_path), &
        int(size(output_path), c_size_t), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_read_error, &
        "equilibrium write overwrote an existing file")
    roundtrip_context = c_null_ptr
    status = equilibrium_create(c_loc(output_path), &
        int(size(output_path), c_size_t), c_loc(roundtrip_context), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok, "round-trip equilibrium create failed")
    status = equilibrium_schema_version(roundtrip_context, &
        c_loc(schema_version), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok .and. schema_version == 1_c_int, &
        "round-trip schema version is wrong")
    status = equilibrium_destroy(c_loc(roundtrip_context), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok, "round-trip destroy failed")

    status = equilibrium_surface_count(context, c_null_ptr, &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null surface-count output was not rejected")
    status = equilibrium_surface_count(context, c_loc(surfaces), &
        c_null_ptr, 1_c_size_t)
    call require(status == status_invalid_argument, &
        "null nonempty error buffer was not rejected")
    status = equilibrium_surface_count(context, c_loc(surfaces), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_ok, "surface count failed")
    call require(surfaces == fixture_ns, "surface count is wrong")

    allocate (s_values(surfaces), d_mercier(surfaces))
    status = mercier_profile_context(context, 32_c_int, 16_c_int, surfaces, &
        c_loc(s_values), c_loc(d_mercier), c_loc(written), &
        c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "context Mercier evaluation failed")
    call require(written == surfaces, "Mercier result size is wrong")
    call require(all(ieee_is_finite(s_values)), "s contains nonfinite values")
    call require(all(ieee_is_finite(d_mercier)), &
        "Mercier profile contains nonfinite values")
    call require(all(s_values(2:) > s_values(:size(s_values) - 1)), &
        "s is not strictly increasing")
    allocate (legacy_s(surfaces), legacy_d(surfaces))
    call mercier_profile_legacy(path, int(size(path), c_int), 32_c_int, &
        16_c_int, int(surfaces, c_int), legacy_surfaces, legacy_s, legacy_d, &
        legacy_status)
    call require(legacy_status == status_ok, "legacy Mercier call failed")
    call require(legacy_surfaces == int(surfaces, c_int), &
        "legacy surface count differs")
    call require(all(legacy_s == s_values), &
        "context changed the legacy radial coordinates")
    call require(all(legacy_d == d_mercier), &
        "context changed the legacy Mercier profile")

    status = mercier_profile_context(context, 32_c_int, 16_c_int, &
        surfaces - 1_c_size_t, c_loc(s_values), c_loc(d_mercier), &
        c_loc(written), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_capacity, "small capacity was not rejected")
    call require(written == surfaces, "required capacity was not returned")
    call require(index(error_text(error_buffer), "capacity") > 0, &
        "capacity error text is missing")

    status = mercier_profile_context(context, 32_c_int, 16_c_int, surfaces, &
        c_null_ptr, c_loc(d_mercier), c_loc(written), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null output was not rejected")
    status = mercier_profile_context(context, 0_c_int, 16_c_int, surfaces, &
        c_loc(s_values), c_loc(d_mercier), c_loc(written), &
        c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "zero angular resolution was not rejected")
    status = mercier_profile_context(context, 32_c_int, 16_c_int, surfaces, &
        c_loc(s_values), c_loc(d_mercier), c_null_ptr, c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null written output was not rejected")

    status = equilibrium_destroy(c_loc(context), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "destroy failed")
    call require(.not. c_associated(context), "destroy did not clear handle")
    status = equilibrium_destroy(c_loc(context), c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_ok, "repeated destroy failed")
    status = equilibrium_surface_count(context, c_loc(surfaces), &
        c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null handle was not rejected")
    status = equilibrium_destroy(c_null_ptr, c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null handle pointer was not rejected")
    status = equilibrium_create(c_loc(path), int(size(path), c_size_t), &
        c_null_ptr, c_loc(error_buffer), int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "null handle output was not rejected")

    call copy_chars("missing.nc", missing_path)
    status = equilibrium_create(c_loc(missing_path), &
        int(size(missing_path), c_size_t), c_loc(context), &
        c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_read_error, "missing file status is wrong")
    call require(.not. c_associated(context), &
        "failed create returned a live handle")
    call require(index(error_text(error_buffer), "read") > 0, &
        "read error text is missing")

    embedded_nul = [character(c_char) :: "a", c_null_char, "b"]
    status = equilibrium_create(c_loc(embedded_nul), &
        int(size(embedded_nul), c_size_t), c_loc(context), &
        c_loc(error_buffer), &
        int(size(error_buffer), c_size_t))
    call require(status == status_invalid_argument, &
        "embedded null was not rejected")
    status = equilibrium_create(c_null_ptr, 1_c_size_t, c_loc(context), &
        c_null_ptr, 0_c_size_t)
    call require(status == status_invalid_argument, &
        "null path without error buffer was not rejected")

    open (unit=13, file=fixture, status="old")
    close (13, status="delete")
    open (unit=13, file=output_file, status="old")
    close (13, status="delete")
    write (*, "(a)") "PASS"

contains

    subroutine copy_chars(source, destination)
        character(len=*), intent(in) :: source
        character(c_char), intent(out) :: destination(:)
        integer :: i

        call require(len(source) == size(destination), "character size mismatch")
        do i = 1, size(destination)
            destination(i) = source(i:i)
        end do
    end subroutine copy_chars

    function error_text(buffer) result(text)
        character(c_char), intent(in) :: buffer(:)
        character(len=size(buffer)) :: text
        integer :: i

        text = ""
        do i = 1, size(buffer)
            if (buffer(i) == c_null_char) exit
            text(i:i) = buffer(i)
        end do
    end function error_text

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_gliss_capi
