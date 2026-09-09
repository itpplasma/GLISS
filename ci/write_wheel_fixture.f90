program write_wheel_fixture
    use cylinder_fixture, only: create_cylinder_fixture
    implicit none

    character(len=4096) :: path
    integer :: status

    if (command_argument_count() /= 1) error stop "expected output path"
    call get_command_argument(1, path, status=status)
    if (status /= 0) error stop "invalid output path"
    call create_cylinder_fixture(trim(path))
end program write_wheel_fixture
