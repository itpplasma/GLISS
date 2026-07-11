program test_terpsichore_normalization
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use terpsichore_normalization, only: &
        map_gliss_export_flip_pol_cell_to_terpsichore, &
        map_gliss_internal_flip_pol_cell_to_terpsichore, &
        terpsichore_normalization_invalid, terpsichore_normalization_ok
    implicit none

    integer, parameter :: field_periods = 3, points = 4
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: exported_phi_slope = -0.5_dp
    real(dp), parameter :: iota = 0.67_dp
    real(dp), parameter :: exported_chi_slope = iota * exported_phi_slope
    real(dp) :: signed_jacobian(points), signed_bjac(points)
    real(dp) :: ftp, fpp, full_volume_radians, full_volume_turns
    real(dp) :: radial_weight
    integer :: info

    signed_jacobian = [-28.0_dp, -30.0_dp, -29.0_dp, -31.0_dp]
    call map_gliss_export_flip_pol_cell_to_terpsichore(field_periods, 4, &
        0.25_dp, 0.5_dp, signed_jacobian, exported_phi_slope, &
        exported_chi_slope, signed_bjac, ftp, fpp, radial_weight, info)
    call require(info == terpsichore_normalization_ok, &
        "valid TERPSICHORE export normalization failed")
    call require(maxval(abs(signed_bjac &
        - field_periods * signed_jacobian / (4.0_dp * pi**2))) < 1.0e-14_dp, &
        "TERPSICHORE BJAC coordinate scale is wrong")
    call require(ftp == exported_phi_slope, &
        "TERPSICHORE toroidal export flux is wrong")
    call require(fpp == exported_chi_slope, &
        "TERPSICHORE poloidal export flux is wrong")
    call require(abs(fpp / ftp - iota) < 2.0e-16_dp, &
        "TERPSICHORE export map changed iota")
    call require(radial_weight == 1.0_dp, &
        "uniform TERPSICHORE radial weight is not unity")
    full_volume_turns = field_periods * sum(signed_jacobian) / points
    full_volume_radians = 4.0_dp * pi**2 * sum(signed_bjac) / points
    call require(abs(full_volume_turns - full_volume_radians) < 2.0e-14_dp, &
        "TERPSICHORE coordinate map changed the full-volume measure")
    call check_qas_values()
    call check_internal_map()
    call check_straight_cylinder()
    call check_nonuniform_weight()
    call check_invalid_inputs()

    write (*, "(a)") "PASS"

contains

    subroutine check_qas_values()
        real(dp), parameter :: qas_jacobian = -29.226284_dp
        real(dp), parameter :: qas_phi_slope = -0.5_dp
        real(dp), parameter :: qas_iota = 0.6739816378256_dp
        real(dp) :: local_bjac(1), local_ftp, local_fpp, local_weight
        integer :: local_info

        call map_gliss_export_flip_pol_cell_to_terpsichore(3, 64, 0.0_dp, &
            1.0_dp / 64.0_dp, [qas_jacobian], qas_phi_slope, &
            qas_iota * qas_phi_slope, local_bjac, local_ftp, local_fpp, &
            local_weight, local_info)
        call require(local_info == terpsichore_normalization_ok, &
            "QAS normalization failed")
        call require(abs(local_bjac(1) + 2.2210336_dp) < 2.0e-4_dp, &
            "QAS BJAC does not reproduce the TERPSICHORE scale")
        call require(local_ftp == qas_phi_slope &
            .and. abs(local_fpp / local_ftp - qas_iota) < 2.0e-16_dp, &
            "QAS flux slopes do not reproduce TERPSICHORE")
    end subroutine check_qas_values

    subroutine check_internal_map()
        real(dp) :: local_bjac(points), local_ftp, local_fpp, local_weight
        integer :: local_info

        call map_gliss_internal_flip_pol_cell_to_terpsichore(field_periods, &
            4, 0.25_dp, 0.5_dp, signed_jacobian, -exported_phi_slope, &
            -exported_chi_slope / field_periods, local_bjac, local_ftp, &
            local_fpp, local_weight, local_info)
        call require(local_info == terpsichore_normalization_ok, &
            "internal TERPSICHORE normalization failed")
        call require(maxval(abs(local_bjac - signed_bjac)) == 0.0_dp &
            .and. local_ftp == ftp .and. local_fpp == fpp &
            .and. local_weight == radial_weight, &
            "export and internal normalization maps disagree")
    end subroutine check_internal_map

    subroutine check_straight_cylinder()
        real(dp), parameter :: radius = 0.7_dp, length = 5.0_dp
        real(dp) :: cylinder_jacobian(1), cylinder_bjac(1)
        real(dp) :: cylinder_ftp, cylinder_fpp, cylinder_weight
        integer :: local_info

        cylinder_jacobian = -pi * radius**2 * length
        call map_gliss_export_flip_pol_cell_to_terpsichore(1, 1, 0.0_dp, &
            1.0_dp, cylinder_jacobian, -pi * radius**2 * 1.4_dp, 0.0_dp, &
            cylinder_bjac, cylinder_ftp, cylinder_fpp, cylinder_weight, &
            local_info)
        call require(local_info == terpsichore_normalization_ok, &
            "straight-cylinder normalization failed")
        call require(abs(4.0_dp * pi**2 * cylinder_bjac(1) &
            + pi * radius**2 * length) < 2.0e-14_dp, &
            "straight-cylinder volume is wrong")
        call require(cylinder_ftp == -pi * radius**2 * 1.4_dp &
            .and. cylinder_fpp == 0.0_dp .and. cylinder_weight == 1.0_dp, &
            "straight-cylinder flux or radial normalization is wrong")
    end subroutine check_straight_cylinder

    subroutine check_nonuniform_weight()
        real(dp) :: local_bjac(points), local_ftp, local_fpp, local_weight
        integer :: local_info

        call map_gliss_export_flip_pol_cell_to_terpsichore(field_periods, 8, &
            0.2_dp, 0.35_dp, signed_jacobian, exported_phi_slope, &
            exported_chi_slope, local_bjac, local_ftp, local_fpp, &
            local_weight, local_info)
        call require(local_info == terpsichore_normalization_ok, &
            "nonuniform TERPSICHORE normalization failed")
        call require(abs(local_weight - 1.2_dp) &
            < 8.0_dp * epsilon(1.0_dp), &
            "nonuniform TERPSICHORE radial weight is wrong")
    end subroutine check_nonuniform_weight

    subroutine check_invalid_inputs()
        real(dp) :: local_bjac(points), local_ftp, local_fpp, local_weight
        real(dp) :: wrong_orientation(points)
        integer :: local_info

        call map_gliss_export_flip_pol_cell_to_terpsichore(0, 4, 0.25_dp, &
            0.5_dp, signed_jacobian, exported_phi_slope, exported_chi_slope, &
            local_bjac, local_ftp, local_fpp, local_weight, local_info)
        call require(local_info == terpsichore_normalization_invalid, &
            "zero field periods were accepted")
        call map_gliss_export_flip_pol_cell_to_terpsichore(field_periods, 4, &
            0.5_dp, 0.25_dp, signed_jacobian, exported_phi_slope, &
            exported_chi_slope, local_bjac, local_ftp, local_fpp, &
            local_weight, local_info)
        call require(local_info == terpsichore_normalization_invalid, &
            "reversed radial cell was accepted")
        call map_gliss_export_flip_pol_cell_to_terpsichore(field_periods, 4, &
            0.25_dp, 0.5_dp, signed_jacobian, 0.0_dp, exported_chi_slope, &
            local_bjac, local_ftp, local_fpp, local_weight, local_info)
        call require(local_info == terpsichore_normalization_invalid, &
            "zero toroidal flux slope was accepted")
        wrong_orientation = -signed_jacobian
        call map_gliss_export_flip_pol_cell_to_terpsichore(field_periods, 4, &
            0.25_dp, 0.5_dp, wrong_orientation, exported_phi_slope, &
            exported_chi_slope, local_bjac, local_ftp, local_fpp, &
            local_weight, local_info)
        call require(local_info == terpsichore_normalization_invalid, &
            "right-handed Jacobian was accepted by the flip-pol map")
    end subroutine check_invalid_inputs

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_normalization
