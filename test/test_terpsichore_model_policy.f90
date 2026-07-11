program test_terpsichore_model_policy
    use, intrinsic :: iso_fortran_env, only: error_unit
    use terpsichore_model_policy, only: decode_terpsichore_model, &
        kinetic_norm_physical, kinetic_norm_reduced, &
        potential_model_kruskal_oberman, potential_model_non_interacting, &
        terpsichore_model_config_t, terpsichore_model_invalid, &
        terpsichore_model_ok
    implicit none

    type(terpsichore_model_config_t) :: config
    integer :: info

    call check_model(0, potential_model_non_interacting, kinetic_norm_reduced)
    call check_model(1, potential_model_kruskal_oberman, kinetic_norm_reduced)
    call check_model(2, potential_model_non_interacting, kinetic_norm_physical)
    call check_model(3, potential_model_kruskal_oberman, kinetic_norm_physical)
    call decode_terpsichore_model(-1, config, info)
    call require(info == terpsichore_model_invalid, &
        "negative MODELK was accepted")
    call decode_terpsichore_model(4, config, info)
    call require(info == terpsichore_model_invalid, &
        "out-of-range MODELK was accepted")

    write (*, "(a)") "PASS"

contains

    subroutine check_model(modelk, potential_model, kinetic_norm)
        integer, intent(in) :: modelk, potential_model, kinetic_norm

        call decode_terpsichore_model(modelk, config, info)
        call require(info == terpsichore_model_ok, "valid MODELK was rejected")
        call require(config%potential_model == potential_model, &
            "MODELK potential model is wrong")
        call require(config%kinetic_norm == kinetic_norm, &
            "MODELK kinetic norm is wrong")
    end subroutine check_model

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_model_policy
