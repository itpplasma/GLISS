module terpsichore_model_policy
    implicit none
    private

    integer, parameter, public :: terpsichore_model_ok = 0
    integer, parameter, public :: terpsichore_model_invalid = -1

    integer, parameter, public :: potential_model_non_interacting = 0
    integer, parameter, public :: potential_model_kruskal_oberman = 1
    integer, parameter, public :: kinetic_norm_reduced = 0
    integer, parameter, public :: kinetic_norm_physical = 1

    type, public :: terpsichore_model_config_t
        integer :: potential_model = -1
        integer :: kinetic_norm = -1
    end type terpsichore_model_config_t

    public :: decode_terpsichore_model

contains

    pure subroutine decode_terpsichore_model(modelk, config, info)
        integer, intent(in) :: modelk
        type(terpsichore_model_config_t), intent(out) :: config
        integer, intent(out) :: info

        config = terpsichore_model_config_t()
        info = terpsichore_model_invalid
        select case (modelk)
        case (0)
            config%potential_model = potential_model_non_interacting
            config%kinetic_norm = kinetic_norm_reduced
        case (1)
            config%potential_model = potential_model_kruskal_oberman
            config%kinetic_norm = kinetic_norm_reduced
        case (2)
            config%potential_model = potential_model_non_interacting
            config%kinetic_norm = kinetic_norm_physical
        case (3)
            config%potential_model = potential_model_kruskal_oberman
            config%kinetic_norm = kinetic_norm_physical
        case default
            return
        end select
        info = terpsichore_model_ok
    end subroutine decode_terpsichore_model

end module terpsichore_model_policy
