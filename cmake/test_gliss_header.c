#include "gliss.h"

#include <stddef.h>
#include <string.h>

int main(void) {
    const char *missing_fort23 = "__gliss_missing_fort23__";
    gliss_equilibrium *equilibrium = NULL;
    gliss_stability_problem *problem = NULL;
    gliss_spectrum_summary summary;
    gliss_energy_terms energy;
    gliss_terpsichore_fixed_boundary_result terpsichore;
    gliss_terpsichore_pseudoplasma_result pseudoplasma;
    gliss_axisymmetric_spectrum_result axisymmetric;
    gliss_cas3d_marginality_result marginality;
    int32_t mode_m[1] = {1};
    int32_t mode_n[1] = {1};
    int32_t envelope_m[1] = {0};
    int32_t envelope_n[1] = {0};
    size_t surfaces = 1;
    int32_t schema_version = -1;
    char error[128];

    memset(&summary, 0, sizeof(summary));
    summary.struct_size = sizeof(summary);
    memset(&energy, 0, sizeof(energy));
    energy.struct_size = sizeof(energy);
    memset(&terpsichore, 0, sizeof(terpsichore));
    memset(&pseudoplasma, 0, sizeof(pseudoplasma));
    memset(&axisymmetric, 0, sizeof(axisymmetric));
    axisymmetric.struct_size = sizeof(axisymmetric);
    axisymmetric.mode_count = 123;
    memset(&marginality, 0, sizeof(marginality));
    marginality.struct_size = sizeof(marginality);
    marginality.mode_count = 456;

    if (gliss_abi_version() != GLISS_ABI_VERSION) {
        return 1;
    }
    if (gliss_equilibrium_destroy(&equilibrium, error, sizeof(error)) !=
        GLISS_STATUS_OK) {
        return 2;
    }
    if (gliss_equilibrium_surface_count(
            equilibrium, &surfaces, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 3;
    }
    if (surfaces != 0 || error[0] == '\0') {
        return 4;
    }
    if (gliss_equilibrium_schema_version(
            equilibrium, &schema_version, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 5;
    }
    if (schema_version != 0 || error[0] == '\0') {
        return 6;
    }
    if (gliss_stability_problem_destroy(&problem, error, sizeof(error)) !=
        GLISS_STATUS_OK) {
        return 7;
    }
    if (summary.struct_size != sizeof(gliss_spectrum_summary)) {
        return 8;
    }
    if (gliss_stability_problem_full_spectrum(
            problem, 1, 0, NULL, NULL, NULL, NULL, 0, NULL, &surfaces,
            &surfaces, error, sizeof(error)) != GLISS_STATUS_INVALID_ARGUMENT) {
        return 9;
    }
    if (gliss_stability_problem_energy(
            problem, 1, 0, NULL, &energy, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 10;
    }
    if (energy.struct_size != sizeof(gliss_energy_terms)) {
        return 11;
    }
    if (gliss_terpsichore_fixed_boundary(
            missing_fort23, strlen(missing_fort23), &terpsichore, error,
            sizeof(error)) != GLISS_STATUS_INVALID_ARGUMENT) {
        return 12;
    }
    if (error[0] == '\0') {
        return 13;
    }
    terpsichore.struct_size = sizeof(terpsichore);
    terpsichore.unknowns = 123;
    if (gliss_terpsichore_fixed_boundary(
            missing_fort23, strlen(missing_fort23), &terpsichore, error,
            sizeof(error)) != GLISS_STATUS_READ_ERROR) {
        return 14;
    }
    if (terpsichore.unknowns != 123 || error[0] == '\0') {
        return 15;
    }
    if (gliss_terpsichore_fixed_boundary(
            missing_fort23, strlen(missing_fort23), NULL, error,
            sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 16;
    }
    if (gliss_terpsichore_pseudoplasma(
            missing_fort23, strlen(missing_fort23), 1, missing_fort23,
            strlen(missing_fort23), &pseudoplasma, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 30;
    }
    pseudoplasma.struct_size = sizeof(pseudoplasma);
    pseudoplasma.unknowns = 123;
    if (gliss_terpsichore_pseudoplasma(
            missing_fort23, strlen(missing_fort23), 0, missing_fort23,
            strlen(missing_fort23), &pseudoplasma, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 31;
    }
    if (pseudoplasma.unknowns != 123 || error[0] == '\0') {
        return 32;
    }
    if (gliss_terpsichore_pseudoplasma(
            missing_fort23, strlen(missing_fort23), 1, missing_fort23,
            strlen(missing_fort23), &pseudoplasma, error, sizeof(error)) !=
        GLISS_STATUS_READ_ERROR) {
        return 33;
    }
    if (pseudoplasma.unknowns != 123 || error[0] == '\0') {
        return 34;
    }
    if (gliss_terpsichore_pseudoplasma(
            missing_fort23, strlen(missing_fort23), 1, missing_fort23,
            strlen(missing_fort23), NULL, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 35;
    }
    if (gliss_axisymmetric_spectrum(
            equilibrium, 1, 8, 1, 1, &axisymmetric, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 17;
    }
    if (axisymmetric.mode_count != 123 || error[0] == '\0') {
        return 18;
    }
    if (gliss_axisymmetric_spectrum(
            equilibrium, 1, 8, 1, 1, NULL, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 19;
    }
    if (gliss_cas3d_marginality(
            equilibrium, 1, mode_m, mode_n, 1, 1, 8, 8, 1, &marginality,
            error, sizeof(error)) != GLISS_STATUS_INVALID_ARGUMENT) {
        return 20;
    }
    if (marginality.mode_count != 456 || error[0] == '\0') {
        return 21;
    }
    if (gliss_cas3d_marginality(
            equilibrium, 1, mode_m, mode_n, 1, 1, 8, 8, 1, NULL, error,
            sizeof(error)) != GLISS_STATUS_INVALID_ARGUMENT) {
        return 22;
    }
    if (gliss_cas3d_phase_envelope(
            equilibrium, 3, 2, 1, envelope_m, envelope_n, 1, 1, 8, 8, 1,
            &marginality, error, sizeof(error)) !=
        GLISS_STATUS_INVALID_ARGUMENT) {
        return 23;
    }
    if (marginality.mode_count != 456 || error[0] == '\0') {
        return 24;
    }
    if (gliss_cas3d_phase_envelope(
            equilibrium, 3, 2, 1, envelope_m, envelope_n, 1, 1, 8, 8, 1,
            NULL, error, sizeof(error)) != GLISS_STATUS_INVALID_ARGUMENT) {
        return 25;
    }
    return 0;
}
