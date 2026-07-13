#include "gliss.h"

#include <stddef.h>
#include <string.h>

int main(void) {
    gliss_equilibrium *equilibrium = NULL;
    gliss_stability_problem *problem = NULL;
    gliss_spectrum_summary summary;
    gliss_energy_terms energy;
    size_t surfaces = 1;
    int32_t schema_version = -1;
    char error[128];

    memset(&summary, 0, sizeof(summary));
    summary.struct_size = sizeof(summary);
    memset(&energy, 0, sizeof(energy));
    energy.struct_size = sizeof(energy);

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
    return 0;
}
