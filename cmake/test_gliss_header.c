#include "gliss.h"

#include <stddef.h>
#include <string.h>

int main(void) {
    gliss_equilibrium *equilibrium = NULL;
    gliss_stability_problem *problem = NULL;
    gliss_spectrum_summary summary;
    size_t surfaces = 1;
    char error[128];

    memset(&summary, 0, sizeof(summary));
    summary.struct_size = sizeof(summary);

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
    if (gliss_stability_problem_destroy(&problem, error, sizeof(error)) !=
        GLISS_STATUS_OK) {
        return 5;
    }
    if (summary.struct_size != sizeof(gliss_spectrum_summary)) {
        return 6;
    }
    return 0;
}
