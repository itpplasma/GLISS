#include "gliss.h"

#include <stddef.h>

int main(void) {
    gliss_equilibrium *equilibrium = NULL;
    size_t surfaces = 1;
    char error[128];

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
    return 0;
}
