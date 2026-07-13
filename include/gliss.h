#ifndef GLISS_H
#define GLISS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GLISS_ABI_VERSION 1

typedef struct gliss_equilibrium gliss_equilibrium;

/* Status values and existing function signatures do not change within one
 * ABI version. Functions that accept an error buffer clear it on success and
 * write a truncated, null-terminated message on failure. The error pointer may
 * be NULL only when error_capacity is zero. */
typedef enum gliss_status {
    GLISS_STATUS_OK = 0,
    GLISS_STATUS_READ_ERROR = 1,
    GLISS_STATUS_COMPUTE_ERROR = 2,
    GLISS_STATUS_CAPACITY = 3,
    GLISS_STATUS_INVALID_ARGUMENT = 4,
    GLISS_STATUS_ALLOCATION_ERROR = 5,
    GLISS_STATUS_INTERNAL_ERROR = 6
} gliss_status;

void gliss_version(char *buffer, int32_t length);
int32_t gliss_abi_version(void);

/* Contexts own all native allocations. They may coexist, but calls using the
 * same context must not overlap. Concurrent creation also requires a
 * thread-safe NetCDF C library. No allocation returned by GLISS is freed by
 * the caller. */
gliss_status gliss_equilibrium_create(
    const char *path,
    size_t path_length,
    gliss_equilibrium **equilibrium,
    char *error,
    size_t error_capacity);

gliss_status gliss_equilibrium_destroy(
    gliss_equilibrium **equilibrium,
    char *error,
    size_t error_capacity);

gliss_status gliss_equilibrium_surface_count(
    const gliss_equilibrium *equilibrium,
    size_t *surface_count,
    char *error,
    size_t error_capacity);

/* s_values and d_mercier are caller-owned contiguous arrays of capacity
 * doubles. On GLISS_STATUS_CAPACITY, written reports the required capacity and
 * neither array is modified. On earlier failures, written is zero. */
gliss_status gliss_mercier_profile_context(
    const gliss_equilibrium *equilibrium,
    int32_t n_theta,
    int32_t n_zeta,
    size_t capacity,
    double *s_values,
    double *d_mercier,
    size_t *written,
    char *error,
    size_t error_capacity);

/* Retained for ABI-v1 compatibility. New code should reuse an equilibrium
 * context instead of loading the same file for every diagnostic. */
void gliss_mercier_profile(
    const char *path,
    int32_t path_length,
    int32_t n_theta,
    int32_t n_zeta,
    int32_t capacity,
    int32_t *surfaces,
    double *s_values,
    double *d_mercier,
    int32_t *status);

#ifdef __cplusplus
}
#endif

#endif
