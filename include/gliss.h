#ifndef GLISS_H
#define GLISS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GLISS_ABI_VERSION 1

typedef struct gliss_equilibrium gliss_equilibrium;
typedef struct gliss_stability_problem gliss_stability_problem;

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

typedef struct gliss_spectrum_summary {
    size_t struct_size;
    int32_t has_chart_metric;
    int32_t has_eigenvector;
    int32_t field_periods;
    int32_t parity_class;
    int32_t radial_quadrature;
    int32_t angular_theta;
    int32_t angular_zeta;
    size_t mode_count;
    size_t unknowns;
    size_t normal_unknowns;
    size_t eta_unknowns;
    size_t mu_unknowns;
    size_t negative_count;
    size_t floor_count;
    double adiabatic_index;
    double density_kg_m3;
    double zero_floor;
    double lowest_eigenvalue;
    double certificate;
    double eigenpair_residual;
    double eigenpair_resolution;
    double inertia_interval;
} gliss_spectrum_summary;

/* A fixed-boundary problem copies and assembles all data it needs, so the
 * equilibrium may be destroyed after successful construction. mode_m and
 * mode_n are mode_count contiguous int32_t values. radial_quadrature is 1 for
 * midpoint or 2 for two-point Gauss quadrature. */
gliss_status gliss_stability_problem_create(
    const gliss_equilibrium *equilibrium,
    double adiabatic_index,
    double density_kg_m3,
    double zero_floor,
    size_t mode_count,
    const int32_t *mode_m,
    const int32_t *mode_n,
    int32_t radial_quadrature,
    gliss_stability_problem **problem,
    char *error,
    size_t error_capacity);

gliss_status gliss_stability_problem_destroy(
    gliss_stability_problem **problem,
    char *error,
    size_t error_capacity);

gliss_status gliss_stability_problem_unknown_count(
    const gliss_stability_problem *problem,
    int32_t parity_class,
    size_t *unknown_count,
    char *error,
    size_t error_capacity);

/* Set summary->struct_size to sizeof(*summary) before calling. eigenvector is
 * caller-owned and uses the documented dynamic component order: fixed-edge
 * normal unknowns, then eta, then mu. It is normalized by x^T M x = 1. On a
 * capacity error, written reports the required count and eigenvector is not
 * modified. If has_eigenvector is zero, written is zero and eigenvector may be
 * NULL. */
gliss_status gliss_stability_problem_solve_class(
    const gliss_stability_problem *problem,
    int32_t parity_class,
    size_t capacity,
    double *eigenvector,
    size_t *written,
    gliss_spectrum_summary *summary,
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
