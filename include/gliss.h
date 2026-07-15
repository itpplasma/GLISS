#ifndef GLISS_H
#define GLISS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GLISS_ABI_VERSION 2

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

/* Legacy exports report schema version 0. GLISS exports report version 1.
 * Direct writes use exclusive creation and never replace an existing path. */
gliss_status gliss_equilibrium_schema_version(
    const gliss_equilibrium *equilibrium,
    int32_t *schema_version,
    char *error,
    size_t error_capacity);

gliss_status gliss_equilibrium_write(
    const gliss_equilibrium *equilibrium,
    const char *path,
    size_t path_length,
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
    int32_t degree;
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

typedef struct gliss_energy_terms {
    size_t struct_size;
    double field_line_bending;
    double magnetic_shear;
    double magnetic_compression;
    double pressure_drive;
    double plasma_compressibility;
    double potential_energy;
    double kinetic_energy;
    double rayleigh_quotient;
    double closure_error;
    double closure_tolerance;
} gliss_energy_terms;

typedef struct gliss_solver_tolerances {
    size_t struct_size;
    double eigenvalue_relative;
    double residual_relative;
    double negative_bracket_relative;
    double negative_bracket_floor;
    int32_t inverse_iteration_limit;
    int32_t bracket_iteration_limit;
} gliss_solver_tolerances;

typedef struct gliss_terpsichore_fixed_boundary_result {
    size_t struct_size;
    size_t unknowns;
    size_t negative_count;
    double eigenvalue;
    double certificate;
    double residual;
    double resolution;
} gliss_terpsichore_fixed_boundary_result;

typedef struct gliss_terpsichore_pseudoplasma_result {
    size_t struct_size;
    size_t unknowns;
    size_t negative_count;
    double eigenvalue;
    double certificate;
    double residual;
    double resolution;
    double growth_rate;
    double reference_eigenvalue;
    double reference_potential;
    double computed_potential;
    double reference_kinetic;
    double computed_kinetic;
    double reference_residual;
    double mode_overlap;
} gliss_terpsichore_pseudoplasma_result;

typedef struct gliss_axisymmetric_spectrum_result {
    size_t struct_size;
    int32_t has_eigenpair;
    int32_t field_periods;
    int32_t toroidal_mode;
    int32_t poloidal_max;
    size_t mode_count;
    size_t radial_surfaces;
    int32_t parity_class;
    int32_t degree;
    size_t negative_count;
    double lowest_eigenvalue;
    double certificate;
    double eigenpair_residual;
    double force_balance_residual;
} gliss_axisymmetric_spectrum_result;

typedef struct gliss_cas3d_marginality_result {
    size_t struct_size;
    int32_t has_eigenpair;
    int32_t field_periods;
    size_t mode_count;
    size_t radial_surfaces;
    int32_t parity_class;
    int32_t degree;
    int32_t angular_theta;
    int32_t angular_zeta;
    size_t negative_count;
    double lowest_eigenvalue;
    double certificate;
    double eigenpair_residual;
    double force_balance_residual;
} gliss_cas3d_marginality_result;

/* Solve the lowest negative eigenpair represented by a TERPSICHORE FORT.23
 * file produced with IVAC=0 and MODELK=0. Set result->struct_size to
 * sizeof(*result). The file's Fourier table, reduced kinetic normalization,
 * and fixed-boundary radial topology are used without reinterpretation. The
 * result is not modified on failure. */
gliss_status gliss_terpsichore_fixed_boundary(
    const char *path,
    size_t path_length,
    gliss_terpsichore_fixed_boundary_result *result,
    char *error,
    size_t error_capacity);

/* Solve the lowest negative eigenpair represented by TERPSICHORE FORT.23 and
 * FORT.24 files produced with IVAC>0 and MODELK=0. vacuum_intervals must match
 * both files. The pressureless pseudo-plasma vacuum is eliminated by its Schur
 * complement. The reference fields report diagnostics of the TERPSICHORE
 * solution stored in FORT.23. Set result->struct_size to sizeof(*result). The
 * result is not modified on failure. */
gliss_status gliss_terpsichore_pseudoplasma(
    const char *matrix_path,
    size_t matrix_path_length,
    int32_t vacuum_intervals,
    const char *vacuum_path,
    size_t vacuum_path_length,
    gliss_terpsichore_pseudoplasma_result *result,
    char *error,
    size_t error_capacity);

/* Evaluate the fixed-boundary, sine-parity axisymmetric family on an existing
 * equilibrium. The native mode table is (0,+n), then (m,-n),(m,+n) through
 * poloidal_max, with the regular-axis powers used by gliss_axisymmetric.
 * degree selects the compatible radial FEEC degree from 1 through 4.
 * solve_eigenpair is 0 for inertia only or 1 for the certified lowest pair.
 * Set result->struct_size to sizeof(*result). The result is unchanged on
 * failure. */
gliss_status gliss_axisymmetric_spectrum(
    const gliss_equilibrium *equilibrium,
    int32_t toroidal_mode,
    int32_t poloidal_max,
    int32_t degree,
    int32_t solve_eigenpair,
    gliss_axisymmetric_spectrum_result *result,
    char *error,
    size_t error_capacity);

/* Evaluate the compatible two-component incompressible functional on an
 * explicit 3-D mode table. Its perpendicular normalization preserves the
 * inertia and marginal boundary but does not define a physical growth rate.
 * The regular-axis factor s^(m/2) is derived from each nonnegative poloidal
 * mode. degree must be between 1 and 4. parity_class must be 1 or 2.
 * solve_eigenpair is 0 for inertia only or 1 for the certified lowest pair.
 * Set result->struct_size to sizeof(*result). The result is unchanged on
 * failure. */
gliss_status gliss_cas3d_marginality(
    const gliss_equilibrium *equilibrium,
    size_t mode_count,
    const int32_t *mode_m,
    const int32_t *mode_n,
    int32_t parity_class,
    int32_t degree,
    int32_t angular_theta,
    int32_t angular_zeta,
    int32_t solve_eigenpair,
    gliss_cas3d_marginality_result *result,
    char *error,
    size_t error_capacity);

/* Evaluate the CAS3D2MN phase-envelope representation. base_m and base_n
 * use the GLISS phase 2*pi*(m*theta - n*zeta/N_T). Envelope modes use
 * 2*pi*(m*theta - n*zeta) on one field period and must begin with (0,0).
 * Each later envelope mode expands to two labeled physical sidebands.
 * Coincident sidebands are assembled once in the physical FEEC space while
 * result mode_count remains the labeled count 2*envelope_count-1.
 * Other arguments and the result contract match gliss_cas3d_marginality. */
gliss_status gliss_cas3d_phase_envelope(
    const gliss_equilibrium *equilibrium,
    int32_t base_m,
    int32_t base_n,
    size_t envelope_count,
    const int32_t *envelope_m,
    const int32_t *envelope_n,
    int32_t parity_class,
    int32_t degree,
    int32_t angular_theta,
    int32_t angular_zeta,
    int32_t solve_eigenpair,
    gliss_cas3d_marginality_result *result,
    char *error,
    size_t error_capacity);

/* A fixed-boundary problem copies and assembles all data it needs, so the
 * equilibrium may be destroyed after successful construction. mode_m and
 * mode_n are mode_count contiguous int32_t values. degree selects the
 * compatible radial FEEC degree from 1 through 4. */
gliss_status gliss_stability_problem_create(
    const gliss_equilibrium *equilibrium,
    double adiabatic_index,
    double density_kg_m3,
    double zero_floor,
    size_t mode_count,
    const int32_t *mode_m,
    const int32_t *mode_n,
    int32_t degree,
    gliss_stability_problem **problem,
    char *error,
    size_t error_capacity);

gliss_status gliss_stability_problem_destroy(
    gliss_stability_problem **problem,
    char *error,
    size_t error_capacity);

/* Replace only the stopping controls of an assembled problem. The historical
 * defaults are 1e-13, 1e-12, 1e-9, 1e-3, 500 and 200 in field order. Set
 * tolerances->struct_size to sizeof(*tolerances). Matrices, discretization,
 * floor-band classification and normalization are unchanged. */
gliss_status gliss_stability_problem_set_solver_tolerances(
    gliss_stability_problem *problem,
    const gliss_solver_tolerances *tolerances,
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

/* Evaluate the five physical contributions to x^T K x, the independent total
 * x^T K x, x^T M x and their quotient for one caller-owned vector in dynamic
 * component order. Set terms->struct_size to sizeof(*terms). vector_count must
 * equal gliss_stability_problem_unknown_count for parity_class. */
gliss_status gliss_stability_problem_energy(
    const gliss_stability_problem *problem,
    int32_t parity_class,
    size_t vector_count,
    const double *vector,
    gliss_energy_terms *terms,
    char *error,
    size_t error_capacity);

/* Apply the reverse derivative of x^T K x / x^T M x with respect to the
 * displacement coefficients in dynamic component order. The primal vector
 * must have positive kinetic norm. gradient_capacity must equal the unknown
 * count. The output is cotangent times 2 (K x - q M x) / (x^T M x), where
 * q is the Rayleigh quotient. No output is modified on failure. */
gliss_status gliss_stability_problem_rayleigh_vjp(
    const gliss_stability_problem *problem,
    int32_t parity_class,
    size_t vector_count,
    const double *vector,
    double cotangent,
    size_t gradient_capacity,
    double *gradient,
    char *error,
    size_t error_capacity);

/* Compute every generalized eigenpair with a dense LAPACK solve. This is an
 * explicit O(unknowns^3)-time, O(unknowns^2)-memory operation; use
 * gliss_stability_problem_solve_class when only the certified lowest pair is
 * required. Eigenvectors contains unknowns contiguous vectors in ascending
 * eigenvalue order, each in dynamic component order and normalized by
 * x^T M x = 1. rayleigh_quotients independently reevaluates x^T K x / x^T M x;
 * residuals and resolutions report the scaled backward error and roundoff
 * resolution for each pair. On a capacity error, both written outputs report
 * the required counts and no data array is modified. eigenvalue_capacity
 * applies to eigenvalues, residuals, resolutions, and rayleigh_quotients.
 * Output arrays must not overlap. */
gliss_status gliss_stability_problem_full_spectrum(
    const gliss_stability_problem *problem,
    int32_t parity_class,
    size_t eigenvalue_capacity,
    double *eigenvalues,
    double *residuals,
    double *resolutions,
    double *rayleigh_quotients,
    size_t eigenvector_capacity,
    double *eigenvectors,
    size_t *eigenvalues_written,
    size_t *eigenvectors_written,
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
