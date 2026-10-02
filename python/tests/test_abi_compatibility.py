"""Native byte-canary checks against frozen earlier C output layouts."""

import ctypes

import pytest

import gliss

pytestmark = pytest.mark.native


class LegacyEnergy(ctypes.Structure):
    # The ABI-3 layout before vacuum_energy was appended.
    _fields_ = [("struct_size", ctypes.c_size_t)] + [
        (name, ctypes.c_double) for name in (
            "field_line_bending", "magnetic_shear", "magnetic_compression",
            "pressure_drive", "plasma_compressibility", "potential_energy",
            "kinetic_energy", "rayleigh_quotient", "closure_error", "closure_tolerance",
        )
    ]


class LegacyMarginality(ctypes.Structure):
    # The ABI-3 layout before quotient_rank/labeled_nullity/block width.
    _fields_ = [
        ("struct_size", ctypes.c_size_t), ("has_eigenpair", ctypes.c_int32),
        ("field_periods", ctypes.c_int32), ("mode_count", ctypes.c_size_t),
        ("radial_surfaces", ctypes.c_size_t), ("parity_class", ctypes.c_int32),
        ("degree", ctypes.c_int32), ("angular_theta", ctypes.c_int32),
        ("angular_zeta", ctypes.c_int32), ("negative_count", ctypes.c_size_t),
        ("lowest_eigenvalue", ctypes.c_double), ("certificate", ctypes.c_double),
        ("eigenpair_residual", ctypes.c_double),
        ("force_balance_residual", ctypes.c_double), ("zero_floor", ctypes.c_double),
    ]


def _guarded(layout):
    class Guarded(ctypes.Structure):
        _fields_ = [("value", layout), ("tail", ctypes.c_ubyte * 64)]

    guarded = Guarded()
    guarded.value.struct_size = ctypes.sizeof(layout)
    for index in range(64):
        guarded.tail[index] = 0xA7
    return guarded


def test_native_size_qualified_legacy_outputs_preserve_canaries(
    native_library, test_data
):
    # A separate CDLL owns these explicit historic signatures, keeping the
    # production Python bindings' current pointer types untouched.
    library = ctypes.CDLL(native_library._name)
    energy = library.gliss_stability_problem_energy
    energy.argtypes = (
        ctypes.c_void_p, ctypes.c_int32, ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double), ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
    )
    energy.restype = ctypes.c_int
    marginality = library.gliss_cas3d_marginality
    marginality.argtypes = (
        ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32), ctypes.c_int32, ctypes.c_int32,
        ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_void_p,
        ctypes.c_void_p, ctypes.c_size_t,
    )
    marginality.restype = ctypes.c_int
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with gliss.StabilityProblem(
            equilibrium, [(1, 1)], degree=1, angular_theta=64,
            angular_zeta=8, radial_cells=2,
        ) as problem:
            vector = problem.solve_class(1).eigenvector
            current = problem.energy(1, vector)
            output = _guarded(LegacyEnergy)
            error = ctypes.create_string_buffer(512)
            status = energy(
                problem._handle, 1, vector.size,
                vector.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(output.value), error, len(error),
            )
            assert status == 0, error.value
            assert bytes(output.tail) == bytes([0xA7]) * 64
            assert output.value.rayleigh_quotient == current.rayleigh_quotient
        modes = (ctypes.c_int32 * 2)(1, 2)
        toroidal = (ctypes.c_int32 * 2)(1, 1)
        output = _guarded(LegacyMarginality)
        status = marginality(
            equilibrium._handle, 2, modes, toroidal, 1, 1, 64, 8, 0,
            ctypes.byref(output.value), error, len(error),
        )
        assert status == 0, error.value
        assert bytes(output.tail) == bytes([0xA7]) * 64
        reference = gliss.cas3d_marginality_inertia(
            equilibrium, [(1, 1), (2, 1)], degree=1, angular_theta=64, angular_zeta=8,
        )
        assert output.value.negative_count == reference.negative_count
        assert output.value.zero_floor == reference.inertia_zero_floor
        # The ABI-2 marginality layout ends before zero_floor. It is too
        # small and must fail without touching its payload or trailing bytes.
        class Abi2Marginality(ctypes.Structure):
            _fields_ = LegacyMarginality._fields_[:-1]

        rejected = _guarded(Abi2Marginality)
        before = bytes(rejected)
        status = marginality(
            equilibrium._handle, 2, modes, toroidal, 1, 1, 64, 8, 0,
            ctypes.byref(rejected.value), error, len(error),
        )
        assert status == 4 and bytes(rejected) == before
