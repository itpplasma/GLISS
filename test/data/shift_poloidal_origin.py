"""Write a copy of a GLISS equilibrium with the poloidal angle origin moved.

usage: python shift_poloidal_origin.py INPUT OUTPUT SHIFT

f(theta + SHIFT) for every stored harmonic c cos(phi) + s sin(phi) with the
phase 2 pi (m theta - n zeta / N_T): c' = c cos(d) + s sin(d) and
s' = s cos(d) - c sin(d), d = 2 pi m SHIFT. The equilibrium is physically
unchanged but stores both parities of every field, so the file is written with
stellarator_symmetry = "False". test/data/solovev_q1.035_shifted.nc is
solovev_q1.035.nc shifted by 3/32 (grid aligned for 32 poloidal points).
"""

import sys

import netCDF4
import numpy as np


def main(source: str, target: str, shift: float) -> None:
    with netCDF4.Dataset(source) as data, netCDF4.Dataset(
        target, "w", format=data.data_model
    ) as out:
        out.setncatts({name: data.getncattr(name) for name in data.ncattrs()})
        out.stellarator_symmetry = "False"
        for name, dimension in data.dimensions.items():
            out.createDimension(name, len(dimension))
        angle = 2.0 * np.pi * np.asarray(data["m"][:], dtype=float) * shift
        fields = sorted(
            {name[:-4] for name in data.variables if name[-4:] in ("_mnc", "_mns")}
        )
        for name, variable in data.variables.items():
            if name[-4:] in ("_mnc", "_mns"):
                continue
            copy = out.createVariable(name, variable.dtype, variable.dimensions)
            copy.setncatts({key: variable.getncattr(key) for key in variable.ncattrs()})
            copy[...] = variable[...]
        for field in fields:
            template = data[f"{field}_mnc"] if f"{field}_mnc" in data.variables else (
                data[f"{field}_mns"]
            )
            zero = np.zeros(template.shape)
            cosine = np.asarray(data[f"{field}_mnc"][:]) if f"{field}_mnc" in data.variables else zero
            sine = np.asarray(data[f"{field}_mns"][:]) if f"{field}_mns" in data.variables else zero
            factor_c = np.cos(angle)[None, :, None]
            factor_s = np.sin(angle)[None, :, None]
            rotated = {
                "mnc": cosine * factor_c + sine * factor_s,
                "mns": sine * factor_c - cosine * factor_s,
            }
            for suffix, values in rotated.items():
                copy = out.createVariable(f"{field}_{suffix}", "d", template.dimensions)
                copy.setncatts({key: template.getncattr(key) for key in template.ncattrs()})
                copy[...] = values


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(eval(sys.argv[3], {}, {})))
