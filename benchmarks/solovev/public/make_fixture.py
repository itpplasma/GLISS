#!/usr/bin/env python3
"""Shrink a patched ``pygvec to-cas3d --stellsym`` export into a test fixture.

The stellarator-symmetric export keeps only the populated parity of every
harmonic field; the fixture drops run-specific attributes and is written in
the compact classic 64-bit NetCDF format.
"""
import argparse

import xarray as xr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export")
    parser.add_argument("fixture")
    args = parser.parse_args()
    data = xr.open_dataset(args.export).load()
    for name in ("statefile", "conversion_time"):
        data.attrs.pop(name, None)
    data.to_netcdf(args.fixture, format="NETCDF3_64BIT")


if __name__ == "__main__":
    main()
