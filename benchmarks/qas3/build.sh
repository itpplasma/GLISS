#!/bin/bash
# Download and build the public reference codes for the QAS3 benchmark.
#
#   STELLOPT VMEC2000 (PrincetonUniversity/STELLOPT, pinned revision)
#   TERPSICHORE 1.2 (public tag; GitHub mirror of gitlab.epfl.ch)
#   booz_xform (PyPI) for gliss.convert_vmec
#
# Debian/Ubuntu packages: gfortran make libopenblas-dev libscalapack-openmpi-dev
# libopenmpi-dev libnetcdf-dev libnetcdff-dev.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HERE/work}
JOBS=${JOBS:-4}
STELLOPT_REV=2f181f0d
TERPSICHORE_URL=${TERPSICHORE_URL:-https://github.com/krystophny/terpsichore.git}
TERPSICHORE_REV=04dcf9a
mkdir -p "$WORK"
cd "$WORK"

if [ ! -x STELLOPT/VMEC2000/Release/xvmec2000 ]; then
    [ -d STELLOPT ] || git clone -q https://github.com/PrincetonUniversity/STELLOPT.git
    git -C STELLOPT checkout -q "$STELLOPT_REV"
    # make_debian.inc without the optional HDF5, DKES and NEO dependencies.
    sed -e 's/LHDF5 = T/LHDF5 = F/' -e 's/LDKES = T/LDKES = F/' \
        -e 's/LNEO  = T/LNEO = F/' -e 's/MYHOME = /MYHOME ?= /' \
        STELLOPT/SHARE/make_debian.inc > STELLOPT/SHARE/make_local.inc
    (cd STELLOPT && MACHINE=local STELLOPT_PATH=$PWD MYHOME=$WORK/bin \
        ./build_all -o release -j "$JOBS" LIBSTELL VMEC2000)
fi

if [ ! -x terpsichore/srcs/tpr_ap.x ]; then
    [ -d terpsichore ] || git clone -q "$TERPSICHORE_URL" terpsichore
    git -C terpsichore checkout -q "$TERPSICHORE_REV"
    make -C terpsichore/srcs -j "$JOBS"
    make -C terpsichore/srcs -f Makefile_vmecv92terps
fi

python3 -m pip install --quiet 'booz-xform>=0.0.9,<0.2' scipy
echo "VMEC:        $WORK/STELLOPT/VMEC2000/Release/xvmec2000"
echo "TERPSICHORE: $WORK/terpsichore/srcs/tpr_ap.x"
