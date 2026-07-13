FROM quay.io/pypa/manylinux_2_28_x86_64@sha256:b04887b645dde99b9e955aeae3ff4da414992d0bd88259f046295b56361c5614

RUN dnf install -y \
        cmake-3.26.5-2.el8 \
        gcc-gfortran-8.5.0-28.el8_10.alma.1 \
        lapack-devel-3.8.0-9.el8_10 \
        netcdf-devel-4.7.0-3.el8 \
        pkgconf-pkg-config-1.4.2-1.el8 \
    && dnf clean all

RUN /opt/python/cp39-cp39/bin/python -m pip install --no-cache-dir \
    auditwheel==6.4.2
