# gliss (Python skeleton)

Python bindings skeleton for [GLISS](https://github.com/itpplasma/GLISS), the
Global Linear Ideal Stability Solver. This package reserves the `gliss` name
on PyPI. Only `gliss.version()` is implemented so far, proving `ctypes`
linkage to the compiled Fortran C API; the Mercier objective binding is a
later milestone.

Build the shared library first:

```sh
cmake -S .. -B ../build -G Ninja
cmake --build ../build --target gliss_c
```

Then:

```python
import gliss
print(gliss.__version__)   # package version, no library needed
print(gliss.version())     # version reported by the compiled library
```

`gliss.version()` loads the shared library via `ctypes`. Set `GLISS_LIB` to
the built `libgliss_c.so` path if it is not already on the loader search
path.
