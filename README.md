# VERSE

Variable-rate Selective Excitation (VERSE) algorithms and helpers in C, MATLAB, and Python, including pypulseq integration.

## Features
- Minimum-time VERSE (mintverse)
- Minimum-SAR VERSE (minsarverse; fixed duration)
- Core C implementation for performance
- MATLAB:
  - Pure MATLAB implementations: `mintverse.m`, `minsarverse.m`
  - C-accelerated MEX wrappers: `mintverse_c.m`, `minsarverse_c.m`
- Python:
  - ctypes wrapper around the C library: `pyverse_c.py`
  - pypulseq helpers: `pyverse_pulseq.py` (including off-center phase utilities)

## Repository layout
- `c/` — C core (`verse.c`, `verse.h`)
- `matlab/` — MATLAB API (pure and MEX wrappers)
- `python/` — Python API (ctypes + pypulseq helpers)
- `build_python_lib.sh` — Build script for Python shared library
- `build_matlab_mex.m` — Build script for MATLAB MEX files
- `archive/` — Legacy/testing implementations

## Requirements
- C compiler (gcc/clang/MSVC)
- MATLAB R2018a+ recommended (MEX compiled with `-R2017b` API)
- Python 3.8+ (NumPy, SciPy, matplotlib, pypulseq)

## Installation

### Python
Install the package (editable install recommended during development):
```bash
pip install -e .
```

Optionally build the C shared library used by the ctypes wrapper:
```bash
bash build_python_lib.sh
```
This produces:
- Linux: `python/libverse.so`
- macOS: `python/libverse.dylib`
- Windows (MSYS/Cygwin): `python/verse.dll`

### MATLAB
Pure MATLAB functions require no build:
- `matlab/mintverse.m`
- `matlab/minsarverse.m`

To build the C-accelerated MEX functions, run inside MATLAB from the repo root:
```matlab
build_matlab_mex
```
This compiles with the `-R2017b` MEX API and outputs:
- `matlab/mintverse_mex.<mexext>`
- `matlab/minsarverse_mex.<mexext>`

Wrappers `mintverse_c.m` and `minsarverse_c.m` call these MEX binaries.

If no compiler is configured, run:
```matlab
mex -setup C
```

## Usage

### Python (pypulseq)
```python
import numpy as np
import pypulseq as pp
from pyverse_pulseq import verse, calculateoffcenterphase

system = pp.Opts()
rf, gz, gzr = pp.make_sinc_pulse(
    flip_angle=np.deg2rad(90), duration=4e-3, slice_thickness=5e-3, return_gz=True, system=system
)

# Minimum-time VERSE
rfv, gv = verse(rf=rf, grad=gz, system=system, type="mintime")

# Minimum-SAR VERSE
rfv, gv = verse(rf=rf, grad=gz, system=system, type="minsar")

# Off-center phase (example: 20 mm)
phase = calculateoffcenterphase(g=gv, offset=20e-3, rf=rf, system=system)
```

### MATLAB
Pure MATLAB:
```matlab
[b1v, gv] = mintverse(b1, g, dt, bmax, gmax, smax); % minimum-time
[b1v, gv] = minsarverse(b1, g, dt, gmax, smax);     % minimum-SAR
```

C-accelerated (after building MEX):
```matlab
[b1v, gv] = mintverse_c(b1, g, dt, bmax, gmax, smax);
[b1v, gv] = minsarverse_c(b1, g, dt, gmax, smax);
```

## Development notes
- Build scripts:
  - Python: `bash build_python_lib.sh`
  - MATLAB: `build_matlab_mex` (inside MATLAB)

## License
See `LICENSE`

## Links
- Source: https://github.com/JosephGWoods/VERSE
- Issues: https://github.com/JosephGWoods/VERSE/issues