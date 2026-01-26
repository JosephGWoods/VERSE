# verse_python.py - Python wrapper for verse.c functions
"""
Python wrapper for VERSE (Variable-rate Selective Excitation) C functions.

This module provides Python interfaces to minsarverse(), mintverse(), and
calculateoffcenterphase() from verse.c using ctypes.

Author: Joseph G. Woods
Date: December 2025
"""

import numpy as np
import ctypes
import os
import platform

# Auto-detect library name based on OS
_system = platform.system()
if _system == 'Darwin':  # macOS
    _libname = 'libverse.dylib'
elif _system == 'Windows':
    _libname = 'verse.dll'
else:  # Linux
    _libname = 'libverse.so'

# Try to load the library from the same directory as this module
_lib_path = os.path.join(os.path.dirname(__file__), _libname)
if not os.path.exists(_lib_path):
    raise FileNotFoundError(
        f"VERSE library not found at {_lib_path}. "
        f"Please compile verse.c first using: bash build_verse.sh"
    )

_lib = ctypes.CDLL(_lib_path)

# Define C function signatures
_lib.calculateoffcenterphase.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # g
    ctypes.c_double,                                                 # offset
    ctypes.c_long,                                                   # nrf
    ctypes.c_long,                                                   # rfup
    ctypes.c_long,                                                   # gup
    ctypes.c_double,                                                 # gamma
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # phase (out)
]
_lib.calculateoffcenterphase.restype = None

_lib.minsarverse.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # br
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # bi
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # g
    ctypes.c_double,                                                 # dt
    ctypes.c_long,                                                   # n
    ctypes.c_double,                                                 # gmax
    ctypes.c_double,                                                 # smax
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # brv (out)
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # biv (out)
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # gv (out)
]
_lib.minsarverse.restype = None

# Note: mintverse allocates memory internally, so we need a different approach
_lib.mintverse.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # br
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # bi
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # g
    ctypes.c_double,                                                 # dt
    ctypes.c_long,                                                   # n
    ctypes.c_double,                                                 # bmax
    ctypes.c_double,                                                 # gmax
    ctypes.c_double,                                                 # smax
    ctypes.c_double,                                                 # emax
    ctypes.POINTER(ctypes.c_long),                                   # nout
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),                 # brv
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),                 # biv
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),                 # gv
]
_lib.mintverse.restype = None


def calculateoffcenterphase(g, offset, nrf, rfup, gup, gamma):
    """
    Calculate off-center slice RF phase waveform for an arbitrary gradient.

    Useful for calculating variable phase variation due to off-center VERSE pulse.

    Parameters
    ----------
    g : ndarray
        Gradient waveform (arbitrary units)
    offset : float
        Slice offset (same units of distance as g)
    nrf : int
        Number of points in RF (allows different gradient and RF raster)
    rfup : int
        RF raster time (µs)
    gup : int
        Gradient raster time (µs)
    gamma : float
        Gyromagnetic ratio (Hz per same unit of field strength as g)

    Returns
    -------
    phase : ndarray
        Slice offset phase waveform (rad), length nrf

    Notes
    -----
    - gup should be a multiple of rfup (not checked)
    - Assumes pulse is symmetrical with respect to setting center phase to 0
    - If gup > rfup, assumes g is constant during each gradient block
    - Example units: g (T/m), offset (m), rfup (µs), gamma (Hz/T)
      or:            g (G/cm), offset (cm), rfup (µs), gamma (Hz/G)
      or:            g (mT/m), offset (mm), rfup (µs), gamma (MHz/T)
    - Output phase is not phase-wrapped
    
    Examples
    --------
    >>> g = np.array([0, 5, 10, 10, 5, 0], dtype=np.float64)  # mT/m
    >>> offset = 20e-3  # 20 mm
    >>> gamma = 42.576e6  # Hz/T = 42.576 MHz/T = 0.042576 MHz/(mT/m)
    >>> phase = calculateoffcenterphase(g, offset, nrf=60, rfup=10, gup=10, gamma=0.042576)
    """
    # Ensure g is C-contiguous float64 array
    g = np.ascontiguousarray(g, dtype=np.float64)

    # Allocate output array
    phase = np.zeros(nrf, dtype=np.float64)

    # Call C function
    _lib.calculateoffcenterphase(g, offset, nrf, rfup, gup, gamma, phase)

    return phase


def minsarverse(br, bi, g, dt, gmax, smax):
    """
    Convert an RF/gradient pair to the minimum-SAR VERSE equivalent pair.
    
    This function:
    1) Converts non-zero B1 to have a flat B1 amplitude
    2) Removes gradient gmax and slew rate violations
    3) Adjusts dts to maintain constant duration
    4) Resamples waveforms uniformly
    
    Parameters
    ----------
    br : ndarray
        Real part of RF waveform (arbitrary units)
    bi : ndarray
        Imaginary part of RF waveform (arbitrary units)
    g : ndarray
        Gradient waveform (arbitrary units)
    dt : float
        Time increment (arbitrary units)
    gmax : float
        Max gradient value (units of g)
    smax : float
        Max gradient slew rate (units of g per unit dt)
    
    Returns
    -------
    brv : ndarray
        Output real B1 pulse array (same length as input)
    biv : ndarray
        Output imaginary B1 pulse array (same length as input)
    gv : ndarray
        Output gradient array (same length as input)
    
    Notes
    -----
    - Input arrays br, bi, g must be same length
    - Points of zero RF or g will be left unchanged
    - No element of g should exceed gmax initially
    
    Examples
    --------
    >>> br = np.array([0, 1, 2, 2, 1, 0], dtype=np.float64)
    >>> bi = np.zeros(6, dtype=np.float64)
    >>> g  = np.array([0, 5, 10, 10, 5, 0], dtype=np.float64)
    >>> brv, biv, gv = minsarverse(br, bi, g, dt=1e-6, gmax=15.0, smax=1e5)
    """
    # Ensure inputs are C-contiguous float64 arrays
    br = np.ascontiguousarray(br, dtype=np.float64)
    bi = np.ascontiguousarray(bi, dtype=np.float64)
    g  = np.ascontiguousarray(g , dtype=np.float64)
    
    n = len(br)
    if len(bi) != n or len(g) != n:
        raise ValueError("br, bi, and g must have the same length")
    
    # Allocate output arrays
    brv = np.zeros(n, dtype=np.float64)
    biv = np.zeros(n, dtype=np.float64)
    gv  = np.zeros(n, dtype=np.float64)
    
    # Call C function
    _lib.minsarverse(br, bi, g, dt, n, gmax, smax, brv, biv, gv)
    
    return brv, biv, gv


def mintverse(br, bi, g, dt, bmax, gmax, smax, emax=-1.0):
    """
    Convert an RF/gradient pair to the minimum-time VERSE equivalent pair.
    
    This function:
    1) Compresses B1/gradient so one is always maximized
    2) Removes gradient slew rate violations
    3) Resamples waveforms uniformly
    
    Parameters
    ----------
    br : ndarray
        Real part of RF waveform (arbitrary units)
    bi : ndarray
        Imaginary part of RF waveform (arbitrary units)
    g : ndarray
        Gradient waveform (arbitrary units)
    dt : float
        Time increment (arbitrary units)
    bmax : float
        Max RF value (units of b)
    gmax : float
        Max gradient value (units of g)
    smax : float
        Max gradient slew rate (units of g per unit dt)
    emax : float, optional
        Max RF energy (units of b*b*dt). Use -1 to not constrain (default)
    
    Returns
    -------
    brv : ndarray
        Output real B1 pulse array (may be shorter than input)
    biv : ndarray
        Output imaginary B1 pulse array (may be shorter than input)
    gv : ndarray
        Output gradient array (may be shorter than input)
    
    Notes
    -----
    - Output arrays may have different length than input
    - Memory for output is allocated by the C function
    - Points of zero RF or g will be left unchanged
    
    Examples
    --------
    >>> br = np.array([0, 1, 2, 2, 1, 0], dtype=np.float64)
    >>> bi = np.zeros(6, dtype=np.float64)
    >>> g  = np.array([0, 5, 10, 10, 5, 0], dtype=np.float64)
    >>> brv, biv, gv = mintverse(br, bi, g, dt=1e-6, bmax=3.0, gmax=15.0, smax=1e5)
    """
    # Ensure inputs are C-contiguous float64 arrays
    br = np.ascontiguousarray(br, dtype=np.float64)
    bi = np.ascontiguousarray(bi, dtype=np.float64)
    g  = np.ascontiguousarray(g , dtype=np.float64)
    
    n = len(br)
    if len(bi) != n or len(g) != n:
        raise ValueError("br, bi, and g must have the same length")
    
    # Prepare output pointers
    nout = ctypes.c_long()
    brv_ptr = ctypes.POINTER(ctypes.c_double)()
    biv_ptr = ctypes.POINTER(ctypes.c_double)()
    gv_ptr  = ctypes.POINTER(ctypes.c_double)()
    
    # Call C function (it allocates memory internally)
    _lib.mintverse(br, bi, g, dt, n, bmax, gmax, smax, emax,
                   ctypes.byref(nout),
                   ctypes.byref(brv_ptr),
                   ctypes.byref(biv_ptr),
                   ctypes.byref(gv_ptr))
    
    # Convert C arrays to NumPy arrays
    nout_val = nout.value
    brv = np.ctypeslib.as_array(brv_ptr, shape=(nout_val,)).copy()
    biv = np.ctypeslib.as_array(biv_ptr, shape=(nout_val,)).copy()
    gv  = np.ctypeslib.as_array(gv_ptr, shape=(nout_val,)).copy()
    
    # Free C-allocated memory
    ctypes.CDLL(None).free(brv_ptr)
    ctypes.CDLL(None).free(biv_ptr)
    ctypes.CDLL(None).free(gv_ptr)
    
    return brv, biv, gv
