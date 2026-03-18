# pyverse_pulseq.py - pypulseq helpers for VERSE
"""
Provides functions to run VERSE on pypulseq RF and gradient objects, including
handling of waveform extraction, alignment, padding, and conversion back to
pypulseq arbitrary RF/grad objects.

Also includes a function to calculate off-center phase variation for a given
gradient and slice offset, which can be useful for including in the RF signal
itself, rather than relying on vendor-specific arbitrary gradient slice-offset
calculations.

Author: Joseph G. Woods
Date: December 2025
"""
from types import SimpleNamespace
from typing import Tuple, Union
from copy import copy

import numpy as np

from pypulseq.opts import Opts
from pypulseq.points_to_waveform import points_to_waveform
from pypulseq.make_arbitrary_grad import make_arbitrary_grad
from pypulseq.make_arbitrary_rf import make_arbitrary_rf


import pyverse_c  # assumes libverse.dylib/verse.dll/libverse.so is alongside pyverse_c.py

def _grad_to_waveform(g, system):
    """Return (waveform, delay, channel) from a pypulseq gradient or ndarray."""

    # Simple ndarray
    if isinstance(g, np.ndarray):
        return np.ascontiguousarray(g, dtype=np.float64), 0.0, 'z'

    # pypulseq arbitrary gradient object
    if getattr(g, 'type', '') == 'grad' and hasattr(g, 'waveform'):
        return np.ascontiguousarray(g.waveform, dtype=np.float64), float(getattr(g, 'delay', 0.0)), getattr(g, 'channel', 'z')

    # pypulseq trapezoid gradient object
    if getattr(g, 'type', '') == 'trap':
        delay = float(getattr(g, 'delay', 0.0))
        rt, ft, rft = g.rise_time, g.flat_time, g.fall_time
        if ft > 0:
            times = np.array([delay,
                              delay + rt,
                              delay + rt + ft,
                              delay + rt + ft + rft])
            amps = np.array([0, g.amplitude, g.amplitude, 0], dtype=np.float64)
        else:
            times = np.array([delay,
                              delay + rt,
                              delay + rt + rft])
            amps = np.array([0, g.amplitude, 0], dtype=np.float64)
        waveform = points_to_waveform(amps, system.grad_raster_time, times)
        return np.ascontiguousarray(waveform, dtype=np.float64), delay, getattr(g, 'channel', 'z')

    raise ValueError("Unsupported gradient input. Provide ndarray, trap, or grad.")


def _rf_to_waveform(rf, system):
    """Return (complex waveform, delay) from pypulseq RF or ndarray."""

    # Simple ndarray
    if isinstance(rf, np.ndarray):
        return np.ascontiguousarray(rf, dtype=np.complex128), 0.0

    # pypulseq RF object
    if getattr(rf, 'type', '') == 'rf' and hasattr(rf, 'signal'):
        waveform = np.ascontiguousarray(rf.signal, dtype=np.complex128)
        delay = float(getattr(rf, 'delay', 0.0))
        system.rf_raster_time = round(rf.shape_dur / len(rf.signal), 9) # Calculate RF raster time from shape_dur and modify system accordingly
        return waveform, delay

    raise ValueError("Unsupported RF input. Provide ndarray or pypulseq RF.")


def _align_and_pad(wf_a, delay_a, wf_b, delay_b, dt, debug_level=0):
    """
    Pad two waveforms (complex or real) to a common start/length.
    Returns (a, b, common_start, pad_a_front, pad_a_back, pad_b_front, pad_b_back)
    """

    # Determine start indices
    start_a = int(round(delay_a / dt))
    start_b = int(round(delay_b / dt))
    common_start = min(start_a, start_b)
    if debug_level > 1:
        print(f"Aligning waveforms: start_a={start_a}, start_b={start_b}, common_start={common_start}")

    # Pad starts of the waveforms
    pad_a_front = start_a - common_start
    pad_b_front = start_b - common_start

    a = np.pad(wf_a, (pad_a_front, 0), mode='constant')
    b = np.pad(wf_b, (pad_b_front, 0), mode='constant')

    # Pad ends of the waveforms
    n = max(len(a), len(b))
    pad_a_back = max(0, n - len(a))
    pad_b_back = max(0, n - len(b))

    a = np.pad(a, (0, pad_a_back), mode='constant')
    b = np.pad(b, (0, pad_b_back), mode='constant')

    return a, b, common_start, pad_a_front, pad_a_back, pad_b_front, pad_b_back


def _get_padded_waveforms(rf, grad, system=None, dt_g=None, dt_rf=None,
                          debug_level=0, return_all=True, interp_to_grad_raster=False):
    """
    Return padded (br, bi, g, dt, common_start, channel, rf_pad_front, rf_pad_back) 
    from pypulseq RF/gradient (or ndarray) inputs.
    rf_pad_front, rf_pad_back: padding counts for the RF waveform
    """

    if system is None:
        system = Opts.default
    else:
        system = copy(system) # avoid modifying user system internally
    if dt_rf is not None:
        system.rf_raster_time = dt_rf # Allows user to pass custom RF raster time
    if dt_g is not None:
        system.grad_raster_time = dt_g # Allows user to pass custom gradient raster time

    rf_wave, rf_delay   = _rf_to_waveform(rf, system)
    g_wave, g_delay, ch = _grad_to_waveform(grad, system)

    # Set shorthands for raster times
    dt_rf = system.rf_raster_time
    dt_g  = system.grad_raster_time

    # Resample RF to grad raster or grad to rf raster if needed
    if abs(dt_g-dt_rf) > 1e-9:
        if debug_level > 1:
            print(f"Raster times differ: dt_rf={dt_rf}, dt_g={dt_g}")

        if interp_to_grad_raster: # Resample (downsample) RF to gradient raster time
            if debug_level > 1:
                print(f"Resampling RF from dt={dt_rf} to dt={dt_g}")
            t_rf = (np.arange(len(rf_wave)) + 0.5) * dt_rf
            t_g  = (np.arange(int(np.ceil((len(rf_wave)*dt_rf)/dt_g))) + 0.5) * dt_g
            rf_wave = np.interp(t_g, t_rf, rf_wave.real) + 1j*np.interp(t_g, t_rf, rf_wave.imag)
            dt = dt_g

        else:  # Resample (upsample) gradient to RF raster time
            if debug_level > 1:
                print(f"Resampling gradient from dt={dt_g} to dt={dt_rf}")
            t_g  = (np.arange(len(g_wave)) + 0.5) * dt_g
            t_rf = (np.arange(int(np.ceil((len(g_wave)*dt_g)/dt_rf))) + 0.5) * dt_rf
            g_wave = np.interp(x=t_rf, xp=t_g, fp=g_wave)
            dt = dt_rf

    else:
        if debug_level > 1:
            print(f"No resampling needed; RF and gradient have the same raster time: dt_rf={dt_rf}, dt_g={dt_g}.")
        dt = dt_g # dt_g and dt_rf are the same

    # Align and pad
    rf_wave, g_wave, common_start, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back = \
        _align_and_pad(rf_wave, rf_delay, g_wave, g_delay, dt, debug_level=debug_level)

    br = np.ascontiguousarray(rf_wave.real, dtype=np.float64)
    bi = np.ascontiguousarray(rf_wave.imag, dtype=np.float64)
    g  = np.ascontiguousarray(g_wave      , dtype=np.float64)

    if return_all:
        return br, bi, g, \
        dt, common_start, ch, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back
    else:
        return br, bi, g


def _unpad_waveforms(rf_waveform, g_waveform, dt, common_start,
                     rf_pad_front, rf_pad_back, g_pad_front, g_pad_back,
                     system=None, debug_level=0, dt_rf=None, dt_g=None, interp_to_grad_raster=False):

    if system is None:
        system = Opts.default
    else:
        system = copy(system) # avoid modifying user system internally
    if dt_rf is not None:
        system.rf_raster_time = dt_rf # Allows user to pass custom RF raster time
    if dt_g is not None:
        system.grad_raster_time = dt_g # Allows user to pass custom gradient raster time

    # Strip RF padding
    if rf_pad_back > 0:
        rf_waveform = rf_waveform[rf_pad_front:-rf_pad_back]
    else:
        rf_waveform = rf_waveform[rf_pad_front:]

    # Strip gradient padding
    if g_pad_back > 0:
        g_waveform = g_waveform[g_pad_front:-g_pad_back]
    else:
        g_waveform = g_waveform[g_pad_front:]

    # Adjust RF delay (should be back on gradient raster time by default)
    rf_delay = (common_start + rf_pad_front) * dt
    rf_delay = round(rf_delay, 9)

    # Adjust gradient delay (should be back on gradient raster time by default)
    g_delay = (common_start + g_pad_front) * dt
    g_delay = round(g_delay, 9)

    if interp_to_grad_raster and abs(dt-system.rf_raster_time) > 1e-9:
        # If we interpolated the rf waveform to the gradient raster, resample it back to the rf raster time
        if debug_level > 1:
            print(f"Resampling RF back from dt={dt} to dt={system.rf_raster_time}")
        dt_rf = system.rf_raster_time # original RF raster time
        t_rf = (np.arange(len(rf_waveform)) + 0.5) * dt # time vector at gradient raster
        t_orig = (np.arange(int(np.ceil((len(rf_waveform)*dt)/dt_rf))) + 0.5) * dt_rf # time vector at original RF raster
        rf_waveform = np.interp(t_orig, t_rf, rf_waveform.real) + 1j*np.interp(t_orig, t_rf, rf_waveform.imag)

    elif not interp_to_grad_raster and abs(dt-system.grad_raster_time) > 1e-9:
        # If we interpolated the gradient waveform to the RF raster, resample it back to the gradient raster time
        if debug_level > 1:
            print(f"Resampling gradient back from dt={dt} to dt={system.grad_raster_time}")
        dt_g = system.grad_raster_time # original gradient raster time
        t_g = (np.arange(len(g_waveform)) + 0.5) * dt # time vector at RF raster
        t_orig = (np.arange(int(np.ceil((len(g_waveform)*dt)/dt_g))) + 0.5) * dt_g # time vector at original gradient raster
        g_waveform = np.interp(t_orig, t_g, g_waveform)

        # We also need to make sure the RF duration is a multiple of the gradient raster time again
        dur_rf = round(len(rf_waveform) * system.rf_raster_time, 9) # RF duration
        dur_rf_on_g_raster = round(np.ceil(dur_rf / system.grad_raster_time) * system.grad_raster_time, 9) # RF duration rounded up to the gradient raster time

        if dur_rf_on_g_raster > dur_rf:
            n_needed = int(round(dur_rf_on_g_raster / system.rf_raster_time))
            # Pad RF waveform front and back to maintain alignment
            if debug_level > 1:
                print(f"Padding RF waveform from {len(rf_waveform)} to {n_needed} samples to match gradient raster.")
            n_pad = n_needed - len(rf_waveform)
            pad_front = n_pad // 2
            pad_back = n_pad - pad_front
            rf_waveform = np.pad(rf_waveform, (pad_front, pad_back), mode='constant')
            # The delay should be the same as before

    return rf_waveform, rf_delay, g_waveform, g_delay


def verse(
    rf: Union[SimpleNamespace, np.ndarray],
    grad: Union[SimpleNamespace, np.ndarray],
    type: str = "mintime",
    max_grad: Union[float, None] = None,
    max_slew: Union[float, None] = None,
    bmax: Union[float, None] = None,
    emax: float = -1.0,
    system: Union[Opts, None] = None,
    debug_level: int = 0,
    interp_to_grad_raster: bool = False
    ) -> Tuple[SimpleNamespace, SimpleNamespace]:
    """
    Run VERSE on pypulseq RF/gradient (or ndarray) inputs.
    Returns arbitrary_rf and arbitrary_grad pypulseq objects.

    Parameters
    ----------
    rf :        SimpleNamespace or numpy ndarray
                Pypulseq RF object or RF waveform.
    grad :      SimpleNamespace or numpy ndarray
                Pypulseq gradient object (trap or grad) or gradient waveform (Hz/m)
    type :      str, optional
                VERSE type to run. Options are: "minsar" or "mintime".
                default: "mintime"
    max_grad :  float, optional
                Maximum gradient strength of accompanying slice select trapezoidal event.
                default: None, which uses system.max_grad.
    max_slew :  float, optional
                Maximum gradient slew rate of accompanying slice select trapezoidal event.
                default: None, which uses system.max_slew.
    bmax :      float, optional
                Max RF amplitude (units of RF, e.g. Hz).
                Only used for mintime VERSE. Ignored for minsar VERSE.
                default: None, which uses system.max_rf or max RF amplitude in input RF.
    emax :      float, optional
                Max RF energy (units of RF^2 * dt).
                Only used for mintime VERSE. Ignored for minsar VERSE.
                default: -1, which does not constrain energy.
    system :    Opts, optional
                System limits object containing raster times and max RF/grad
                default: None, which uses Opts.default.
    debug_level : int, optional
                Debug verbosity level (default: 0). Higher values print more debug info.
    interp_to_grad_raster : bool, optional
                True: interpolates the RF waveform to the gradient raster time before running
                VERSE, and then back to the original RF raster time after.
                If False, interpolates the gradient to the RF raster time and then back instead.
                default: False (interpolate gradient to RF raster).

    Returns
    -------
    rf_out : SimpleNamespace
        Pypulseq arbitrary RF object containing the VERSEd RF waveform and parameters.
    grad_out : SimpleNamespace
        Pypulseq arbitrary gradient object containing the VERSEd gradient waveform and parameters.
    """

    if system is None:
        system = Opts.default
    else:
        system = copy(system) # avoid modifying user system internally
    if max_grad is None:
        max_grad = system.max_grad
    if max_slew is None:
        max_slew = system.max_slew
    if type.lower() == "mintime" and (bmax is None or bmax < 0):
        if getattr(system, "max_rf", None) is not None and system.max_rf > 0:
            bmax = system.max_rf
        elif hasattr(rf, 'signal'):
            bmax = np.max(np.abs(rf.signal))
        elif isinstance(rf, np.ndarray):
            bmax = np.max(np.abs(rf))
        else:
            raise ValueError("Cannot determine bmax for mintime VERSE. Please provide bmax parameter.")

    # Error checking
    if type.lower() not in ["minsar", "mintime"]:
        raise ValueError(f"Unsupported VERSE type: {type}. Supported types are 'minsar' and 'mintime'.")
    if max_grad <= 0:
        raise ValueError("max_grad must be greater than 0.")
    if max_slew <= 0:
        raise ValueError("max_slew must be greater than 0.")

    # Extract and pad waveforms to match lengths
    br, bi, g, dt, common_start, ch, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back =_get_padded_waveforms(
        rf, grad, system=system, debug_level=debug_level, return_all=True, interp_to_grad_raster=interp_to_grad_raster)

    # Call C VERSE functions
    if type.lower() == "minsar":
        brv, biv, gv = pyverse_c.minsarverse(br, bi, g, dt, max_grad, max_slew)
    elif type.lower() == "mintime":
        brv, biv, gv = pyverse_c.mintverse(br, bi, g, dt, bmax, max_grad, max_slew, emax)
    else:
        raise ValueError(f"Unsupported VERSE type: {type}. Supported types are 'minsar' and 'mintime'.")

    # Strip padding
    rfv_waveform, rfv_delay, gv_waveform, gv_delay = _unpad_waveforms(
        (brv + 1j*biv), gv, dt, common_start,
        rf_pad_front, rf_pad_back, g_pad_front, g_pad_back,
        system=system, debug_level=debug_level, interp_to_grad_raster=interp_to_grad_raster)

    # Build pypulseq arbitrary RF/grad
    rf_out = make_arbitrary_rf(
        signal            = rfv_waveform,
        flip_angle        = 0, # Not used due to no_signal_scaling=True, but pypulseq requires a value
        dwell             = system.rf_raster_time,
        delay             = rfv_delay,
        no_signal_scaling = True,
        freq_offset       = rf.freq_offset if hasattr(rf, 'freq_offset') else 0,
        phase_offset      = rf.phase_offset if hasattr(rf, 'phase_offset') else 0,
        use               = rf.use if hasattr(rf, 'use') else '',
        system            = system
    )
    grad_out = make_arbitrary_grad(
        channel  = ch,
        waveform = gv_waveform,
        delay    = gv_delay,
        system   = system
    )

    return rf_out, grad_out


def apply_offcenter_phase(
    rf: Union[SimpleNamespace, np.ndarray],
    grad: Union[SimpleNamespace, np.ndarray],
    offset: float,
    gamma: float = 1,
    system: Union[Opts, None] = None,
    dt_g: Union[float, None] = None,
    dt_rf: Union[float, None] = None,
    debug_level: int = 0
    ) -> SimpleNamespace:
    """
    Calculate the off-center slice phase waveform needed for an arbitrary gradient shape,
    add it to the RF waveform and return the modified RF as a pypulseq RF object.

    Useful for calculating variable phase variation due to off-center VERSE pulse.

    Parameters
    ----------
    rf : Pypulseq RF object or ndarray
        Pypulseq RF object or complex ndarray of RF waveform (arbitrary units).
    grad : gradient object or ndarray
        Pypulseq gradient object (trap or grad) or ndarray of gradient waveform
        (arbitrary units, but expect Hz/m by default).
    offset : float
        Slice offset ((same units of distance as grad, e.g. m if gradient is in Hz/m).
    gamma : float, optional
        Gyromagnetic ratio (Hz per same unit of field strength as g. default: 1 assuming Hz/Hz)
    system : pypulseq.Opts, optional
        System limits object containing raster times (default: None, which uses Opts.default)
    dt_g : float, optional
        Custom gradient raster time in seconds (overrides system.grad_raster_time)
    dt_rf : float, optional
        Custom RF raster time in seconds (overrides system.rf_raster_time)
    debug_level : int, optional
        Debug verbosity level (default: 0)

    Returns
    -------
    rf_out : RF object
        Modified Pypulseq RF object with off-center phase applied

    Notes
    -----
    - Assumes gradient raster time is a multiple of RF raster time
    - Assumes pulse is symmetrical with respect to setting center phase to 0

    Examples
    --------
    >>> import pypulseq as pp
    >>> system = pp.Opts()
    >>> g = pp.make_trapezoid(channel='z', amplitude=10e-3, duration=4e-3, system=system)
    >>> rf = pp.make_sinc_pulse(flip_angle=np.pi/2, duration=4e-3, system=system)
    >>> offset = 20e-3  # 20 mm
    >>> rf_out = apply_offcenter_phase(rf, g, offset, system=system)
    """

    if system is None:
        system = Opts.default

    # Extract, pad and interpolate waveforms to match lengths
    br, bi, g, dt, common_start, ch, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back =_get_padded_waveforms(
        rf, grad, system=system, dt_g=dt_g, dt_rf=dt_rf, debug_level=debug_level, return_all=True, interp_to_grad_raster=False)

    # Determine nrf for input to C function
    nrf = br.shape[0]

    # Convert the now common raster time to microseconds for C function
    rfup = int(round(dt * 1e6)) # µs
    gup  = int(round(dt * 1e6)) # µs

    # Call C function
    rf_phase = pyverse_c.calculateoffcenterphase(
        g, offset, nrf, rfup, gup, gamma
    )

    # Add off-center phase waveform to the original RF waveform
    rf_waveform_out = (br + 1j*bi) * np.exp(1j*rf_phase)

    # Strip padding and undo interpolation
    rf_waveform_out, rf_delay, _, _ = _unpad_waveforms(
        rf_waveform_out, g, dt, common_start,
        rf_pad_front, rf_pad_back, g_pad_front, g_pad_back,
        system=system, debug_level=debug_level, interp_to_grad_raster=False)

    # Build pypulseq arbitrary RF/grad
    rf_out = make_arbitrary_rf(
        signal            = rf_waveform_out,
        flip_angle        = 0, # Not used due to no_signal_scaling=True, but pypulseq requires a value
        dwell             = system.rf_raster_time,
        delay             = rf_delay,
        no_signal_scaling = True,
        freq_offset       = rf.freq_offset if hasattr(rf, 'freq_offset') else 0,
        phase_offset      = rf.phase_offset if hasattr(rf, 'phase_offset') else 0,
        use               = rf.use if hasattr(rf, 'use') else '',
        system            = system
    )

    return rf_out
