# verse_pypulseq.py - pypulseq helpers for VERSE
import numpy as np
import pypulseq as pp
import warnings

from pypulseq.points_to_waveform import points_to_waveform
from pypulseq.make_arbitrary_grad import make_arbitrary_grad
from pypulseq.make_arbitrary_rf import make_arbitrary_rf
from copy import copy

import pyverse_c  # assumes libverse.dylib/verse.dll/libverse.so is alongside pyverse_c.py

def _grad_to_waveform(g, system):
    """Return (waveform, delay, channel) from a pypulseq gradient or ndarray."""

    # Simple ndarray
    if isinstance(g, np.ndarray):
        return np.ascontiguousarray(g, dtype=np.float64), 0.0, 'z'

    # pypulseq arbitrary gradient object
    if getattr(g, 'type', '') == 'grad' and hasattr(g, 'waveform'):
        ch = getattr(g, 'channel', 'z')
        return np.ascontiguousarray(g.waveform, dtype=np.float64), float(getattr(g, 'delay', 0.0)), ch

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
        wf = points_to_waveform(amps, system.grad_raster_time, times)
        return np.ascontiguousarray(wf, dtype=np.float64), delay, getattr(g, 'channel', 'z')

    raise ValueError("Unsupported gradient input. Provide ndarray, trap, or grad.")


def _rf_to_waveform(rf, system):
    """Return (complex waveform, delay) from pypulseq RF or ndarray."""

    # Simple ndarray
    if isinstance(rf, np.ndarray):
        return np.ascontiguousarray(rf, dtype=np.complex128), 0.0

    # pypulseq RF object
    if getattr(rf, 'type', '') == 'rf' and hasattr(rf, 'signal'):
        sig = np.ascontiguousarray(rf.signal, dtype=np.complex128)
        delay = float(getattr(rf, 'delay', 0.0))
        system.rf_raster_time = round(rf.shape_dur / len(rf.signal), 9)
        return sig, delay

    raise ValueError("Unsupported RF input. Provide ndarray or pypulseq RF.")


def _align_and_pad(wf_a, delay_a, wf_b, delay_b, dt, debugLevel=0):
    """
    Pad two waveforms (complex or real) to a common start/length.
    Returns (a, b, common_start, pad_a_front, pad_a_back, pad_b_front, pad_b_back)
    """

    # Determine start indices
    start_a = int(round(delay_a / dt))
    start_b = int(round(delay_b / dt))
    common_start = min(start_a, start_b)
    if debugLevel > 1:
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
                          debugLevel=0, return_all=True, interp_to_grad_raster=True):
    """
    Return padded (br, bi, g, dt, common_start, channel, rf_pad_front, rf_pad_back) 
    from pypulseq RF/gradient (or ndarray) inputs.
    rf_pad_front, rf_pad_back: padding counts for the RF waveform
    """

    if system is None:
        system = pp.Opts.default
    if dt_rf is not None:
        system.rf_raster_time = dt_rf # Allows user to pass custom RF raster time
    if dt_g is not None:
        system.grad_raster_time = dt_g # Allows user to pass custom gradient raster time

    rf_wave, rf_delay = _rf_to_waveform(rf, system)
    g_wave, g_delay, ch = _grad_to_waveform(grad, system)

    # Set shorthands for raster times
    dt_rf = system.rf_raster_time
    dt_g = system.grad_raster_time

    # Resample RF to grad raster or grad to rf raster if needed
    if abs(dt_g-dt_rf) > 1e-9:
        if debugLevel > 1:
            print(f"Raster times differ: dt_rf={dt_rf}, dt_g={dt_g}")

        if interp_to_grad_raster: # Resample (downsample) RF to gradient raster time
            if debugLevel > 1:
                print(f"Resampling RF from dt={dt_rf} to dt={dt_g}")
            t_rf = (np.arange(len(rf_wave)) + 0.5) * dt_rf
            t_g  = (np.arange(int(np.ceil((len(rf_wave)*dt_rf)/dt_g))) + 0.5) * dt_g
            rf_wave = np.interp(t_g, t_rf, rf_wave.real) + 1j*np.interp(t_g, t_rf, rf_wave.imag)
            rf_delay = 0.0  # now aligned to grad grid
            dt = dt_g

        else:  # Resample (upsample) gradient to RF raster time
            if debugLevel > 1:
                print(f"Resampling gradient from dt={dt_g} to dt={dt_rf}")
            t_g  = (np.arange(len(g_wave)) + 0.5) * dt_g
            t_rf = (np.arange(int(np.ceil((len(g_wave)*dt_g)/dt_rf))) + 0.5) * dt_rf
            g_wave = np.interp(t_rf, t_g, g_wave)
            g_delay = 0.0  # now aligned to rf grid
            dt = dt_rf

    else:
        if debugLevel > 1:
            print(f"No resampling needed; RF and gradient have the same raster time: dt_rf={dt_rf}, dt_g={dt_g}.")
        dt = dt_g # dt_g and dt_rf are the same

    # Align and pad
    rf_wave, g_wave, common_start, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back = \
        _align_and_pad(rf_wave, rf_delay, g_wave, g_delay, dt, debugLevel=debugLevel)

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
                     system=None, debugLevel=0, dt_rf=None, dt_g=None, interp_to_grad_raster=True):

    if system is None:
        system = pp.Opts.default
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
        if debugLevel > 1:
            print(f"Resampling RF back from dt={dt} to dt={system.rf_raster_time}")
        dt_rf = system.rf_raster_time # original RF raster time
        t_rf = (np.arange(len(rf_waveform)) + 0.5) * dt # time vector at gradient raster
        t_orig = (np.arange(int(np.ceil((len(rf_waveform)*dt)/dt_rf))) + 0.5) * dt_rf # time vector at original RF raster
        rf_waveform = np.interp(t_orig, t_rf, rf_waveform.real) + 1j*np.interp(t_orig, t_rf, rf_waveform.imag)

    elif not interp_to_grad_raster and abs(dt-system.grad_raster_time) > 1e-9:
        # If we interpolated the gradient waveform to the RF raster, resample it back to the gradient raster time
        if debugLevel > 1:
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
            if debugLevel > 1:
                print(f"Padding RF waveform from {len(rf_waveform)} to {n_needed} samples to match gradient raster.")
            n_pad = n_needed - len(rf_waveform)
            pad_front = n_pad // 2
            pad_back = n_pad - pad_front
            rf_waveform = np.pad(rf_waveform, (pad_front, pad_back), mode='constant')
            # The delay should be the same as before

    return rf_waveform, rf_delay, g_waveform, g_delay


def verse(rf, grad, type="mintime", max_grad=None, max_slew=None, bmax=None, emax=-1.0, system=None, debugLevel=0, interp_to_grad_raster=True):
    """
    Run VERSE on pypulseq RF/gradient (or ndarray) inputs.
    Returns (arbitrary_rf, arbitrary_grad) pypulseq objects.

    type: "minsar" or "mintime"
    bmax: max RF amplitude for mintime (Hz). If None, uses system.max_rf or max RF amplitude in input RF. (mintime only)
    emax: max RF slew rate for mintime (Hz/s). If -1, uses system.max_rf_slew. (mintime only)
    """

    if system is None:
        system = pp.Opts.default
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
            bmax = np.max(np.sqrt(rf.signal.real**2 + rf.signal.imag**2))
        elif isinstance(rf, np.ndarray):
            bmax = np.max(np.sqrt(rf.real**2 + rf.imag**2))
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
        rf, grad, system, debugLevel=debugLevel, return_all=True, interp_to_grad_raster=interp_to_grad_raster)

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
        system=system, debugLevel=debugLevel, interp_to_grad_raster=interp_to_grad_raster)

    # Build pypulseq arbitrary RF/grad
    rf_out = make_arbitrary_rf(
        signal=rfv_waveform,
        flip_angle=None,
        dwell=system.rf_raster_time,
        delay=rfv_delay,
        no_signal_scaling=True,
        freq_offset=rf.freq_offset if hasattr(rf, 'freq_offset') else 0,
        phase_offset=rf.phase_offset if hasattr(rf, 'phase_offset') else 0,
        use=rf.use if hasattr(rf, 'use') else '',
        system=system,
    )
    grad_out = make_arbitrary_grad(
        channel=ch,
        waveform=gv_waveform,
        delay=gv_delay,
        system=system,
    )
    return rf_out, grad_out


def calculateoffcenterphase(grad, offset, rf=None, nrf=None, gamma=42.576e6, system=None, debugLevel=0):
    """
    Calculate off-center slice RF phase waveform for a pypulseq gradient.

    Useful for calculating variable phase variation due to off-center VERSE pulse.

    Parameters
    ----------
    grad : gradient object or ndarray
        Pypulseq gradient object (trap or grad) or ndarray of gradient waveform (Hz/m)
    offset : float
        Slice offset in meters (m)
    rf : RF object, optional
        Pypulseq RF object to determine nrf automatically. If None, nrf is used.
        If nrf is also None, number of gradient points is used.
    nrf : int, optional
        Number of points in RF waveform. If None, derived from rf object.
        If rf is also None, number of gradient points is used.
    gamma : float, optional
        Gyromagnetic ratio in Hz/T (default: 42.576e6 for 1H)
    system : pypulseq.Opts, optional
        System limits object containing raster times
    debugLevel : int, optional
        Debug verbosity level (default: 0)

    Returns
    -------
    phase : ndarray
        Slice offset phase waveform in radians, length nrf

    Notes
    -----
    - Assumes gradient raster time is a multiple of RF raster time
    - Assumes pulse is symmetrical with respect to setting center phase to 0
    - Default units: gradient (Hz/m), offset (m), gamma (2π), phase (radians)
    - Output phase is not phase-wrapped

    Examples
    --------
    >>> import pypulseq as pp
    >>> system = pp.Opts.default
    >>> g = pp.make_trapezoid(channel='z', amplitude=10e-3, duration=4e-3, system=system)
    >>> rf = pp.make_sinc_pulse(flip_angle=np.pi/2, duration=4e-3, system=system)
    >>> offset = 20e-3  # 20 mm
    >>> phase = calculateoffcenterphase(g, offset, rf=rf, system=system)
    """

    # TODO: Need to pad array to correctly calculate phase

    if system is None:
        system = pp.Opts.default

    # Convert gradient to waveform
    g_wave, g_delay, _ = _grad_to_waveform(grad, system)

    # Determine nrf
    if nrf is None:
        if rf is not None:
            if hasattr(rf, 'signal'):
                nrf = len(rf.signal)
            else:
                nrf = len(g_wave)
                warnings.warn("RF object provided has no 'signal' attribute; using number of gradient points for nrf.", RuntimeWarning)
        else:
            nrf = len(g_wave)
            warnings.warn("RF object not provided; using number of gradient points for nrf.", RuntimeWarning)

    if debugLevel > 1:
        print(f"Gradient waveform: {len(g_wave)} points, delay={g_delay}s")
        print(f"RF points: {nrf}")
        print(f"Offset: {offset*1e3:.2f} mm")

    # Convert raster times to microseconds for C function
    rfup = int(round(system.rf_raster_time * 1e6))  # µs
    gup  = int(round(system.grad_raster_time * 1e6))  # µs

    if debugLevel > 1:
        print(f"RF raster: {rfup} µs")
        print(f"Gradient raster: {gup} µs")
        print(f"Gamma: {gamma/1e6:.4f} MHz/T")

    # Call C function
    phase = pyverse_c.calculateoffcenterphase(
        g_wave, offset, nrf, rfup, gup, gamma
    )

    return phase
