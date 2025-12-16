# verse_pypulseq.py - pypulseq helpers for VERSE
import numpy as np
import pypulseq as pp

from pypulseq.points_to_waveform import points_to_waveform
from pypulseq.make_arbitrary_grad import make_arbitrary_grad
from pypulseq.make_arbitrary_rf import make_arbitrary_rf

import verse_python as verse  # assumes libverse.dylib/verse.dll/libverse.so is alongside verse_python.py


def _grad_to_waveform(g, system):
    """Return (waveform, delay, channel) from a pypulseq gradient or ndarray."""
    if isinstance(g, np.ndarray):
        return np.ascontiguousarray(g, dtype=np.float64), 0.0, 'z'
    # Arbitrary gradient
    if getattr(g, 'type', '') == 'grad' and hasattr(g, 'waveform'):
        ch = getattr(g, 'channel', 'z')
        return np.ascontiguousarray(g.waveform, dtype=np.float64), float(getattr(g, 'delay', 0.0)), ch
    # Trapezoid
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
    if isinstance(rf, np.ndarray):
        return np.ascontiguousarray(rf, dtype=np.complex128), 0.0, system.rf_raster_time
    if getattr(rf, 'type', '') == 'rf' and hasattr(rf, 'signal'):
        sig = np.ascontiguousarray(rf.signal, dtype=np.complex128)
        delay = float(getattr(rf, 'delay', 0.0))
        dt_rf = round(rf.shape_dur / len(rf.signal), 9)
        return sig, delay, dt_rf
    raise ValueError("Unsupported RF input. Provide ndarray or pypulseq RF.")


def _align_and_pad(wf_a, delay_a, wf_b, delay_b, dt, debugLevel=0):
    """
    Pad two waveforms (complex or real) to a common start/length.
    Returns (a, b, common_start, pad_a_front, pad_a_back, pad_b_front, pad_b_back)
    """
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


def _get_padded_waveforms(rf, grad, system=None, dt_g=None, dt_rf=None, debugLevel=0, return_all=True):
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

    dt_g  = system.grad_raster_time

    rf_wave, rf_delay, dt_rf = _rf_to_waveform(rf, system)
    g_wave, g_delay, ch = _grad_to_waveform(grad, system)

    # Resample RF to grad raster if needed
    if not np.isclose(dt_rf, dt_g):
        if debugLevel > 1:
            print(f"Resampling RF from dt={dt_rf} to dt={dt_g}")
        t_rf = (np.arange(len(rf_wave)) + 0.5) * dt_rf
        t_g  = (np.arange(int(np.ceil((len(rf_wave)*dt_rf)/dt_g))) + 0.5) * dt_g
        rf_wave = np.interp(t_g, t_rf, rf_wave.real) + 1j*np.interp(t_g, t_rf, rf_wave.imag)
        rf_delay = 0.0  # now aligned to grad grid
    dt = dt_g

    # Align and pad
    rf_wave, g_wave, common_start, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back = _align_and_pad(rf_wave, rf_delay, g_wave, g_delay, dt, debugLevel=debugLevel)

    br = np.ascontiguousarray(rf_wave.real, dtype=np.float64)
    bi = np.ascontiguousarray(rf_wave.imag, dtype=np.float64)
    g  = np.ascontiguousarray(g_wave      , dtype=np.float64)

    if return_all:
        return br, bi, g, dt, common_start, ch, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back
    else:
        return br, bi, g
    

def _unpad_waveforms(rf_waveform, g_waveform, dt, common_start, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back):

    # Strip RF padding
    if rf_pad_back > 0:
        rf_waveform = rf_waveform[rf_pad_front:-rf_pad_back]
    else:
        rf_waveform = rf_waveform[rf_pad_front:]

    # Adjust RF delay
    rf_delay = (common_start + rf_pad_front) * dt
    rf_delay = round(rf_delay, 9)

    # Strip gradient padding
    if g_pad_back > 0:
        g_waveform = g_waveform[g_pad_front:-g_pad_back]
    else:
        g_waveform = g_waveform[g_pad_front:]
    
    # Adjust gradient delay
    g_delay = (common_start + g_pad_front) * dt
    g_delay = round(g_delay, 9)

    return rf_waveform, rf_delay, g_waveform, g_delay


def minsarverse(rf, grad, max_grad=None, max_slew=None, system=None, debugLevel=0):
    """
    Run minsarverse on pypulseq RF/gradient (or ndarray) inputs.
    Returns (arbitrary_rf, arbitrary_grad) pypulseq objects.
    """

    if system is None:
        system = pp.Opts.default
    if max_grad is None:
        max_grad = system.max_grad
    if max_slew is None:
        max_slew = system.max_slew

    # Extract and pad waveforms to match lengths
    br, bi, g, dt, common_start, ch, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back =_get_padded_waveforms(
        rf, grad, system, debugLevel=debugLevel)

    # Run minsarverse
    brv, biv, gv = verse.minsarverse(br, bi, g, dt, max_grad, max_slew)

    # Strip padding
    rfv_waveform, rfv_delay, gv_waveform, gv_delay = _unpad_waveforms(
        (brv + 1j*biv), gv, dt, common_start, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back)

    # Build pypulseq arbitrary RF/grad
    rf_out = make_arbitrary_rf(
        signal=rfv_waveform,
        flip_angle=None,
        dwell=dt,
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
        delay=gv_delay, # avoid floating point issues
        system=system,
    )
    return rf_out, grad_out


def mintverse(rf, grad, bmax=np.inf, emax=-1.0, max_grad=None, max_slew=None, system=None, debugLevel=0):
    """
    Run mintverse on pypulseq RF/gradient (or ndarray) inputs.
    Returns (arbitrary_rf, arbitrary_grad) pypulseq objects.
    """

    if system is None:
        system = pp.Opts.default
    if max_grad is None:
        max_grad = system.max_grad
    if max_slew is None:
        max_slew = system.max_slew

    # Extract and pad waveforms to match lengths
    br, bi, g, dt, common_start, ch, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back =_get_padded_waveforms(
        rf, grad, system, debugLevel=debugLevel)

    # Run mintverse
    brv, biv, gv = verse.mintverse(br, bi, g, dt, bmax, max_grad, max_slew, emax)

    # Strip padding
    rfv_waveform, rfv_delay, gv_waveform, gv_delay = _unpad_waveforms(
        (brv + 1j*biv), gv, dt, common_start, rf_pad_front, rf_pad_back, g_pad_front, g_pad_back)

    # Build pypulseq arbitrary RF/grad
    rf_out = make_arbitrary_rf(
        signal=rfv_waveform,
        flip_angle=None,
        dwell=dt,
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
        delay=gv_delay, # avoid floating point issues
        system=system,
    )
    return rf_out, grad_out