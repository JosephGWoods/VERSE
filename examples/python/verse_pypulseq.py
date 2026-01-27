# test_verse.py - Test script for VERSE python wrapper
"""
Test script to verify the VERSE python wrapper is working correctly.
Run this after building the library to ensure everything is set up properly.
"""

import pypulseq as pp
import numpy as np
import pyverse_pulseq as ppverse

def example_verse(type="minsar"):
    """
    Example functions to demonstrate minsarverse and mintverse function integration with pypulseq.
    Generates a test RF pulse and gradient waveform, applies VERSE, and prints summary statistics.
    """

    # interp_to_grad_raster:
    #   True means interpolate RF to gradient raster time
    #   False means interpolate gradient to RF raster time
    interp_to_grad_raster = False

    system = pp.Opts(max_grad=80, grad_unit="mT/m", max_slew=200, slew_unit="T/m/s",
                        rf_ringdown_time=30e-6, rf_dead_time=100e-6, rf_raster_time=1e-6, grad_raster_time=10e-6)

    seq = pp.Sequence(system)

    # Create test RF and gradient waveforms
    rf, gz, gzr = pp.make_sinc_pulse(
        flip_angle      = 90*np.pi/180,
        duration        = 2.56e-3,
        slice_thickness = 5.0e-3,
        time_bw_product = 10,
        delay           = system.rf_dead_time,
        system          = system,
        return_gz       = True
    )

    if type == "minsar":
        print("\nRunning minsar VERSE...")
        rfv, gzv = ppverse.verse(rf, gz, type="minsar", system=system, interp_to_grad_raster=interp_to_grad_raster)
        print(f"minsar VERSE completed successfully")

    elif type == "mintime":
        print("\nRunning mintime VERSE...")
        bmax = np.max(np.sqrt(rf.signal.real**2 + rf.signal.imag**2)) # Max B1 amplitude of input RF
        rfv, gzv = ppverse.verse(rf, gz, type="mintime", bmax=bmax, emax=-1.0, system=system, interp_to_grad_raster=interp_to_grad_raster)
        print(f"mintime VERSE completed successfully")

    # Print summary statistics
    print(f"Input RF length: {len(rf.signal)}, Output RF length: {len(rfv.signal)}")
    print(f"Input RF duration: {len(rf.signal)*system.rf_raster_time*1e3:.2f} ms, Output RF duration: {len(rfv.signal)*system.rf_raster_time*1e3:.2f} ms")
    print(f"Input max gradient: {gz.amplitude*1e3/system.gamma:.2f} mT/m (limit: {system.max_grad*1e3/system.gamma} mT/m)")
    print(f"Output max gradient: {np.max(np.abs(gzv.waveform))*1e3/system.gamma:.2f} mT/m (limit: {system.max_grad*1e3/system.gamma} mT/m)")
    print(f"Input RF peak amplitude: {np.max(np.sqrt(rf.signal.real**2 + rf.signal.imag**2))*1e6/system.gamma:.2f} µT")
    print(f"Output RF peak amplitude: {np.max(np.sqrt(rfv.signal.real**2 + rfv.signal.imag**2))*1e6/system.gamma:.2f}  µT")
    print(f"Input RF energy: {1e12 * np.sum(rf.signal.real**2 + rf.signal.imag**2) * system.rf_raster_time / system.gamma :.6f} µT^2*s")
    print(f"Output RF energy: {1e12 * np.sum(rfv.signal.real**2 + rfv.signal.imag**2) * system.rf_raster_time / system.gamma:.6f} µT^2*s")

    # Add to sequence and plot using pypulseq
    seq.add_block(rf, gz)
    seq.add_block(gzr)
    seq.add_block(rfv, gzv)
    seq.add_block(gzr)
    seq.plot()

    # Extract waveforms for manual plotting overlaying the waveforms
    br, bi, g    = ppverse._get_padded_waveforms(rf,  gz,  system=system, return_all=False, interp_to_grad_raster=interp_to_grad_raster)
    brv, biv, gv = ppverse._get_padded_waveforms(rfv, gzv, system=system, return_all=False, interp_to_grad_raster=interp_to_grad_raster)

    # Plot input (br, bi, gz) vs results (brv, biv, gzv)
    import matplotlib.pyplot as plt
    t = np.arange(len(br)) * system.rf_raster_time * 1e3  # in ms
    tv = np.arange(len(brv)) * system.rf_raster_time * 1e3  # in ms
    plt.figure(figsize=(12, 8))
    plt.subplot(3, 1, 1)
    plt.plot(t, br, label='Input B1+ real', linestyle='--')
    plt.plot(tv, brv, label='Output B1+ real')
    plt.title('RF Real Part')
    plt.xlabel('Time (ms)')
    plt.ylabel('Amplitude (Hz)')
    plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(t, bi, label='Input B1+ imaginary', linestyle='--')
    plt.plot(tv, biv, label='Output B1+ imaginary')
    plt.title('RF Imaginary Part')
    plt.xlabel('Time (ms)')
    plt.ylabel('Amplitude (Hz)')
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(t, g, label='Input Gz', linestyle='--')
    plt.plot(tv, gv, label='Output Gz')
    plt.title('Gradient Z')
    plt.xlabel('Time (ms)')
    plt.ylabel('Amplitude (Hz/m)')
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    print("="*60)
    print("VERSE pypulseq example script")
    print("="*60)

    example_verse(type="minsar")
    example_verse(type="mintime")