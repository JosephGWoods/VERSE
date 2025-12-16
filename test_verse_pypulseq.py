# test_verse.py - Test script for VERSE python wrapper
"""
Test script to verify the VERSE pyhton wrapper is working correctly.
Run this after building the library to ensure everything is set up properly.
"""

import pypulseq as pp
import numpy as np
import sys
import os

# Add current directory to path if needed
sys.path.insert(0, os.path.dirname(__file__))

try:
    import verse_pypulseq as verse
    print("Successfully imported verse_pypulseq")
except ImportError as e:
    print(f"Failed to import verse_pypulseq: {e}")
    print("\nMake sure you have:")
    print("1. Compiled verse.c using: bash build_verse.sh")
    print("2. The library file (libverse.dylib) is in the same directory as verse_pypulseq.py")
    sys.exit(1)


def test_verse(type, debugLevel=2):
    """Test verse functions"""
    print("\n--- Testing verse ---")

    system = pp.Opts(max_grad=30, grad_unit="mT/m", max_slew=150, slew_unit="T/m/s",
                     rf_ringdown_time=30e-6, rf_dead_time=100e-6, rf_raster_time=1e-6, grad_raster_time=10e-6)

    # Create test RF and gradient waveforms
    rf, gz, gzr = pp.make_sinc_pulse(
        flip_angle=20*np.pi/180,
        duration=0.5e-3,
        dwell=system.grad_raster_time,
        slice_thickness=8e-3,
        time_bw_product=4,
        system=system,
        return_gz=True
        )
    
    # Extract waveforms
    br, bi, g = verse._get_padded_waveforms(rf, gz, system=system, debugLevel=debugLevel, return_all=False)
    
    try:
        if type == "minsarverse":
            print("Running minsarverse...")
            rfv, gzv = verse.minsarverse(rf, gz, system=system, debugLevel=debugLevel)
            print(f"minsarverse completed successfully")

        elif type == "mintverse":
            print("Running mintverse...")
            bmax = np.max(np.sqrt(br**2 + bi**2))
            rfv, gzv = verse.mintverse(rf, gz, bmax=bmax, system=system, debugLevel=debugLevel)
            print(f"mintverse completed successfully")

        brv, biv, gv = verse._get_padded_waveforms(rfv, gzv, system=system, debugLevel=debugLevel, return_all=False)

        # Print summary statistics
        print(f"Input length: {len(br)}, Output length: {len(brv)}")
        print(f"Input duration: {len(br)*system.grad_raster_time*1e3:.2f} ms, Output duration: {len(brv)*system.grad_raster_time*1e3:.2f} ms")
        print(f"Input max gradient: {np.max(np.abs(g))*1e3/system.gamma:.2f} mT/m (limit: {system.max_grad*1e3/system.gamma} mT/m)")
        print(f"Output max gradient: {np.max(np.abs(gv))*1e3/system.gamma:.2f} mT/m (limit: {system.max_grad*1e3/system.gamma} mT/m)")
        print(f"Input RF peak amplitude: {np.max(np.sqrt(br**2 + bi**2))*1e6/system.gamma:.2f} µT")
        print(f"Output RF peak amplitude: {np.max(np.sqrt(brv**2 + biv**2))*1e6/system.gamma:.2f}  µT")
        print(f"Input RF energy: {np.sum(br**2 + bi**2) * system.grad_raster_time / system.gamma :.6f} T*s")
        print(f"Output RF energy: {np.sum(brv**2 + biv**2) * system.grad_raster_time / system.gamma:.6f} T*s")

        # Plot input (br, bi, gz) vs results (brv, biv, gzv)
        import matplotlib.pyplot as plt
        t = np.arange(len(br)) * system.grad_raster_time * 1e3  # in ms
        tv = np.arange(len(brv)) * system.grad_raster_time * 1e3  # in ms
        plt.figure(figsize=(12, 8))
        plt.subplot(3, 1, 1)
        plt.plot(t, br, label='Input B1+ real', linestyle='--')
        plt.plot(tv, brv, label='Output B1+ real')
        plt.title('RF Real Part')
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.subplot(3, 1, 2)
        plt.plot(t, bi, label='Input B1+ imaginary', linestyle='--')
        plt.plot(tv, biv, label='Output B1+ imaginary')
        plt.title('RF Imaginary Part')
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude')
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

        return True
    except Exception as e:
        print(f"minsarverse failed: {e}")
        return False



if __name__ == "__main__":
    print("="*60)
    print("VERSE pypulseq wrapper test script")
    print("="*60)
    
    results = []
    results.append(test_verse("minsarverse"))
    results.append(test_verse("mintverse"))
    
    print("\n" + "="*60)
    if all(results):
        print("All tests passed!")
    else:
        print("Some tests failed. Please check the error messages above.")
        sys.exit(1)