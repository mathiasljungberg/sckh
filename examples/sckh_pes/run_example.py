#!/usr/bin/env python
"""Example: SCKH spectrum calculation for a PES model system.

This example demonstrates:
1. Loading configuration from YAML
2. Running classical trajectory dynamics
3. Computing X-ray emission spectrum via FFT
4. Comparing 'standard' vs 'fortran' compatibility modes

Usage:
    cd examples/sckh_pes
    python run_example.py
"""

from pathlib import Path

import python_scripts.dynamics_1d as dyn


def main():
    # Load configuration
    config = dyn.load_full_config("dynamics.yaml")

    print("=" * 60)
    print("SCKH Spectrum Calculation Example")
    print("=" * 60)
    print(f"\nConfiguration:")
    print(f"  Reduced mass: {config.dynamics.mu} amu")
    print(f"  Grid: {config.dynamics.grid.npoints} points, "
          f"dx={config.dynamics.grid.dx} Å")
    print(f"  Time: {config.dynamics.time.nsteps} steps, "
          f"dt={config.dynamics.time.dt} fs")
    print(f"  Sampling: {config.dynamics.sampling.npoints_x} x "
          f"{config.dynamics.sampling.npoints_mom} = "
          f"{config.dynamics.sampling.npoints_x * config.dynamics.sampling.npoints_mom} trajectories")
    print(f"  Broadening (FWHM): {config.spectrum.gamma_fwhm} eV")
    print(f"  Final states: {len(config.spectrum.pes_final_list)}")
    print(f"  Compatibility mode: {config.spectrum.compatibility_mode}")

    # Run dynamics
    print("\n" + "-" * 60)
    print("Running classical trajectory dynamics...")
    runner = dyn.trajectory.DynamicsRunner(config)
    result_dyn = runner.run(verbose=False)
    print(f"  Generated {len(result_dyn.trajectories)} trajectories")

    # Compute spectrum
    print("\nComputing spectrum via FFT...")
    calculator = dyn.SpectrumCalculator(config)
    calculator.load_surfaces()
    result_sp = calculator.compute_spectrum(result_dyn.trajectories)

    print(f"  Mean transition energy: {result_sp.E_mean:.2f} eV")
    print(f"  Frequency range: {result_sp.omega[0]:.2f} to {result_sp.omega[-1]:.2f} eV")

    # Find peak
    import numpy as np
    peak_idx = np.argmax(result_sp.sigma_tot)
    print(f"  Peak position: {result_sp.omega[peak_idx]:.2f} eV")
    print(f"  Peak intensity: {result_sp.sigma_tot[peak_idx]:.4f}")

    # Save results
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)
    runner.save_results(result_dyn, output_dir=str(output_dir))
    calculator.save_results(result_sp, output_dir)
    print(f"\nResults saved to {output_dir}/")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
