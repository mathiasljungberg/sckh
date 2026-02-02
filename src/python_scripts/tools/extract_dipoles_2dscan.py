#!/usr/bin/env python3
"""
Extract transition dipole components from .2dscan files.

The .2dscan file contains ground/excited state energies and transition data
(energy, μx, μy, μz) for each grid point. This script extracts the individual
dipole components and writes them to separate .dip.dat files.

Output format: x  y  μx  μy  μz
"""

import argparse
import numpy as np
from pathlib import Path
import sys


def parse_header(lines, idx=0):
    """Parse the 4-line header of a .2dscan file.

    Returns:
        nx, ny: Grid dimensions
        dx: Step size (same for x and y)
        x0, y0: Initial position
        idx: Updated line index
    """
    nx = int(lines[idx].strip())
    ny = int(lines[idx + 1].strip())
    dx = float(lines[idx + 2].strip())
    x0_y0 = lines[idx + 3].split()
    x0 = float(x0_y0[0])
    y0 = float(x0_y0[1])
    return nx, ny, dx, x0, y0, idx + 4


def parse_2dscan(filepath):
    """Parse a .2dscan file and extract all data.

    Args:
        filepath: Path to the .2dscan file

    Returns:
        nx, ny: Grid dimensions
        dx: Step size
        x0, y0: Initial position
        ground_energies: Array of ground state energies [ny, nx]
        excited_energies: Array of excited state energies [ny, nx]
        binding_energies: Array of 1s binding energies [ny, nx]
        transition_energies: Array of transition energies [ny, nx, n_trans]
        dipoles: Array of dipole components [ny, nx, n_trans, 3] (μx, μy, μz)
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Parse header
    nx, ny, dx, x0, y0, idx = parse_header(lines)

    # Initialize arrays - we'll determine n_trans from first grid point
    n_trans = None
    ground_energies = np.zeros((ny, nx))
    excited_energies = np.zeros((ny, nx))
    binding_energies = np.zeros((ny, nx))
    transition_energies = None
    dipoles = None

    # Parse data blocks (X-fast = Fortran ordering)
    # Outer loop: Y (slow), Inner loop: X (fast)
    for j in range(ny):
        for i in range(nx):
            # Read 3 energy lines
            ground_energies[j, i] = float(lines[idx].strip())
            excited_energies[j, i] = float(lines[idx + 1].strip())
            binding_energies[j, i] = float(lines[idx + 2].strip())
            idx += 3

            # Read XES line to get number of transitions
            xes_line = lines[idx].split()
            assert xes_line[0] == 'XES', f"Expected 'XES' at line {idx+1}, got {xes_line[0]}"
            n_trans_current = int(xes_line[1])
            idx += 1

            # Initialize arrays on first grid point
            if n_trans is None:
                n_trans = n_trans_current
                transition_energies = np.zeros((ny, nx, n_trans))
                dipoles = np.zeros((ny, nx, n_trans, 3))
            else:
                assert n_trans_current == n_trans, \
                    f"Inconsistent number of transitions: {n_trans_current} vs {n_trans}"

            # Read transition data
            for k in range(n_trans):
                trans_data = lines[idx].split()
                transition_energies[j, i, k] = float(trans_data[0])
                dipoles[j, i, k, 0] = float(trans_data[1])  # μx
                dipoles[j, i, k, 1] = float(trans_data[2])  # μy
                dipoles[j, i, k, 2] = float(trans_data[3])  # μz
                idx += 1

    return (nx, ny, dx, x0, y0, ground_energies, excited_energies,
            binding_energies, transition_energies, dipoles)


def write_dipole_files(nx, ny, dx, x0, y0, dipoles, output_dir, prefix):
    """Write dipole component files for each transition.

    Args:
        nx, ny: Grid dimensions
        dx: Step size
        x0, y0: Initial position
        dipoles: Array of dipole components [ny, nx, n_trans, 3]
        output_dir: Directory for output files
        prefix: File prefix (e.g., 'tmpVH')
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_trans = dipoles.shape[2]

    # Compute grid coordinates
    x = x0 + np.arange(nx) * dx
    y = y0 + np.arange(ny) * dx

    for k in range(n_trans):
        filename = output_dir / f"{prefix}{k+1:03d}.dip.dat"
        with open(filename, 'w') as f:
            # Write in Fortran order (X-fast)
            for j in range(ny):
                for i in range(nx):
                    f.write(f"  {x[i]:11.6f}  {y[j]:11.6f}  "
                            f"{dipoles[j, i, k, 0]:13.8f}  "
                            f"{dipoles[j, i, k, 1]:13.8f}  "
                            f"{dipoles[j, i, k, 2]:13.8f}\n")
        print(f"Written: {filename}")


def verify_with_osc_files(nx, ny, dx, x0, y0, transition_energies, dipoles,
                          osc_dir, prefix, tolerance=1e-4):
    """Verify extracted dipoles by comparing computed OSC with existing files.

    The oscillator strength formula is: (μx² + μy² + μz²) × E³

    Args:
        nx, ny: Grid dimensions
        dx: Step size
        x0, y0: Initial position
        transition_energies: Array [ny, nx, n_trans]
        dipoles: Array [ny, nx, n_trans, 3]
        osc_dir: Directory containing .osc.dat files
        prefix: File prefix
        tolerance: Relative tolerance for comparison

    Returns:
        True if verification passes, False otherwise
    """
    osc_dir = Path(osc_dir)
    n_trans = dipoles.shape[2]

    # Compute grid coordinates
    x = x0 + np.arange(nx) * dx
    y = y0 + np.arange(ny) * dx

    all_passed = True

    for k in range(n_trans):
        osc_file = osc_dir / f"{prefix}{k+1:03d}.osc.dat"
        if not osc_file.exists():
            print(f"Warning: {osc_file} not found, skipping verification")
            continue

        # Read existing OSC file
        osc_data = np.loadtxt(osc_file)

        # Compute OSC from dipoles: (μx² + μy² + μz²) × E³
        computed_osc = []
        for j in range(ny):
            for i in range(nx):
                mu_sq = (dipoles[j, i, k, 0]**2 +
                         dipoles[j, i, k, 1]**2 +
                         dipoles[j, i, k, 2]**2)
                E = transition_energies[j, i, k]
                osc = mu_sq * E**3
                computed_osc.append([x[i], y[j], osc])

        computed_osc = np.array(computed_osc)

        # Compare coordinates
        coord_diff_x = np.max(np.abs(computed_osc[:, 0] - osc_data[:, 0]))
        coord_diff_y = np.max(np.abs(computed_osc[:, 1] - osc_data[:, 1]))

        if coord_diff_x > 1e-6 or coord_diff_y > 1e-6:
            print(f"ERROR: Coordinate mismatch in transition {k+1}")
            print(f"  Max X diff: {coord_diff_x}")
            print(f"  Max Y diff: {coord_diff_y}")
            all_passed = False
            continue

        # Compare OSC values (relative error)
        # Avoid division by zero for very small values
        mask = np.abs(osc_data[:, 2]) > 1e-10
        if np.any(mask):
            rel_err = np.abs(computed_osc[mask, 2] - osc_data[mask, 2]) / np.abs(osc_data[mask, 2])
            max_rel_err = np.max(rel_err)
            mean_rel_err = np.mean(rel_err)
        else:
            max_rel_err = 0.0
            mean_rel_err = 0.0

        # Also check absolute error for small values
        abs_err = np.abs(computed_osc[:, 2] - osc_data[:, 2])
        max_abs_err = np.max(abs_err)

        if max_rel_err > tolerance and max_abs_err > 1e-6:
            print(f"ERROR: OSC mismatch in transition {k+1}")
            print(f"  Max relative error: {max_rel_err:.6e}")
            print(f"  Mean relative error: {mean_rel_err:.6e}")
            print(f"  Max absolute error: {max_abs_err:.6e}")
            all_passed = False
        else:
            print(f"Transition {k+1:3d}: OK (max rel err: {max_rel_err:.2e}, "
                  f"max abs err: {max_abs_err:.2e})")

    return all_passed


def main():
    parser = argparse.ArgumentParser(
        description='Extract dipole components from .2dscan files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python extract_dipoles_2dscan.py hdl_high.2dscan --verify
  python extract_dipoles_2dscan.py hdl_high.2dscan --output-dir dipoles --prefix tmpVH
        """
    )

    parser.add_argument('input_file', help='Input .2dscan file')
    parser.add_argument('--output-dir', '-o', default='.',
                        help='Directory for output files (default: current)')
    parser.add_argument('--prefix', '-p', default='tmpVH',
                        help='File prefix (default: tmpVH)')
    parser.add_argument('--verify', '-v', action='store_true',
                        help='Compare with existing .osc.dat files')
    parser.add_argument('--verify-only', action='store_true',
                        help='Only verify, do not write dipole files')
    parser.add_argument('--tolerance', '-t', type=float, default=1e-4,
                        help='Relative tolerance for verification (default: 1e-4)')

    args = parser.parse_args()

    input_file = Path(args.input_file)
    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    print(f"Parsing {input_file}...")
    (nx, ny, dx, x0, y0, ground_energies, excited_energies,
     binding_energies, transition_energies, dipoles) = parse_2dscan(input_file)

    n_trans = dipoles.shape[2]
    print(f"Grid: {nx} x {ny}, step: {dx}")
    print(f"Origin: ({x0}, {y0})")
    print(f"Transitions: {n_trans}")

    if args.verify or args.verify_only:
        print("\nVerifying against existing .osc.dat files...")
        osc_dir = input_file.parent
        passed = verify_with_osc_files(nx, ny, dx, x0, y0, transition_energies,
                                       dipoles, osc_dir, args.prefix, args.tolerance)
        if not passed:
            print("\nVerification FAILED")
            if args.verify_only:
                sys.exit(1)
        else:
            print("\nVerification PASSED")

    if not args.verify_only:
        print(f"\nWriting dipole files to {args.output_dir}...")
        write_dipole_files(nx, ny, dx, x0, y0, dipoles, args.output_dir, args.prefix)
        print(f"\nWrote {n_trans} dipole files")


if __name__ == '__main__':
    main()
