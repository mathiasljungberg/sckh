# Python Tools

## extract_dipoles_2dscan.py

Extract transition dipole components (μx, μy, μz) from `.2dscan` files.

### Purpose

The `.2dscan` files contain transition data including energies and dipole moments for each grid point. The existing `.osc.dat` files only contain the combined oscillator strength `(μx² + μy² + μz²) × E³`. This tool extracts the individual dipole components needed for polarization-resolved spectrum calculations.

### Usage

```bash
# Basic usage - extract dipoles and write .dip.dat files
python extract_dipoles_2dscan.py hdl_high.2dscan

# Extract and verify against existing .osc.dat files
python extract_dipoles_2dscan.py hdl_high.2dscan --verify

# Only verify without writing files
python extract_dipoles_2dscan.py hdl_high.2dscan --verify-only

# Custom output directory and file prefix
python extract_dipoles_2dscan.py hdl_high.2dscan --output-dir dipoles --prefix tmpVH

# Set verification tolerance
python extract_dipoles_2dscan.py hdl_high.2dscan --verify --tolerance 1e-3
```

### Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--output-dir` | `-o` | `.` | Directory for output files |
| `--prefix` | `-p` | `tmpVH` | File prefix (e.g., tmpVH for valence hole states) |
| `--verify` | `-v` | off | Compare computed OSC with existing .osc.dat files |
| `--verify-only` | | off | Only verify, do not write dipole files |
| `--tolerance` | `-t` | `1e-4` | Relative tolerance for verification |

### Input Format

The `.2dscan` file has a 4-line header followed by data blocks:

```
34                           # Grid size X (nx)
34                           # Grid size Y (ny)
0.050000                     # Step size (dx = dy)
-0.350000    -0.350000       # Initial position (x0, y0)
```

Each grid point contains:
```
-608.86379346                # Ground state energy
-589.08766396                # Excited state energy
538.13617035                 # 1s electron binding energy
XES     50                   # Number of transitions
  19.46120   0.00040   0.00070   0.00040   # E_trans, μx, μy, μz
  ...                                       # (N transition lines)
```

Grid ordering is X-fast (Fortran style).

### Output Format

Creates one `.dip.dat` file per transition:

```
# tmpVH001.dip.dat
x          y          μx            μy            μz
-0.350000  -0.350000   0.00040000   0.00070000   0.00040000
-0.300000  -0.350000   0.00050000   0.00120000   0.00060000
...
```

### Example

```bash
cd /path/to/data
python /home/mathias/projects/sckh/src/python_scripts/tools/extract_dipoles_2dscan.py \
    hdl_high.2dscan --verify

# Output:
# Parsing hdl_high.2dscan...
# Grid: 34 x 34, step: 0.05
# Origin: (-0.35, -0.35)
# Transitions: 50
#
# Verifying against existing .osc.dat files...
# Transition   1: OK (max rel err: 3.25e-07, max abs err: 5.00e-09)
# ...
# Verification PASSED
#
# Writing dipole files to ....
# Written: tmpVH001.dip.dat
# ...
# Wrote 50 dipole files
```
