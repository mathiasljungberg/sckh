"""Physical constants matching Fortran m_constants.F90 for numerical compatibility."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Constants:
    """Physical constants in SI units.

    Values match those in src/m_constants.F90 for numerical consistency
    when comparing Python and Fortran results.
    """

    # Bohr radius in meters
    bohr: float = 0.52917725e-10

    # Pi
    pi: float = 3.14159265358979

    # Planck constant over 2*pi (J*s)
    hbar: float = 1.05457168e-34

    # Atomic mass unit (kg)
    u: float = 1.66053873e-27

    # Speed of light (m/s)
    c: float = 2.99792458e8

    # Hartree energy (J)
    hartree: float = 4.35974418e-18

    # Planck constant (J*s)
    h: float = 6.6260688e-34

    # cm^-1 to J^-1 conversion
    cm: float = 5.03411759319722e22

    # Electron volt (J)
    eV: float = 1.60217646e-19

    # Hartree to eV conversion
    hartree2eV: float = 27.211396132

    # Boltzmann constant (J/K)
    k_b: float = 1.3806503e-23

    # Atomic time unit (s)
    autime: float = 2.418884326505e-17


# Global instance for easy access
CONST = Constants()
