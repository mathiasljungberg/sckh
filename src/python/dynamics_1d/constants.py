"""Physical constants for dynamics calculations.

Uses scipy.constants for up-to-date CODATA values.
Fortran-compatible values available via CONST_FORTRAN for validation.
"""

from dataclasses import dataclass

from scipy import constants as sp


@dataclass(frozen=True)
class Constants:
    """Physical constants in SI units from scipy.constants."""

    # Bohr radius (m)
    bohr: float = sp.physical_constants["Bohr radius"][0]

    # Pi
    pi: float = sp.pi

    # Planck constant over 2*pi (J*s)
    hbar: float = sp.hbar

    # Atomic mass unit (kg)
    u: float = sp.atomic_mass

    # Speed of light (m/s)
    c: float = sp.c

    # Hartree energy (J)
    hartree: float = sp.physical_constants["Hartree energy"][0]

    # Planck constant (J*s)
    h: float = sp.h

    # Electron volt (J)
    eV: float = sp.eV

    # Hartree to eV conversion
    hartree2eV: float = sp.physical_constants["Hartree energy in eV"][0]

    # Boltzmann constant (J/K)
    k_b: float = sp.k

    # Atomic time unit (s)
    autime: float = sp.physical_constants["atomic unit of time"][0]

    # cm^-1 per Joule (to convert J to wavenumbers)
    @property
    def cm(self) -> float:
        return 1.0 / (100.0 * self.h * self.c)


@dataclass(frozen=True)
class ConstantsFortran:
    """Physical constants matching Fortran m_constants.F90.

    Use these for numerical validation against Fortran output.
    """

    bohr: float = 0.52917725e-10
    pi: float = 3.14159265358979
    hbar: float = 1.05457168e-34
    u: float = 1.66053873e-27
    c: float = 2.99792458e8
    hartree: float = 4.35974418e-18
    h: float = 6.6260688e-34
    cm: float = 5.03411759319722e22  # cm^-1 per Joule
    eV: float = 1.60217646e-19
    hartree2eV: float = 27.211396132
    k_b: float = 1.3806503e-23
    autime: float = 2.418884326505e-17


# Default: use scipy constants
CONST = Constants()

# For Fortran validation
CONST_FORTRAN = ConstantsFortran()
