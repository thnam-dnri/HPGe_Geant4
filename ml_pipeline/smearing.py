from __future__ import annotations

from typing import Optional

import numpy as np


def fwhm_keV(energy_keV: np.ndarray | float, a2: float = 10.26, b: float = 0.00168) -> np.ndarray:
    """HPGe energy resolution FWHM model in keV.

    FWHM(E) = sqrt(a^2 + b * E), with a^2 (keV^2) and b (keV^-1).
    Defaults follow the project plan.
    """
    E = np.asarray(energy_keV, dtype=float)
    return np.sqrt(np.maximum(0.0, a2 + b * E))


def gaussian_smear(
    energies_keV: np.ndarray,
    seed: Optional[int] = None,
    a2: float = 10.26,
    b: float = 0.00168,
) -> np.ndarray:
    """Apply Gaussian smearing with energy-dependent sigma derived from FWHM(E).

    Parameters:
    - energies_keV: input energies to smear (e.g., Edep_keV)
    - seed: optional RNG seed for reproducibility
    - a2, b: resolution parameters
    """
    rng = np.random.default_rng(seed)
    E = np.asarray(energies_keV, dtype=float)
    sigma = fwhm_keV(E, a2=a2, b=b) / 2.355
    # Guard against zero/NaN sigma
    sigma = np.where(np.isfinite(sigma), sigma, 0.0)
    return E + rng.normal(loc=0.0, scale=sigma, size=E.shape)

