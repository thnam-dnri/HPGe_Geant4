from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np


def _global_stats(spectrum: np.ndarray, energy_axis_keV: np.ndarray) -> Tuple[List[float], List[str]]:
    c = spectrum.astype(float)
    total = float(c.sum())
    if total <= 0:
        return [0.0, 0.0, 0.0, 0.0], ["total_counts", "mean_energy_keV", "std_energy_keV", "entropy"]
    mean_e = float((c * energy_axis_keV).sum() / total)
    var_e = float((c * (energy_axis_keV - mean_e) ** 2).sum() / total)
    std_e = float(np.sqrt(max(0.0, var_e)))
    p = c / total
    # Entropy in nats
    eps = 1e-12
    entropy = float(-(p * np.log(p + eps)).sum())
    return [total, mean_e, std_e, entropy], ["total_counts", "mean_energy_keV", "std_energy_keV", "entropy"]


def _energy_fractions(spectrum: np.ndarray, energy_axis_keV: np.ndarray) -> Tuple[List[float], List[str]]:
    c = spectrum.astype(float)
    total = float(c.sum())
    if total <= 0:
        return [0.0, 0.0, 0.0], ["frac_low_0_300", "frac_mid_300_1000", "frac_high_1000_3000"]
    e = energy_axis_keV
    low = float(c[(e >= 0.0) & (e < 300.0)].sum() / total)
    mid = float(c[(e >= 300.0) & (e < 1000.0)].sum() / total)
    high = float(c[(e >= 1000.0)].sum() / total)
    return [low, mid, high], ["frac_low_0_300", "frac_mid_300_1000", "frac_high_1000_3000"]


def _topk_peak_features(peaks: Dict[str, np.ndarray], K: int = 6) -> Tuple[List[float], List[str]]:
    energies = peaks.get("energies_keV", np.array([], dtype=float))
    areas = peaks.get("areas", np.array([], dtype=float))
    widths = peaks.get("fwhm_est_keV", np.array([], dtype=float))
    if energies.size == 0:
        feats: List[float] = []
        names: List[str] = []
        for k in range(K):
            feats.extend([0.0, 0.0, 0.0])
            names.extend([f"peak{k+1}_energy_keV", f"peak{k+1}_area_frac", f"peak{k+1}_fwhm_keV"])
        # area ratios
        feats.extend([0.0, 0.0, 0.0])
        names.extend(["area_ratio_p1_p2", "area_ratio_p1_p3", "area_ratio_p2_p3"])
        return feats, names

    # Sort peaks by descending area
    order = np.argsort(areas)[::-1]
    energies = energies[order]
    areas = areas[order]
    widths = widths[order]

    total_area = float(areas.sum()) if areas.sum() > 0 else 1.0
    feats = []
    names = []
    for k in range(K):
        if k < areas.size:
            feats.extend([float(energies[k]), float(areas[k] / total_area), float(widths[k])])
        else:
            feats.extend([0.0, 0.0, 0.0])
        names.extend([f"peak{k+1}_energy_keV", f"peak{k+1}_area_frac", f"peak{k+1}_fwhm_keV"])

    # Area ratio features for first three peaks (robust and compact)
    def safe_ratio(a: float, b: float) -> float:
        return float(a / b) if b > 0 else 0.0

    a1 = float(areas[0]) if areas.size > 0 else 0.0
    a2 = float(areas[1]) if areas.size > 1 else 0.0
    a3 = float(areas[2]) if areas.size > 2 else 0.0
    feats.extend([safe_ratio(a1, a2), safe_ratio(a1, a3), safe_ratio(a2, a3)])
    names.extend(["area_ratio_p1_p2", "area_ratio_p1_p3", "area_ratio_p2_p3"])
    return feats, names


def compute_spectral_features(
    spectrum: np.ndarray,
    energy_axis_keV: np.ndarray,
    peaks: Dict[str, np.ndarray],
    K_peaks: int = 6,
    extra_match_ratio: float | None = None,
) -> Tuple[np.ndarray, List[str]]:
    """Compute Phase 2 spectral features and names.

    - spectrum: counts per bin
    - energy_axis_keV: bin centers
    - peaks: peak dict from peak_finder
    - K_peaks: maximum number of peaks to encode
    - extra_match_ratio: optional fraction of matched truth lines to append
    """
    gvals, gnames = _global_stats(spectrum, energy_axis_keV)
    evals, enames = _energy_fractions(spectrum, energy_axis_keV)
    pvals, pnames = _topk_peak_features(peaks, K=K_peaks)

    feats = gvals + evals + pvals
    names = gnames + enames + pnames

    if extra_match_ratio is not None:
        feats.append(float(extra_match_ratio))
        names.append("truth_match_ratio")

    return np.asarray(feats, dtype=float), names

