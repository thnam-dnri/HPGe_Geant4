from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from scipy.signal import find_peaks, savgol_filter

from ..smearing import fwhm_keV


def _bin_width_keV(axis_keV: np.ndarray) -> float:
    if axis_keV.size < 2:
        return 1.0
    # Use median in case of minor irregularities
    return float(np.median(np.diff(axis_keV)))


def _local_centroid(energy_axis: np.ndarray, counts: np.ndarray, idx_center: int, half_window_bins: int) -> float:
    i0 = max(0, idx_center - half_window_bins)
    i1 = min(counts.size, idx_center + half_window_bins + 1)
    c = counts[i0:i1].astype(float)
    e = energy_axis[i0:i1].astype(float)
    s = c.sum()
    if s <= 0:
        return float(energy_axis[idx_center])
    return float((c * e).sum() / s)


def _local_area(counts: np.ndarray, idx_center: int, half_window_bins: int, baseline: float = 0.0) -> float:
    i0 = max(0, idx_center - half_window_bins)
    i1 = min(counts.size, idx_center + half_window_bins + 1)
    # Simple baseline subtraction (floor at 0)
    vals = np.maximum(0.0, counts[i0:i1].astype(float) - baseline)
    return float(vals.sum())


def find_photopeaks(
    spectrum: np.ndarray,
    energy_axis_keV: np.ndarray,
    prominence_rel: float = 5.0,
    smooth_window: int = 11,
    smooth_polyorder: int = 2,
    max_peaks: int | None = None,
) -> Dict[str, np.ndarray]:
    """Detect photopeaks in a 1D spectrum.

    Parameters:
    - spectrum: 1D counts per energy bin
    - energy_axis_keV: bin centers in keV
    - prominence_rel: threshold factor relative to baseline noise (median absolute deviation)
    - smooth_window: Savitzky-Golay filter window (odd)
    - smooth_polyorder: polynomial order for Savitzky-Golay
    - max_peaks: optional cap on number of peaks (kept by descending prominence)

    Returns a dict with keys:
      indices (int), energies_keV, centroids_keV, areas, fwhm_est_keV, prominences
    """
    if spectrum.size != energy_axis_keV.size:
        raise ValueError("Spectrum and energy axis must have the same length")

    spec = spectrum.astype(float)

    # Smooth to reduce high-frequency noise; ensure window is valid
    w = max(5, smooth_window)
    if w % 2 == 0:
        w += 1
    if w > spec.size:
        w = spec.size - (1 - spec.size % 2)  # make odd and <= size
    smoothed = savgol_filter(spec, window_length=w, polyorder=min(smooth_polyorder, w - 1)) if spec.size >= w else spec

    # Estimate noise scale using MAD of first differences
    diffs = np.diff(smoothed)
    mad = np.median(np.abs(diffs - np.median(diffs))) if diffs.size > 0 else 0.0
    noise = mad if mad > 0 else max(1.0, np.std(smoothed) * 0.1)
    prominence = max(1.0, prominence_rel * noise)

    # Find peaks with minimum prominence and reasonable width
    peaks, props = find_peaks(smoothed, prominence=prominence)
    prominences = props.get("prominences", np.zeros_like(peaks, dtype=float))

    if peaks.size == 0:
        empty = np.array([], dtype=float)
        return {
            "indices": peaks.astype(int),
            "energies_keV": empty,
            "centroids_keV": empty,
            "areas": empty,
            "fwhm_est_keV": empty,
            "prominences": empty,
        }

    # Keep top-N by prominence if requested
    order = np.argsort(prominences)[::-1]
    if max_peaks is not None and max_peaks > 0:
        order = order[:max_peaks]
    peaks = peaks[order]
    prominences = prominences[order]

    # Estimate FWHM at each peak energy via model
    bin_width = _bin_width_keV(energy_axis_keV)
    peak_energies = energy_axis_keV[peaks]
    fwhm_vals = fwhm_keV(peak_energies)
    half_windows = np.clip(np.rint(fwhm_vals / bin_width).astype(int), 1, max(1, spectrum.size // 20))

    # Estimate baseline as median of spectrum (simple global baseline)
    baseline = float(np.median(smoothed))

    centroids = np.empty_like(peak_energies, dtype=float)
    areas = np.empty_like(peak_energies, dtype=float)
    for i, (p, hw) in enumerate(zip(peaks, half_windows)):
        centroids[i] = _local_centroid(energy_axis_keV, smoothed, int(p), int(hw))
        areas[i] = _local_area(smoothed, int(p), int(hw), baseline=baseline)

    return {
        "indices": peaks.astype(int),
        "energies_keV": peak_energies.astype(float),
        "centroids_keV": centroids.astype(float),
        "areas": areas.astype(float),
        "fwhm_est_keV": fwhm_vals.astype(float),
        "prominences": prominences.astype(float),
    }

