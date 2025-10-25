from __future__ import annotations

import argparse
import json
from typing import Dict, List, Tuple

import h5py
import numpy as np

from .features.peak_finder import find_photopeaks
from .features.spectral_features import compute_spectral_features
from .smearing import fwhm_keV


def _truth_match_ratio(
    peaks_energies: np.ndarray,
    truth_e: np.ndarray,
    axis_keV: np.ndarray,
    fwhm_func=fwhm_keV,
    min_intensity: float = 0.0,
) -> float:
    if truth_e is None or truth_e.size == 0 or peaks_energies.size == 0:
        return 0.0
    # Allow matches within 2*FWHM(E_truth)
    matched = 0
    used = np.zeros(peaks_energies.size, dtype=bool)
    for e_true in truth_e:
        fwhm = float(fwhm_func(e_true))
        tol = 2.0 * fwhm
        # find nearest unused peak
        idx = np.argmin(np.abs(peaks_energies - e_true))
        if not used[idx] and abs(float(peaks_energies[idx]) - float(e_true)) <= tol:
            matched += 1
            used[idx] = True
    denom = max(1, min(len(truth_e), len(peaks_energies)))
    return float(matched) / float(denom)


def build_features(hdf5_path: str, feature_set: str = "phase2", K_peaks: int = 6) -> None:
    with h5py.File(hdf5_path, "a") as f:
        if "energy_axis" not in f:
            raise RuntimeError("HDF5 missing energy_axis dataset")
        axis_keV = f["energy_axis"][:]
        if axis_keV.ndim != 1:
            raise RuntimeError("energy_axis must be 1D")

        # Ensure features group
        feat_root = f.require_group("features")
        feat_grp = feat_root.require_group(feature_set)

        feature_names: List[str] | None = None

        # Iterate isotopes in spectra
        if "spectra" not in f:
            raise RuntimeError("HDF5 missing spectra group")
        spectra_grp = f["spectra"]
        truth_grp = f.get("truth_lines")

        for iso in spectra_grp.keys():
            specs = spectra_grp[iso][:]
            # Optional truth lines per isotope
            truth_e = None
            if truth_grp is not None and iso in truth_grp:
                truth_e = truth_grp[iso]["E_gamma_keV"][:]

            rows: List[np.ndarray] = []
            for spec in specs:
                peaks = find_photopeaks(spec, axis_keV)
                match_ratio = _truth_match_ratio(peaks.get("energies_keV", np.array([])), truth_e, axis_keV) if truth_e is not None else None
                feats, names = compute_spectral_features(spec, axis_keV, peaks, K_peaks=K_peaks, extra_match_ratio=match_ratio)
                if feature_names is None:
                    feature_names = names
                rows.append(feats)

            data = np.stack(rows, axis=0) if rows else np.zeros((0, len(feature_names or [])), dtype=float)
            # write per-isotope dataset
            if iso in feat_grp:
                del feat_grp[iso]
            feat_grp.create_dataset(iso, data=data, compression="gzip", compression_opts=6, shuffle=True)

        # Write metadata
        meta = f.require_group("features_meta")
        meta_set = meta.require_group(feature_set)
        if "feature_names" in meta_set:
            del meta_set["feature_names"]
        names_arr = np.asarray(feature_names or [], dtype=h5py.string_dtype(encoding="utf-8"))
        meta_set.create_dataset("feature_names", data=names_arr)
        meta_set.attrs["K_peaks"] = int(K_peaks)

    print(f"Features written to group: features/{feature_set}")


def main():
    parser = argparse.ArgumentParser(description="Compute Phase 2 features and write into HDF5")
    parser.add_argument("--dataset", type=str, default="data/processed/hpge_dataset.h5")
    parser.add_argument("--feature-set", type=str, default="phase2")
    parser.add_argument("--k-peaks", type=int, default=6)
    args = parser.parse_args()

    build_features(args.dataset, args.feature_set, args.k_peaks)


if __name__ == "__main__":
    main()

