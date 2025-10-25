from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import uproot


def collect_root_files(root_dir: str) -> List[str]:
    """Return a list of ROOT files under root_dir (non-recursive)."""
    if not os.path.isdir(root_dir):
        return []
    return [
        os.path.join(root_dir, f)
        for f in os.listdir(root_dir)
        if f.lower().endswith(".root") and os.path.isfile(os.path.join(root_dir, f))
    ]


def _safe_first(values):
    """Return the first element of an uproot array column, decoding bytes if needed."""
    if values is None:
        return None
    if isinstance(values, (list, tuple, np.ndarray)):
        if len(values) == 0:
            return None
        val = values[0]
    else:
        val = values
    if isinstance(val, (bytes, bytearray)):
        return val.decode("utf-8", errors="ignore")
    return val


def read_runinfo(path: str) -> Dict[str, object]:
    """Read RunInfo ntuple into a flat dict of first-row values.

    Falls back to empty dict if tree not present.
    """
    info: Dict[str, object] = {}
    with uproot.open(path) as f:
        if "RunInfo" not in f:
            return info
        try:
            arrs = f["RunInfo"].arrays(library="np")  # type: ignore[assignment]
        except Exception:
            return info
        for k, v in arrs.items():
            info[k] = _safe_first(v)
    return info


def read_energy_axis(path: str, default_min: float = 0.0, default_max: float = 3000.0, default_bins: int = 8192) -> np.ndarray:
    """Load energy bin centers from the Axes_energy_centers ntuple, or return default axis.

    Returns a 1D numpy array of bin centers in keV.
    """
    with uproot.open(path) as f:
        if "Axes_energy_centers" in f and "energy_keV" in f["Axes_energy_centers"]:
            try:
                centers = f["Axes_energy_centers"]["energy_keV"].array(library="np")  # type: ignore[index]
                if centers is not None and len(centers) > 0:
                    return np.asarray(centers, dtype=float)
            except Exception:
                pass
    # Fallback: uniform axis
    return np.linspace(default_min, default_max, default_bins, endpoint=False, dtype=float)


def read_events(path: str) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Read Events ntuple columns Edep_keV and Etrue_keV.

    Returns (Edep_keV, Etrue_keV_or_None).
    """
    with uproot.open(path) as f:
        if "Events" not in f:
            return np.empty(0, dtype=float), None
        tree = f["Events"]
        cols = []
        if "Edep_keV" in tree:
            cols.append("Edep_keV")
        if "Etrue_keV" in tree:
            cols.append("Etrue_keV")
        arrays = tree.arrays(cols, library="np")
        e_dep = np.asarray(arrays.get("Edep_keV", np.empty(0)), dtype=float)
        e_true = arrays.get("Etrue_keV")
        if e_true is not None:
            e_true = np.asarray(e_true, dtype=float)
        return e_dep, e_true  # type: ignore[return-value]


def read_truth_lines(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read Truth_gamma_lines ntuple (E_gamma_keV, I_per_decay) if present.

    Returns (energies_keV, intensities) as 1D numpy arrays (float64). Empty arrays if missing.
    """
    try:
        with uproot.open(path) as f:
            if "Truth_gamma_lines" not in f:
                return np.array([], dtype=float), np.array([], dtype=float)
            tree = f["Truth_gamma_lines"]
            cols = []
            if "E_gamma_keV" in tree:
                cols.append("E_gamma_keV")
            if "I_per_decay" in tree:
                cols.append("I_per_decay")
            if not cols:
                return np.array([], dtype=float), np.array([], dtype=float)
            arrays = tree.arrays(cols, library="np")
            e = np.asarray(arrays.get("E_gamma_keV", np.array([], dtype=float)), dtype=float)
            i = np.asarray(arrays.get("I_per_decay", np.array([], dtype=float)), dtype=float)
            if e.shape != i.shape:
                # pad/truncate to the shorter length
                n = min(e.size, i.size)
                e, i = e[:n], i[:n]
            return e, i
    except Exception:
        return np.array([], dtype=float), np.array([], dtype=float)


def infer_isotope_label(path: str, runinfo: Optional[Dict[str, object]] = None) -> str:
    """Determine isotope label from RunInfo["Nuclide"] or filename prefix before first underscore.

    Normalizes isotope names like "Cs137"; allows fallback to "Unknown".
    """
    if runinfo and isinstance(runinfo.get("Nuclide"), str):
        label = str(runinfo.get("Nuclide"))
        if label:
            return label

    base = os.path.basename(path)
    m = re.match(r"([A-Za-z]{1,2}\d{1,3})[_\-.].*", base)
    if m:
        return m.group(1)
    # Last resort: filename without extension
    stem = os.path.splitext(base)[0]
    if stem:
        return stem
    return "Unknown"


def axis_to_edges(centers: np.ndarray) -> np.ndarray:
    """Compute bin edges from bin centers assuming near-uniform spacing.

    For the first/last edge, extrapolate by half-step.
    """
    centers = np.asarray(centers, dtype=float)
    if centers.size < 2:
        # Single bin: create +/- 0.5 keV edges around center
        c = centers[0] if centers.size == 1 else 0.0
        return np.array([c - 0.5, c + 0.5], dtype=float)
    steps = np.diff(centers)
    step = np.median(steps)
    # Handle possible non-uniform: build edges via midpoints
    mids = centers[:-1] + 0.5 * np.diff(centers)
    first = centers[0] - 0.5 * step
    last = centers[-1] + 0.5 * step
    edges = np.concatenate(([first], mids, [last])).astype(float)
    return edges
