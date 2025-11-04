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
    """Backward-compatible reader for (E_gamma_keV, I_per_decay)."""
    e, i, _, _ = read_truth_lines_detailed(path)
    return e, i


def read_truth_lines_detailed(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read Truth_gamma_lines with emitter info if present.

    Returns (E_gamma_keV, I_per_decay, Emitter[str], EmitterHalfLife_s)
    as arrays of equal length. If emitter fields are missing, returns
    empty strings and zeros for those arrays.
    """
    try:
        with uproot.open(path) as f:
            if "Truth_gamma_lines" not in f:
                z = np.array([], dtype=float)
                return z, z, np.array([], dtype="S1"), z
            tree = f["Truth_gamma_lines"]
            cols = []
            have_emit = False
            have_hl = False
            if "E_gamma_keV" in tree:
                cols.append("E_gamma_keV")
            if "I_per_decay" in tree:
                cols.append("I_per_decay")
            if "Emitter" in tree:
                cols.append("Emitter")
                have_emit = True
            if "EmitterHalfLife_s" in tree:
                cols.append("EmitterHalfLife_s")
                have_hl = True
            if not cols:
                z = np.array([], dtype=float)
                return z, z, np.array([], dtype="S1"), z
            arrays = tree.arrays(cols, library="np")
            e = np.asarray(arrays.get("E_gamma_keV", np.array([], dtype=float)), dtype=float)
            i = np.asarray(arrays.get("I_per_decay", np.array([], dtype=float)), dtype=float)
            if e.shape != i.shape:
                n = min(e.size, i.size)
                e, i = e[:n], i[:n]
            if have_emit:
                em = arrays.get("Emitter")
                if em is None:
                    emit = np.array([b"" for _ in range(e.size)], dtype="S10")
                else:
                    emit = np.asarray(em)
                    if emit.shape[0] != e.shape[0]:
                        n = min(emit.shape[0], e.shape[0])
                        emit = emit[:n]
                        e, i = e[:n], i[:n]
            else:
                emit = np.array([b"" for _ in range(e.size)], dtype="S10")
            if have_hl:
                hl = arrays.get("EmitterHalfLife_s")
                half = np.asarray(hl if hl is not None else np.zeros_like(e), dtype=float)
                if half.shape[0] != e.shape[0]:
                    n = min(half.shape[0], e.shape[0])
                    half = half[:n]
                    e, i, emit = e[:n], i[:n], emit[:n]
            else:
                half = np.zeros_like(e)
            return e, i, emit, half
    except Exception:
        z = np.array([], dtype=float)
        return z, z, np.array([], dtype="S1"), z


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
