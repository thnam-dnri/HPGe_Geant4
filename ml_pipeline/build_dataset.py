from __future__ import annotations

import argparse
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np

from . import io as io_utils
from .smearing import gaussian_smear


def _histogram_spectrum(energies_keV: np.ndarray, centers_keV: np.ndarray) -> np.ndarray:
    edges = io_utils.axis_to_edges(centers_keV)
    hist, _ = np.histogram(energies_keV, bins=edges)
    return hist.astype(np.int64, copy=False)


def build_dataset(
    root_dir: str,
    output_path: str,
    isotopes_filter: Optional[List[str]] = None,
    smear: bool = True,
    use_true: bool = False,
    seed: Optional[int] = 42,
    default_axis: Tuple[float, float, int] = (0.0, 3000.0, 8192),
) -> None:
    """Aggregate spectra from ROOT files into an HDF5 dataset grouped by isotope.

    - root_dir: directory containing ROOT files to read
    - output_path: HDF5 file to write
    - isotopes_filter: optional whitelist of isotope labels
    - smear: apply Gaussian smearing using the project resolution model
    - use_true: if True, smear Etrue_keV instead of Edep_keV (if available)
    - seed: RNG seed for deterministic smearing
    - default_axis: (emin, emax, nbins) if file lacks Axes_energy_centers
    """
    files = io_utils.collect_root_files(root_dir)
    if not files:
        raise FileNotFoundError(f"No ROOT files found in: {root_dir}")

    # Determine a common energy axis
    centers = None
    for p in files:
        centers = io_utils.read_energy_axis(p, *default_axis)
        if centers is not None and centers.size > 0:
            break
    if centers is None or centers.size == 0:
        emin, emax, nbins = default_axis
        centers = np.linspace(emin, emax, nbins, endpoint=False, dtype=float)

    # Aggregate spectra by isotope
    spectra_by_iso: Dict[str, List[np.ndarray]] = defaultdict(list)
    file_counts: Dict[str, int] = defaultdict(int)
    truth_by_iso: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = defaultdict(list)

    for path in files:
        runinfo = io_utils.read_runinfo(path)
        label = io_utils.infer_isotope_label(path, runinfo)
        if isotopes_filter and label not in isotopes_filter:
            continue

        e_dep, e_true = io_utils.read_events(path)
        if e_dep.size == 0 and (e_true is None or e_true.size == 0):
            continue

        energies = e_true if (use_true and e_true is not None and e_true.size > 0) else e_dep
        # Guard: keep only finite and positive energies
        energies = energies[np.isfinite(energies) & (energies > 0.0)]
        if smear and energies.size > 0:
            energies = gaussian_smear(energies, seed=seed)

        spec = _histogram_spectrum(energies, centers)
        spectra_by_iso[label].append(spec)
        file_counts[label] += 1

        # Truth lines (optional)
        t_e, t_i = io_utils.read_truth_lines(path)
        if t_e.size > 0 and t_i.size > 0:
            truth_by_iso[label].append((t_e, t_i))

    if not spectra_by_iso:
        raise RuntimeError("No spectra were produced. Check ROOT files and filters.")

    # Ensure directory exists
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    # Write HDF5
    with h5py.File(output_path, "w") as f:
        spectra_grp = f.create_group("spectra")
        for iso, spectra_list in spectra_by_iso.items():
            if len(spectra_list) == 0:
                continue
            data = np.stack(spectra_list, axis=0)
            spectra_grp.create_dataset(
                iso,
                data=data,
                compression="gzip",
                compression_opts=6,
                shuffle=True,
            )

        # Truth lines (if available)
        if len(truth_by_iso) > 0:
            truth_grp = f.create_group("truth_lines")
            for iso, lists in truth_by_iso.items():
                if not lists:
                    continue
                # Concatenate and de-duplicate close energies (0.01 keV tolerance)
                e_all = np.concatenate([e for e, _ in lists])
                i_all = np.concatenate([i for _, i in lists])
                if e_all.size == 0:
                    continue
                order = np.argsort(e_all)
                e_all = e_all[order]
                i_all = i_all[order]
                # Merge close energies by taking max intensity
                merged_e = [e_all[0]]
                merged_i = [i_all[0]]
                tol = 0.01
                for e_val, i_val in zip(e_all[1:], i_all[1:]):
                    if abs(e_val - merged_e[-1]) <= tol:
                        merged_i[-1] = max(merged_i[-1], float(i_val))
                    else:
                        merged_e.append(float(e_val))
                        merged_i.append(float(i_val))
                g = truth_grp.create_group(iso)
                g.create_dataset("E_gamma_keV", data=np.asarray(merged_e, dtype=np.float32))
                g.create_dataset("I_per_decay", data=np.asarray(merged_i, dtype=np.float32))

        # Metadata and axes
        f.create_dataset("energy_axis", data=centers.astype(np.float32))
        f.attrs["n_isotopes"] = len(spectra_by_iso)
        f.attrs["n_bins"] = int(centers.size)
        f.attrs["energy_range_keV"] = [float(centers[0]), float(centers[-1])]
        f.attrs["smeared"] = bool(smear)
        f.attrs["use_true"] = bool(use_true)
        f.attrs["file_counts_json"] = str(dict(file_counts))

    # Summary on stdout
    print("Wrote:", output_path)
    total_spectra = sum(len(v) for v in spectra_by_iso.values())
    print(f"Isotopes: {len(spectra_by_iso)}, total spectra: {total_spectra}, bins: {centers.size}")
    for iso, spectra_list in spectra_by_iso.items():
        print(f"  {iso}: {len(spectra_list)} spectra")


def main():
    parser = argparse.ArgumentParser(description="Build HDF5 dataset from HPGe ROOT files")
    parser.add_argument("--root-dir", type=str, default="build/", help="Directory with ROOT files")
    parser.add_argument("--output", type=str, default="data/processed/hpge_dataset.h5", help="Output HDF5 path")
    parser.add_argument("--isotopes", nargs="*", help="Optional list of isotope labels to include")
    parser.add_argument("--no-smear", action="store_true", help="Disable Gaussian smearing")
    parser.add_argument("--use-true", action="store_true", help="Use Etrue_keV instead of Edep_keV if available")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for smearing")

    args = parser.parse_args()
    build_dataset(
        root_dir=args.root_dir,
        output_path=args.output,
        isotopes_filter=args.isotopes,
        smear=not args.no_smear,
        use_true=args.use_true,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
