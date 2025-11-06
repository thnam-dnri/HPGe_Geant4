# HPGe Single Detector Simulation - Project Context

## Project Overview
Geant4-based simulation of a High-Purity Germanium (HPGe) detector for gamma-ray spectroscopy. This is Phase 1 of the HPGe ORTEC detector simulation project.

## Environment
- **Geant4 Version**: 11.3.2 (with multithreading support)
- **Geant4 Installation**: /home/nam/geant4-install/
- **CMake Version**: 3.28.3
- **Compiler**: GCC 13.3.0
- **Platform**: Linux 6.14.0-33-generic
- **C++ Standard**: C++17 (required by Geant4 11.x)

## Project Structure
```
/home/nam/Dropbox/HPGe_ORTEC/Single_Det/
├── CMakeLists.txt          # Build configuration
├── HPGeSingle.cc           # Main application
├── include/                # Header files
│   ├── DetectorConstruction.hh
│   ├── EventAction.hh
│   ├── PhysicsList.hh
│   ├── PrimaryGeneratorAction.hh
│   ├── RunAction.hh
│   ├── Run.hh
│   └── SteppingAction.hh
├── src/                    # Source files
│   ├── DetectorConstruction.cc
│   ├── EventAction.cc
│   ├── PhysicsList.cc
│   ├── PrimaryGeneratorAction.cc
│   ├── RunAction.cc
│   ├── Run.cc
│   └── SteppingAction.cc
├── build/                  # Build directory
├── *.mac                   # Geant4 macro files
└── .claude/docs/          # Project documentation
```

## Progress Update - 2025-10-24

- Implemented Phase 1 shield geometry (graded-Z):
  - Added 10 cm lead shield (G4_Pb) with 1 mm copper side liner (G4_Cu)
  - Added top lead cover (10 cm) with a 1 mm copper top liner (cap)
  - Added bottom lead cover (10 cm) with a 1 mm copper bottom liner (cap)
  - Introduced cylindrical nested volumes and visualization attributes
  - Ensured detector volumes are daughters of an air-filled shield cavity to avoid overlaps
- Verified build succeeds with Geant4 11.3.2
- Added quick geometry macro `geom_check.mac` for batch startup and geometry probe
- Updated source to a 2 cm × 1 mm disk source (uniform over area and thickness)
- Added CLI arg `--src-gap <value[mm|cm|m]>` for surface-to-surface distance between disk source and detector front window; `0cm` means just touching without overlap

### Current State
- **Build Status**: Compiles and links cleanly with shield geometry
- **Executable**: build/HPGeSingle present; runtime OK when `LD_LIBRARY_PATH` includes Geant4 libs
- **Geometry**: Lead + copper + cavity present; HPGe assembly placed inside shield cavity
- **World**: Expanded to 80 cm box to accommodate top/bottom covers
- **Known Warnings**: None new introduced by shield changes (compiler warnings remain enabled)

### Files Modified
- include/DetectorConstruction.hh — Added shield parameters/materials and logical volumes (added top cover volumes)
- src/DetectorConstruction.cc — Implemented lead + copper side liner and cavity; added top and bottom lead covers with copper liners; updated mothers for HPGe housing/vacuum; added vis attributes
- src/PrimaryGeneratorAction.cc — Implemented disk source sampling (2 cm diameter, 1 mm thickness)
- .claude/docs/context.md — Updated with latest progress
- Added: geom_check.mac — simple geometry startup/exit macro

### Build Instructions
```bash
# Clean build (recommended after CMakeLists.txt changes)
rm -rf build && mkdir build

# Configure
cmake -S . -B build

# Build
cmake --build build -j4

# Run (ensure Geant4 libs on runtime path)
LD_LIBRARY_PATH=/home/nam/geant4-install/lib:$LD_LIBRARY_PATH ./build/HPGeSingle [macro]

# Example (batch, geometry probe)
LD_LIBRARY_PATH=/home/nam/geant4-install/lib:$LD_LIBRARY_PATH ./build/HPGeSingle geom_check.mac

# Example (set disk source gap to 5 mm and run macro)
LD_LIBRARY_PATH=/home/nam/geant4-install/lib:$LD_LIBRARY_PATH ./build/HPGeSingle --src-gap 5mm run.mac
```

### Issues
- Runtime library path: needed to export `LD_LIBRARY_PATH` to locate Geant4 shared libs (e.g., `libG4ptl.so.3`)
- `/geometry/test/run` command not available in current app state; visualization route recommended for overlap checks

### TODO
 - Add `shielded_run.mac` with source/physics verbosity and optional EM deexcitation verification
 - Optionally adjust PrimaryGeneratorAction to position source inside shield (or control via macro)
 - Visualize geometry (`init_vis.mac`, `gui.mac`) to confirm shield alignment and absence of overlaps
 - Optional: region-specific production cuts for shield volumes if performance becomes an issue
 - Record and compare spectra with and without copper liner (Phase 3 validation)

## Progress Update - 2025-10-25

- Added isotope-driven primary generation using JSON files under `./isotope_data/`
  - New CLI: `-isotope <Symbol>` (e.g., `-isotope Cs137`) selects the source isotope
  - Follows decay chain until a stable daughter is reached
  - Emits only gamma lines from each chosen decay branch; non-gamma modes (beta/alpha) are not generated
  - Singles mode: emits one gamma per Geant4 event (prevents sum peaks) using an internal queue
  - Tracks counters: `Ndecays` (parent decays sampled) and `NgammaPrimaries` (gammas emitted)
  - Provides last-true-energy for QA and truth gamma lines aggregated along the chain
  - Falls back to prior behavior (RAINIER or Co-60 test) when `-isotope` is not given

- ROOT output per ML spec (via Geant4 analysis manager):
  - File name: `<Nuclide>_<DetID>_<GapMM>_<runid>.root` (if `-isotope` used), else `gamma_spectrum.root`
  - Ntuple `Events`: columns `Edep_keV`, `Etrue_keV`
  - Ntuple `RunInfo`: strings (Nuclide, Generator, G4Version, PhysicsLists, DetID, Geometry, Binning, SourceJSON, SourceHash) and ints (Ndecays, NgammaPrimaries)
  - Ntuple `Axes_energy_centers`: column `energy_keV` with one row per bin center (0–3000 keV, 8192 bins)
  - Ntuple `Truth_gamma_lines`: columns `E_gamma_keV`, `I_per_decay` (aggregated along decay chain)
  - Note: uses ntuples instead of ROOT directories/TObjString for compatibility without linking ROOT directly

### Current State
- Build passes; `./build/HPGeSingle -isotope Cs137 run.mac` runs and uses `isotope_data/Cs137.json`
- Gamma sampling per branch uses `absolute_intensity` as an emission probability (or mean multiplicity if >1)
- Chain traversal honors `branching_ratio` to select daughter; stops at stable isotope

### Files Modified
- HPGeSingle.cc — Parse `-isotope` CLI and pass to generator
- HPGeSingle.cc — Added `--det-id` CLI, passed to RunAction for file naming
- include/PrimaryGeneratorAction.hh — Added isotope symbol setter and generation hook
- src/PrimaryGeneratorAction.cc — Implemented singles gamma emission (queue), counters, truth lines
- include/IsotopeDecay.hh — New minimal JSON loader (no external deps)
- src/IsotopeDecay.cc — Parser implementation for `decay_modes[].gammas[]`
- include/RunAction.hh, src/RunAction.cc — New multi-ntuple ROOT layout (Events, RunInfo, Axes, Truth), dynamic filename
- src/EventAction.cc — Fill Events ntuple (Edep_keV, Etrue_keV)

### Dependencies
- No new external libraries; minimal JSON parsing implemented in-house

### Issues
- Coincidence/cascade correlations (`coincident_gammas`) are not modeled; each listed gamma is sampled independently by intensity
- Minimal JSON parser assumes the current file structure; significant schema changes may require parser updates
 - `Axes_energy_centers` stored as one value per row (vector column not used); Python should read and stack the rows
 - `G4Version` stored as "unknown" placeholder (can be improved if needed)

### TODO
- Optionally model gamma-gamma coincidences using `coincident_gammas` correlations for cascades (e.g., Co-60 1173–1332 keV pair)
- Add CLI to control activity or per-event multiplicity behavior
- Add a quick validation macro demonstrating several isotopes (Cs-137, Co-60) and expected peak lists
 - Optionally switch to direct ROOT (TFile/TTree) if available to mirror directory layout exactly

### Last Updated
2025-10-25 — Added `-isotope` gamma-chain singles generator and ML-oriented ROOT output; build verified

## Progress Update - 2025-10-29

- Updated visualization styling around the HPGe crystal
  - Restored the original semi-transparent solids for the lead shield, copper liner, housing, window stack, and dead layers
  - Outer aluminum cup now renders in wireframe so the crystal interface stays visible without hiding surrounding parts
- Corrected bore hole orientation so it now extends from the rear (far) face of the germanium crystal
- Added optional interactive mode for visualization macros
  - Running `./build/HPGeSingle init_vis.mac` now keeps the UI session alive (mirrors launching without a macro)
  - New `--interactive`/`--ui` flag can force the UI to stay open after any macro executes

### Current State
- Visualization macros should now display the crystal inside the shield without hiding external geometry
- Build/run behavior unchanged from prior update; no new warnings introduced

### Files Modified
- HPGeSingle.cc — Allow `init_vis.mac` (or `--interactive`) to hand off to a `G4UIExecutive` session after executing macros, keeping visualization open
- src/DetectorConstruction.cc — Restored shield/housing transparency, made only the outer aluminum cup wireframe, and repositioned the crystal bore hole to exit the back face

### Dependencies
- No changes

### Issues
- Interactive visualization not re-run in this session; manual confirmation still recommended

### TODO
- Launch `./build/HPGeSingle` with `init_vis.mac` (or `--interactive`) to verify the updated cup wireframe looks correct in your environment

### Last Updated
2025-10-29 — Adjusted detector visualization to keep the germanium crystal visible during rendering

## Progress Update - 2025-10-30

- Removed legacy RAINIER integration and standardized on isotope JSON primaries
  - `PrimaryGeneratorAction` now defaults to Co-60 singles sourced from `isotope_data` and no longer parses external spectrum files
  - CLI defaults to the same Co-60 dataset while `--isotope` or a bare nuclide argument can override it
  - Entry point banner and logs now reflect the isotope-driven workflow
- ROOT output now targets the `training_data/` directory; the folder is created automatically if missing
- Added a GCC 8 compatibility link step (`stdc++fs`) in `CMakeLists.txt` so the project builds on the remote saho-a host (CMake 3.20.6, GCC 8.3.1, Geant4 11.3.2)
- Switched all `target_link_libraries` calls in CMake to keyword form to satisfy older CMake policies on saho-a
- Updated run metadata strings and repository guidelines to drop obsolete RAINIER references

### Current State
- Executable builds cleanly (`cmake --build build`) with the simplified generator
- Simulations require a valid `isotope_data/<Nuclide>.json`; default Co-60 path confirmed during initialization

### Files Modified
- HPGeSingle.cc — removed `--rainier` handling, defaulted to Co-60 isotope, refreshed CLI messaging
- include/PrimaryGeneratorAction.hh, src/PrimaryGeneratorAction.cc — excised RAINIER data structures, reset logic, and added isotope setter reset hooks
- src/RunAction.cc — metadata generator label now reflects isotope JSON singles
- AGENTS.md — testing guideline now records isotope selections instead of RAINIER inputs

### Dependencies
- No changes

### Issues
- None observed; ensure required isotope JSON files remain in `./isotope_data/`

### TODO
- Consider exposing CLI to select default isotope list via configuration file if more flexibility is needed

### Last Updated
2025-10-30 — Removed RAINIER pathways and enforced isotope-based primary generation

## Research Notes
- 2025-10-24 Lead Shield Architecture: 10 cm Pb + 1 mm Cu liner recommended for soil analysis (15x background reduction, 50% Pb X-ray suppression)
  Key insight: Graded shield design suppresses 85 keV lead fluorescence using copper liner; nested cylindrical geometry optimal for implementation
  See: lead-shield-plan.md for full architectural plan with material definitions, geometry design, physics configuration, and 3-phase implementation roadmap
- 2025-10-25 ML Data Pipeline Architecture: Python pipeline using uproot (5-10x faster I/O), HDF5 storage (3-4x compression), Gaussian energy smearing (FWHM = √(10.26 + 0.00168·E) keV)
  Key insight: Add time-stamped events, spatial coordinates, and activity parameters to Geant4 output for temporal analysis, position-sensitive efficiency, and activity quantification ML tasks
  See: ml-data-pipeline-plan.md for complete architecture with ROOT processing, energy smearing, feature engineering (peak finding, spectral features), 4 ML tasks (isotope ID, activity quant, calibration QC, anomaly detection), dataset schemas (HDF5/Parquet), 6-week implementation roadmap, and recommended Geant4 output additions

### Dependencies
- Geant4 11.0 or higher (with ui_all and vis_all components)
- OpenGL (for visualization)
- X11 (for GUI)
- Threads (for multithreading support)
- No new dependencies introduced in the latest update

### Last Updated
2025-10-24 - Implemented lead + copper shield geometry; build verified; documented runtime path requirement

## Progress Update - 2025-10-25 (ML Pipeline Phase 1)

- Created Python virtual environment at `.venv` and installed ML/data packages
  - numpy, h5py, uproot, awkward, pandas, scikit-learn
- Implemented Phase 1 of ML data pipeline:
  - `ml_pipeline/io.py` — uproot-based readers for `Events`, `RunInfo`, and `Axes_energy_centers`; isotope label inference; axis utilities
  - `ml_pipeline/smearing.py` — HPGe resolution model and Gaussian smearing using FWHM(E) = sqrt(10.26 + 0.00168·E)
  - `ml_pipeline/build_dataset.py` — CLI to aggregate ROOT files into HDF5 grouped by isotope with optional smearing
  - `requirements.txt` — pinned core Python deps for reproducibility

### Current State
- Virtual environment active at `.venv` with required packages installed
- Dataset builder runs as: `./.venv/bin/python -m ml_pipeline.build_dataset --root-dir build --output data/processed/hpge_dataset.h5`
- Uses per-file energy axis when available, else defaults to 0–3000 keV with 8192 bins
- Groups spectra by isotope label from `RunInfo.Nuclide` or filename prefix

### Files Modified / Added
- Added: `ml_pipeline/__init__.py`
- Added: `ml_pipeline/io.py`
- Added: `ml_pipeline/smearing.py`
- Added: `ml_pipeline/build_dataset.py`
- Added: `requirements.txt`
- Updated: `.claude/docs/context.md`

### Dependencies
- New Python runtime deps (installed into `.venv`): numpy, h5py, uproot, awkward, pandas, scikit-learn
- No changes to C++/Geant4 dependencies

### Issues
- None observed during installation; relies on internet access for pip during initial setup
- If ROOT files lack `Axes_energy_centers`, a default uniform axis is used

### TODO
- Phase 2: Feature engineering and optional training scripts (e.g., RandomForest classifier) per plan
- Add CLI switches for custom energy ranges/binning and per-isotope filtering presets
- Optionally persist per-file `RunInfo` metadata into the HDF5 for traceability

### Last Updated
2025-10-25 — Implemented ML pipeline Phase 1 (I/O, smearing, dataset builder); venv configured

## Progress Update - 2025-10-25 (ML Pipeline Phase 2)

- Extended dataset builder to embed truth gamma lines:
  - `truth_lines/<isotope>/E_gamma_keV`, `truth_lines/<isotope>/I_per_decay` aggregated across files with 0.01 keV de-duplication
- Implemented feature extraction modules:
  - `ml_pipeline/features/peak_finder.py` — Savitzky–Golay smoothing + prominence-based peak detection; centroid, area, FWHM estimates
  - `ml_pipeline/features/spectral_features.py` — global stats, energy band fractions, top-K peak features, area ratios; optional truth-match ratio
- Added CLI to compute and persist features:
  - `ml_pipeline/build_features.py` writes features under `features/phase2/<isotope>` with metadata in `features_meta/phase2/feature_names`
- Updated `requirements.txt` to include `scipy`

### Current State
- Phase 2 features can be generated from the HDF5 dataset built in Phase 1
- Truth lines are available in HDF5 when present in ROOT inputs and are used to compute a simple truth-match ratio feature

### Usage
- Build spectra HDF5 (Phase 1):
  - `./.venv/bin/python -m ml_pipeline.build_dataset --root-dir build --output data/processed/hpge_dataset.h5`
- Compute Phase 2 features and write to HDF5:
  - `./.venv/bin/python -m ml_pipeline.build_features --dataset data/processed/hpge_dataset.h5 --feature-set phase2 --k-peaks 6`

### Files Modified / Added
- Updated: `ml_pipeline/build_dataset.py` — writes aggregated truth lines
- Added: `ml_pipeline/features/__init__.py`
- Added: `ml_pipeline/features/peak_finder.py`
- Added: `ml_pipeline/features/spectral_features.py`
- Added: `ml_pipeline/build_features.py`
- Updated: `requirements.txt`
- Updated: `.claude/docs/context.md`

### Issues
- Peak baseline subtraction is simple (median-based); can be refined with local sideband fits
- Truth-match ratio uses nearest-peak within 2×FWHM heuristic; intensity weighting not yet applied

### TODO
- Phase 3: Baseline models and validation notebooks per plan
- Optional: store per-file RunInfo into HDF5 for better traceability
- Optional: flatten per-isotope features into a single X/y table for scikit-learn convenience

### Last Updated
2025-10-25 — Implemented ML pipeline Phase 2 (truth lines, peak features, feature CLI)

## Progress Update - 2025-10-30 (Lead Shield Toggle)

- **Completed**: Added a `constexpr bool kEnableLeadShield` switch at the top of `src/DetectorConstruction.cc` that gates construction/visualization of the entire lead + copper shield stack (ring, covers, liners, cavity) and falls back to placing the detector directly in the world when disabled.
- **Current State**: Shield geometry is enabled by default; flip the constant to `false` before rebuilding to run an unshielded configuration. Build/tests not rerun after this refactor.
- **Files Modified**: `src/DetectorConstruction.cc`
- **Dependencies**: None added or removed.
- **Issues**: Remember to rebuild and rerun macros after changing the toggle to ensure geometry overlaps are still absent.
- **TODO**: Consider adding a runtime CLI or macro-driven control if frequent toggling is needed without recompilation.
- **Last Updated**: 2025-10-30 — Added compile-time switch for shield geometry.

## Progress Update - 2025-11-04 (Isotope emission mode)

- Completed: Added runtime control for isotope gamma emission: `--isotope-mode parent-only|full-chain`. Parent-only emits only the selected isotope’s immediate gamma lines (no daughter tracking). Full-chain preserves prior behavior.
- Default: Parent-only, to better match sources like Am-241 (e.g., 59.5 keV line; long-lived daughters suppressed).
- Files Modified: `include/PrimaryGeneratorAction.hh`, `src/PrimaryGeneratorAction.cc`, `HPGeSingle.cc`.
- Current State: Truth gamma lines now respect the selected mode: parent-only aggregates parent lines only; full-chain traverses to stability.
- Issues: For isotopes where main peaks originate from short-lived daughters (e.g., Cs-137 → Ba-137m), use `--isotope-mode full-chain` to include those lines.
- TODO: Optionally add a prompt-only mode with daughter half-life cutoff for more realism without full chains.
- Last Updated: 2025-11-04 — Introduced `--isotope-mode` flag and parent-only default.

## Progress Update - 2025-11-04 (Emitter-aware truth export)

- Completed (C++): Truth ntuple now includes emitter nuclide and half-life for each gamma line.
  - Added `IsotopeInfo.half_life_seconds` in loader and parsed from JSON.
  - Added `PrimaryGeneratorAction::GetTruthGammaLinesDetailed()` accumulating `(E, I, emitter, hl_s)` and respecting `--isotope-mode`.
  - Extended `Truth_gamma_lines` tree with `Emitter` and `EmitterHalfLife_s` columns.
- Completed (Python): ML pipeline reads detailed truth and writes to HDF5 under `truth_lines/<iso>/`:
  - Datasets: `E_gamma_keV`, `I_per_decay`, `Emitter` (bytes), `EmitterHalfLife_s`.
  - For near-duplicate energies, keeps the record with higher intensity and its emitter label.
- Files Modified: `include/IsotopeDecay.hh`, `src/IsotopeDecay.cc`, `include/PrimaryGeneratorAction.hh`, `src/PrimaryGeneratorAction.cc`, `src/RunAction.cc`, `ml_pipeline/io.py`, `ml_pipeline/build_dataset.py`.
- TODO: Optional prompt-only mode with daughter half-life cutoff; consider adding an `IsotopeMode` string to RunInfo for traceability.
- Last Updated: 2025-11-04 — Emitter-aware truth plumbed through ROOT → HDF5.

## Progress Update - 2025-11-04 (Isotope speed optimization)

- Completed: Eliminated slow per-event JSON loads and retry loops when using `--isotope`.
  - Added in-memory cache to `IsotopeDataLoader` to parse each `isotope_data/<Nuclide>.json` only once.
  - Precompute a "singles" sampling distribution from truth gamma lines and sample one gamma per event using `std::discrete_distribution`.
  - Default event generation now performs zero JSON I/O and no decay retries during the event loop.
- Current State: Simulations with low-intensity parents run at similar speed to Co-60; full-chain remains available and is now precomputed once per run as well.
- Files Modified:
  - `include/IsotopeDecay.hh` — added cache members and mutex.
  - `src/IsotopeDecay.cc` — implemented cache lookup/store in `Load()`.
  - `include/PrimaryGeneratorAction.hh` — added singles distribution members and `RebuildSinglesDistribution()`; ensured mode switch invalidates singles.
  - `src/PrimaryGeneratorAction.cc` — switched `GeneratePrimaries` to use precomputed singles distribution; preserved legacy path as fallback.
- Dependencies: No new external libraries; uses standard `<unordered_map>`, `<mutex>`, and `<random>`.
- Issues: `Ndecays` counter semantics now reflect one sampled gamma per event rather than simulated parent decay attempts; metadata remains compatible but values change (closer to NgammaPrimaries).
- TODO: Optionally expose a flag to force legacy decay-queue behavior for comparison/validation; optionally record `IsotopeMode` into `RunInfo`.
- Last Updated: 2025-11-04 — Cached loader + precomputed singles sampler for fast `--isotope` runs.

## Progress Update - 2025-11-05 (Run progress verbosity)

- Completed: Added periodic run progress output every 100,000 events.
  - Set `G4RunManager::SetPrintProgress(100000)` after constructing the run manager in `HPGeSingle.cc`.
- Current State: In batch runs (e.g., `run.mac` with 1,000,000 events), the console prints a progress line at each 100k-event boundary. No impact on analysis or output files.
- Files Modified: `HPGeSingle.cc`.
- Dependencies: No changes.
- Issues: None observed.
- TODO: Optionally add a CLI flag (e.g., `--print-progress <N>`) to make the interval configurable.
- Last Updated: 2025-11-05 — Enabled periodic progress printing.
## Progress Update - 2025-11-06 (Gamma-less isotope guard)

- Completed: End run and suppress output when the selected isotope produces no gamma lines under the current `--isotope-mode`.
  - In `HPGeSingle.cc`, after setting the isotope and mode, the app queries truth gamma lines; if empty, it exits immediately without executing macros or creating any ROOT file.
  - In `RunAction::BeginOfRunAction`, added a second guard that aborts the run and returns before `OpenFile()` if truth gamma lines are empty (defensive in case a `beamOn` is triggered anyway).
- Current State: For gamma-less isotopes, no `training_data/*.root` is created; for others, behavior unchanged.
- Files Modified: `HPGeSingle.cc`, `include/RunAction.hh`, `src/RunAction.cc`.
- Dependencies: No changes.
- Issues: Build not executed in this session; please rebuild to validate.
- TODO: Optionally hint users to try `--isotope-mode full-chain` when parent-only yields no lines.
- Last Updated: 2025-11-06 — Guard added to prevent output for gamma-less isotopes.
