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

## Research Notes
- 2025-10-24 Lead Shield Architecture: 10 cm Pb + 1 mm Cu liner recommended for soil analysis (15x background reduction, 50% Pb X-ray suppression)
  Key insight: Graded shield design suppresses 85 keV lead fluorescence using copper liner; nested cylindrical geometry optimal for implementation
  See: lead-shield-plan.md for full architectural plan with material definitions, geometry design, physics configuration, and 3-phase implementation roadmap

### Dependencies
- Geant4 11.0 or higher (with ui_all and vis_all components)
- OpenGL (for visualization)
- X11 (for GUI)
- Threads (for multithreading support)
- No new dependencies introduced in the latest update

### Last Updated
2025-10-24 - Implemented lead + copper shield geometry; build verified; documented runtime path requirement
