# Repository Guidelines

## Project Structure & Module Organization
- `HPGeSingle.cc` wires together detector geometry, physics, primary generator, and user actions.
- `include/` holds user-defined Geant4 classes (`DetectorConstruction.hh`, `RunAction.hh`, etc.).
- `src/` provides the implementations for the headers, plus the custom `Run` analysis class.
- Macro scripts (`run.mac`, `init_vis.mac`, `vis.mac`, `gui.mac`) live at the repository root and are copied into `build/` during configuration.
- `build/` is the default out-of-source build directory created by CMake; avoid hand-editing files there.

## Build, Test, and Development Commands
- `cmake -S . -B build` configures the project against your Geant4 installation (requires Geant4 11.x and C++17).
- `cmake --build build` compiles `HPGeSingle` with warnings enabled (`-Wall -Wextra -pedantic`).
- `./build/HPGeSingle run.mac` runs a quiet 1 000-event batch simulation and writes the per-event spectrum to `gamma_spectrum.root`.
- `./build/HPGeSingle` (no macro) launches the interactive UI; use `init_vis.mac` and `gui.mac` for visualization.

## Coding Style & Naming Conventions
- Use C++17, 4-space indentation, and keep brace placement consistent with existing files (opening brace on the same line).
- Prefer Geant4 typedefs (`G4double`, `G4int`, `G4ThreeVector`) and constants from `G4SystemOfUnits`.
- Header/implementation pairs follow `ClassName.hh` and `ClassName.cc`; new classes should mirror this pattern in `include/` and `src/`.
- Keep logging via `G4cout`/`G4endl`; reserve `std::cout` for non-Geant4 contexts.

## Testing Guidelines
- No unit-test framework is integrated; validation is performed by running defined macros.
- For geometry or physics changes, run `./build/HPGeSingle run.mac`, review console summaries, and inspect the `EventData` ntuple inside `gamma_spectrum.root`.
- For visualization changes, start the interactive session and execute `/control/execute init_vis.mac` to confirm geometry renders correctly.
- Record the isotope configuration used for each test macro (e.g., `./build/HPGeSingle --isotope Cs137 run.mac`) so reviewers can reproduce results.

## Commit & Pull Request Guidelines
- Write concise, imperative commit titles (e.g., `Add Li-doped contact material`) and include context in the body when behavior changes.
- Reference related issues or task IDs in the body (`Refs: #42`) and summarize test evidence (`Test: ./build/HPGeSingle run.mac`).
- Pull requests should outline the motivation, list major code paths touched, attach relevant macro outputs or spectra, and call out any new dependencies or runtime switches.

## Context Management
**CRITICAL: Context is PROJECT-SPECIFIC based on current working directory**

### Context Location
- Always use `./.claude/docs/` relative to current working directory
- Each project maintains its own `.claude/` directory structure

## Workflow
### 1. Project Initialization (start of conversation)
- ALWAYS run `pwd` to confirm current working directory
- Verify `.claude/docs/` exists: `ls -la .claude/docs/`
- If missing, create: `mkdir -p .claude/docs/`
- Read `./.claude/docs/context.md` to understand project state
- If context.md missing, ask user for project overview and create it

### 2. Before Implementation
- Check if relevant `./.claude/docs/*-plan.md` exists
- Read ALL existing plans in `./.claude/docs/` for full picture
- Understand the full research before implementing

### 4. After Implementation
Update `./.claude/docs/context.md` with:
- **Completed**: Specific features/functions implemented
- **Current State**: Working/broken/partial/tested
- **Files Modified**: List of changed files
- **Dependencies**: New libraries/packages added
- **Issues**: Any problems encountered
- **TODO**: Remaining tasks from the plan
- **Last Updated**: Date and brief description
