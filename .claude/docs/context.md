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

## Progress Update - 2025-10-23

### Completed
- Fixed CMake cache conflict from previous macOS installation path
- Updated CMakeLists.txt with modern best practices for Geant4 11.x:
  - Explicitly set C++17 standard requirement
  - Added project metadata (VERSION, DESCRIPTION, LANGUAGES)
  - Minimum Geant4 version check (11.0)
  - Modern target-based configuration using target_include_directories
  - Added compiler warnings (-Wall -Wextra -pedantic) for better code quality
  - Status messages showing Geant4 version and multithreading support

### Current State
- **Build Status**: Successfully compiling and linking
- **Executable**: /home/nam/Dropbox/HPGe_ORTEC/Single_Det/build/HPGeSingle (329KB)
- **Known Warnings**:
  - Unused variables in DetectorConstruction.cc:122-125, 147
  - These are non-critical and can be addressed in future cleanup

### Files Modified
- CMakeLists.txt - Updated with modern CMake practices
- .claude/docs/context.md - Created project documentation

### Build Instructions
```bash
# Clean build (recommended after CMakeLists.txt changes)
rm -rf build && mkdir build

# Configure
cmake -S . -B build

# Build
cmake --build build -j4

# Run
./build/HPGeSingle [options]
```

### Issues Resolved
- **CMake cache mismatch**: Removed stale CMakeCache.txt from old path (/Users/namtran/OneDrive/DNRI/GEANT4/HPGe_ORTEC/SingleDec/)
- **Build in source directory**: Cleaned CMake files accidentally generated in project root

### TODO
- Address unused variable warnings in DetectorConstruction.cc
- Consider enabling multithreading in the main application (currently available but not utilized)
- Add unit tests
- Implement visualization macros

### Dependencies
- Geant4 11.0 or higher (with ui_all and vis_all components)
- OpenGL (for visualization)
- X11 (for GUI)
- Threads (for multithreading support)

### Last Updated
2025-10-23 - Fixed CMake configuration and updated build system for Geant4 11.3.2
