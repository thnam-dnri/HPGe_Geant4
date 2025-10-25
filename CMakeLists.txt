# ==============================================================================
# CMakeLists.txt for HPGe Simulation - Phase 1
# ==============================================================================

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# Project definition
project(HPGeSingle
    VERSION 1.0
    DESCRIPTION "HPGe Detector Simulation with Geant4"
    LANGUAGES CXX
)

# Require C++17 standard (required by Geant4 11.x)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find required packages
find_package(Geant4 11.0 REQUIRED ui_all vis_all)

# Report Geant4 version and configuration
message(STATUS "Geant4 version: ${Geant4_VERSION}")
message(STATUS "Geant4 multithreading: ${Geant4_multithreaded_FOUND}")

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable
add_executable(HPGeSingle HPGeSingle.cc ${sources} ${headers})

# Modern target-based configuration
target_include_directories(HPGeSingle
    PRIVATE ${PROJECT_SOURCE_DIR}/include
)

# Link Geant4 libraries
target_link_libraries(HPGeSingle ${Geant4_LIBRARIES})

# Provide absolute isotope data directory to the code for robust loading
target_compile_definitions(HPGeSingle PRIVATE ISOTOPE_DATA_DIR="${PROJECT_SOURCE_DIR}/isotope_data")

# Set compiler warnings for better code quality
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(HPGeSingle PRIVATE -Wall -Wextra -pedantic)
endif()

# Copy all scripts to the build directory
set(GAMMA_CASCADE_SCRIPTS
    init_vis.mac
    vis.mac
    run.mac
    gui.mac
    )

foreach(_script ${GAMMA_CASCADE_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY
        )
endforeach()

# For internal Geant4 use - but has no effect if you build this standalone
add_custom_target(HPGeSingleExe DEPENDS HPGeSingle)

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS HPGeSingle DESTINATION bin)
