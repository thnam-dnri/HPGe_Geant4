"""
ML data pipeline utilities for HPGe Geant4 output.

Modules:
- io: ROOT readers and metadata utilities
- smearing: Energy resolution models and Gaussian smearing
- build_dataset: CLI to convert ROOT files to HDF5 datasets
"""

__all__ = [
    "io",
    "smearing",
]

