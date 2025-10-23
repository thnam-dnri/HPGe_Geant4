// ==============================================================================
// DetectorConstruction.hh - Phase 2 Realistic HPGe Detector
// ==============================================================================

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // Get methods for analysis
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    
private:
    // Detector parameters
    static const G4double fSourceDetectorDistance;  // 10 cm
    static const G4double fWorldSize;              // 50 cm

    // Materials
    void DefineMaterials();
    G4Material* fWorldMaterial;
    G4Material* fGermanium;
    G4Material* fAluminum;
    G4Material* fVacuum;
    G4Material* fMylar;
    G4Material* fLithium;
    G4Material* fBoron;

    // Volumes
    G4LogicalVolume* fWorldLV;
    G4LogicalVolume* fScoringVolume;  // Points to Ge crystal
    G4VPhysicalVolume* fWorldPV;

    // Construction methods
    G4VPhysicalVolume* DefineVolumes();
};

#endif
