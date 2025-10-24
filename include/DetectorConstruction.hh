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
    G4double GetDetectorSurfaceZ() const { return fDetectorSurfaceZ; }
    
private:
    // Detector parameters
    static const G4double fSourceDetectorDistance;  // 10 cm
    static const G4double fWorldSize;              // 80 cm

    // Shield parameters
    static const G4double fShieldInnerRadius;      // 50 mm cavity radius
    static const G4double fShieldThickness;        // 100 mm Pb thickness
    static const G4double fShieldOuterRadius;      // 150 mm outer radius
    static const G4double fShieldHeight;           // 350 mm height
    static const G4double fCopperLinerThickness;   // 1 mm Cu liner
    static const G4double fShieldCenterZ;          // Shield center Z (relative to world)

    // Materials
    void DefineMaterials();
    G4Material* fWorldMaterial;
    G4Material* fGermanium;
    G4Material* fAluminum;
    G4Material* fVacuum;
    G4Material* fMylar;
    G4Material* fLithium;
    G4Material* fBoron;
    G4Material* fLead;
    G4Material* fCopper;

    // Volumes
    G4LogicalVolume* fWorldLV;
    G4LogicalVolume* fScoringVolume;  // Points to Ge crystal
    G4VPhysicalVolume* fWorldPV;

    // Shield volumes
    G4LogicalVolume* fLeadShieldLV;
    G4LogicalVolume* fCopperLinerLV;
    G4LogicalVolume* fShieldCavityLV;
    G4LogicalVolume* fLeadTopLV;
    G4LogicalVolume* fCopperTopLV;
    G4LogicalVolume* fLeadBottomLV;
    G4LogicalVolume* fCopperBottomLV;

    // Construction methods
    G4VPhysicalVolume* DefineVolumes();

    // Cached geometry query points
    G4double fDetectorSurfaceZ {0.0};
};

#endif
