// ==============================================================================
// DetectorConstruction.cc - Phase 2 Realistic HPGe Detector
// ==============================================================================

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <algorithm>

// Static member definitions
const G4double DetectorConstruction::fWorldSize = 80.0*cm; // expanded to accommodate top/bottom covers
const G4double DetectorConstruction::fSourceDetectorDistance = 10.0*cm;

// Lead shield geometry parameters
const G4double DetectorConstruction::fShieldInnerRadius = 50.0*mm;      // cavity radius
const G4double DetectorConstruction::fShieldThickness   = 100.0*mm;     // Pb thickness
const G4double DetectorConstruction::fShieldOuterRadius = 150.0*mm;     // inner + thickness
const G4double DetectorConstruction::fShieldHeight      = 350.0*mm;     // total height
const G4double DetectorConstruction::fCopperLinerThickness = 1.0*mm;    // Cu liner
const G4double DetectorConstruction::fShieldCenterZ     = 0.0*mm;       // center at world origin

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fWorldMaterial(nullptr), fGermanium(nullptr), fAluminum(nullptr), 
  fVacuum(nullptr), fMylar(nullptr), fLithium(nullptr), fBoron(nullptr),
  fLead(nullptr), fCopper(nullptr),
  fWorldLV(nullptr), fScoringVolume(nullptr), fWorldPV(nullptr),
  fLeadShieldLV(nullptr), fCopperLinerLV(nullptr), fShieldCavityLV(nullptr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();
    
    // Air for world
    fWorldMaterial = nist->FindOrBuildMaterial("G4_AIR");
    
    // Vacuum for gaps and housing interior
    fVacuum = nist->FindOrBuildMaterial("G4_Galactic");
    
    // Aluminum for housing and windows
    fAluminum = nist->FindOrBuildMaterial("G4_Al");
    
    // High purity germanium
    fGermanium = nist->FindOrBuildMaterial("G4_Ge");
    
    // Mylar for window layer
    fMylar = nist->FindOrBuildMaterial("G4_MYLAR");
    
    // Create realistic dead layer materials (Ge with dopants)
    // Lithium-doped germanium (n+ contact, outer)
    fLithium = new G4Material("Li_doped_Ge", 5.32*g/cm3, 2);
    fLithium->AddMaterial(fGermanium, 99.9*perCent);
    fLithium->AddMaterial(nist->FindOrBuildMaterial("G4_Li"), 0.1*perCent);
    
    // Boron-doped germanium (p+ contact, inner)
    fBoron = new G4Material("B_doped_Ge", 5.32*g/cm3, 2);
    fBoron->AddMaterial(fGermanium, 99.9*perCent);
    fBoron->AddMaterial(nist->FindOrBuildMaterial("G4_B"), 0.1*perCent);

    // Shield materials
    fLead   = nist->FindOrBuildMaterial("G4_Pb");
    fCopper = nist->FindOrBuildMaterial("G4_Cu");

    G4cout << "Materials defined for realistic HPGe detector" << G4endl;
    G4cout << "Dead layers use doped germanium rather than pure dopants" << G4endl;
    G4cout << "Shield materials: Pb density=" << fLead->GetDensity()/(g/cm3)
           << " g/cm^3, Cu density=" << fCopper->GetDensity()/(g/cm3) << " g/cm^3" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    // World volume
    G4Box* worldS = new G4Box("World", fWorldSize/2, fWorldSize/2, fWorldSize/2);
    fWorldLV = new G4LogicalVolume(worldS, fWorldMaterial, "World");
    fWorldPV = new G4PVPlacement(0, G4ThreeVector(), fWorldLV, "World", 0, false, 0, true);

    // ---------------------------------------------------------------------
    // Lead shield + copper liner + cavity (nested cylindrical volumes)
    // ---------------------------------------------------------------------
    // 1) Lead shield (ring)
    G4double leadInnerR = fShieldInnerRadius;
    G4double leadOuterR = fShieldOuterRadius;
    G4double leadHeight = fShieldHeight;

    G4Tubs* leadShieldS = new G4Tubs("LeadShield",
                                     leadInnerR, leadOuterR,
                                     leadHeight/2,
                                     0*deg, 360*deg);
    fLeadShieldLV = new G4LogicalVolume(leadShieldS, fLead, "LeadShield");
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, fShieldCenterZ),
                      fLeadShieldLV,
                      "LeadShield",
                      fWorldLV,
                      false, 0, true);

    // 2) Copper liner (thin ring, slightly shorter to avoid coincident surfaces)
    G4double cuOuterR = leadInnerR;
    G4double cuInnerR = std::max(0.0, (cuOuterR - fCopperLinerThickness));
    G4double cuHeight = std::max(0.0, (leadHeight - 2.0*mm));

    G4Tubs* copperLinerS = new G4Tubs("CopperLiner",
                                      cuInnerR, cuOuterR,
                                      cuHeight/2,
                                      0*deg, 360*deg);
    fCopperLinerLV = new G4LogicalVolume(copperLinerS, fCopper, "CopperLiner");
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, fShieldCenterZ),
                      fCopperLinerLV,
                      "CopperLiner",
                      fWorldLV,
                      false, 0, true);

    // 3) Shield cavity (air-filled space for detector and source)
    G4double cavityR = cuInnerR;
    G4double cavityH = leadHeight; // full height cavity
    G4Tubs* cavityS = new G4Tubs("ShieldCavity",
                                 0, cavityR,
                                 cavityH/2,
                                 0*deg, 360*deg);
    fShieldCavityLV = new G4LogicalVolume(cavityS, fWorldMaterial, "ShieldCavity");
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, fShieldCenterZ),
                      fShieldCavityLV,
                      "ShieldCavity",
                      fWorldLV,
                      false, 0, true);

    // 4) Top lead cover (same thickness as lead walls), sits above the shield
    G4double topThick = fShieldThickness;
    G4Tubs* leadTopS = new G4Tubs("LeadTopCover",
                                  0, leadOuterR,
                                  topThick/2,
                                  0*deg, 360*deg);
    fLeadTopLV = new G4LogicalVolume(leadTopS, fLead, "LeadTopCover");
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, fShieldCenterZ + (leadHeight/2 + topThick/2)),
                      fLeadTopLV,
                      "LeadTopCover",
                      fWorldLV,
                      false, 0, true);

    // 5) Top copper liner (thin), mounted inside cavity just below the top plane
    G4double cuTopThick = fCopperLinerThickness;
    G4Tubs* cuTopS = new G4Tubs("CopperTopLiner",
                                0, cavityR,
                                cuTopThick/2,
                                0*deg, 360*deg);
    fCopperTopLV = new G4LogicalVolume(cuTopS, fCopper, "CopperTopLiner");
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, fShieldCenterZ + (cavityH/2 - cuTopThick/2)),
                      fCopperTopLV,
                      "CopperTopLiner",
                      fShieldCavityLV,
                      false, 0, true);

    // 6) Bottom lead cover (same thickness), sits below the shield
    G4Tubs* leadBottomS = new G4Tubs("LeadBottomCover",
                                     0, leadOuterR,
                                     topThick/2,
                                     0*deg, 360*deg);
    fLeadBottomLV = new G4LogicalVolume(leadBottomS, fLead, "LeadBottomCover");
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, fShieldCenterZ - (leadHeight/2 + topThick/2)),
                      fLeadBottomLV,
                      "LeadBottomCover",
                      fWorldLV,
                      false, 0, true);

    // 7) Bottom copper liner (thin), mounted inside cavity near bottom plane
    G4double cuBottomThick = fCopperLinerThickness;
    G4Tubs* cuBottomS = new G4Tubs("CopperBottomLiner",
                                   0, cavityR,
                                   cuBottomThick/2,
                                   0*deg, 360*deg);
    fCopperBottomLV = new G4LogicalVolume(cuBottomS, fCopper, "CopperBottomLiner");
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, fShieldCenterZ - (cavityH/2 - cuBottomThick/2)),
                      fCopperBottomLV,
                      "CopperBottomLiner",
                      fShieldCavityLV,
                      false, 0, true);

    // Detector parameters (convert to G4 units)
    G4double housingInnerDiam = 64.46*mm;
    G4double housingOuterDiam = 67.0*mm;
    G4double housingLength = 76.0*mm;
    
    G4double windowDiam = 64.46*mm;
    G4double alWindowThick = 1.27*mm;
    G4double mylarThick = 0.025*mm;
    G4double alFoilThick = 0.025*mm;
    G4double alCupThick = 0.5*mm;
    
    G4double geCrystalDiam = 57.6*mm;
    G4double geCrystalLength = 66.8*mm;
    G4double boreHoleDiam = 10.5*mm;
    G4double boreHoleDepth = 53.5*mm;
    
    G4double liDeadThick = 0.7*mm;
    G4double bDeadThick = 0.3*micrometer;
    
    G4double windowGap = 3.0*mm;
    
    // Calculate positions
    G4double housingZ = fSourceDetectorDistance + housingLength/2;
    
    G4double windowStartZ = fSourceDetectorDistance;
    fDetectorSurfaceZ = windowStartZ; // cache for primary generator positioning
    
    G4double totalWindowThick = alWindowThick + mylarThick + alFoilThick + alCupThick;
    G4double geZ = windowStartZ + totalWindowThick + windowGap + geCrystalLength/2;
    
    G4cout << "\n=== HPGe Detector Geometry ===" << G4endl;
    G4cout << "Source position: (0, 0, 0)" << G4endl;
    G4cout << "Window front face at z = " << windowStartZ/mm << " mm" << G4endl;
    G4cout << "Ge crystal center at z = " << geZ/mm << " mm" << G4endl;
    G4cout << "Housing center at z = " << housingZ/mm << " mm" << G4endl;

    // 1. Aluminum Housing  
    G4Tubs* housingSolid = new G4Tubs("Housing", 
                                     housingInnerDiam/2, housingOuterDiam/2, 
                                     housingLength/2, 0*deg, 360*deg);
    G4LogicalVolume* housingLV = new G4LogicalVolume(housingSolid, fAluminum, "Housing");
    new G4PVPlacement(0, G4ThreeVector(0, 0, housingZ), housingLV, "Housing", fShieldCavityLV, false, 0, true);

    // 1b. Vacuum inside housing (mother volume for detector components)
    G4Tubs* vacuumSolid = new G4Tubs("VacuumInside", 0, housingInnerDiam/2, 
                                    housingLength/2, 0*deg, 360*deg);
    G4LogicalVolume* vacuumLV = new G4LogicalVolume(vacuumSolid, fVacuum, "VacuumInside");
    new G4PVPlacement(0, G4ThreeVector(0, 0, housingZ), vacuumLV, "VacuumInside", fShieldCavityLV, false, 0, true);

    // Calculate relative positions within housing
    G4double relativeGeZ = geZ - housingZ;  // Ge position relative to housing center
    G4double relativeWindowStartZ = windowStartZ - housingZ;  // Window start relative to housing center
    
    // 2. Front Window Assembly (inside vacuum volume)
    // a. Outer Aluminum Window
    G4Tubs* alWindowSolid = new G4Tubs("AlWindow", 0, windowDiam/2, alWindowThick/2, 0*deg, 360*deg);
    G4LogicalVolume* alWindowLV = new G4LogicalVolume(alWindowSolid, fAluminum, "AlWindow");
    new G4PVPlacement(0, G4ThreeVector(0, 0, relativeWindowStartZ + alWindowThick/2), alWindowLV, "AlWindow", vacuumLV, false, 0, true);

    // b. Mylar Layer
    G4Tubs* mylarSolid = new G4Tubs("Mylar", 0, windowDiam/2, mylarThick/2, 0*deg, 360*deg);
    G4LogicalVolume* mylarLV = new G4LogicalVolume(mylarSolid, fMylar, "Mylar");
    new G4PVPlacement(0, G4ThreeVector(0, 0, relativeWindowStartZ + alWindowThick + mylarThick/2), mylarLV, "Mylar", vacuumLV, false, 0, true);

    // c. Aluminum Foil
    G4Tubs* alFoilSolid = new G4Tubs("AlFoil", 0, windowDiam/2, alFoilThick/2, 0*deg, 360*deg);
    G4LogicalVolume* alFoilLV = new G4LogicalVolume(alFoilSolid, fAluminum, "AlFoil");
    new G4PVPlacement(0, G4ThreeVector(0, 0, relativeWindowStartZ + alWindowThick + mylarThick + alFoilThick/2), alFoilLV, "AlFoil", vacuumLV, false, 0, true);

    // d. Aluminum Cup
    G4Tubs* alCupSolid = new G4Tubs("AlCup", 0, windowDiam/2, alCupThick/2, 0*deg, 360*deg);
    G4LogicalVolume* alCupLV = new G4LogicalVolume(alCupSolid, fAluminum, "AlCup");
    new G4PVPlacement(0, G4ThreeVector(0, 0, relativeWindowStartZ + totalWindowThick - alCupThick/2), alCupLV, "AlCup", vacuumLV, false, 0, true);

    // 3. Germanium Crystal with bore hole (inside vacuum volume)
    G4Tubs* geCrystalOuter = new G4Tubs("GeCrystalOuter", 0, geCrystalDiam/2, geCrystalLength/2, 0*deg, 360*deg);
    G4Tubs* boreHole = new G4Tubs("BoreHole", 0, boreHoleDiam/2, boreHoleDepth/2, 0*deg, 360*deg);
    
    // Position bore hole from front face
    G4double boreHoleZ = -geCrystalLength/2 + boreHoleDepth/2;
    G4SubtractionSolid* geCrystalSolid = new G4SubtractionSolid("GeCrystal", geCrystalOuter, boreHole, 
                                                               0, G4ThreeVector(0, 0, boreHoleZ));
    
    G4LogicalVolume* geCrystalLV = new G4LogicalVolume(geCrystalSolid, fGermanium, "GeCrystal");
    new G4PVPlacement(0, G4ThreeVector(0, 0, relativeGeZ), geCrystalLV, "GeCrystal", vacuumLV, false, 0, true);

    // 4. Lithium Dead Layer (outer contact) - inside vacuum volume
    G4double liInnerDiam = geCrystalDiam;
    G4double liOuterDiam = geCrystalDiam + 2*liDeadThick;
    G4Tubs* liDeadSolid = new G4Tubs("LiDead", liInnerDiam/2, liOuterDiam/2, geCrystalLength/2, 0*deg, 360*deg);
    G4LogicalVolume* liDeadLV = new G4LogicalVolume(liDeadSolid, fLithium, "LiDead");
    new G4PVPlacement(0, G4ThreeVector(0, 0, relativeGeZ), liDeadLV, "LiDead", vacuumLV, false, 0, true);

    // 5. Boron Dead Layer (inner contact) - inside vacuum volume
    G4double bInnerDiam = boreHoleDiam - 2*bDeadThick;
    G4double bOuterDiam = boreHoleDiam;
    G4Tubs* bDeadSolid = new G4Tubs("BDead", bInnerDiam/2, bOuterDiam/2, boreHoleDepth/2, 0*deg, 360*deg);
    G4LogicalVolume* bDeadLV = new G4LogicalVolume(bDeadSolid, fBoron, "BDead");
    new G4PVPlacement(0, G4ThreeVector(0, 0, relativeGeZ + boreHoleZ), bDeadLV, "BDead", vacuumLV, false, 0, true);

    // Set scoring volume to germanium crystal
    fScoringVolume = geCrystalLV;

    // Visualization attributes
    fWorldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Lead shield - dark gray, semi-transparent
    G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.6));
    leadVis->SetForceSolid(true);
    fLeadShieldLV->SetVisAttributes(leadVis);
    fLeadTopLV->SetVisAttributes(leadVis);
    fLeadBottomLV->SetVisAttributes(leadVis);

    // Copper liner - orange/brown, semi-transparent
    G4VisAttributes* copperVis = new G4VisAttributes(G4Colour(0.8, 0.5, 0.2, 0.5));
    copperVis->SetForceSolid(true);
    fCopperLinerLV->SetVisAttributes(copperVis);
    fCopperTopLV->SetVisAttributes(copperVis);
    fCopperBottomLV->SetVisAttributes(copperVis);

    // Shield cavity - invisible (to see detector inside)
    fShieldCavityLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Housing - dark gray
    G4VisAttributes* housingVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.8));
    housingVis->SetForceSolid(true);
    housingLV->SetVisAttributes(housingVis);

    // Windows - light gray (aluminum), purple (mylar)
    G4VisAttributes* alVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.6));
    alVis->SetForceSolid(true);
    alWindowLV->SetVisAttributes(alVis);
    alFoilLV->SetVisAttributes(alVis);
    alCupLV->SetVisAttributes(alVis);

    G4VisAttributes* mylarVis = new G4VisAttributes(G4Colour(0.8, 0.2, 0.8, 0.6));
    mylarVis->SetForceSolid(true);
    mylarLV->SetVisAttributes(mylarVis);

    // Germanium crystal - cyan
    G4VisAttributes* geVis = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.8));
    geVis->SetForceSolid(true);
    geCrystalLV->SetVisAttributes(geVis);

    // Dead layers - red (Li), blue (B)
    G4VisAttributes* liVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
    liVis->SetForceSolid(true);
    liDeadLV->SetVisAttributes(liVis);

    G4VisAttributes* bVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.7));
    bVis->SetForceSolid(true);
    bDeadLV->SetVisAttributes(bVis);

    G4cout << "Realistic HPGe detector constructed successfully" << G4endl;
    G4cout << "Active Ge volume: " << geCrystalDiam/mm << " mm diameter, " 
           << geCrystalLength/mm << " mm length" << G4endl;
    G4cout << "Bore hole: " << boreHoleDiam/mm << " mm diameter, " 
           << boreHoleDepth/mm << " mm deep" << G4endl;

    return fWorldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    // Create sensitive detector for germanium crystal only
    G4MultiFunctionalDetector* detector = new G4MultiFunctionalDetector("HPGe");
    G4SDManager::GetSDMpointer()->AddNewDetector(detector);

    G4VPrimitiveScorer* primitive = new G4PSEnergyDeposit("EnergyDeposit");
    detector->RegisterPrimitive(primitive);

    SetSensitiveDetector("GeCrystal", detector);

    G4cout << "Sensitive detector attached to Ge crystal only" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
