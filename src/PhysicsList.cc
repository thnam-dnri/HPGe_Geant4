// ==============================================================================
// PhysicsList.cc - Implementation
// ==============================================================================

#include "PhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
: G4VModularPhysicsList()
{
    SetVerboseLevel(2);

    // Default physics for gamma spectroscopy
    
    // Electromagnetic physics (high precision for gamma spectroscopy)
    RegisterPhysics(new G4EmStandardPhysics_option4());
    
    // Extra EM processes (gamma conversion, etc.)
    RegisterPhysics(new G4EmExtraPhysics());
    
    // Decay physics
    RegisterPhysics(new G4DecayPhysics());
    
    // Radioactive decay physics
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    
    // Ion physics (for completeness)
    RegisterPhysics(new G4IonPhysics());
    
    G4cout << "Physics list initialized for gamma spectroscopy" << G4endl;
    G4cout << "Using G4EmStandardPhysics_option4 for high precision EM" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
    // Set production cuts optimized for gamma spectroscopy
    
    // Default cut values
    SetCutsWithDefault();
    
    // Lower cuts for better accuracy in detector
    // These are production thresholds, not tracking cuts
    SetCutValue(0.1*mm, "gamma");      // Photons
    SetCutValue(0.1*mm, "e-");         // Electrons  
    SetCutValue(0.1*mm, "e+");         // Positrons
    SetCutValue(0.1*mm, "proton");     // Protons
    
    // Even lower cuts in germanium detector for precision
    G4Region* detectorRegion = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld");
    if (detectorRegion) {
        G4ProductionCuts* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(0.01*mm, "gamma");
        cuts->SetProductionCut(0.01*mm, "e-");  
        cuts->SetProductionCut(0.01*mm, "e+");
        cuts->SetProductionCut(0.01*mm, "proton");
        detectorRegion->SetProductionCuts(cuts);
    }
    
    if (verboseLevel > 0) {
        DumpCutValuesTable();
    }
}