// ==============================================================================
// RunAction.hh/cc - Phase 1 Run-level analysis
// ==============================================================================

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include "G4AnalysisManager.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();

    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddEnergyDeposit(G4double edep);
    void SetDetID(const G4String& det) { fDetID = det; }

private:
    G4Accumulable<G4double> fEnergyDeposit;
    G4AnalysisManager* fAnalysisManager;

    // Output controls
    G4String fDetID {"HPGe1"};
    std::vector<G4double> fEnergyAxisCenters_keV; // for Axes ntuple
    G4bool fIncludeTrueEnergy {true};
    G4bool fOutputEnabled {true};
};
#endif
