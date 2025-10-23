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

private:
    G4Accumulable<G4double> fEnergyDeposit;
    G4AnalysisManager* fAnalysisManager;
};
#endif
