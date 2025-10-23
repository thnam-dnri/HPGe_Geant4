
// ==============================================================================
// EventAction.hh/cc - Phase 1 Event-level analysis
// ==============================================================================

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEnergyDeposit(G4double edep) { fEnergyDeposit += edep; }

private:
    RunAction* fRunAction;
    G4double fEnergyDeposit;
};
#endif
