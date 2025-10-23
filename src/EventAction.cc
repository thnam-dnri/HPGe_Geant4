
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// EventAction.cc - Implementation

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEnergyDeposit(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    fEnergyDeposit = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Add energy deposit to run accumulation
    fRunAction->AddEnergyDeposit(fEnergyDeposit);
    
    // Also update the Run class for energy spectrum and record in ROOT
    if (fEnergyDeposit > 0.) {
        Run* currentRun = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
        if (currentRun) {
            currentRun->AddEnergySpectrum(fEnergyDeposit);
        }

        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillNtupleIColumn(0, event ? event->GetEventID() : -1);
        analysisManager->FillNtupleDColumn(1, fEnergyDeposit / keV);
        analysisManager->AddNtupleRow();
    }
}
