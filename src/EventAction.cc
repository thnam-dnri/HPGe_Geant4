
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// EventAction.cc - Implementation

#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

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
        // Ntuple 0: Events
        analysisManager->FillNtupleDColumn(0, 0, fEnergyDeposit / keV);
        // Optional true energy from generator (keV)
        const auto* gen = static_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
        if (gen) {
            analysisManager->FillNtupleDColumn(0, 1, gen->GetLastTrueEnergyKeV());
        }
        analysisManager->AddNtupleRow(0);
    }
}
