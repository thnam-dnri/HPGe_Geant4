//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// RunAction.cc - Implementation

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction(),
  fEnergyDeposit("EnergyDeposit", 0.),
  fAnalysisManager(G4AnalysisManager::Instance())
{
    fAnalysisManager->SetVerboseLevel(1);
    fAnalysisManager->SetFileName("gamma_spectrum");
    fAnalysisManager->SetDefaultFileType("root");
    fAnalysisManager->CreateNtuple("EventData", "Per-event energy deposit");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->CreateNtupleDColumn("edep_keV");
    fAnalysisManager->FinishNtuple();

    // Register accumulable to the accumulable manager
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Register(fEnergyDeposit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
    return new Run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    fAnalysisManager->OpenFile();

    // reset accumulables to their initial values
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Reset();

    G4cout << "\n-------- Starting Run --------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;

    // Merge accumulables
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Merge();

    // Compute dose = total energy deposit in a run and its variance
    G4double energyDeposit = fEnergyDeposit.GetValue();
    G4double energyDeposit2 = fEnergyDeposit.GetValue()*fEnergyDeposit.GetValue();

    G4double rms = energyDeposit2 - energyDeposit*energyDeposit/nofEvents;
    if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

    const DetectorConstruction* detectorConstruction
     = static_cast<const DetectorConstruction*>
       (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
    G4double dose = energyDeposit/mass;
    G4double rmsDose = rms/mass;

    // Run conditions
    const PrimaryGeneratorAction* generatorAction
     = static_cast<const PrimaryGeneratorAction*>
       (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    G4String runCondition;
    if (generatorAction)
    {
      const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
      runCondition += particleGun->GetParticleDefinition()->GetParticleName();
      runCondition += " of ";
      G4double particleEnergy = particleGun->GetParticleEnergy();
      runCondition += G4BestUnit(particleEnergy,"Energy");
    }

    // Print
    if (IsMaster()) {
        G4cout
         << G4endl
         << "-------- End of Global Run --------"
         << G4endl
         << " The run consists of " << nofEvents << " events"
         << G4endl
         << " Cumulative energy deposit: "
         << G4BestUnit(energyDeposit,"Energy") << " rms = "
         << G4BestUnit(rmsDose,"Energy")
         << G4endl
         << " Dose in scoring volume : "
         << G4BestUnit(dose,"Dose") << " rms = "
         << G4BestUnit(rmsDose,"Dose")
         << G4endl
         << "------------------------------------"
         << G4endl
         << G4endl;
    }
    // Print final results and write spectrum file
    Run* localRun = (Run*)run;
    if (IsMaster()) {
        localRun->PrintResults();
    }

    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEnergyDeposit(G4double edep)
{
    fEnergyDeposit  += edep;
}
