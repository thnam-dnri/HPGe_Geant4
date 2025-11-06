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
#include <ctime>
#include <sstream>
#include <filesystem>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction(),
  fEnergyDeposit("EnergyDeposit", 0.),
  fAnalysisManager(G4AnalysisManager::Instance())
{
    fAnalysisManager->SetVerboseLevel(1);
    fAnalysisManager->SetDefaultFileType("root");

    // Create ntuple 0: Events (singles list-mode)
    fAnalysisManager->CreateNtuple("Events", "Singles events");
    fAnalysisManager->CreateNtupleDColumn("Edep_keV"); // col 0
    if (fIncludeTrueEnergy) {
        fAnalysisManager->CreateNtupleDColumn("Etrue_keV"); // col 1
    }
    fAnalysisManager->FinishNtuple(); // id 0

    // Create ntuple 1: RunInfo (one-row metadata)
    fAnalysisManager->CreateNtuple("RunInfo", "Run metadata");
    fAnalysisManager->CreateNtupleSColumn("Nuclide");
    fAnalysisManager->CreateNtupleSColumn("Generator");
    fAnalysisManager->CreateNtupleSColumn("G4Version");
    fAnalysisManager->CreateNtupleSColumn("PhysicsLists");
    fAnalysisManager->CreateNtupleSColumn("DetID");
    fAnalysisManager->CreateNtupleSColumn("Geometry");
    fAnalysisManager->CreateNtupleSColumn("Binning");
    fAnalysisManager->CreateNtupleSColumn("SourceJSON");
    fAnalysisManager->CreateNtupleSColumn("SourceHash");
    fAnalysisManager->CreateNtupleIColumn("Ndecays");
    fAnalysisManager->CreateNtupleIColumn("NgammaPrimaries");
    fAnalysisManager->FinishNtuple(); // id 1

    // Create ntuple 2: Axes (one-row centers)
    fAnalysisManager->CreateNtuple("Axes_energy_centers", "Energy axis centers");
    // We'll attach vector at BeginOfRunAction
    // Use a placeholder scalar to keep structure; we will replace with vector column at fill time if supported
    // fAnalysisManager->CreateNtupleDColumn("energy_keV", fEnergyAxisCenters_keV);
    fAnalysisManager->CreateNtupleDColumn("energy_keV");
    fAnalysisManager->FinishNtuple(); // id 2

    // Create ntuple 3: Truth (gamma lines)
    fAnalysisManager->CreateNtuple("Truth_gamma_lines", "Per-line truth");
    fAnalysisManager->CreateNtupleDColumn("E_gamma_keV");
    fAnalysisManager->CreateNtupleDColumn("I_per_decay");
    // Emitter-aware fields
    fAnalysisManager->CreateNtupleSColumn("Emitter");
    fAnalysisManager->CreateNtupleDColumn("EmitterHalfLife_s");
    fAnalysisManager->FinishNtuple(); // id 3

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

    // Check if current isotope configuration actually produces gamma lines
    const auto* gen = static_cast<const PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    fOutputEnabled = true;
    if (gen) {
        const auto& truthCheck = const_cast<PrimaryGeneratorAction*>(gen)->GetTruthGammaLines();
        if (truthCheck.empty()) {
            G4cout << "[RunAction] Selected isotope has no gamma lines; aborting run and suppressing output file." << G4endl;
            fOutputEnabled = false;
            G4RunManager::GetRunManager()->AbortRun(true);
            return; // do not open files or fill any ntuples
        }
    }

    // Compose output file name per spec if possible
    G4String outfile = "gamma_spectrum"; // fallback
    if (gen && !gen->GetIsotopeSymbol().empty()) {
        // Build: <Nuclide>_<DetID>_<GapMM>_<runid>
        const G4String nuclide = gen->GetIsotopeSymbol().c_str();
        const G4double gapmm = gen->GetSourceSurfaceGap() / CLHEP::mm;
        const int gapInt = static_cast<int>(std::round(gapmm));
        // Simple run id from time
        std::time_t t = std::time(nullptr);
        char buf[32];
        std::strftime(buf, sizeof(buf), "%Y%m%d%H%M%S", std::localtime(&t));
        outfile = nuclide + "_" + fDetID + "_" + std::to_string(gapInt) + "mm_" + buf;
    }
    const std::filesystem::path outputDir{"training_data"};
    std::error_code dirErr;
    std::filesystem::create_directories(outputDir, dirErr);
    if (dirErr) {
        G4cout << "[RunAction] Warning: unable to ensure training_data directory exists ("
               << dirErr.message() << "). Files will be created relative to current path." << G4endl;
    }
    const std::filesystem::path outputPath = outputDir / outfile.c_str();
    fAnalysisManager->SetFileName(outputPath.string());
    fAnalysisManager->OpenFile();

    // reset accumulables to their initial values
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Reset();

    G4cout << "\n-------- Starting Run --------" << G4endl;

    // Prepare energy axis centers (consistent across runs)
    // Default: 0..3000 keV, 8192 bins
    const int n_bins = 8192;
    const double e_min_keV = 0.0;
    const double e_max_keV = 3000.0;
    fEnergyAxisCenters_keV.clear();
    fEnergyAxisCenters_keV.reserve(n_bins);
    const double bin_width = (e_max_keV - e_min_keV) / n_bins;
    for (int i = 0; i < n_bins; ++i) {
        double center = e_min_keV + (i + 0.5) * bin_width;
        fEnergyAxisCenters_keV.push_back(center);
    }
    // Fill ntuple 2: Axes_energy_centers (as a single row with repeated fills)
    for (double e : fEnergyAxisCenters_keV) {
        fAnalysisManager->FillNtupleDColumn(2, 0, e);
        fAnalysisManager->AddNtupleRow(2);
    }

    // Fill ntuple 3: Truth lines from generator (with emitter symbol and half-life if available)
    if (gen) {
        const auto& dtruth = const_cast<PrimaryGeneratorAction*>(gen)->GetTruthGammaLinesDetailed();
        if (!dtruth.empty()) {
            for (const auto& t : dtruth) {
                fAnalysisManager->FillNtupleDColumn(3, 0, t.energy_keV);
                fAnalysisManager->FillNtupleDColumn(3, 1, t.intensity_per_decay);
                fAnalysisManager->FillNtupleSColumn(3, 2, t.emitter_symbol.c_str());
                fAnalysisManager->FillNtupleDColumn(3, 3, t.emitter_half_life_s);
                fAnalysisManager->AddNtupleRow(3);
            }
        } else {
            const auto& truth = const_cast<PrimaryGeneratorAction*>(gen)->GetTruthGammaLines();
            for (const auto& p : truth) {
                fAnalysisManager->FillNtupleDColumn(3, 0, p.first);
                fAnalysisManager->FillNtupleDColumn(3, 1, p.second);
                fAnalysisManager->FillNtupleSColumn(3, 2, "");
                fAnalysisManager->FillNtupleDColumn(3, 3, 0.0);
                fAnalysisManager->AddNtupleRow(3);
            }
        }
    }
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

    // Fill RunInfo (ntuple id 1), one row
    const auto* gen2 = static_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    G4String nuclide = (gen2 && !gen2->GetIsotopeSymbol().empty()) ? G4String(gen2->GetIsotopeSymbol().c_str()) : G4String("N/A");
    G4String generator = gen2 ? G4String("IsotopeJSON gamma singles; chain-to-stable") : G4String("Unconfigured");
    G4String g4ver = G4String("unknown");
    G4String phys = G4String("G4EmStandardPhysics_option4 + G4DecayPhysics + G4RadioactiveDecayPhysics");
    G4double gapmm = gen2 ? gen2->GetSourceSurfaceGap() / CLHEP::mm : 0.0;
    std::ostringstream geom;
    geom << "{\"Gap_mm\":" << gapmm << "}";
    std::ostringstream binjson;
    binjson << "{\"e_min_keV\":0,\"e_max_keV\":3000,\"n_bins\":8192,\"grid\":\"centers\"}";
    G4String srcjson = nuclide != "N/A" ? G4String(("isotope_data/" + nuclide + ".json").c_str()) : G4String("N/A");
    G4String srchash = "N/A"; // hash omitted in this implementation
    G4int ndec = gen2 ? static_cast<G4int>(gen2->GetNDecays()) : 0;
    G4int ngam = gen2 ? static_cast<G4int>(gen2->GetNGammaPrimaries()) : 0;

    if (fOutputEnabled) {
        fAnalysisManager->FillNtupleSColumn(1, 0, nuclide);
        fAnalysisManager->FillNtupleSColumn(1, 1, generator);
        fAnalysisManager->FillNtupleSColumn(1, 2, g4ver);
        fAnalysisManager->FillNtupleSColumn(1, 3, phys);
        fAnalysisManager->FillNtupleSColumn(1, 4, fDetID);
        fAnalysisManager->FillNtupleSColumn(1, 5, geom.str());
        fAnalysisManager->FillNtupleSColumn(1, 6, binjson.str());
        fAnalysisManager->FillNtupleSColumn(1, 7, srcjson);
        fAnalysisManager->FillNtupleSColumn(1, 8, srchash);
        fAnalysisManager->FillNtupleIColumn(1, 9, ndec);
        fAnalysisManager->FillNtupleIColumn(1,10, ngam);
        fAnalysisManager->AddNtupleRow(1);

        fAnalysisManager->Write();
        fAnalysisManager->CloseFile();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEnergyDeposit(G4double edep)
{
    fEnergyDeposit  += edep;
}
