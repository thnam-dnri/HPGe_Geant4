
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Run.cc - Implementation

#include "Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
: G4Run(),
  fTotalEnergyDeposit(0.),
  fTotalEvents(0)
{
    // Initialize the histogram with zeros for all bins
//    for (G4int i = 0; i < fNbins; i++) {
//        fEnergyHistogram[i] = 0;
//    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
    const Run* localRun = static_cast<const Run*>(run);
    
    // Merge histograms
    for (const auto& bin : localRun->fEnergyHistogram) {
        fEnergyHistogram[bin.first] += bin.second;
    }
    
    fTotalEnergyDeposit += localRun->fTotalEnergyDeposit;
    fTotalEvents += localRun->fTotalEvents;
    
    G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEnergySpectrum(G4double energy)
{
    G4int bin = EnergyToBin(energy);
    if (bin >= 0 && bin < fNbins) {
        fEnergyHistogram[bin]++;
    }
    fTotalEnergyDeposit += energy;
    fTotalEvents++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int Run::EnergyToBin(G4double energy) const
{
    // Convert energy from MeV to keV and bin it (1 keV per bin)
    G4double energy_keV = energy / keV;
    if (energy_keV < 0 || energy_keV >= fEmax * 1000) {
        return -1; // Out of range
    }
    return static_cast<G4int>(energy_keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::PrintResults() const
{
    G4cout << "\n========== Energy Spectrum Results ==========" << G4endl;
    G4cout << "Total events processed: " << fTotalEvents << G4endl;
    G4cout << "Total energy deposited: " << G4BestUnit(fTotalEnergyDeposit, "Energy") << G4endl;
    
    // Print significant peaks (>10 counts)
    G4cout << "\nSignificant peaks (>10 counts):" << G4endl;
    G4cout << "Energy [keV]\tCounts" << G4endl;
    for (const auto& bin : fEnergyHistogram) {
        if (bin.second > 10) {
            G4cout << bin.first << "\t\t" << bin.second << G4endl;
        }
    }
    
    G4cout << "\nPer-event energy deposits recorded in gamma_spectrum.root"
           << " (ntuple EventData with columns eventID and edep_keV)" << G4endl;
    
    G4cout << "=============================================\n" << G4endl;
}
