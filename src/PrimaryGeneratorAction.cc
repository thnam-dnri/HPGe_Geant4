
// ==============================================================================
// PrimaryGeneratorAction.cc - Implementation  
// ==============================================================================

#include "PrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4Event.hh"
#include <fstream>
#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(const std::string& rainierFile)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fRAINIERFile(rainierFile),
  fRandomGenerator(std::random_device{}())
{
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    // Default particle type and properties
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName = "gamma";
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
    
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fParticleGun->SetParticleEnergy(1.*MeV);

    // Load gamma data
    if (!fRAINIERFile.empty()) {
        LoadRAINIERData();
    } else {
        LoadTestData();
        G4cout << "Using Co-60 test data (1173 keV and 1332 keV gammas)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Sample a gamma ray from the data
    GammaData gamma = SampleGamma();
    
    // Set gamma energy
    fParticleGun->SetParticleEnergy(gamma.energy * MeV);
    
    // Set source position (point source for now)
    G4ThreeVector sourcePos = SampleSourcePosition();
    fParticleGun->SetParticlePosition(sourcePos);
    
    // Set isotropic direction
    G4ThreeVector direction = SampleDirection();
    fParticleGun->SetParticleMomentumDirection(direction);
    
    // DEBUG: Print every 1000th event
    if (anEvent->GetEventID() % 1000 == 0) {
        G4cout << "Event " << anEvent->GetEventID() 
               << ": Generated " << gamma.energy << " MeV gamma at " 
               << sourcePos << " with direction " << direction << G4endl;
    }
    
    // Generate the primary particle
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::LoadRAINIERData()
{
    // Simple ASCII format reader for now
    // Expected format: energy(MeV) intensity
    // TODO: Implement ROOT file reader for actual RAINIER output
    
    std::ifstream file(fRAINIERFile);
    if (!file.is_open()) {
        G4cout << "Warning: Could not open RAINIER file " << fRAINIERFile << G4endl;
        G4cout << "Using test data instead." << G4endl;
        LoadTestData();
        return;
    }
    
    fGammaData.clear();
    double energy, intensity;
    
    G4cout << "Loading RAINIER gamma data from " << fRAINIERFile << G4endl;
    
    while (file >> energy >> intensity) {
        if (energy > 0.001 && intensity > 0.001) { // Filter low energy/intensity
            GammaData gamma;
            gamma.energy = energy;
            gamma.intensity = intensity;
            fGammaData.push_back(gamma);
        }
    }
    
    file.close();
    
    if (fGammaData.empty()) {
        G4cout << "No valid gamma data found in file. Using test data." << G4endl;
        LoadTestData();
    } else {
        G4cout << "Loaded " << fGammaData.size() << " gamma transitions" << G4endl;
        
        // Print first few gammas for verification
        G4cout << "First few gamma energies:" << G4endl;
        for (size_t i = 0; i < std::min(size_t(5), fGammaData.size()); i++) {
            G4cout << "  " << fGammaData[i].energy << " MeV (intensity: " 
                   << fGammaData[i].intensity << ")" << G4endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::LoadTestData()
{
    // Co-60 test data: 1173.2 keV and 1332.5 keV gammas
    fGammaData.clear();
    
    GammaData gamma1;
    gamma1.energy = 1.1732;    // MeV
    gamma1.intensity = 99.85;  // %
    
    GammaData gamma2;
    gamma2.energy = 1.3325;    // MeV  
    gamma2.intensity = 99.98;  // %
    
    fGammaData.push_back(gamma1);
    fGammaData.push_back(gamma2);
    
    G4cout << "Loaded Co-60 test data:" << G4endl;
    G4cout << "  1173.2 keV gamma (99.85% intensity)" << G4endl;
    G4cout << "  1332.5 keV gamma (99.98% intensity)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GammaData PrimaryGeneratorAction::SampleGamma()
{
    if (fGammaData.empty()) {
        LoadTestData();
    }
    
    // Simple random sampling based on intensity
    // Create cumulative distribution
    std::vector<double> cumulative;
    double sum = 0.0;
    
    for (const auto& gamma : fGammaData) {
        sum += gamma.intensity;
        cumulative.push_back(sum);
    }
    
    // Sample using cumulative distribution
    std::uniform_real_distribution<double> dist(0.0, sum);
    double r = dist(fRandomGenerator);
    
    for (size_t i = 0; i < cumulative.size(); i++) {
        if (r <= cumulative[i]) {
            return fGammaData[i];
        }
    }
    
    // Fallback to last gamma
    return fGammaData.back();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleSourcePosition()
{
    // Point source at origin for now
    // TODO: Add extended source options
    return G4ThreeVector(0., 0., 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleDirection()
{
    // Uniformly sample within a 120° cone aimed along +Z toward the detector
    static const G4double coneHalfAngle = 60.0 * deg;
    const double cosThetaMax = std::cos(coneHalfAngle);

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const double cosTheta = cosThetaMax + (1.0 - cosThetaMax) * dist(fRandomGenerator);
    const double phi = 2.0 * pi * dist(fRandomGenerator);
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    
    return G4ThreeVector(sinTheta * std::cos(phi),
                         sinTheta * std::sin(phi),
                         cosTheta);
}
