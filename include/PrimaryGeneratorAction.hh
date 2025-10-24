// ==============================================================================
// PrimaryGeneratorAction.hh/cc - Phase 1 RAINIER Data Reader & Generator
// ==============================================================================

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <string>
#include <vector>
#include <random>

class DetectorConstruction; // forward declaration

class G4ParticleGun;
class G4Event;

// Simple structure for gamma data from RAINIER
struct GammaData {
    double energy;    // MeV
    double intensity; // relative intensity
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(const std::string& rainierFile = "");
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    // Configuration hooks
    void SetDetectorConstruction(const DetectorConstruction* det) { fDetConstruction = det; }
    void SetSourceSurfaceGap(G4double gap) { fSourceSurfaceGap = gap; }

private:
    G4ParticleGun* fParticleGun;
    std::string fRAINIERFile;
    std::vector<GammaData> fGammaData;
    std::mt19937 fRandomGenerator;

    // Methods for RAINIER data handling
    void LoadRAINIERData();
    void LoadTestData();  // Co-60 test data if no RAINIER file
    GammaData SampleGamma();
    
    // Position and direction sampling
    G4ThreeVector SampleSourcePosition();
    G4ThreeVector SampleDirection();

    // Geometry coupling
    const DetectorConstruction* fDetConstruction {nullptr};
    G4double fSourceSurfaceGap {0.0}; // surface-to-surface distance to detector front window
};
#endif
