// ==============================================================================
// PrimaryGeneratorAction.hh/cc - Isotope-driven gamma generator
// ==============================================================================

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <string>
#include <vector>
#include <random>
#include <deque>
#include <utility>

class DetectorConstruction; // forward declaration

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    // Configuration hooks
    void SetDetectorConstruction(const DetectorConstruction* det) { fDetConstruction = det; }
    void SetSourceSurfaceGap(G4double gap) { fSourceSurfaceGap = gap; }
    void SetIsotopeSymbol(const std::string& s);
    const std::string& GetIsotopeSymbol() const { return fIsotopeSymbol; }
    G4double GetSourceSurfaceGap() const { return fSourceSurfaceGap; }

    // Emission mode: parent-only (no daughter tracking) or full-chain
    enum class DecayEmissionMode { ParentOnly = 0, FullChain = 1 };
    void SetDecayMode(DecayEmissionMode m) { fDecayMode = m; fTruthReady = false; }
    DecayEmissionMode GetDecayMode() const { return fDecayMode; }

    // Last primary true energy (keV) for QA/optional output
    double GetLastTrueEnergyKeV() const { return fLastTrueEnergy_keV; }

    // Counters for normalization
    unsigned long long GetNDecays() const { return fNdecays; }
    unsigned long long GetNGammaPrimaries() const { return fNgammaPrimaries; }

    // Truth gamma lines aggregated to stability: (E_gamma_keV, I_per_decay)
    const std::vector<std::pair<double,double>>& GetTruthGammaLines();
    struct TruthLineDetailed {
        double energy_keV {0.0};
        double intensity_per_decay {0.0};
        std::string emitter_symbol; // emitter nuclide symbol
        double emitter_half_life_s {0.0};
    };
    // Emitter-aware truth lines
    const std::vector<TruthLineDetailed>& GetTruthGammaLinesDetailed();

private:
    G4ParticleGun* fParticleGun;
    std::string fIsotopeSymbol;
    std::mt19937 fRandomGenerator;

    // Position and direction sampling
    G4ThreeVector SampleSourcePosition();
    G4ThreeVector SampleDirection();

    // Geometry coupling
    const DetectorConstruction* fDetConstruction {nullptr};
    G4double fSourceSurfaceGap {0.0}; // surface-to-surface distance to detector front window

    // Isotope-based generation helpers
    void GenerateFromIsotope_Singles(G4Event*);
    bool PrepareNextDecayGammas();
    void AccumulateTruthLines(const std::string& iso, double weight);
    void AccumulateTruthLinesDetailed(const std::string& iso, double weight);
    void AccumulateTruthLinesParentOnly(const std::string& iso);

    // Queue of gamma energies (keV) to emit, one per Geant4 event
    std::deque<double> fGammaQueue_keV;
    double fLastTrueEnergy_keV {0.0};
    unsigned long long fNdecays {0};
    unsigned long long fNgammaPrimaries {0};
    bool fTruthReady {false};
    std::vector<std::pair<double,double>> fTruthLines; // (E_gamma_keV, I_per_decay)
    std::vector<TruthLineDetailed> fTruthLinesDetailed;

    DecayEmissionMode fDecayMode {DecayEmissionMode::ParentOnly};
};
#endif
