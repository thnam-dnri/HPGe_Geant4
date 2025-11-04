
// ==============================================================================
// PrimaryGeneratorAction.cc - Implementation  
// ==============================================================================

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "IsotopeDecay.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4Event.hh"
#include <algorithm>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
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

    SetIsotopeSymbol("Co60");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetIsotopeSymbol(const std::string& s)
{
    fIsotopeSymbol = s;
    fGammaQueue_keV.clear();
    fLastTrueEnergy_keV = 0.0;
    fNdecays = 0;
    fNgammaPrimaries = 0;
    fTruthReady = false;
    fTruthLines.clear();
    fSinglesReady = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Fast path: precomputed singles sampling (one gamma per event)
    if (!fSinglesReady) {
        RebuildSinglesDistribution();
    }

    if (fSinglesReady && !fSinglesEnergies_keV.empty() && !fSinglesWeights.empty()) {
        // Sample one gamma according to intensity weights
        const size_t idx = fSinglesDist(fRandomGenerator);
        fLastTrueEnergy_keV = fSinglesEnergies_keV[idx];
        ++fNgammaPrimaries;
        ++fNdecays; // treat one sample as one decay for counters

        fParticleGun->SetParticleEnergy(fLastTrueEnergy_keV * keV);
        fParticleGun->SetParticlePosition(SampleSourcePosition());
        fParticleGun->SetParticleMomentumDirection(SampleDirection());
        fParticleGun->GeneratePrimaryVertex(anEvent);
        return;
    }

    // Fallback to legacy queue-based decay sampler if singles not available
    GenerateFromIsotope_Singles(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleSourcePosition()
{
    // Disk source: 2 cm diameter (R=10 mm), 1 mm thickness
    static const G4double R = 10.0 * mm;         // radius
    static const G4double T = 1.0 * mm;          // thickness

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Sample uniformly over disk area: r = R * sqrt(u), phi = 2*pi*v
    const double u = dist(fRandomGenerator);
    const double v = dist(fRandomGenerator);
    const double r = R * std::sqrt(u);
    const double phi = 2.0 * pi * v;
    const double x = r * std::cos(phi);
    const double y = r * std::sin(phi);

    // Sample uniformly along thickness around z=0
    const double w = dist(fRandomGenerator);
    const double zLocal = (w - 0.5) * T; // range [-T/2, +T/2]

    // Place disk in front of detector window, at surface-to-surface gap
    // Detector surface Z comes from DetectorConstruction (window front face)
    G4double surfaceZ = 0.0;
    if (fDetConstruction) {
        surfaceZ = fDetConstruction->GetDetectorSurfaceZ();
    }
    // Center of disk sits at: surfaceZ - (gap + T/2)
    const G4double centerZ = surfaceZ - (fSourceSurfaceGap + 0.5*T);
    const G4double z = centerZ + zLocal;

    return G4ThreeVector(x, y, z);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateFromIsotope_Singles(G4Event* anEvent)
{
    // Ensure there is at least one gamma queued; if not, sample decays until there is
    int tries = 0;
    while (fGammaQueue_keV.empty() && tries++ < 256) {
        PrepareNextDecayGammas();
    }

    // If still empty, nothing to shoot
    if (fGammaQueue_keV.empty()) {
        if (anEvent->GetEventID() % 1000 == 0) {
            G4cout << "[Isotope] No gamma available to emit in event " << anEvent->GetEventID() << G4endl;
        }
        // Emit a dummy very low energy gamma to keep event loop consistent
        fLastTrueEnergy_keV = 0.0;
        fParticleGun->SetParticleEnergy(1e-6 * keV);
        fParticleGun->SetParticlePosition(SampleSourcePosition());
        fParticleGun->SetParticleMomentumDirection(SampleDirection());
        fParticleGun->GeneratePrimaryVertex(anEvent);
        return;
    }

    // Pop one gamma and emit
    fLastTrueEnergy_keV = fGammaQueue_keV.front();
    fGammaQueue_keV.pop_front();
    ++fNgammaPrimaries;

    fParticleGun->SetParticleEnergy(fLastTrueEnergy_keV * keV);
    fParticleGun->SetParticlePosition(SampleSourcePosition());
    fParticleGun->SetParticleMomentumDirection(SampleDirection());
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

bool PrimaryGeneratorAction::PrepareNextDecayGammas()
{
    static IsotopeDataLoader loader;
    std::string current = fIsotopeSymbol;
    if (current.empty()) return false;
    ++fNdecays;

    bool anyGamma = false;
    int steps = 0;
    const int maxSteps = 64;
    while (!current.empty() && steps++ < maxSteps) {
        IsotopeInfo iso;
        if (!loader.Load(current, iso)) {
            G4cout << "[Isotope] Warning: could not load data for '" << current
                   << "' from isotope_data. Stopping chain for this decay." << G4endl;
            break;
        }
        if (iso.is_stable || iso.modes.empty()) {
            break;
        }
        // Sample a decay mode by branching ratio
        double totalBR = 0.0;
        for (const auto& m : iso.modes) totalBR += m.branching_ratio;
        if (totalBR <= 0.0) break;
        std::uniform_real_distribution<double> dist(0.0, totalBR);
        const double r = dist(fRandomGenerator);
        const DecayMode* chosen = nullptr;
        double acc = 0.0;
        for (const auto& m : iso.modes) {
            acc += m.branching_ratio;
            if (r <= acc) { chosen = &m; break; }
        }
        if (!chosen) chosen = &iso.modes.back();

        // Sample gamma emissions for this branch
        for (const auto& gl : chosen->gammas) {
            double mean = gl.absolute_intensity;
            if (mean <= 0.0) continue;
            int copies = 0;
            if (mean < 1.0) {
                std::bernoulli_distribution b(mean);
                copies = b(fRandomGenerator) ? 1 : 0;
            } else {
                copies = static_cast<int>(mean);
                double frac = mean - copies;
                if (frac > 1e-12) {
                    std::bernoulli_distribution b(frac);
                    if (b(fRandomGenerator)) ++copies;
                }
            }
            for (int i = 0; i < copies; ++i) {
                fGammaQueue_keV.push_back(gl.energy_keV);
                anyGamma = true;
            }
        }
        // In ParentOnly mode, do not traverse into daughters
        if (fDecayMode == DecayEmissionMode::ParentOnly) {
            break;
        }
        current = chosen->daughter;
    }
    return anyGamma;
}

const std::vector<std::pair<double,double>>& PrimaryGeneratorAction::GetTruthGammaLines()
{
    if (fTruthReady) return fTruthLines;
    fTruthLines.clear();
    if (!fIsotopeSymbol.empty()) {
        if (fDecayMode == DecayEmissionMode::ParentOnly) {
            // Parent-only truth aggregation: immediate gammas only
            // Weighted by parent branching ratios
            static IsotopeDataLoader loader;
            IsotopeInfo iso;
            if (loader.Load(fIsotopeSymbol, iso)) {
                for (const auto& m : iso.modes) {
                    const double w = m.branching_ratio;
                    if (w <= 0.0) continue;
                    for (const auto& gl : m.gammas) {
                        if (gl.absolute_intensity > 0.0) {
                            fTruthLines.emplace_back(gl.energy_keV, w * gl.absolute_intensity);
                        }
                    }
                }
            }
        } else {
            AccumulateTruthLines(fIsotopeSymbol, 1.0);
        }
    }
    fTruthReady = true;
    return fTruthLines;
}

void PrimaryGeneratorAction::AccumulateTruthLines(const std::string& isoSymbol, double weight)
{
    static IsotopeDataLoader loader;
    IsotopeInfo iso;
    if (!loader.Load(isoSymbol, iso)) return;
    if (iso.is_stable || iso.modes.empty()) return;
    for (const auto& m : iso.modes) {
        const double w = weight * m.branching_ratio;
        if (w <= 0.0) continue;
        for (const auto& gl : m.gammas) {
            if (gl.absolute_intensity > 0.0) {
                fTruthLines.emplace_back(gl.energy_keV, w * gl.absolute_intensity);
            }
        }
        if (!m.daughter.empty()) {
            AccumulateTruthLines(m.daughter, w);
        }
    }
}

const std::vector<PrimaryGeneratorAction::TruthLineDetailed>& PrimaryGeneratorAction::GetTruthGammaLinesDetailed()
{
    // Recompute detailed truth lines each call to ensure consistency with mode
    fTruthLinesDetailed.clear();
    if (fIsotopeSymbol.empty()) return fTruthLinesDetailed;
    AccumulateTruthLinesDetailed(fIsotopeSymbol, 1.0);
    return fTruthLinesDetailed;
}

void PrimaryGeneratorAction::RebuildSinglesDistribution()
{
    fSinglesEnergies_keV.clear();
    fSinglesWeights.clear();

    // Build from truth lines respecting current decay mode
    const auto& lines = GetTruthGammaLines();
    double totalW = 0.0;
    for (const auto& p : lines) {
        const double e = p.first;
        const double w = p.second;
        if (e > 0.0 && w > 0.0) {
            fSinglesEnergies_keV.push_back(e);
            fSinglesWeights.push_back(w);
            totalW += w;
        }
    }

    if (!fSinglesWeights.empty() && totalW > 0.0) {
        fSinglesDist = std::discrete_distribution<size_t>(fSinglesWeights.begin(), fSinglesWeights.end());
        fSinglesReady = true;
    } else {
        // No valid lines available; mark ready but empty to trigger fallback
        fSinglesReady = true;
    }
}

void PrimaryGeneratorAction::AccumulateTruthLinesDetailed(const std::string& isoSymbol, double weight)
{
    static IsotopeDataLoader loader;
    IsotopeInfo iso;
    if (!loader.Load(isoSymbol, iso)) return;
    if (iso.modes.empty()) return;
    for (const auto& m : iso.modes) {
        const double w = weight * m.branching_ratio;
        if (w <= 0.0) continue;
        for (const auto& gl : m.gammas) {
            if (gl.absolute_intensity > 0.0) {
                TruthLineDetailed tl;
                tl.energy_keV = gl.energy_keV;
                tl.intensity_per_decay = w * gl.absolute_intensity;
                tl.emitter_symbol = iso.symbol;
                tl.emitter_half_life_s = iso.half_life_seconds;
                fTruthLinesDetailed.push_back(std::move(tl));
            }
        }
        // Recurse into daughters only if FullChain is selected
        if (!m.daughter.empty() && fDecayMode == DecayEmissionMode::FullChain) {
            AccumulateTruthLinesDetailed(m.daughter, w);
        }
    }
}
