
// ==============================================================================
// Run.hh/cc - Custom Run class for analysis
// ==============================================================================

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "globals.hh"
#include <map>

class Run : public G4Run
{
public:
    Run();
    virtual ~Run();

    virtual void Merge(const G4Run*);
    void AddEnergySpectrum(G4double energy);
    void PrintResults() const;

private:
    std::map<G4int, G4int> fEnergyHistogram;  // Simple energy histogram
    G4double fTotalEnergyDeposit;
    G4int fTotalEvents;
    
    static constexpr G4int fNbins = 5000;  // Energy bins (1 keV per bin up to 3 MeV)
    static constexpr G4double fEmax = 5.0; // MeV
    
    G4int EnergyToBin(G4double energy) const;
};
#endif
