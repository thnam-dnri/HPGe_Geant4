// ==============================================================================
// IsotopeDecay.hh - Minimal JSON reader for isotope gamma decay chains
// Parses files in ./isotope_data/<Symbol>.json
// Focus: decay_modes[].{mode,daughter,branching_ratio,gammas[]}
// ==============================================================================

#ifndef IsotopeDecay_h
#define IsotopeDecay_h 1

#include <string>
#include <vector>
#include <unordered_map>
#include <mutex>

struct GammaLine {
    double energy_keV {0.0};
    double absolute_intensity {0.0}; // probability per decay (0..1), may exceed 1 for avg multiplicity
    std::string id;                  // optional correlation id
};

struct DecayMode {
    std::string mode;         // e.g., "it", "b-", "ec", "a"
    std::string daughter;     // e.g., "Ba137m"
    double branching_ratio {0.0};
    std::vector<GammaLine> gammas; // gamma emissions for this branch
};

struct IsotopeInfo {
    std::string symbol;       // e.g., "Cs137"
    bool is_stable {false};
    double half_life_seconds {0.0};
    std::vector<DecayMode> modes;
};

class IsotopeDataLoader {
public:
    // Loads isotope JSON from ./isotope_data/<symbol>.json
    // Returns true on success; fills out info
    bool Load(const std::string& symbol, IsotopeInfo& info) const;

private:
    static bool ExtractString(const std::string& s, size_t start, const std::string& key, std::string& out);
    static bool ExtractNumber(const std::string& s, size_t start, const std::string& key, double& out);
    static bool FindArrayRange(const std::string& s, size_t keyPos, char open, char close, size_t& begin, size_t& end);
    static bool FindObjectRange(const std::string& s, size_t from, size_t& objBegin, size_t& objEnd);

    // Simple in-memory cache to avoid re-reading/parsing JSON repeatedly.
    // Thread-safety: guarded by fMutex; scope is the lifetime of the loader instance.
    mutable std::unordered_map<std::string, IsotopeInfo> fCache;
    mutable std::mutex fMutex;
};

#endif
