// ==============================================================================
// IsotopeDecay.cc - Minimal JSON reader for isotope gamma decay chains
// ==============================================================================

#include "IsotopeDecay.hh"

#include <fstream>
#include <sstream>
#include <cctype>

using std::string;

static inline bool isNumberChar(char c) {
    return std::isdigit(static_cast<unsigned char>(c)) || c=='+' || c=='-' || c=='.' || c=='e' || c=='E';
}

bool IsotopeDataLoader::FindArrayRange(const string& s, size_t keyPos, char open, char close, size_t& begin, size_t& end) {
    size_t pos = s.find(open, keyPos);
    if (pos == string::npos) return false;
    int depth = 0;
    for (size_t i = pos; i < s.size(); ++i) {
        if (s[i] == open) depth++;
        else if (s[i] == close) {
            depth--;
            if (depth == 0) { begin = pos; end = i; return true; }
        }
    }
    return false;
}

bool IsotopeDataLoader::FindObjectRange(const string& s, size_t from, size_t& objBegin, size_t& objEnd) {
    size_t pos = s.find('{', from);
    if (pos == string::npos) return false;
    int depth = 0;
    for (size_t i = pos; i < s.size(); ++i) {
        if (s[i] == '{') depth++;
        else if (s[i] == '}') {
            depth--;
            if (depth == 0) { objBegin = pos; objEnd = i; return true; }
        }
    }
    return false;
}

bool IsotopeDataLoader::ExtractString(const string& s, size_t start, const string& key, string& out) {
    size_t k = s.find("\"" + key + "\"", start);
    if (k == string::npos) return false;
    size_t colon = s.find(':', k);
    if (colon == string::npos) return false;
    size_t q1 = s.find('"', colon + 1);
    if (q1 == string::npos) return false;
    size_t q2 = s.find('"', q1 + 1);
    if (q2 == string::npos) return false;
    out = s.substr(q1 + 1, q2 - (q1 + 1));
    return true;
}

bool IsotopeDataLoader::ExtractNumber(const string& s, size_t start, const string& key, double& out) {
    size_t k = s.find("\"" + key + "\"", start);
    if (k == string::npos) return false;
    size_t colon = s.find(':', k);
    if (colon == string::npos) return false;
    size_t i = colon + 1;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = i;
    while (j < s.size() && isNumberChar(s[j])) ++j;
    if (j == i) return false;
    try {
        out = std::stod(s.substr(i, j - i));
        return true;
    } catch (...) {
        return false;
    }
}

bool IsotopeDataLoader::Load(const string& symbol, IsotopeInfo& info) const {
    info = IsotopeInfo{};
    info.symbol = symbol;

    // Check cache first
    {
        std::lock_guard<std::mutex> lock(fMutex);
        auto it = fCache.find(symbol);
        if (it != fCache.end()) {
            info = it->second;
            return true;
        }
    }

    // Resolve candidate paths
    std::vector<string> candidates;
    candidates.push_back(string("isotope_data/") + symbol + ".json");
    candidates.push_back(string("./isotope_data/") + symbol + ".json");
    candidates.push_back(string("../isotope_data/") + symbol + ".json");
#ifdef ISOTOPE_DATA_DIR
    candidates.push_back(string(ISOTOPE_DATA_DIR) + "/" + symbol + ".json");
#endif

    std::ifstream in;
    string chosen;
    for (const auto& p : candidates) {
        in.open(p);
        if (in.is_open()) { chosen = p; break; }
        in.clear();
    }
    if (!in.is_open()) {
        return false;
    }
    std::ostringstream ss; ss << in.rdbuf();
    const string content = ss.str();

    // Detect stability and read half-life seconds if provided
    // Stable examples have: "unit": "stable" and empty decay_modes
    size_t halfLifePos = content.find("\"half_life\"");
    if (halfLifePos != string::npos) {
        string unit;
        if (ExtractString(content, halfLifePos, "unit", unit)) {
            if (unit == "stable") {
                info.is_stable = true;
            }
        }
        double seconds = 0.0;
        if (ExtractNumber(content, halfLifePos, "seconds", seconds)) {
            info.half_life_seconds = seconds;
        }
    }

    // Parse decay modes array
    size_t modesKey = content.find("\"decay_modes\"");
    if (modesKey == string::npos) {
        // If no key, treat as stable
        return true;
    }
    size_t arrBegin = 0, arrEnd = 0;
    if (!FindArrayRange(content, modesKey, '[', ']', arrBegin, arrEnd)) {
        return true; // assume stable / no modes
    }
    const string modesArray = content.substr(arrBegin + 1, arrEnd - arrBegin - 1);

    size_t p = 0;
    while (true) {
        size_t objB = 0, objE = 0;
        if (!FindObjectRange(modesArray, p, objB, objE)) break;
        const string modeObj = modesArray.substr(objB, objE - objB + 1);
        p = objE + 1;

        DecayMode mode;
        ExtractString(modeObj, 0, "mode", mode.mode);
        ExtractString(modeObj, 0, "daughter", mode.daughter);
        ExtractNumber(modeObj, 0, "branching_ratio", mode.branching_ratio);

        // Parse gammas for this mode (if present)
        size_t gKey = modeObj.find("\"gammas\"");
        if (gKey != string::npos) {
            size_t gB = 0, gE = 0;
            if (FindArrayRange(modeObj, gKey, '[', ']', gB, gE)) {
                const string gArray = modeObj.substr(gB + 1, gE - gB - 1);
                size_t q = 0;
                while (true) {
                    size_t gObjB = 0, gObjE = 0;
                    if (!FindObjectRange(gArray, q, gObjB, gObjE)) break;
                    const string gObj = gArray.substr(gObjB, gObjE - gObjB + 1);
                    q = gObjE + 1;

                    GammaLine gl;
                    ExtractNumber(gObj, 0, "energy_keV", gl.energy_keV);
                    ExtractNumber(gObj, 0, "absolute_intensity", gl.absolute_intensity);
                    ExtractString(gObj, 0, "id", gl.id);
                    if (gl.energy_keV > 0.0 && gl.absolute_intensity >= 0.0) {
                        mode.gammas.push_back(gl);
                    }
                }
            }
        }

        info.modes.push_back(std::move(mode));
    }

    // If there are modes, not stable
    if (!info.modes.empty()) info.is_stable = false;

    // Store in cache
    {
        std::lock_guard<std::mutex> lock(fMutex);
        fCache[symbol] = info;
    }
    return true;
}
