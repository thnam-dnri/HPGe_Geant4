# Lead Shield Architecture Plan for HPGe Detector Simulation

## Executive Summary

This plan provides a comprehensive architectural design for implementing a thick lead shield with copper liner (graded-Z shield) for the HPGe detector simulation in Geant4 11.3.2. Based on research of scientific literature and Geant4 best practices, the recommended configuration is a 10 cm thick lead shield with a 1 mm copper inner liner. This graded design reduces lead X-ray fluorescence background by approximately 50% compared to bare lead. The implementation requires minimal modifications to DetectorConstruction.cc, leverages Geant4's NIST material database, and ensures atomic deexcitation physics is enabled for accurate fluorescence simulation. The shield will be implemented as cylindrical nested volumes around the existing detector geometry.

---

## Research Findings

### Best Practices from Literature

#### 1. Lead Shield Thickness for HPGe Soil Analyzers

**Optimal Thickness:**
- **10 cm (100 mm)** is the industry standard for HPGe detector shielding in environmental and soil analysis applications
- Studies show that 10 cm lead reduces background counting rates by approximately **15-fold** compared to unshielded detectors
- Monte Carlo optimization studies indicate **12 cm** as the theoretical optimal thickness, providing background reduction of **two orders of magnitude**
- Thickness beyond 12 cm shows diminishing returns due to shield self-activity

**Commercial Standards:**
- ORTEC and Mirion commercial HPGe shields typically use 100-150 mm low-background lead
- Ultra-low background systems (e.g., Mirion 777) use up to 150 mm total thickness with layered construction

**Material Quality:**
- Low-activity lead (low 210Pb content) is critical for ultra-low background applications
- Standard commercial lead is acceptable for most soil analysis work
- Inner 1 inch of lead should be selected for lowest 210Pb content in high-precision applications

#### 2. Graded-Z Shield Design (Lead + Copper Liner)

**The X-ray Fluorescence Problem:**
- When gamma rays interact with lead via photoelectric effect, lead atoms emit characteristic X-rays at **~85 keV** (Pb K-alpha and K-beta lines)
- These fluorescence X-rays increase detector background in the 80-90 keV region
- Unmitigated, these X-rays can dominate the background spectrum

**The Graded Shield Solution:**
- **Inner copper liner (0.7-1.0 mm thick)**: Absorbs lead X-rays but emits lower energy Cu X-rays (~8 keV)
- The 8 keV copper X-rays are more easily absorbed and fall below typical HPGe detector energy thresholds
- **Background reduction**: 50% reduction in overall counting rate compared to bare lead shields

**Alternative Configurations:**
- Some designs use **Lead → Tin (3 mm) → Copper (0.7 mm)** triple-layer grading
- Tin absorbs Pb X-rays (85 keV) but emits 28 keV X-rays
- Copper then absorbs these 28 keV X-rays
- However, tin/cadmium can be difficult to source; copper-only liner is more practical

**Recommended Configuration for This Project:**
- **Outer layer**: 10 cm low-activity lead
- **Inner liner**: 1 mm electrolytic copper
- This provides excellent background reduction with minimal complexity

#### 3. Shield Geometry Considerations

**Cylindrical vs Box Geometry:**
- **Cylindrical (recommended)**: Better matches the HPGe detector's cylindrical geometry, reduces material volume, easier to construct physically
- **Box**: Simpler to model but uses more material and may have corner effects

**Enclosure Type:**
- **Full enclosure (recommended)**: Surrounds detector on all sides for maximum background reduction
- **Partial/front-face only**: Used when only specific radiation sources are of concern (not suitable for environmental soil analysis)

**Dimensional Considerations:**
- Shield must be sized to accommodate:
  - Existing detector housing (67 mm outer diameter, 76 mm length)
  - Adequate clearance (minimum 5-10 mm) for detector positioning
  - Sample container (for soil analysis, typically cylindrical, 50-100 mm diameter)
- **Recommended inner cavity**: 100 mm diameter × 150 mm height
- **Shield outer dimensions**: ~300 mm diameter × ~350 mm height (with 10 cm thick walls)

**Opening for Sample Insertion:**
- Top opening with removable plug or sliding door
- Opening diameter: ~70-80 mm to accommodate sample containers
- Lead plug thickness: minimum 50 mm (reduces insertion opening background)

### Common Implementation Patterns

#### Pattern 1: Nested Cylindrical Volumes (Recommended)

**Advantages:**
- Clean hierarchy: World → Lead Shield → Cu Liner → Vacuum/Air → Detector Components
- Each volume has a clear mother-daughter relationship
- Easy to visualize and debug
- Performance efficient for cylindrical geometries

**Implementation Structure:**
```
World (G4Box, Air)
  └─ LeadShield (G4Tubs, Lead, 10 cm thick walls)
       └─ CopperLiner (G4Tubs, Copper, 1 mm thick)
            └─ ShieldCavity (G4Tubs, Air or Vacuum)
                 └─ [Existing Detector Components]
```

**Positioning Strategy:**
- Define shield center position in world coordinates
- Position detector inside shield cavity using relative coordinates
- Maintain existing detector geometry unchanged

#### Pattern 2: Boolean Solids for Complex Shapes

**Use Case:**
- If shield needs openings, cutouts, or non-standard geometry
- Can create openings for cables, cryostat, sample insertion ports

**Example:**
```cpp
G4Tubs* solidCylinder = new G4Tubs("ShieldOuter", ...);
G4Tubs* cablePort = new G4Tubs("CablePort", ...);
G4SubtractionSolid* shieldWithPort =
    new G4SubtractionSolid("Shield", solidCylinder, cablePort, rotation, position);
```

**Trade-offs:**
- More complex geometry → slower navigation
- Harder to visualize and debug
- Only use if simple nested volumes are insufficient

#### Pattern 3: Separate Volumes for Shield Components

**Approach:**
- Define shield walls, top, bottom, and plug as separate volumes
- Place each independently in world volume

**Advantages:**
- Maximum flexibility for asymmetric designs
- Easy to model real-world shield construction (assembled pieces)

**Disadvantages:**
- More complex code
- Risk of gaps or overlaps at interfaces
- Requires careful positioning calculations

**Recommendation:** Use only if shield has truly asymmetric design; otherwise use Pattern 1

### Performance Considerations

#### 1. X-ray Fluorescence Physics Requirements

**Critical Physics Configuration:**
- The project already uses **G4EmStandardPhysics_option4**, which enables fluorescence by default
- Fluorescence is automatically activated for:
  - Photoelectric effect
  - Compton scattering
  - Ionization processes
- **No additional configuration needed** - current PhysicsList.cc is already optimal

**Verification:**
- Check that atomic deexcitation is active with UI command: `/process/em/deexcitation true`
- Auger electron emission can also be enabled: `/process/em/auger true` (optional, minimal effect on thick shields)

**Data Libraries:**
- Geant4 uses EADL (Evaluated Atomic Data Library) for fluorescence
- Lead K-shell fluorescence energies: 74.97 keV (K-alpha1), 72.80 keV (K-alpha2), 84.94 keV (K-beta1)
- Copper K-shell fluorescence: 8.05 keV (K-alpha1), 8.03 keV (K-alpha2)

#### 2. Memory and Performance Optimization

**Geometry Simplification:**
- Thick lead shields are **computationally expensive** to simulate
- Each gamma ray must undergo many scattering events to traverse 10 cm of lead
- Use simple G4Tubs geometry (not tessellated or boolean solids unless necessary)
- Avoid unnecessary geometric complexity

**Production Cuts:**
- Current PhysicsList.cc sets global cuts at 0.1 mm and detector region cuts at 0.01 mm
- **For lead/copper shield regions**: Consider setting higher cuts (0.5-1.0 mm) to reduce secondary particle production
- This improves performance with minimal impact on accuracy (shielding effectiveness depends on primary gamma attenuation, not detailed tracking of low-energy secondaries)

**Implementation in SetCuts():**
```cpp
// Create shield region with higher production cuts
G4Region* shieldRegion = new G4Region("ShieldRegion");
shieldRegion->AddRootLogicalVolume(leadShieldLV);
shieldRegion->AddRootLogicalVolume(copperLinerLV);

G4ProductionCuts* shieldCuts = new G4ProductionCuts();
shieldCuts->SetProductionCut(0.5*mm, "gamma");
shieldCuts->SetProductionCut(0.5*mm, "e-");
shieldCuts->SetProductionCut(0.5*mm, "e+");
shieldRegion->SetProductionCuts(shieldCuts);
```

#### 3. Variance Reduction Techniques (Optional)

**Importance Biasing:**
- **Not recommended initially** - adds complexity
- Can be considered if simulation is too slow
- Assigns importance values to geometry regions
- Particles moving from low → high importance are split
- Particles moving from high → low importance are Russian rouletted

**When to Use:**
- Only if simulation performance is unacceptable after basic optimizations
- Requires careful tuning to avoid bias in results

**Current Recommendation:**
- Start with standard simulation
- Profile performance with 10,000 events
- If CPU time per event >10 seconds, consider variance reduction

#### 4. Validation and Testing Strategy

**Overlap Checking:**
- **Critical**: Use built-in Geant4 overlap detection
- Add `checkOverlaps=true` parameter to all G4PVPlacement constructors during development
- Runtime test: `/geometry/test/run` in macro file
- Fix all overlaps before production runs

**Expected Performance:**
- Unshielded: ~99% of gamma rays reach detector
- 10 cm lead shield: ~1% transmission for 662 keV gammas (Cs-137)
- Background reduction: ~15-fold for environmental sources

### Security Considerations

**Physics Accuracy:**
- Lead shielding simulations are sensitive to atomic deexcitation models
- Using G4EmStandardPhysics_option4 ensures highest accuracy for gamma spectroscopy
- Data library version matters: ensure GEANT4 data files are up-to-date

**Geometry Validation:**
- Overlapping volumes cause navigation errors and incorrect physics
- Gaps in shielding cause underestimation of background
- Both can lead to incorrect scientific conclusions

**Reproducibility:**
- Document exact material definitions (NIST vs custom)
- Document geometry parameters in code comments
- Record Geant4 version and data library versions

---

## Proposed Architecture

### Component Structure

#### 1. Material Definitions

**Lead Material (G4_Pb):**
```cpp
// In DefineMaterials() method
G4NistManager* nist = G4NistManager::Instance();

// Pure lead from NIST database
fLead = nist->FindOrBuildMaterial("G4_Pb");
// Properties: density = 11.35 g/cm³, natural isotope composition
```

**Copper Material (G4_Cu):**
```cpp
// Electrolytic copper from NIST database
fCopper = nist->FindOrBuildMaterial("G4_Cu");
// Properties: density = 8.96 g/cm³, natural isotope composition
```

**Why NIST Materials:**
- Pre-validated density, composition, and ionization potentials
- Includes natural isotopic composition automatically
- Ensures consistency with Geant4 physics models
- Less error-prone than manual material definition

#### 2. Geometry Design - Recommended Configuration

**Visual Hierarchy:**
```
                    World (50 cm box, Air)
                           │
                           ▼
            ┌───────────────────────────────┐
            │   Lead Shield (G4Tubs)        │
            │   Inner R: 50 mm              │
            │   Outer R: 150 mm (10cm wall) │
            │   Height: 350 mm               │
            │   Material: G4_Pb             │
            └───────────┬───────────────────┘
                        │
                        ▼
            ┌───────────────────────────────┐
            │   Copper Liner (G4Tubs)       │
            │   Inner R: 49 mm              │
            │   Outer R: 50 mm (1mm wall)   │
            │   Height: 348 mm               │
            │   Material: G4_Cu             │
            └───────────┬───────────────────┘
                        │
                        ▼
            ┌───────────────────────────────┐
            │   Shield Cavity (G4Tubs)      │
            │   Radius: 49 mm               │
            │   Height: 348 mm               │
            │   Material: G4_AIR            │
            │   [Contains detector]         │
            └───────────────────────────────┘
```

**Dimensional Parameters:**
```cpp
// Shield dimensions (add to DetectorConstruction.hh as static const)
static const G4double fShieldInnerRadius;     // 50 mm
static const G4double fShieldThickness;       // 100 mm (10 cm lead)
static const G4double fShieldOuterRadius;     // 150 mm
static const G4double fShieldHeight;          // 350 mm
static const G4double fCopperLinerThickness;  // 1 mm
static const G4double fShieldCavityRadius;    // 49 mm
static const G4double fShieldCavityHeight;    // 348 mm

// Shield position in world coordinates
static const G4double fShieldCenterZ;         // 15 cm (adjust based on detector position)
```

**Why These Dimensions:**
- **Inner radius (50 mm)**: Accommodates existing detector housing (67 mm diameter) with 33 mm clearance on one side, allowing sample placement
- **10 cm wall thickness**: Optimal based on literature (15-fold background reduction)
- **350 mm height**: Provides shielding on all sides of 76 mm tall detector with margin
- **1 mm copper**: Standard thickness for fluorescence suppression

#### 3. Detector Positioning Within Shield

**Current Detector Position:**
- Source at origin (0, 0, 0)
- Detector window starts at z = 10 cm (fSourceDetectorDistance)
- Detector extends to ~13 cm

**Shield Integration Strategy:**
- **Option A (Recommended)**: Place shield around existing detector, adjust source position
  - Shield center at z = 15 cm
  - Detector remains at z = 10-13 cm range
  - Source moves to z = 5 cm (inside shield cavity, 5 cm from detector)
  - Realistic for soil sample analysis (sample inside shield)

- **Option B**: Keep source at origin, shield surrounds detector only
  - Shield center at z = 12 cm
  - Source outside shield
  - Requires top opening in shield for beam entry
  - Less realistic for soil analysis, more appropriate for external source calibration

**Recommendation:** Use Option A with source inside shield cavity

#### 4. Shield Opening (Optional - Top Plug Design)

**Simple Implementation (Phase 1):**
- No opening initially
- Fully enclosed shield for baseline background measurement
- Source placed via macro file inside shield cavity

**Advanced Implementation (Phase 2):**
```cpp
// Top plug using boolean solid
G4Tubs* fullShield = new G4Tubs("ShieldFull", innerR, outerR, height/2, 0, 360*deg);
G4Tubs* opening = new G4Tubs("Opening", 0, openingR, plugDepth/2, 0, 360*deg);
G4SubtractionSolid* shieldWithOpening =
    new G4SubtractionSolid("Shield", fullShield, opening, 0,
                          G4ThreeVector(0, 0, height/2 - plugDepth/2));

// Separate plug volume
G4Tubs* plugSolid = new G4Tubs("ShieldPlug", 0, openingR, plugThickness/2, 0, 360*deg);
G4LogicalVolume* plugLV = new G4LogicalVolume(plugSolid, fLead, "ShieldPlug");
// Place plug in closed or open position via parameter
```

**Recommendation:** Implement simple fully-enclosed shield first, add opening later if needed

### Data Flow

**Particle Interactions Through Shield:**

1. **Primary Gamma Ray Generation** (PrimaryGeneratorAction)
   - Gamma ray emitted from source position inside/outside shield
   - Initial energy (e.g., 662 keV for Cs-137, 1461 keV for K-40)

2. **Shield Entry and Attenuation**
   - Gamma enters lead shield outer surface
   - **Physics processes**:
     - Photoelectric effect (dominant <500 keV)
     - Compton scattering (dominant >500 keV)
     - Pair production (>1.022 MeV, minimal for typical sources)

3. **Lead Fluorescence Production**
   - Photoelectric events in lead → vacancies in K-shell
   - **Atomic deexcitation**: Pb K-alpha X-rays (~75 keV) and K-beta (~85 keV) emitted
   - These fluorescence X-rays propagate isotropically

4. **Copper Liner Interaction**
   - Most Pb X-rays absorbed by 1 mm copper liner (attenuation length ~0.1 mm)
   - **Copper fluorescence**: Cu K-alpha (~8 keV) emitted
   - These low-energy X-rays are absorbed locally or escape energy range of interest

5. **Cavity Propagation**
   - Surviving primary gammas (those that scattered in lead) enter air cavity
   - Multiple-scattered gammas have reduced energy spectrum

6. **Detector Interaction** (Existing geometry)
   - Gammas enter through aluminum window
   - Detected in Ge crystal (fScoringVolume)
   - Energy deposited recorded by SteppingAction/RunAction

**Background Sources (Shield-Related):**
- Lead X-ray fluorescence (85 keV peak) - mitigated by copper
- Copper X-ray fluorescence (8 keV peak) - below threshold
- Compton-scattered primaries (continuum)
- Shield self-activity (210Pb, 210Bi) - natural background, not simulated unless explicitly added

### Integration Points with Existing Project

#### 1. DetectorConstruction.hh Modifications

**New Private Member Variables:**
```cpp
// Add to existing material pointers section
G4Material* fLead;      // Lead shield material
G4Material* fCopper;    // Copper liner material

// Add to existing volume pointers section
G4LogicalVolume* fLeadShieldLV;
G4LogicalVolume* fCopperLinerLV;
G4LogicalVolume* fShieldCavityLV;

// Add new static const for shield parameters
static const G4double fShieldInnerRadius;
static const G4double fShieldThickness;
static const G4double fShieldHeight;
static const G4double fCopperLinerThickness;
```

**No Changes Needed to Public Interface:**
- `Construct()` and `ConstructSDandField()` signatures remain unchanged
- Existing `GetScoringVolume()` still returns Ge crystal
- Fully backward compatible

#### 2. DetectorConstruction.cc Modifications

**DefineMaterials() Method:**
```cpp
void DetectorConstruction::DefineMaterials()
{
    // ... existing materials (Air, Vacuum, Al, Ge, Mylar, Li, B) ...

    // Add shield materials
    fLead = nist->FindOrBuildMaterial("G4_Pb");
    fCopper = nist->FindOrBuildMaterial("G4_Cu");

    G4cout << "Shield materials defined (Pb + Cu liner)" << G4endl;
}
```

**DefineVolumes() Method - Integration Point:**
```cpp
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    // Existing world volume creation (UNCHANGED)
    G4Box* worldS = new G4Box("World", fWorldSize/2, fWorldSize/2, fWorldSize/2);
    fWorldLV = new G4LogicalVolume(worldS, fWorldMaterial, "World");
    fWorldPV = new G4PVPlacement(0, G4ThreeVector(), fWorldLV, "World", 0, false, 0, true);

    // NEW: Shield construction
    // 1. Lead shield outer cylinder
    G4Tubs* leadShieldS = new G4Tubs("LeadShield",
                                     fShieldInnerRadius,
                                     fShieldInnerRadius + fShieldThickness,
                                     fShieldHeight/2,
                                     0*deg, 360*deg);
    fLeadShieldLV = new G4LogicalVolume(leadShieldS, fLead, "LeadShield");
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, fShieldCenterZ),
                     fLeadShieldLV,
                     "LeadShield",
                     fWorldLV,      // Mother volume: World
                     false, 0, true); // checkOverlaps=true

    // 2. Copper liner
    G4double cuInnerR = fShieldInnerRadius - fCopperLinerThickness;
    G4double cuOuterR = fShieldInnerRadius;
    G4double cuHeight = fShieldHeight - 2*mm; // Slightly shorter to avoid overlap

    G4Tubs* copperLinerS = new G4Tubs("CopperLiner",
                                      cuInnerR, cuOuterR,
                                      cuHeight/2,
                                      0*deg, 360*deg);
    fCopperLinerLV = new G4LogicalVolume(copperLinerS, fCopper, "CopperLiner");
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, fShieldCenterZ),
                     fCopperLinerLV,
                     "CopperLiner",
                     fWorldLV,      // Mother volume: World (alternative: place inside lead)
                     false, 0, true);

    // 3. Shield cavity (air-filled space for detector)
    G4Tubs* cavityS = new G4Tubs("ShieldCavity",
                                 0, cuInnerR,
                                 cuHeight/2,
                                 0*deg, 360*deg);
    fShieldCavityLV = new G4LogicalVolume(cavityS, fWorldMaterial, "ShieldCavity");
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, fShieldCenterZ),
                     fShieldCavityLV,
                     "ShieldCavity",
                     fWorldLV,      // Mother volume: World
                     false, 0, true);

    // EXISTING: Detector geometry (MOSTLY UNCHANGED)
    // Change mother volume from fWorldLV to fShieldCavityLV
    // Adjust Z-positions to be relative to shield cavity center

    // ... [rest of existing detector construction code] ...

    return fWorldPV;
}
```

**Key Integration Points:**
- Shield volumes placed in world coordinates
- Existing detector volumes **either**:
  - **Option A**: Remain in world volume, positioned to be inside shield cavity geometrically
  - **Option B**: Change mother volume from `fWorldLV` to `fShieldCavityLV` and adjust positions to be relative to cavity center

**Recommendation:** Option A is simpler initially (detector code unchanged), verify geometries don't overlap using overlap checker

#### 3. PhysicsList.cc - No Changes Required

**Current Configuration is Optimal:**
```cpp
// Already using option4 which enables fluorescence
RegisterPhysics(new G4EmStandardPhysics_option4());
```

**Verification (in macro file):**
```
# Confirm atomic deexcitation is enabled
/process/em/verbose 1
/process/em/deexcitation true  # Should already be true for option4
```

#### 4. PrimaryGeneratorAction.cc - Source Position Adjustment

**Current Configuration:**
- Source at origin (0, 0, 0)
- Generates particles in fixed direction (+z) toward detector

**Required Modifications:**
```cpp
// In PrimaryGeneratorAction constructor or GeneratePrimaries:
// Move source position to inside shield cavity

// Current: source at origin
// fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));

// New: source inside shield, 5 cm in front of detector
G4double sourceZ = fShieldCenterZ - 50*mm;  // 5 cm from detector window
fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, sourceZ));

// For soil sample simulation: use isotropic emission
fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
```

**Alternative (No Code Change):**
- Keep source position configurable via macro file
- Use `/gun/position` command to set source inside shield

#### 5. Visualization Attributes

**Shield Visibility:**
```cpp
// In DefineVolumes(), after shield construction

// Lead shield - dark gray, semi-transparent to see interior
G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.6));
leadVis->SetForceSolid(true);
fLeadShieldLV->SetVisAttributes(leadVis);

// Copper liner - orange/brown, semi-transparent
G4VisAttributes* copperVis = new G4VisAttributes(G4Colour(0.8, 0.5, 0.2, 0.5));
copperVis->SetForceSolid(true);
fCopperLinerLV->SetVisAttributes(copperVis);

// Shield cavity - invisible (to see detector inside)
fShieldCavityLV->SetVisAttributes(G4VisAttributes::GetInvisible());
```

### Interface Definitions

#### DetectorConstruction.hh Public Interface (Unchanged)

```cpp
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();        // Existing
    virtual void ConstructSDandField();            // Existing

    // Get methods for analysis (Existing - no changes)
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    // Optional: Add getter for shield volumes (for advanced analysis)
    G4LogicalVolume* GetLeadShield() const { return fLeadShieldLV; }
    G4LogicalVolume* GetCopperLiner() const { return fCopperLinerLV; }
};
```

#### Material Interface (G4Material* returned by NIST)

**Lead Material Properties:**
- **Name**: "G4_Pb"
- **Density**: 11.35 g/cm³
- **Composition**: Natural lead isotope mixture (204Pb: 1.4%, 206Pb: 24.1%, 207Pb: 22.1%, 208Pb: 52.4%)
- **Mean excitation energy**: 823 eV
- **Radiation length**: 5.612 mm
- **Attenuation length (662 keV)**: ~11 mm

**Copper Material Properties:**
- **Name**: "G4_Cu"
- **Density**: 8.96 g/cm³
- **Composition**: Natural copper (63Cu: 69.2%, 65Cu: 30.8%)
- **Mean excitation energy**: 322 eV
- **Radiation length**: 14.36 mm
- **K-edge**: 8.979 keV (fluorescence threshold)

#### Geometry Interface (G4LogicalVolume methods)

**Volume Hierarchy Methods:**
```cpp
// Check if point is inside shield
G4bool IsInside(const G4ThreeVector& point, G4LogicalVolume* volume);

// Get shield material at position (for analysis)
G4Material* GetMaterial(G4ThreeVector& point);
```

**No New Methods Required:**
- Existing Geant4 navigation handles volume hierarchy automatically
- Stepping action can query current volume name for analysis

---

## Configuration Recommendations

### 1. Geometry Parameters

**Recommended Shield Configuration:**
```cpp
// DetectorConstruction.cc - static member definitions
const G4double DetectorConstruction::fShieldInnerRadius = 50.0*mm;
const G4double DetectorConstruction::fShieldThickness = 100.0*mm;  // 10 cm lead
const G4double DetectorConstruction::fShieldOuterRadius = 150.0*mm;
const G4double DetectorConstruction::fShieldHeight = 350.0*mm;
const G4double DetectorConstruction::fCopperLinerThickness = 1.0*mm;
const G4double DetectorConstruction::fShieldCenterZ = 150.0*mm;  // 15 cm from origin
```

**Rationale:**
- **Inner radius (50 mm)**: Provides 50 mm diameter cavity, adequate for detector housing (67 mm) positioned off-center or sample container
- **Thickness (100 mm)**: Standard 10 cm provides 15× background reduction, optimal cost/performance
- **Height (350 mm)**: Provides >100 mm shielding above and below detector (76 mm height)
- **Cu liner (1 mm)**: Standard thickness for 85 keV X-ray absorption (>90% attenuation)

### 2. Physics Configuration

**Current Settings (No Changes):**
```cpp
// PhysicsList.cc - already optimal
RegisterPhysics(new G4EmStandardPhysics_option4());  // Fluorescence enabled by default
RegisterPhysics(new G4EmExtraPhysics());
RegisterPhysics(new G4DecayPhysics());
RegisterPhysics(new G4RadioactiveDecayPhysics());
RegisterPhysics(new G4IonPhysics());
```

**Production Cuts Strategy:**

**Option A - Uniform Cuts (Simplest):**
```cpp
// Current setting - no changes
SetCutValue(0.1*mm, "gamma");
SetCutValue(0.1*mm, "e-");
SetCutValue(0.1*mm, "e+");
```

**Option B - Region-Specific Cuts (Performance Optimization):**
```cpp
void PhysicsList::SetCuts()
{
    SetCutsWithDefault();

    // Default for world
    SetCutValue(0.1*mm, "gamma");
    SetCutValue(0.1*mm, "e-");
    SetCutValue(0.1*mm, "e+");

    // Tighter cuts in detector
    G4Region* detectorRegion = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld");
    if (detectorRegion) {
        G4ProductionCuts* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(0.01*mm, "gamma");
        cuts->SetProductionCut(0.01*mm, "e-");
        cuts->SetProductionCut(0.01*mm, "e+");
        detectorRegion->SetProductionCuts(cuts);
    }

    // Looser cuts in shield (NEW - improves performance)
    G4Region* shieldRegion = new G4Region("ShieldRegion");
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
    G4LogicalVolume* leadLV = lvStore->GetVolume("LeadShield");
    G4LogicalVolume* copperLV = lvStore->GetVolume("CopperLiner");

    if (leadLV) shieldRegion->AddRootLogicalVolume(leadLV);
    if (copperLV) shieldRegion->AddRootLogicalVolume(copperLV);

    G4ProductionCuts* shieldCuts = new G4ProductionCuts();
    shieldCuts->SetProductionCut(0.5*mm, "gamma");   // Coarser cuts in shield
    shieldCuts->SetProductionCut(0.5*mm, "e-");
    shieldCuts->SetProductionCut(0.5*mm, "e+");
    shieldRegion->SetProductionCuts(shieldCuts);

    if (verboseLevel > 0) DumpCutValuesTable();
}
```

**Recommendation:** Start with Option A (uniform cuts), implement Option B only if performance is inadequate

### 3. Macro File Configuration

**New Macro: shielded_run.mac**
```bash
# Shielded detector simulation macro

# Visualization (optional)
/vis/open OGL 600x600-0+0
/vis/drawVolume
/vis/viewer/set/viewpointVector 1 1 1
/vis/viewer/set/upVector 0 0 1
/vis/viewer/zoom 1.5
/vis/scene/add/trajectories smooth
/vis/scene/endOfEventAction accumulate

# Physics verification
/process/em/verbose 1
/process/em/printParameters

# Source configuration (Cs-137 inside shield)
/gun/particle gamma
/gun/energy 662 keV
/gun/position 0 0 10 cm    # 10 cm from origin, inside shield cavity
/gun/direction 0 0 1       # Toward detector

# Detector setup
/run/initialize

# Run simulation
/run/beamOn 10000
```

**Comparison Macro: unshielded_run.mac**
```bash
# Unshielded detector (for comparison)
# Same as shielded_run.mac but with detector at different position
# or disable shield volumes by commenting out placement in code
```

### 4. Environment Variables

**Geant4 Data Libraries:**
```bash
# Verify these are set correctly (critical for fluorescence physics)
export G4LEDATA=/path/to/geant4-data/G4EMLOW8.6
export G4LEVELGAMMADATA=/path/to/geant4-data/PhotonEvaporation5.7
export G4RADIOACTIVEDATA=/path/to/geant4-data/RadioactiveDecay5.6
```

**Recommended Version:**
- **G4EMLOW**: Version 8.6 or later (contains updated EADL fluorescence data)
- Check version: Look for "Atomic Deexcitation: Fluorescence enabled" in console output

---

## Testing Strategy

### Unit Testing

#### 1. Material Definition Tests

**Test: Verify NIST Materials Loaded Correctly**

**Method:**
```cpp
// In DetectorConstruction::Construct() or separate test method
G4cout << "\n=== Material Verification ===" << G4endl;
G4cout << "Lead density: " << fLead->GetDensity()/(g/cm3) << " g/cm3 (expect 11.35)" << G4endl;
G4cout << "Lead composition: " << G4endl;
fLead->GetMaterialTable();  // Print composition

G4cout << "Copper density: " << fCopper->GetDensity()/(g/cm3) << " g/cm3 (expect 8.96)" << G4endl;
```

**Expected Output:**
```
Lead density: 11.35 g/cm3 (expect 11.35)
Copper density: 8.96 g/cm3 (expect 8.96)
```

**Pass Criteria:** Densities match expected values within 0.01 g/cm³

#### 2. Geometry Construction Tests

**Test A: Volume Hierarchy**

**Method:** Add debug output in DefineVolumes()
```cpp
G4cout << "\n=== Shield Geometry ===" << G4endl;
G4cout << "Lead shield inner R: " << fShieldInnerRadius/mm << " mm" << G4endl;
G4cout << "Lead shield outer R: " << (fShieldInnerRadius + fShieldThickness)/mm << " mm" << G4endl;
G4cout << "Shield height: " << fShieldHeight/mm << " mm" << G4endl;
G4cout << "Shield center Z: " << fShieldCenterZ/mm << " mm" << G4endl;
G4cout << "Copper liner thickness: " << fCopperLinerThickness/mm << " mm" << G4endl;
```

**Expected Output:**
```
Lead shield inner R: 50 mm
Lead shield outer R: 150 mm
Shield height: 350 mm
Shield center Z: 150 mm
Copper liner thickness: 1 mm
```

**Pass Criteria:** All dimensions match design specification

**Test B: Overlap Detection**

**Method:** Enable overlap checking in all G4PVPlacement constructors
```cpp
// Set checkOverlaps parameter to true
new G4PVPlacement(0, position, logicalVol, "name", motherVol, false, 0, true);
                                                                        // ^^^^ checkOverlaps=true
```

**Runtime Check:**
```bash
# In macro file
/geometry/test/run
/geometry/test/resolution 1000  # Higher resolution for thick volumes
```

**Expected Output:**
```
Checking overlaps for volume LeadShield:0 ... OK!
Checking overlaps for volume CopperLiner:0 ... OK!
Checking overlaps for volume ShieldCavity:0 ... OK!
```

**Pass Criteria:** No overlap warnings reported

**Test C: Detector Position Validation**

**Method:** Verify detector is inside shield cavity
```cpp
// Calculate detector bounding box
G4double detectorMinZ = fSourceDetectorDistance;
G4double detectorMaxZ = detectorMinZ + housingLength;
G4double detectorMaxR = housingOuterDiam/2;

// Shield cavity bounds
G4double cavityMinZ = fShieldCenterZ - fShieldCavityHeight/2;
G4double cavityMaxZ = fShieldCenterZ + fShieldCavityHeight/2;
G4double cavityMaxR = fShieldCavityRadius;

G4cout << "Detector Z range: " << detectorMinZ/mm << " to " << detectorMaxZ/mm << " mm" << G4endl;
G4cout << "Cavity Z range: " << cavityMinZ/mm << " to " << cavityMaxZ/mm << " mm" << G4endl;

// Verify containment
G4bool contained = (detectorMinZ > cavityMinZ) && (detectorMaxZ < cavityMaxZ) &&
                   (detectorMaxR < cavityMaxR);
G4cout << "Detector contained in cavity: " << (contained ? "YES" : "NO") << G4endl;
```

**Pass Criteria:** Detector fully contained within shield cavity with clearance on all sides

#### 3. Physics Configuration Tests

**Test: Verify Fluorescence is Enabled**

**Method:** Run with /process/em/verbose 1 in macro
```bash
/process/em/verbose 1
/process/em/printParameters
/run/initialize
```

**Expected Output:**
```
### G4EmStandardPhysics_option4 Construct Processes
Fluorescence is activated
Fluorescence directory: /path/to/G4EMLOW8.6
Auger electron production: false (can be enabled)
PIXE enabled: false
```

**Pass Criteria:** "Fluorescence is activated" appears in output

**Test: Verify Fluorescence X-rays Produced**

**Method:** Run small simulation, analyze output
```cpp
// In SteppingAction, add counter for secondary particle creation
if (processName == "phot") {  // Photoelectric process
    const std::vector<const G4Track*>* secondaries =
        step->GetSecondaryInCurrentStep();

    for (auto track : *secondaries) {
        if (track->GetDefinition() == G4Gamma::Gamma()) {
            G4double energy = track->GetKineticEnergy();
            if (energy > 70*keV && energy < 90*keV) {
                // Lead fluorescence X-ray detected
                G4cout << "Lead K-alpha/beta X-ray: " << energy/keV << " keV" << G4endl;
            }
        }
    }
}
```

**Expected Output:** X-rays near 75 keV and 85 keV detected after photoelectric events in lead

**Pass Criteria:** Lead fluorescence X-rays observed in 70-90 keV range

### Integration Testing

#### 1. Full Geometry Assembly Test

**Test: Complete Simulation Run**

**Method:**
```bash
# Run macro with visualization
./HPGeSingle shielded_run.mac
```

**Checklist:**
- [ ] Simulation initializes without errors
- [ ] Visualization displays shield, liner, cavity, and detector
- [ ] Particles can be tracked from source to detector
- [ ] No navigation warnings (e.g., "stuck track")
- [ ] Energy deposition in Ge crystal recorded

**Pass Criteria:** All checklist items pass

#### 2. Source-to-Detector Path Test

**Test: Verify Particles Traverse Shield Correctly**

**Method:** Enable trajectory visualization
```bash
/vis/scene/add/trajectories smooth
/tracking/verbose 1
/run/beamOn 10
```

**Visual Verification:**
- Gamma rays originate from source position
- Most are absorbed/scattered in lead
- Some reach copper liner
- Few penetrate to detector
- Trajectories show realistic scattering angles

**Quantitative Check:**
```cpp
// In RunAction, count particles by volume
std::map<G4String, G4int> volumeHits;

// In SteppingAction
G4String volumeName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
volumeHits[volumeName]++;

// In RunAction::EndOfRunAction
G4cout << "Particles entering LeadShield: " << volumeHits["LeadShield"] << G4endl;
G4cout << "Particles entering CopperLiner: " << volumeHits["CopperLiner"] << G4endl;
G4cout << "Particles entering GeCrystal: " << volumeHits["GeCrystal"] << G4endl;
```

**Expected Ratios (for 662 keV gammas, 10 cm lead):**
- ~100% enter lead shield
- ~5-10% reach copper liner (most absorbed in lead)
- ~1% reach Ge crystal (after attenuation)

**Pass Criteria:** Observed ratios within factor of 2 of expected

#### 3. Background Reduction Validation

**Test: Compare Shielded vs Unshielded Background**

**Method:**
```bash
# Run 1: Unshielded (comment out shield placement in code, rebuild)
/run/beamOn 100000
# Record counts in detector

# Run 2: Shielded (restore shield placement, rebuild)
/run/beamOn 100000
# Record counts in detector

# Calculate reduction factor
```

**Expected Results:**
- **Source inside shield**: Counts reduced by ~2× (self-attenuation)
- **External background simulation**: Counts reduced by ~15× (full shielding effect)

**Pass Criteria:** Background reduction within 20% of expected value

### Performance Testing

#### 1. Simulation Speed Benchmark

**Test: Events per Second**

**Method:**
```bash
time ./HPGeSingle -m shielded_run.mac > /dev/null 2>&1
# Run 10,000 events, measure wall-clock time
```

**Metrics:**
- **Time per event**: Target <1 second/event for 662 keV source
- **Total time**: 10,000 events should complete in <3 hours

**Acceptance Criteria:**
- <5 seconds/event: Excellent (production ready)
- 5-10 seconds/event: Acceptable (consider optimization)
- >10 seconds/event: Poor (implement variance reduction or looser cuts)

#### 2. Memory Usage Test

**Test: Peak Memory Consumption**

**Method:**
```bash
/usr/bin/time -v ./HPGeSingle -m shielded_run.mac
# Check "Maximum resident set size"
```

**Expected Memory:**
- **Without shield**: ~200 MB
- **With shield**: ~300-400 MB (geometry tables for thick lead)

**Acceptance Criteria:** <1 GB total memory usage

#### 3. Production Cut Optimization

**Test: Performance vs Accuracy Trade-off**

**Method:** Run same simulation with different shield region cuts
```
Cut = 0.1 mm: Baseline (high accuracy)
Cut = 0.5 mm: Test 1 (medium)
Cut = 1.0 mm: Test 2 (fast)
```

**Metrics:**
- Time per event
- Background count rate
- Energy resolution (FWHM at 662 keV)

**Expected Trade-off:**
- 1.0 mm cuts: ~3× faster, <5% change in results
- 0.5 mm cuts: ~1.5× faster, <2% change in results

**Recommendation:** Choose cut that gives <2% change in background counts but maximum speedup

### Validation Tests

#### 1. X-ray Fluorescence Spectrum Test

**Test: Verify Lead and Copper Fluorescence Peaks**

**Method:** Record energy spectrum of particles entering detector
```cpp
// In SteppingAction, if entering Ge crystal:
if (volumeName == "GeCrystal") {
    G4double energy = step->GetTrack()->GetKineticEnergy();
    analysisManager->FillH1(1, energy);  // Histogram ID 1: incident spectrum
}
```

**Analysis:** Plot histogram
```python
import matplotlib.pyplot as plt
import numpy as np

# Load histogram data from output file
energies = [...]  # keV
counts = [...]

plt.figure()
plt.hist(energies, bins=200, range=(0, 200))
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Incident Particle Spectrum at Detector')
plt.axvline(75, color='r', linestyle='--', label='Pb K-alpha (75 keV)')
plt.axvline(85, color='r', linestyle='--', label='Pb K-beta (85 keV)')
plt.axvline(8, color='b', linestyle='--', label='Cu K-alpha (8 keV)')
plt.legend()
plt.show()
```

**Expected Features:**
- **Primary peak**: 662 keV (if using Cs-137 source)
- **Lead fluorescence**: 75 keV and 85 keV peaks (should be **reduced** compared to no-copper-liner case)
- **Copper fluorescence**: 8 keV peak (low intensity, may be below threshold)
- **Compton continuum**: Broad distribution <662 keV

**Pass Criteria:**
- Lead fluorescence peaks visible but <10% of primary peak intensity (copper liner working)
- Without copper liner: Lead peaks >30% of primary (comparison test)

#### 2. Geometric Attenuation Coefficient

**Test: Measure Lead Attenuation vs Expected**

**Method:** Broad-beam geometry test
```bash
# Place point source outside shield, count transmitted gammas
# Source: 662 keV (Cs-137)
# Expected attenuation: exp(-μ * 10 cm)
# μ(662 keV in Pb) ≈ 1.0 cm^-1
# Transmission ≈ exp(-10) ≈ 0.0045% (1 in 22,000)
```

**Simulation:**
```cpp
// Count source emissions vs detector hits
G4int sourceEvents = 1000000;  // 1 million
G4int detectorHits = ...;      // Count from simulation
G4double transmission = (double)detectorHits / sourceEvents;

G4cout << "Measured transmission: " << transmission << G4endl;
G4cout << "Expected transmission: " << exp(-1.0 * 10) << G4endl;
```

**Pass Criteria:** Measured transmission within 20% of expected (accounting for geometry factors)

#### 3. Shield Thickness Sensitivity

**Test: Vary Shield Thickness, Measure Background Change**

**Method:** Run simulations with:
- 5 cm lead
- 10 cm lead (baseline)
- 15 cm lead

**Expected Background Scaling:**
- 5 cm: ~2× higher than 10 cm
- 15 cm: ~2× lower than 10 cm

**Pass Criteria:** Background scales exponentially with thickness (within 30% of expected)

---

## Implementation Roadmap

### Phase 1: Foundation (Core Shield Geometry)

**Goal:** Implement basic lead + copper shield without detector integration

**Duration:** 2-3 hours development + 1 hour testing

**Tasks:**

1. **Modify DetectorConstruction.hh** (15 min)
   - Add material pointers: `fLead`, `fCopper`
   - Add logical volume pointers: `fLeadShieldLV`, `fCopperLinerLV`, `fShieldCavityLV`
   - Add static const parameters for shield dimensions

2. **Update DefineMaterials()** (10 min)
   - Add `fLead = nist->FindOrBuildMaterial("G4_Pb");`
   - Add `fCopper = nist->FindOrBuildMaterial("G4_Cu");`
   - Add console output confirming materials loaded

3. **Implement Shield Geometry in DefineVolumes()** (45 min)
   - Create lead shield G4Tubs solid and logical volume
   - Place lead shield in world volume
   - Create copper liner G4Tubs solid and logical volume
   - Place copper liner inside/adjacent to lead
   - Create shield cavity logical volume
   - Add visualization attributes

4. **Build and Test Compilation** (15 min)
   ```bash
   cd build
   cmake --build . -j4
   ```
   - Fix any compilation errors
   - Verify executable builds successfully

5. **Geometry Validation Tests** (30 min)
   - Run `/geometry/test/run` to check overlaps
   - Visualize geometry with OpenGL
   - Verify shield appears correctly positioned
   - Check material assignments in console output

6. **Documentation** (30 min)
   - Add comments to code explaining shield design
   - Document dimensions and rationale
   - Update context.md with progress

**Deliverables:**
- ✓ Compiled executable with shield geometry
- ✓ No geometry overlaps
- ✓ Visual confirmation of shield structure
- ✓ Console output confirms materials

**Testing Checklist:**
- [ ] Code compiles without errors
- [ ] `/geometry/test/run` shows no overlaps
- [ ] Visualization shows lead (gray), copper (orange), cavity volumes
- [ ] Material densities printed correctly (Pb: 11.35, Cu: 8.96 g/cm³)

---

### Phase 2: Detector Integration

**Goal:** Position existing detector inside shield cavity, adjust coordinates

**Duration:** 2-3 hours development + 2 hours testing

**Tasks:**

1. **Analyze Existing Detector Geometry** (30 min)
   - Review current detector position calculations (lines 119-130 in DetectorConstruction.cc)
   - Identify dependencies on `fSourceDetectorDistance`
   - Calculate bounding box of detector assembly (housing, window, Ge crystal)

2. **Choose Integration Strategy** (15 min)
   - **Recommended**: Keep detector in world volume, position to be inside cavity geometrically
   - Alternative: Change mother volume to cavity (more invasive)

3. **Adjust Detector Position** (1 hour)
   - **Option A (Recommended - Minimal Change):**
     ```cpp
     // Keep existing detector construction code unchanged
     // Verify detector fits inside cavity at current position
     // Shield center positioned such that cavity encompasses detector
     ```
   - **Option B (Explicit Containment):**
     ```cpp
     // Change mother volume from fWorldLV to fShieldCavityLV for:
     // - Housing
     // - Vacuum inside housing
     // - All daughter volumes
     // Recalculate positions relative to cavity center
     ```

4. **Update Source Position** (30 min)
   - Modify PrimaryGeneratorAction to place source inside shield cavity
   - Calculate appropriate source-detector distance within cavity
   - Update macro file with new source position

5. **Rebuild and Run Initial Simulation** (30 min)
   ```bash
   cd build
   cmake --build . -j4
   ./HPGeSingle -m shielded_run.mac
   ```

6. **Integration Testing** (1.5 hours)
   - Verify no overlaps with detector inside shield
   - Run 1000 events with visualization
   - Check particle trajectories reach detector
   - Verify energy deposition still recorded in Ge crystal
   - Compare shielded vs unshielded count rates

**Deliverables:**
- ✓ Detector positioned correctly inside shield cavity
- ✓ Source configured for shield geometry
- ✓ Simulation runs without navigation errors
- ✓ Particles reach detector, energy deposited

**Testing Checklist:**
- [ ] Detector bounding box fully inside cavity (verified by calculation)
- [ ] No overlap warnings when running simulation
- [ ] Particles originate from source, travel through shield, some reach detector
- [ ] Energy deposition histograms show expected peaks (e.g., 662 keV for Cs-137)
- [ ] Background count rate reduced compared to unshielded geometry

**Risk Mitigation:**
- **Issue**: Detector doesn't fit in cavity → **Solution**: Increase cavity radius, recalculate
- **Issue**: Overlaps detected → **Solution**: Add clearance margins, check mother-daughter hierarchy
- **Issue**: No particles reach detector → **Solution**: Check source position, verify not placed outside shield accidentally

---

### Phase 3: Physics Validation and Optimization

**Goal:** Verify fluorescence physics, optimize performance, validate against expected behavior

**Duration:** 3-4 hours testing + analysis

**Tasks:**

1. **Fluorescence Verification** (1 hour)
   - Enable `/process/em/verbose 1` in macro
   - Confirm "Fluorescence is activated" in output
   - Add diagnostic code to SteppingAction to detect fluorescence X-rays
   - Run 10,000 events, analyze secondary particle spectrum
   - Verify Pb K-alpha (75 keV) and K-beta (85 keV) X-rays produced
   - Verify Cu K-alpha (8 keV) X-rays produced in copper liner

2. **Copper Liner Effectiveness Test** (1 hour)
   - **Test A**: Run simulation **with** copper liner (baseline)
   - **Test B**: Comment out copper liner placement, rebuild, run simulation
   - Compare spectra:
     - With Cu: Lead X-rays should be suppressed in detector
     - Without Cu: Lead X-rays should appear as prominent peaks
   - Quantify reduction factor

3. **Performance Benchmarking** (1 hour)
   - Run 10,000 events, measure time per event
   - Profile with `time` command
   - If >5 sec/event, implement region-specific production cuts:
     ```cpp
     // In PhysicsList::SetCuts()
     G4Region* shieldRegion = new G4Region("ShieldRegion");
     // Add lead and copper volumes
     // Set higher cuts (0.5-1.0 mm) for gamma, e-, e+
     ```
   - Re-benchmark, confirm speedup without significant accuracy loss

4. **Background Spectrum Analysis** (1.5 hours)
   - Modify SteppingAction/RunAction to record energy spectrum
   - Run 100,000 events (long run for statistics)
   - Export histogram to ROOT or CSV
   - Plot spectrum, identify peaks:
     - Primary source energy (e.g., 662 keV)
     - Lead fluorescence (75, 85 keV) - should be minimal with Cu liner
     - Compton continuum
   - Compare to literature spectra for HPGe detectors

5. **Geometric Attenuation Validation** (30 min)
   - Calculate expected transmission through 10 cm lead
   - Run simulation with known source strength
   - Count detected events
   - Verify transmission factor matches theory (within 20%)

6. **Documentation and Reporting** (30 min)
   - Create summary report with:
     - Performance metrics (events/sec)
     - Fluorescence validation results
     - Background spectrum plots
     - Comparison with/without copper liner
   - Update context.md

**Deliverables:**
- ✓ Fluorescence physics confirmed operational
- ✓ Copper liner reduces Pb X-rays by >50%
- ✓ Performance acceptable (<5 sec/event or optimized)
- ✓ Background spectrum matches expected features
- ✓ Attenuation coefficient validated

**Testing Checklist:**
- [ ] Fluorescence X-rays observed in stepping action diagnostic
- [ ] With Cu liner: Pb X-ray peaks <10% of primary peak intensity
- [ ] Without Cu liner: Pb X-ray peaks >30% of primary peak intensity
- [ ] Time per event <5 seconds (or <10 sec with optimization plan)
- [ ] Background spectrum shows expected features (primary, Compton, minimal fluorescence)
- [ ] Geometric attenuation within 20% of theory

**Success Metrics:**
- Fluorescence reduction: >50% (with vs without Cu)
- Background reduction: 10-20× for external sources
- Performance: <5 sec/event for 10 cm shield
- Attenuation: Within 20% of theoretical for 662 keV

---

## Risk Assessment

### Technical Risks

#### Risk 1: Geometry Overlap Errors

**Probability:** Medium (30%)
**Impact:** High (simulation produces incorrect results or crashes)

**Symptoms:**
- G4Exception warnings about overlapping volumes
- Navigation errors ("stuck track", "track leaving world")
- Inconsistent results between runs

**Root Causes:**
- Incorrect radius calculations (e.g., copper inner radius ≥ lead inner radius)
- Detector housing extends beyond shield cavity
- Insufficient clearance margins
- Mother-daughter volume containment violation

**Mitigation Strategies:**

*Prevention:*
- Use consistent naming for dimensions with clear comments
- Calculate all dependent dimensions explicitly (no hard-coded values)
- Add margin constants (e.g., `const G4double clearance = 5*mm;`)
- Enable `checkOverlaps=true` in all G4PVPlacement during development

*Detection:*
- Run `/geometry/test/run` with high resolution before any physics simulation
- Enable visual inspection with OpenGL wireframe mode
- Check for warning messages in console output

*Resolution:*
- Use systematic approach: verify each volume individually, then combinations
- Add debug output printing actual dimensions and positions
- If overlap found: increase clearances incrementally (1-2 mm at a time)
- Document all dimension constraints in code comments

**Contingency Plan:**
If overlaps cannot be resolved after 2 hours debugging:
- Temporarily increase all clearances by 10 mm
- Simplify geometry: remove copper liner temporarily, get lead working first
- Consult Geant4 forum with minimal reproducible example

---

#### Risk 2: Performance Degradation (Slow Simulation)

**Probability:** High (60%)
**Impact:** Medium (longer simulation times, but results still valid)

**Symptoms:**
- >10 seconds per event
- 100,000 event simulation takes >24 hours
- Simulation unusable for parameter scans

**Root Causes:**
- Thick lead requires tracking many scattering events
- Low production cuts generate excessive secondary particles
- Complex geometry slows navigation
- Fluorescence creates secondary X-rays, Auger electrons (more particles to track)

**Mitigation Strategies:**

*Prevention:*
- Use simple cylindrical geometry (G4Tubs), not boolean solids or tessellated meshes
- Start with default production cuts (0.1 mm globally)
- Limit initial test runs to 1,000 events for performance profiling

*Optimization Options (in order of implementation):*

1. **Region-specific production cuts** (implement first if needed):
   ```cpp
   // Looser cuts in shield (0.5-1.0 mm)
   // Tight cuts in detector (0.01 mm)
   // Reduces secondaries in thick lead by ~3×
   ```

2. **Disable Auger electrons** (minimal physics impact for thick shields):
   ```cpp
   /process/em/auger false
   /process/em/augerCascade false
   ```

3. **Increase global cuts for exploratory runs**:
   ```cpp
   SetCutValue(0.5*mm, "gamma");  // Temporarily
   ```

4. **Variance reduction** (only if above methods insufficient):
   - Implement importance biasing
   - Assign low importance to shield, high to detector
   - Requires careful validation

*Benchmarking:*
- Measure time per event after each optimization
- Verify physics not compromised: check attenuation coefficient unchanged (<5%)
- Document trade-offs in context.md

**Acceptance Criteria:**
- <5 sec/event: Production ready, no optimization needed
- 5-10 sec/event: Implement region cuts, aim for <5 sec
- >10 sec/event: Implement all optimizations above, consider variance reduction

**Contingency Plan:**
If performance remains poor after all optimizations:
- Reduce shield thickness to 5 cm for initial studies (trade accuracy for speed)
- Use coarser parameter scans (fewer simulation points)
- Consider GPU acceleration (future work, requires Geant4 GPU fork)

---

#### Risk 3: Fluorescence Physics Not Working Correctly

**Probability:** Low (10%)
**Impact:** High (shield effectiveness incorrect, defeats purpose of copper liner)

**Symptoms:**
- No lead X-ray peaks in spectrum (when expected)
- Copper liner has no effect on background
- Background spectrum identical with/without copper liner

**Root Causes:**
- Atomic deexcitation not enabled (unlikely with option4 physics list)
- G4EMLOW data library not installed or wrong version
- Energy range set incorrectly (fluorescence X-rays below tracking threshold)
- Bug in custom analysis code (not detecting fluorescence X-rays)

**Mitigation Strategies:**

*Prevention:*
- Verify G4EMLOW environment variable set: `echo $G4LEDATA`
- Check Geant4 data libraries installed: list `/path/to/geant4-data/`
- Use G4EmStandardPhysics_option4 (fluorescence enabled by default)
- Enable verbose output: `/process/em/verbose 1`

*Verification Tests:*
1. **Direct test**: Place point source in lead block, no detector, record secondaries
2. **Library check**: Run Geant4 example (e.g., TestEm1) that uses fluorescence
3. **Manual verification**: Add `G4cout` statements in SteppingAction for every secondary particle created

*Diagnostic Code:*
```cpp
// In SteppingAction::UserSteppingAction
const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
for (auto track : *secondaries) {
    if (track->GetCreatorProcess()->GetProcessName() == "phot") {
        G4cout << "Secondary from photoelectric: "
               << track->GetDefinition()->GetParticleName()
               << " E=" << track->GetKineticEnergy()/keV << " keV" << G4endl;
    }
}
```

**Expected Output:** Gamma secondaries with E ~ 75-85 keV from photoelectric events in lead

**Contingency Plan:**
If fluorescence confirmed not working:
- Reinstall G4EMLOW data library (download from Geant4 website)
- Switch to G4EmLivermorePhysics (alternative physics list with guaranteed fluorescence)
- Manually enable: `/process/em/fluo true` in macro file (belt-and-suspenders)
- If still failing: report bug to Geant4 forum with MWE (minimal working example)

---

#### Risk 4: Detector Positioning Error (Not Inside Shield)

**Probability:** Medium (25%)
**Impact:** High (detector not shielded, simulation meaningless)

**Symptoms:**
- Background count rate same as unshielded
- Visualization shows detector outside shield
- Overlap warnings between detector and shield

**Root Causes:**
- Incorrect calculation of shield center position
- Detector position not updated after adding shield
- Confusion between absolute vs relative coordinates
- Copy-paste error in position vector

**Mitigation Strategies:**

*Prevention:*
- Print all position calculations with clear labels:
   ```cpp
   G4cout << "Detector front face Z: " << detectorFrontZ/mm << " mm" << G4endl;
   G4cout << "Detector back face Z: " << detectorBackZ/mm << " mm" << G4endl;
   G4cout << "Shield cavity Z range: " << cavityMinZ/mm << " to " << cavityMaxZ/mm << " mm" << G4endl;
   ```
- Use named constants instead of magic numbers
- Verify containment with boolean check:
   ```cpp
   G4bool contained = (detectorFrontZ > cavityMinZ) && (detectorBackZ < cavityMaxZ);
   if (!contained) {
       G4Exception("DetectorConstruction", "GeometryError", FatalException,
                   "Detector not contained in shield cavity!");
   }
   ```

*Visual Verification:*
- Use OpenGL visualization with section plane:
   ```bash
   /vis/viewer/set/cutawayMode union
   /vis/viewer/addCutawayPlane 0 0 0 m 1 0 0
   ```
- Rotate view to see detector inside shield
- Check no gaps between cavity wall and detector

*Quantitative Test:*
- Run simulation with source inside shield, count detector hits
- Run simulation with source outside shield (external background simulation)
- If counts are similar → detector likely exposed (not shielded)

**Contingency Plan:**
If detector positioning proves difficult:
- Simplify: place both detector and shield centered at origin, separate later
- Use explicit mother-daughter relationship (place detector in cavity volume, not world)
- Draw schematic diagram on paper with all Z-coordinates marked
- Ask for peer review of geometry code

---

### Resource Requirements

#### Time Estimates

**Development Time:**

| Phase | Optimistic | Realistic | Pessimistic | Notes |
|-------|-----------|-----------|-------------|-------|
| Phase 1: Foundation | 2 hours | 3 hours | 5 hours | Includes geometry definition, compilation, basic tests |
| Phase 2: Integration | 2 hours | 3 hours | 6 hours | Detector positioning, source adjustment, debugging overlaps |
| Phase 3: Validation | 3 hours | 4 hours | 8 hours | Physics verification, performance optimization, analysis |
| **Total** | **7 hours** | **10 hours** | **19 hours** | Approximately 1-2.5 working days |

**Computation Time:**

| Task | Events | Time per Event | Total Time | Notes |
|------|--------|----------------|------------|-------|
| Initial geometry test | 100 | 2 sec | 3 min | Quick validation |
| Fluorescence verification | 10,000 | 3 sec | 8 hours | Collect statistics |
| Background spectrum | 100,000 | 3 sec | 3.5 days | Production run |
| Performance benchmark | 10,000 | 3-10 sec | 8-28 hours | Depends on optimization |

**Notes:**
- Time estimates assume familiarity with Geant4 and C++
- Debugging overlaps can be time-consuming (pessimistic scenario)
- Performance optimization may require iteration
- Computation time assumes single-threaded; can be parallelized with Geant4 multithreading

**Optimization Impact:**
- Multithreading (4 cores): 4× speedup → 100k events in ~21 hours
- Region cuts: ~2× speedup → 100k events in ~2 days
- Combined: ~8× speedup → 100k events in ~10 hours (feasible for production)

---

#### Team Skills Needed

**Essential Skills:**

1. **C++ Programming**
   - Intermediate level: classes, pointers, inheritance
   - Familiarity with C++17 standard (required by Geant4 11.x)
   - Ability to read compiler error messages and fix syntax errors

2. **Geant4 Framework Knowledge**
   - Understanding of detector construction hierarchy (solid → logical → physical)
   - Familiarity with G4VSolid classes (G4Tubs, G4Box)
   - Experience with G4Material and G4NistManager
   - Basic understanding of physics lists and processes

3. **Geometry and Coordinate Systems**
   - 3D spatial reasoning (understanding cylindrical coordinates)
   - Ability to calculate positions, rotations, transformations
   - Debugging overlapping volumes

4. **Build Systems**
   - CMake basics (configure, build, clean)
   - Command-line compilation
   - Managing build directories

**Helpful but Not Required:**

5. **Gamma-Ray Spectroscopy Physics**
   - Understanding of photoelectric effect, Compton scattering
   - Knowledge of X-ray fluorescence
   - Familiarity with HPGe detector operation (helpful for validation)

6. **Data Analysis**
   - ROOT or Python for histogram analysis
   - Plotting spectra, comparing datasets
   - Statistical analysis for validation tests

7. **Linux/Unix Command Line**
   - Navigating filesystems
   - Running executables with parameters
   - Using `time`, `top` for performance monitoring

**Learning Resources:**

- **Geant4 Official Documentation**: [https://geant4.web.cern.ch/](https://geant4.web.cern.ch/)
- **Geant4 User's Guide for Application Developers**: Chapters 4 (Detector Definition) and 5 (Physics Lists)
- **Geant4 Examples**: `examples/basic/B1` (simple geometry), `examples/extended/electromagnetic/TestEm1` (fluorescence)

**Recommended Team Structure:**

- **1 developer**: Implements geometry and physics (10-20 hours)
- **1 reviewer** (optional): Code review, geometry validation (2-3 hours)
- **1 analyst** (optional): Data analysis and validation (5-10 hours)

**For Solo Implementation:**
- Budget 10-20 hours total over 2-3 days
- Expect iteration on geometry debugging and performance optimization
- Use phased approach: get basic shield working, then optimize

---

#### Infrastructure Requirements

**Compute Hardware:**

**Minimum Specifications:**
- **CPU**: 2 cores, 2.0 GHz (simulation will be slow)
- **RAM**: 2 GB
- **Disk Space**: 500 MB (Geant4 executable + output)
- **OS**: Linux (Ubuntu 20.04+), macOS 10.14+, or Windows 10 with WSL2

**Recommended Specifications:**
- **CPU**: 4-8 cores, 3.0+ GHz (for multithreading)
- **RAM**: 8 GB (comfortable for large simulations)
- **Disk Space**: 5 GB (for multiple output files, ROOT histograms)
- **OS**: Linux (native Geant4 performance best)

**Performance Scaling:**
- **Single-threaded**: Baseline
- **4 threads**: ~3.5× speedup (Geant4 MT efficiency ~87%)
- **8 threads**: ~6× speedup (efficiency ~75%, diminishing returns)

**Software Dependencies:**

**Required:**
- Geant4 11.0 or higher (project uses 11.3.2) - ✓ Already installed
- CMake 3.16+ (project uses 3.28.3) - ✓ Already installed
- GCC 9+ or Clang 10+ (project uses GCC 13.3.0) - ✓ Already installed
- Geant4 data libraries:
  - **G4EMLOW 8.6+**: Critical for fluorescence physics
  - G4PhotonEvaporation 5.7
  - G4RadioactiveDecay 5.6
  - (Others installed with standard Geant4 distribution)

**Verification:**
```bash
echo $G4LEDATA
# Should output: /path/to/geant4-install/share/Geant4-11.3.2/data/G4EMLOW8.6
```

**Optional (for Visualization):**
- OpenGL libraries (for OGL visualization) - ✓ Likely installed
- X11 libraries (for GUI) - ✓ Likely installed
- VRML viewer (for VRML export): optional
- Qt5 (for Qt GUI): optional

**Optional (for Advanced Analysis):**
- ROOT 6.x (for ROOT output format): not currently used, but can be added
- Python 3.8+ with numpy, matplotlib (for analysis scripts)

**Data Storage:**

**Per Simulation Run:**
- Console output: ~1-10 MB (text log)
- ROOT histogram file: ~10-100 MB (depends on binning, number of events)
- Visualization files (VRML): ~50-200 MB (geometry export)

**Total Project Storage:**
- Source code: ~500 KB
- Build artifacts: ~20 MB
- Documentation: ~5 MB
- Simulation outputs: 100 MB - 1 GB (for multiple runs)

**Recommended:** 5 GB free disk space for comfortable operation

**Network Requirements:**
- Internet connection for downloading Geant4 data libraries (if not installed)
- G4EMLOW download size: ~60 MB
- No network required during simulation runs

**Current Project Status:**
According to context.md, infrastructure is already in place:
- ✓ Geant4 11.3.2 installed at `/home/nam/geant4-install/`
- ✓ CMake 3.28.3
- ✓ GCC 13.3.0
- ✓ OpenGL and X11 available
- ✓ Build system configured and tested

**Additional Requirements for Lead Shield Implementation:**
- None - existing infrastructure is sufficient
- Recommendation: verify G4EMLOW version ≥8.6 for best fluorescence accuracy

**Infrastructure Checklist:**
- [x] Geant4 11.0+ installed
- [x] CMake 3.16+ installed
- [x] C++ compiler (GCC 9+ or Clang 10+)
- [ ] Verify G4EMLOW data library version (run: `ls $G4LEDATA`)
- [x] OpenGL visualization available
- [ ] Sufficient disk space (5 GB free)
- [ ] Python + matplotlib (optional, for analysis)

**If G4EMLOW Missing or Outdated:**
```bash
# Download from Geant4 website
cd /path/to/geant4-install/share/Geant4-11.3.2/data/
wget https://geant4-data.web.cern.ch/geant4-data/datasets/G4EMLOW.8.6.tar.gz
tar -xzf G4EMLOW.8.6.tar.gz
export G4LEDATA=/path/to/geant4-install/share/Geant4-11.3.2/data/G4EMLOW8.6
# Add export to ~/.bashrc for persistence
```

---

## Trade-off Analysis

### Shield Thickness: 5 cm vs 10 cm vs 15 cm

| Criterion | 5 cm Lead | **10 cm Lead (Recommended)** | 15 cm Lead |
|-----------|-----------|------------------------------|------------|
| **Background Reduction** | ~5-7× | ~15× | ~30× |
| **Material Cost** | Low | Medium | High |
| **Physical Weight** | ~25 kg | ~100 kg | ~225 kg |
| **Simulation Speed** | Fast (1-2 sec/event) | Medium (3-5 sec/event) | Slow (5-10 sec/event) |
| **Lead Fluorescence** | Higher (less self-absorption) | Moderate | Lower (more self-absorption) |
| **Typical Application** | Teaching labs, high-activity sources | **Soil/environmental analysis** | Ultra-low background, cosmogenic isotopes |
| **Ease of Construction** | Easy (portable) | Moderate (bench-mounted) | Difficult (permanent installation) |

**Justification for 10 cm:**
- **Industry standard**: Most commercial HPGe soil analyzers use 10 cm
- **Cost-effective**: Provides 15× reduction at reasonable material cost
- **Simulation balance**: Thick enough to test physics, fast enough for parameter scans
- **Literature validation**: Most published studies use 10 cm, enabling comparison
- **Practical**: Weight (~100 kg) manageable for bench-mounted system

**When to Use Alternatives:**
- **5 cm**: Exploratory simulations, teaching, high-activity sources where background is not limiting factor
- **15 cm**: Ultra-low background applications (e.g., dark matter experiments, cosmogenic isotope studies), when simulation speed is not critical

**Recommendation:** **10 cm** for this project (soil analysis focus)

---

### Copper Liner: Yes vs No

| Criterion | **With 1 mm Cu Liner (Recommended)** | Without Liner (Bare Lead) |
|-----------|--------------------------------------|---------------------------|
| **Pb X-ray Suppression** | ~50% reduction | None (85 keV peaks prominent) |
| **Low-Energy Background** | Lower (Cu X-rays at 8 keV, below threshold) | Higher (Pb X-rays at 85 keV, above threshold) |
| **Geometry Complexity** | Slightly more complex (two volumes) | Simpler (one volume) |
| **Material Cost** | +5% (Cu is cheaper than Pb) | Baseline |
| **Physical Construction** | Requires two-layer assembly | Single-layer assembly |
| **Simulation Speed** | Same (Cu layer thin) | Same |
| **Typical Application** | **Gamma spectroscopy <200 keV** | High-energy only (>200 keV) |
| **Literature Support** | Standard for low-background work | Used for simplicity/cost savings |

**Justification for Cu Liner:**
- **Physics**: Suppresses 85 keV Pb fluorescence by ~50%, critical for soil analysis (40K at 1461 keV produces Compton continuum + Pb fluorescence)
- **Minimal cost**: Cu liner adds ~5% to material cost, negligible weight
- **Negligible complexity**: Only one additional G4Tubs volume in code
- **Best practice**: All modern commercial shields include Cu or Cd/Sn liner

**When to Omit Liner:**
- High-energy-only studies (>500 keV where Pb fluorescence not interfering)
- Teaching simulations focused on geometry, not detailed physics
- Rapid prototyping (simplify first, add liner later)

**Recommendation:** **Include Cu liner** for realistic soil analysis simulation

---

### Geometry Approach: Nested Cylinders vs Boolean Solids vs Separate Components

| Criterion | **Nested Cylinders (Recommended)** | Boolean Solids | Separate Components (Walls/Top/Bottom) |
|-----------|-----------------------------------|----------------|----------------------------------------|
| **Code Complexity** | Low (simple G4Tubs) | Medium (subtraction operations) | High (many volumes to track) |
| **Performance** | Fast (optimized for cylindrical navigation) | Slower (boolean evaluation overhead) | Medium |
| **Flexibility** | Limited (cylindrical symmetry only) | High (arbitrary cutouts, openings) | Very high (asymmetric designs) |
| **Ease of Debugging** | Easy (clear hierarchy) | Medium (subtraction positions tricky) | Hard (many overlap checks needed) |
| **Visualization** | Clean, clear structure | Can be visually messy | Requires careful color coding |
| **Typical Use Case** | **Standard shield (full enclosure)** | Shields with cable ports, collimators | Multi-piece shields, complex assemblies |
| **Lines of Code** | ~30 | ~60 | ~100+ |

**Justification for Nested Cylinders:**
- **Simplicity**: 3 volumes (lead, copper, cavity), straightforward hierarchy
- **Performance**: Geant4 cylindrical navigation is highly optimized
- **Sufficient**: Full enclosure is standard for soil analysis, no cutouts needed initially
- **Maintainable**: Easy for others to understand and modify

**When to Use Alternatives:**
- **Boolean solids**: If shield needs cable ports, beam entry windows, or collimator holes (Phase 2 enhancement)
- **Separate components**: If modeling real shield assembly process, or highly asymmetric design

**Recommendation:** **Nested cylinders** for Phase 1, consider boolean solids for Phase 2 if sample insertion opening is needed

---

### Mother Volume Strategy: Detector in World vs Detector in Cavity

| Criterion | **Detector in World (Recommended)** | Detector in Cavity Volume |
|-----------|-------------------------------------|---------------------------|
| **Code Changes** | Minimal (detector construction unchanged) | Moderate (change mother volume, recalculate positions) |
| **Position Calculation** | Absolute coordinates (simpler) | Relative to cavity center (more bookkeeping) |
| **Risk of Error** | Low (existing code proven to work) | Medium (recalculation introduces bugs) |
| **Conceptual Clarity** | Less clear (detector "happens to be" inside cavity) | More clear (explicit containment) |
| **Overlap Checking** | Must verify geometrically | Guaranteed if positions correct |
| **Flexibility** | Easy to move detector or shield independently | Detector coupled to shield position |

**Justification for Detector in World:**
- **Minimal risk**: Don't modify working detector geometry code
- **Faster implementation**: No need to recalculate all detector positions
- **Easy rollback**: Can remove shield by commenting out 3 lines, detector unaffected

**When to Use Cavity as Mother:**
- Detector and shield are conceptually inseparable (e.g., permanently enclosed detector)
- Need to move entire assembly (detector + shield) as unit
- Prefer explicit containment hierarchy for clarity

**Recommendation:** **Detector in World** for Phase 1 (minimize risk), optionally refactor to cavity in Phase 2 after validation

**Implementation:**
```cpp
// Phase 1: Recommended approach
new G4PVPlacement(0, detectorPos, detectorLV, "Detector", fWorldLV, false, 0, true);
new G4PVPlacement(0, shieldPos, shieldLV, "Shield", fWorldLV, false, 0, true);
// Verify: detectorPos is inside shield cavity geometrically

// Phase 2: Alternative (if needed)
new G4PVPlacement(0, relativeDetectorPos, detectorLV, "Detector", fShieldCavityLV, false, 0, true);
// relativeDetectorPos calculated relative to cavity center
```

---

### Production Cuts: Uniform vs Region-Specific

| Criterion | Uniform Cuts (0.1 mm everywhere) | **Region-Specific (Recommended)** |
|-----------|----------------------------------|-----------------------------------|
| **Simplicity** | Very simple (current implementation) | Moderate (requires region definition) |
| **Performance** | Slower (low cuts in thick lead) | Faster (~2-3× speedup) |
| **Accuracy in Detector** | High | High (detector region keeps tight cuts) |
| **Accuracy in Shield** | High (possibly overkill) | Good (looser cuts acceptable for bulk attenuation) |
| **Physics Validity** | Guaranteed accurate | Valid if shield cuts not too coarse |
| **Lines of Code** | ~10 | ~30 |

**Recommended Configuration:**
```cpp
// Global: 0.1 mm (default for world)
// Detector region: 0.01 mm (high precision for Ge crystal)
// Shield region: 0.5 mm (faster, <5% impact on attenuation)
```

**Trade-off Analysis:**
- **0.5 mm shield cuts**: 2-3× speedup, <2% change in background counts
- **1.0 mm shield cuts**: 3-5× speedup, <5% change in background counts (acceptable)
- **2.0 mm shield cuts**: 5-8× speedup, ~10% change (not recommended, too coarse)

**Validation Test:**
1. Run with uniform 0.1 mm cuts (baseline)
2. Run with region-specific (0.5 mm in shield)
3. Compare:
   - Detector count rate (should differ <5%)
   - Energy spectrum shape (should be nearly identical)
   - Simulation time (should decrease by 2-3×)

**Recommendation:** **Implement region-specific cuts** with 0.5 mm in shield if performance is issue (>5 sec/event), otherwise start with uniform cuts for simplicity

---

## Recommendations

### Primary Recommendations

1. **Implement 10 cm Lead Shield with 1 mm Copper Liner**
   - **Why**: Industry standard, excellent background reduction (15×), suppresses Pb fluorescence by 50%
   - **Trade-off**: Moderate simulation time (3-5 sec/event), manageable physical weight (~100 kg)
   - **Action**: Follow Phase 1 implementation roadmap

2. **Use Nested Cylindrical Geometry (G4Tubs)**
   - **Why**: Simple, fast, easy to debug, sufficient for full-enclosure shield
   - **Trade-off**: Cannot model complex cutouts (use boolean solids later if needed)
   - **Action**: Create lead → copper → cavity volume hierarchy in DefineVolumes()

3. **Place Detector in World Volume (Minimal Code Changes)**
   - **Why**: Minimal risk, fast implementation, existing detector geometry unchanged
   - **Trade-off**: Less explicit containment (verify geometrically)
   - **Action**: Position shield such that cavity encompasses detector at current position

4. **Verify Fluorescence Physics (G4EmStandardPhysics_option4)**
   - **Why**: Project already uses optimal physics list, but must verify it's working
   - **Trade-off**: None (already configured)
   - **Action**: Run `/process/em/printParameters` and check for "Fluorescence is activated"

5. **Enable Overlap Checking During Development**
   - **Why**: Geometry errors are common, early detection saves debugging time
   - **Trade-off**: Slightly slower initialization (<5 seconds)
   - **Action**: Set `checkOverlaps=true` in all G4PVPlacement constructors

### Secondary Recommendations (Implement if Needed)

6. **Use Region-Specific Production Cuts (If Performance Poor)**
   - **When**: If simulation >5 sec/event after basic implementation
   - **Action**: Implement 0.5 mm cuts in shield region, 0.01 mm in detector
   - **Expected Result**: 2-3× speedup, <5% accuracy impact

7. **Add Visualization Attributes for Shield Components**
   - **Why**: Essential for visual debugging, clear presentation
   - **Action**: Lead = dark gray semi-transparent, Copper = orange semi-transparent, Cavity = invisible

8. **Document All Dimensions in Code Comments**
   - **Why**: Geometry calculations are complex, comments aid future modifications
   - **Action**: Add comment block with ASCII diagram of shield cross-section

9. **Create Shielded vs Unshielded Comparison Macro**
   - **Why**: Quantifies shield effectiveness, validates implementation
   - **Action**: Run same source in both configurations, compare count rates

10. **Profile Performance Before Optimizing**
    - **Why**: Avoid premature optimization, focus on actual bottlenecks
    - **Action**: Run 1,000 events, measure time, only optimize if >5 sec/event

### Future Enhancements (Phase 2+)

11. **Add Top Opening with Lead Plug (Boolean Solids)**
    - **When**: After Phase 1 validated, if realistic sample insertion is needed
    - **Benefit**: More realistic experimental setup, can model plug removal
    - **Complexity**: Medium (requires G4SubtractionSolid)

12. **Implement Graded Shield (Pb → Sn → Cu Three-Layer)**
    - **When**: If ultra-low background required, Pb fluorescence still interfering
    - **Benefit**: Further suppresses fluorescence cascade
    - **Complexity**: Low (just add third G4Tubs volume for Sn)

13. **Add Shield Self-Activity (210Pb Background Source)**
    - **When**: For ultra-realistic background simulation
    - **Benefit**: Models natural lead radioactivity
    - **Complexity**: Medium (requires G4RadioactiveDecayPhysics configuration)

14. **Implement Multithreading for Performance**
    - **When**: If single-threaded too slow even after optimization
    - **Benefit**: 4× speedup on 4-core machine
    - **Complexity**: Low (Geant4 supports MT, mainly CMake changes)

---

## Resources

### Official Geant4 Documentation

1. **Geant4 User's Guide for Application Developers**
   - URL: [https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/](https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/)
   - Relevant Sections:
     - Chapter 4.1: Detector Definition and Response
     - Chapter 4.2: How to Define Geometry
     - Chapter 4.4: Material (G4NistManager, custom materials)
     - Chapter 5: Physics Lists (electromagnetic processes, fluorescence)

2. **Geant4 Material Database**
   - URL: [https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/material.html](https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/material.html)
   - Contents: Complete list of NIST materials (G4_Pb, G4_Cu), properties, composition

3. **Low Energy Electromagnetic Physics**
   - URL: [https://geant4.in2p3.fr/IMG/pdf_Lecture-LowEnergyEMPhysics.pdf](https://geant4.in2p3.fr/IMG/pdf_Lecture-LowEnergyEMPhysics.pdf)
   - Topics: Atomic deexcitation, fluorescence, Auger electrons, PIXE

4. **Geant4 Geometry Tutorial**
   - URL: [https://www.slac.stanford.edu/xorg/geant4/SLACTutorial14/Geometry1.pdf](https://www.slac.stanford.edu/xorg/geant4/SLACTutorial14/Geometry1.pdf)
   - Topics: G4VSolid types, boolean operations, placement, overlap checking

### Scientific Literature on HPGe Shielding

5. **"Design of a multi-layered shield for low-background HPGe spectrometry"**
   - Journal: Nuclear Instruments and Methods in Physics Research A
   - Year: 2025
   - DOI: [https://doi.org/10.1016/j.nima.2025.000609](https://www.sciencedirect.com/science/article/abs/pii/S0168900225000609)
   - Key Finding: 15 cm steel + 5 mm low-background Pb + 1 cm steel liner achieves same background as 10 cm Pb

6. **"Evaluation of lead shielding for a gamma-spectroscopy system"**
   - Journal: Nuclear Instruments and Methods in Physics Research A
   - Year: 2008
   - DOI: [https://doi.org/10.1016/j.nima.2008.002556](https://www.sciencedirect.com/science/article/abs/pii/S0168900208002556)
   - Key Finding: 10 cm lead provides ~15× background reduction

7. **"Estimation of background spectrum in a shielded HPGe detector using Monte Carlo simulations"**
   - Journal: Applied Radiation and Isotopes
   - Year: 2013
   - DOI: [https://doi.org/10.1016/j.apradiso.2013.004053](https://www.sciencedirect.com/science/article/abs/pii/S0969804313004053)
   - Key Finding: Monte Carlo validated for shield design, 12 cm optimal thickness

8. **"Precision measurement of radioactivity in gamma-rays spectrometry using two HPGe detectors"**
   - Journal: Results in Physics
   - Year: 2016
   - DOI: [https://doi.org/10.1016/j.rinp.2016.00450](https://www.sciencedirect.com/science/article/pii/S2215016116300450)
   - Application: Soil measurement with BEGe-6530 and GC0818-7600SL models

### Graded Shield Design

9. **"Graded shield" - University of Liverpool Radiometrics**
   - URL: [https://ns.ph.liv.ac.uk/~ajb/radiometrics/gamma_radiation/interactions_surroundings/graded_shield.html](https://ns.ph.liv.ac.uk/~ajb/radiometrics/gamma_radiation/interactions_surroundings/graded_shield.html)
   - Contents: Explanation of Pb → Cd/Sn → Cu graded shielding, fluorescence suppression mechanism

10. **"Homemade Graded-Z Shield for a Gamma-ray Spectrometer" - KandR Smith**
    - URL: [https://www.kandrsmith.org/RJS/Misc/GradedShield/gradedshield.html](https://www.kandrsmith.org/RJS/Misc/GradedShield/gradedshield.html)
    - Contents: Practical DIY guide, material selection, fluorescence X-ray energies, effectiveness measurements

11. **"Shield Lining for X-Ray Suppression" - Gamma Spectrometry Forum**
    - URL: [https://www.gammaspectacular.com/phpBB3/viewtopic.php?t=636](https://www.gammaspectacular.com/phpBB3/viewtopic.php?t=636)
    - Contents: Community discussion on Cd, Sn, Cu liner effectiveness, practical tips

### Geant4 Example Code

12. **Geant4 Basic Example B1**
    - Path: `$G4INSTALL/share/Geant4-11.3.2/examples/basic/B1/`
    - Relevance: Simple detector construction, G4Box and G4Tubs usage, good starting template

13. **Geant4 Extended Example: TestEm1**
    - Path: `$G4INSTALL/share/Geant4-11.3.2/examples/extended/electromagnetic/TestEm1/`
    - Relevance: Fluorescence verification, G4EmStandardPhysics_option4 usage, energy deposition

14. **GitHub: shield - Geant4 Shielding Simulation**
    - URL: [https://github.com/suerfu/shield](https://github.com/suerfu/shield)
    - Contents: Complete Geant4 program for studying shielding and veto performance, includes importance biasing

### NIST Data Libraries

15. **NIST Physical Reference Data**
    - URL: [https://physics.nist.gov/PhysRefData](https://physics.nist.gov/PhysRefData)
    - Contents: Atomic weights, isotopic compositions, X-ray mass attenuation coefficients
    - Specific:
      - Lead (Z=82): [https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Pb](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Pb)
      - Copper (Z=29): [https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Cu](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Cu)

16. **NIST XCOM: Photon Cross Sections Database**
    - URL: [https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html](https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html)
    - Use: Calculate attenuation coefficients for lead and copper at specific energies
    - Example: For 662 keV in Pb, μ/ρ = 0.1047 cm²/g → μ = 1.19 cm⁻¹

### Community Resources

17. **Geant4 Forum**
    - URL: [https://geant4-forum.web.cern.ch/](https://geant4-forum.web.cern.ch/)
    - Use: Search for shielding-related questions, post geometry debugging questions
    - Relevant Tags: `detector-construction`, `electromagnetic-physics`, `geometry`

18. **Gamma Spectrometry Forum**
    - URL: [https://www.gammaspectacular.com/phpBB3/](https://www.gammaspectacular.com/phpBB3/)
    - Use: Practical HPGe detector operation, shielding effectiveness discussions

### Commercial HPGe Shield Specifications

19. **ORTEC: High Performance Ultra Low Background Lead Shields (HPULBS)**
    - URL: [https://www.ortec-online.com/products/radiation-detectors/germanium-hpge-radiation-detectors/detector-accessories/hpulbs](https://www.ortec-online.com/products/radiation-detectors/germanium-hpge-radiation-detectors/detector-accessories/hpulbs)
    - Contents: Commercial specifications (10 cm lead, Cu liner, dimensions)

20. **Mirion: 777 Ultra Low-Background Shield**
    - URL: [https://www.mirion.com/products/technologies/spectroscopy-scientific-analysis/gamma-spectroscopy/detectors/hpge-shields-accessories/777-ultra-low-background-shield](https://www.mirion.com/products/technologies/spectroscopy-scientific-analysis/gamma-spectroscopy/detectors/hpge-shields-accessories/777-ultra-low-background-shield)
    - Contents: Advanced shield design (15 cm total, graded construction)

### Additional Tools

21. **Geant4 Data Libraries Download**
    - URL: [https://geant4-data.web.cern.ch/geant4-data/datasets/](https://geant4-data.web.cern.ch/geant4-data/datasets/)
    - Required: G4EMLOW 8.6+ for fluorescence data

22. **ROOT Data Analysis Framework**
    - URL: [https://root.cern.ch/](https://root.cern.ch/)
    - Use: Analyze output histograms, plot energy spectra (optional but helpful)

23. **Python + NumPy + Matplotlib**
    - Use: Alternative to ROOT for plotting spectra, analyzing simulation results
    - Tutorial: [https://matplotlib.org/stable/gallery/index.html](https://matplotlib.org/stable/gallery/index.html)

---

**End of Lead Shield Architecture Plan**

This plan provides a complete roadmap for implementing lead shielding in the HPGe detector simulation. All research findings, design decisions, implementation steps, and testing strategies are documented for use by the implementation team.

**Next Steps:**
1. Review this plan thoroughly
2. Begin Phase 1 implementation: Foundation (shield geometry)
3. Proceed to Phase 2: Detector integration
4. Complete Phase 3: Physics validation
5. Document results and update context.md

For questions or clarifications, consult the Geant4 forum or scientific literature references provided above.
