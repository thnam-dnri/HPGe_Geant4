// ==============================================================================
// GammaCascade.cc - Phase 1 Main Application
// Isotope-driven HPGe gamma cascade simulation
// ==============================================================================

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#include <string>
#include <sstream>
#include <cctype>
#include <algorithm>
#include "G4SystemOfUnits.hh"

int main(int argc, char** argv)
{
    // Print welcome message
    G4cout << "\n========================================" << G4endl;
    G4cout << "  Gamma Cascade Simulation - Phase 1" << G4endl;
    G4cout << "  Isotope-driven gamma singles generator" << G4endl;
    G4cout << "========================================\n" << G4endl;

    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager();

    // Set mandatory initialization classes
    
    // Detector construction
    DetectorConstruction* detConstruction = new DetectorConstruction();
    runManager->SetUserInitialization(detConstruction);

    // Physics list
    PhysicsList* physicsList = new PhysicsList();
    runManager->SetUserInitialization(physicsList);

    // Parse arguments: --src-gap <value[mm|cm|m]>, --isotope <Symbol>, --det-id <ID>, <macro.mac>
    std::string isotopeSymbol = "Co60";
    bool isotopeExplicit = false;
    std::string detID = "HPGe1";
    G4String macroFile;
    G4double srcGap = 0.0; // default 0 distance (touching surfaces)
    bool forceInteractive = false;

    auto parseLength = [](const std::string& s)->G4double {
        // Accept forms: 10mm, 1cm, 0.5m, or plain number = mm
        std::string num, unit;
        for (char c : s) {
            if (std::isdigit(c) || c=='.' || c=='e' || c=='E' || c=='+' || c=='-') num.push_back(c);
            else unit.push_back(c);
        }
        double val = 0.0;
        std::istringstream(num) >> val;
        if (unit == "mm" || unit.empty()) return val * mm;
        if (unit == "cm") return val * cm;
        if (unit == "m")  return val * m;
        // Fallback: treat as mm
        return val * mm;
    };

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--isotope" || arg == "-isotope") {
            if (i+1 < argc) {
                isotopeSymbol = argv[++i];
                isotopeExplicit = true;
            }
        } else if (arg.rfind("--isotope=",0)==0 || arg.rfind("-isotope=",0)==0) {
            size_t eq = arg.find('=');
            isotopeSymbol = arg.substr(eq+1);
            isotopeExplicit = true;
        } else if (arg == "--src-gap" || arg == "-g") {
            if (i+1 < argc) srcGap = parseLength(argv[++i]);
        } else if (arg.rfind("--src-gap=",0)==0) {
            srcGap = parseLength(arg.substr(10));
        } else if (arg == "--det-id") {
            if (i+1 < argc) detID = argv[++i];
        } else if (arg.rfind("--det-id=",0)==0) {
            detID = arg.substr(10);
        } else if (arg.size() > 4 && arg.substr(arg.size()-4) == ".mac") {
            macroFile = arg.c_str();
        } else if (arg == "--interactive" || arg == "--ui") {
            forceInteractive = true;
        } else if (!arg.empty() && arg[0] != '-') {
            isotopeSymbol = arg;
            isotopeExplicit = true;
        }
    }

    G4cout << "Using isotope data: " << isotopeSymbol << " (from ./isotope_data)" << G4endl;
    if (!isotopeExplicit) {
        G4cout << "  (default selection; override with --isotope <Symbol>)" << G4endl;
    }

    PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction();
    primaryGenerator->SetIsotopeSymbol(isotopeSymbol);
    primaryGenerator->SetDetectorConstruction(detConstruction);
    primaryGenerator->SetSourceSurfaceGap(srcGap);
    runManager->SetUserAction(primaryGenerator);

    // Set user action classes
    RunAction* runAction = new RunAction();
    runAction->SetDetID(detID.c_str());
    runManager->SetUserAction(runAction);

    EventAction* eventAction = new EventAction(runAction);
    runManager->SetUserAction(eventAction);

    SteppingAction* steppingAction = new SteppingAction(eventAction);
    runManager->SetUserAction(steppingAction);

    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    if (!macroFile.empty()) {
        // batch mode (optionally keep UI session for visualization macros)
        std::string macroName = macroFile;
        size_t slash = macroName.find_last_of("/\\");
        if (slash != std::string::npos) {
            macroName = macroName.substr(slash + 1);
        }
        std::string macroLower = macroName;
        std::transform(macroLower.begin(), macroLower.end(), macroLower.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        bool macroRequestsVisSession = (macroLower == "init_vis.mac");
        bool startUISession = forceInteractive || macroRequestsVisSession;

        G4UIExecutive* ui = nullptr;
        if (startUISession) {
            ui = new G4UIExecutive(argc, argv);
        }

        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);

        if (ui && ui->IsGUI() && macroRequestsVisSession) {
            UImanager->ApplyCommand("/control/execute gui.mac");
        }

        if (ui) {
            ui->SessionStart();
            delete ui;
            ui = nullptr;
        }
    } else {
        // interactive mode
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        
        // Execute initialization macro
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        if (ui->IsGUI()) {
            UImanager->ApplyCommand("/control/execute gui.mac");
        }
        
        // Start session
        ui->SessionStart();
        delete ui;
    }

    // Clean up
    delete visManager;
    delete runManager;

    G4cout << "\nSimulation completed successfully!" << G4endl;
    return 0;
}
