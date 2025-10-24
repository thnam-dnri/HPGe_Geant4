// ==============================================================================
// GammaCascade.cc - Phase 1 Main Application
// Simple Geant4 gamma cascade simulation using RAINIER data
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
#include "G4SystemOfUnits.hh"

int main(int argc, char** argv)
{
    // Print welcome message
    G4cout << "\n========================================" << G4endl;
    G4cout << "  Gamma Cascade Simulation - Phase 1" << G4endl;
    G4cout << "  RAINIER + Geant4 Integration" << G4endl;
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

    // Parse arguments: --rainier <file>, --src-gap <value[mm|cm|m]>, <macro.mac>
    std::string rainierFile;
    G4String macroFile;
    G4double srcGap = 0.0; // default 0 distance (touching surfaces)

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
        if (arg == "--rainier" || arg == "-r") {
            if (i+1 < argc) rainierFile = argv[++i];
        } else if (arg.rfind("--rainier=",0)==0) {
            rainierFile = arg.substr(10);
        } else if (arg == "--src-gap" || arg == "-g") {
            if (i+1 < argc) srcGap = parseLength(argv[++i]);
        } else if (arg.rfind("--src-gap=",0)==0) {
            srcGap = parseLength(arg.substr(10));
        } else if (arg.size() > 4 && arg.substr(arg.size()-4) == ".mac") {
            macroFile = arg.c_str();
        } else if (rainierFile.empty()) {
            // Backward-compat: first bare arg as RAINIER file if not a macro
            rainierFile = arg;
        }
    }

    if (!rainierFile.empty()) {
        G4cout << "Using RAINIER input file: " << rainierFile << G4endl;
    } else {
        G4cout << "No RAINIER file specified. Using test Co-60 cascade." << G4endl;
    }
    
    PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction(rainierFile);
    primaryGenerator->SetDetectorConstruction(detConstruction);
    primaryGenerator->SetSourceSurfaceGap(srcGap);
    runManager->SetUserAction(primaryGenerator);

    // Set user action classes
    RunAction* runAction = new RunAction();
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
        // batch mode
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);
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
