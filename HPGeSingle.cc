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

    // Primary generator action
    std::string rainierFile = "";
    if (argc > 1) {
        rainierFile = std::string(argv[1]);
        G4cout << "Using RAINIER input file: " << rainierFile << G4endl;
    } else {
        G4cout << "No RAINIER file specified. Using test Co-60 cascade." << G4endl;
    }
    
    PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction(rainierFile);
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
    if (argc != 1) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[argc-1];
    if (fileName.find(".mac") != G4String::npos) {
            UImanager->ApplyCommand(command + fileName);
        } else {
            // Run with default parameters
            UImanager->ApplyCommand("/run/initialize");
            UImanager->ApplyCommand("/run/beamOn 10000");
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
