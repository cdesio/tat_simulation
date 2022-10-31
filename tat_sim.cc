#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4ScoringManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "QGSP_BERT.hh" //physicslist for hadronic processes

#include "construction.hh"
#include "physics.hh"
#include "action.hh"

int main(int argc, char **argv)
{
    G4RunManager *runManager = new G4RunManager();
    G4ScoringManager::GetScoringManager();
    // initialize detector
    runManager->SetUserInitialization(new MyDetectorConstruction());

    // initialize physicslist
    runManager->SetUserInitialization(new MyPhysicsList());

    // initialize action
    runManager->SetUserInitialization(new MyActionInitialization());

    // for atmosphere
    // G4VModularPhysicsList *physics = new QGSP_BERT();
    // physics->RegisterPhysics(new G4DecayPhysics());
    // runManager->SetUserInitialization(physics);

    runManager->Initialize();
    // set ui to 0, only to create ui if in interactive mode (no macro)
    G4UIExecutive *ui = 0;

    if (argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
    }
    G4VisManager *visManager = new G4VisExecutive();
    visManager->Initialize();

    // Get ui manager pointer
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (ui)
    {
        // use ui manager to set commands
        // UImanager->ApplyCommand("/vis/open OGL");                         // open visualizer
        // UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1"); // standard view point
        // UImanager->ApplyCommand("/vis/drawVolume");                       // draw volume
        // UImanager->ApplyCommand("/vis/viewer/set/autorefresh true");      // update every event
        // UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");    // draw trajectory
        // UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
    }
    else // if more than 1 argument (macro mode)
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
    return 0;
}