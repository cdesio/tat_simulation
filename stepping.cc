#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
    fEventAction = eventAction; // gain access to the object we have created
    G4String filename("PSfile.bin");
    PSfile.open(filename, std::ios::out | std::ios::binary);
}

MySteppingAction::~MySteppingAction()
{
}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
    G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4LogicalVolume *fScoringDetector = detectorConstruction->GetScoringDetector();
    // G4LogicalVolume *fScoringTarget = detectorConstruction->GetScoringTarget();

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    G4double edepStep = step->GetTotalEnergyDeposit();
    if (edepStep <= 0)
        return;
    G4ThreeVector step_pos = step->GetPreStepPoint()->GetPosition();

    G4int trackID = step->GetTrack()->GetTrackID();
    G4double time = step->GetPreStepPoint()->GetGlobalTime();

    const G4VProcess *process = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4String processname = process->GetProcessName();

    const G4VTouchable *touchable = step->GetPreStepPoint()->GetTouchable();
    G4VPhysicalVolume *physVol = touchable->GetVolume();

    G4int copyno = touchable->GetCopyNumber();
    G4String particle_name = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4double energy = step->GetPreStepPoint()->GetKineticEnergy();

    if (volume == fScoringDetector)
    {

        // G4cout << "STEPPING - copyNo: " << copyno << G4endl;

        man->FillNtupleDColumn(3, 0, edepStep);
        man->FillNtupleIColumn(3, 1, copyno);
        man->FillNtupleDColumn(3, 2, step_pos[0]);
        man->FillNtupleDColumn(3, 3, step_pos[1]);
        man->FillNtupleDColumn(3, 4, step_pos[2]);
        man->FillNtupleIColumn(3, 5, trackID);
        man->FillNtupleIColumn(3, 6, time);
        man->FillNtupleSColumn(3, 7, processname);
        man->FillNtupleSColumn(3, 8, particle_name);
        man->FillNtupleDColumn(3, 9, energy);
        man->AddNtupleRow(3);
    }

    fEventAction->AddEdep(edepStep);

    // saving to PSfile

    // alphas that lose energy inside the logicDetector
    if (((particle_name == "alpha") || (particle_name == "e-") || (particle_name == "gamma")) && (volume == fScoringDetector) && (edepStep > 0))
    {
        G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
        G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentumDirection();

        G4int particle_ID;

        if (particle_name == "alpha")
            particle_ID = 0;
        if (particle_name == "gamma")
            particle_ID = 1;
        if (particle_name == "e-")
            particle_ID = 2;

        float output[9];
        output[0] = step_pos[0] / um;
        output[1] = step_pos[1] / um;
        output[2] = step_pos[2] / um;
        output[3] = momentum[0];
        output[4] = momentum[1];
        output[5] = momentum[2];
        output[6] = energy / MeV;
        output[7] = eventID;
        output[8] = particle_ID;
        // G4cout << "evt ID: " << eventID << " saving " << particle_name << " as " << particle_ID << G4endl;
        PSfile.write((char *)&output, sizeof(output));
    }
    // alphas entering the logicDetector from outside
    // if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "alpha") && (volume->GetName() == "logicWorld"))
    // {
    //     if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "logicDetector")
    //     {
    //         G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    //         G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentumDirection();
    //         G4double energy = step->GetPreStepPoint()->GetKineticEnergy();

    //         G4int particle_ID;

    // if (particle_name == "alpha")
    //     particle_ID = 0;
    // if (particle_name == "gamma")
    //     particle_ID = 1;
    // if (particle_name == "e-")
    //     particle_ID = 2;

    // float output[9];
    // output[0] = step_pos[0] / um;
    // output[1] = step_pos[1] / um;
    // output[2] = step_pos[2] / um;
    // output[3] = momentum[0];
    // output[4] = momentum[1];
    // output[5] = momentum[2];
    // output[6] = energy / MeV;
    // output[7] = eventID;
    // output[8] = particle_ID;
    //         G4cout << "event id: " << eventID << G4endl;
    //         PSfile.write((char *)&output, sizeof(output));
    //     }
    // }
}