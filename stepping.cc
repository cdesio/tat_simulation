#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
    fEventAction = eventAction; // gain access to the object we have created
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
    // G4cout << time * ns << G4endl;
    //  G4int parentID = step->GetTrack()->GetParentID();
    //  G4int stepNo = step->GetTrack()->GetCurrentStepNumber();
    const G4VProcess *process = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4String processname = process->GetProcessName();

    const G4VTouchable *touchable = step->GetPreStepPoint()->GetTouchable();
    G4VPhysicalVolume *physVol = touchable->GetVolume();

    G4int copyno = touchable->GetCopyNumber();

    if (volume == fScoringDetector)
    {

        // G4cout << "STEPPING - copyNo: " << copyno << G4endl;

        man->FillNtupleDColumn(3, 0, edepStep);
        man->FillNtupleIColumn(3, 1, copyno);
        man->FillNtupleDColumn(3, 2, step_pos[0]);
        man->FillNtupleDColumn(3, 3, step_pos[1]);
        man->FillNtupleDColumn(3, 4, step_pos[2]);
        // man->FillNtupleDColumn(3, 5, posDetector[0]);
        // man->FillNtupleDColumn(3, 6, posDetector[1]);
        // man->FillNtupleDColumn(3, 7, posDetector[2]);

        man->FillNtupleIColumn(3, 5, trackID);
        man->FillNtupleIColumn(3, 6, time);
        man->FillNtupleSColumn(3, 7, processname);
        man->AddNtupleRow(3);
    }

    fEventAction->AddEdep(edepStep);
}