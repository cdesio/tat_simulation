#include "detector.hh"
#include <math.h>

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
}

MySensitiveDetector::~MySensitiveDetector()
{
}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    // access track of the particle which enters sensitive volume
    G4Track *track = aStep->GetTrack();
    // track->SetTrackStatus(fStopAndKill); // kill track when photon enters the detector and do not propagate
    // G4double stepl = track->GetStepLength();

    // G4cout << "step length: " << stepl << G4endl;

    const G4ParticleDefinition *particle = track->GetParticleDefinition();

    G4String name = particle->GetParticleName();
    G4int pid = particle->GetPDGEncoding();

    G4int trackID = track->GetTrackID();
    G4int parentID = track->GetParentID();
    // G4cout << "Pid: " << pid << ", name: " << name << ", track: " << track->GetTrackID() << ", parent: " << par_id << G4endl;
    //   G4int Z = particle->GetAtomicNumber();
    //   G4int A = particle->GetAtomicMass();
    //   G4double charge = particle->GetPDGCharge();
    G4double energy = track->GetKineticEnergy();
    G4double tot_energy = track->GetTotalEnergy();
    // G4cout << "K Energy: " << tot_energy << G4endl;
    // if (tot_energy != energy)
    // {
    //     G4cout << "tot Energy: " << tot_energy << G4endl;
    // };
    //  start and end of track
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    // access position of the photon when it enters
    G4ThreeVector posParticle = preStepPoint->GetPosition();
    // G4cout << "Photon position: " << posPhoton << G4endl;

    G4ThreeVector momParticle = preStepPoint->GetMomentum();

    // get
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

    G4int copyNo = touchable->GetCopyNumber(); // get id of the sensitive voxel
    G4String volName = touchable->GetVolume()->GetName();

    //  access position of sensitive detector

    G4VPhysicalVolume *physVol = touchable->GetVolume();

    G4ThreeVector posDetector = physVol->GetTranslation();
    G4RotationMatrix rotDetector = physVol->GetObjectRotationValue();
    // G4cout << "inverted: " << rotDetector.invert() << G4endl;
    // G4double theta = acos(rotDetector[0][0]);

    // G4Rotate3D rotZ(theta, G4ThreeVector(0, 0, 1));

    // G4RotationMatrix *rot = physVol->GetRotation();

    G4ThreeVector localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(posParticle);

    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    G4AnalysisManager *man = G4AnalysisManager::Instance();
    // fill ntuple, first index is ntuple id
    // mc
    man->FillNtupleIColumn(0, 0, evt);
    man->FillNtupleDColumn(0, 1, posParticle[0]);
    man->FillNtupleDColumn(0, 2, posParticle[1]);
    man->FillNtupleDColumn(0, 3, posParticle[2]);
    man->FillNtupleIColumn(0, 4, pid);
    man->FillNtupleSColumn(0, 5, name);
    man->FillNtupleDColumn(0, 6, energy);
    man->FillNtupleIColumn(0, 7, trackID);
    man->FillNtupleIColumn(0, 8, parentID);
    // man->FillNtupleDColumn(0, 4, wlen);
    man->AddNtupleRow(0);
    // hit
    man->FillNtupleIColumn(1, 0, evt);
    man->FillNtupleDColumn(1, 1, posDetector[0]);
    man->FillNtupleDColumn(1, 2, posDetector[1]);
    man->FillNtupleDColumn(1, 3, posDetector[2]);
    man->FillNtupleIColumn(1, 4, copyNo);
    man->AddNtupleRow(1);

    return true;
}
