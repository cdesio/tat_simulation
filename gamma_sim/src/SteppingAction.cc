//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
#include "SteppingAction.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <map>
#include "globals.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "CommandLineParser.hh"
#include "EventAction.hh"
using namespace G4DNAPARSER;

// using MapOfDelayedLists =
// std::map<double, std::map<int, G4TrackList*> >;
// using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(/*DetectorConstruction* fpDet*/)
    : G4UserSteppingAction(), fpEventAction(0)
// , fpDetector(fpDet)
{
  fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);

  G4String fileName{"PSfile.bin"};
  if ((command = parser->GetCommandIfActive("-out")) != 0)

    if (command->GetOption().empty() == false)
    {
      fileName = command->GetOption() + ".bin";
    }

  PSfile.open(fileName, std::ios::out | std::ios::binary);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{
  PSfile.close();
}

// void SteppingAction::Initialize()
// {
//   fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SteppingAction::UserSteppingAction(const G4Step *step)
{

  G4double dE = step->GetTotalEnergyDeposit();
  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  
  const G4String &volumeName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  //G4cout << "particle name: " << particleName<< " in " << volumeName  << ", dE: " << dE << G4endl;
    if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "gamma") && (volumeName == "chromatinSegment"))
  {
    G4DNAPARSER::CommandLineParser *parser = G4DNAPARSER::CommandLineParser::GetParser();
    // G4DNAPARSER::Command *command(0);
    if (parser->GetCommandIfActive("-out") == 0)
      return;
    
    fpEventAction->AddPathLength(step->GetStepLength());
    //G4cout << "pl: " << step->GetStepLength()/nm << G4endl;
  }

  if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() != "e-") && (volumeName == "chromatinSegment") /*&& (dE > 0*/)
  {
    G4DNAPARSER::CommandLineParser *parser = G4DNAPARSER::CommandLineParser::GetParser();
    G4DNAPARSER::Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
      return;

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    // G4ThreeVector prePoint = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
    // G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);
    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4int particleID = particleMap[particleName];
    G4int stepID = step->GetTrack()->GetCurrentStepNumber();
    auto particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    G4double steplength = step->GetStepLength();
    //G4cout << "particle: " << particleName << " in " << volumeName << G4endl;

    analysisManager->FillNtupleIColumn(0, 0, eventID);
    analysisManager->FillNtupleDColumn(0, 1, dE / eV);
    analysisManager->FillNtupleDColumn(0, 2, postPoint.x() / nanometer);
    analysisManager->FillNtupleDColumn(0, 3, postPoint.y() / nanometer);
    analysisManager->FillNtupleDColumn(0, 4, postPoint.z() / nanometer);
    analysisManager->FillNtupleIColumn(0, 5, particleID);
    analysisManager->FillNtupleSColumn(0, 6, particleName);
    analysisManager->FillNtupleIColumn(0, 7, stepID);
    analysisManager->FillNtupleDColumn(0, 8, particleEnergy/MeV);
    analysisManager->FillNtupleDColumn(0, 9, steplength);
    analysisManager->AddNtupleRow(0);

    fpEventAction->AddEdep(dE);
  }

  if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-") && (volumeName == "waterBox")) // e- from outside the tracking volume
  {
    if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "TrackingVol") // step ends on boundary to tracking volume
    {
    
    G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector worldMomentum = step->GetPostStepPoint()->GetMomentumDirection();

    auto particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    G4int copyNo = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
    G4double time = step->GetPreStepPoint()->GetGlobalTime();
    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4int particleID = particleMap[particleName];
    G4cout << "PS1 particleID: " << particleID << " in " << volumeName << G4endl;
    

    float output[12];
    output[0] = worldPos.x() / mm;
    output[1] = worldPos.y() / mm;
    output[2] = worldPos.z() / mm;
    output[3] = worldMomentum.x();
    output[4] = worldMomentum.y();
    output[5] = worldMomentum.z();
    output[6] = particleEnergy / MeV;
    output[7] = eventID;
    output[8] = particleID;
    output[9] = copyNo;
    output[10] = time;
    output[11] = 2; 

    PSfile.write((char *)&output, sizeof(output));
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

    fpEventAction->AddSecondary();
    }
  }

  if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-") && ((volumeName == "TrackingVol") || (volumeName == "chromatinSegment"))) // e- from interactions inside the tracking volume
  {

    G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector worldMomentum = step->GetPostStepPoint()->GetMomentumDirection();

    auto particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    G4int copyNo = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
    G4double time = step->GetPreStepPoint()->GetGlobalTime();
    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4int particleID = particleMap[particleName];
    G4cout << "PS2: particleID: " << particleID << " in " << volumeName << G4endl;
    
    float output[12];
    output[0] = worldPos.x() / mm;
    output[1] = worldPos.y() / mm;
    output[2] = worldPos.z() / mm;
    output[3] = worldMomentum.x();
    output[4] = worldMomentum.y();
    output[5] = worldMomentum.z();
    output[6] = particleEnergy / MeV;
    output[7] = eventID;
    output[8] = particleID;
    output[9] = copyNo;
    output[10] = time;
    output[11] = 2; 

    PSfile.write((char *)&output, sizeof(output));
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

    fpEventAction->AddSecondary();
  }
}
