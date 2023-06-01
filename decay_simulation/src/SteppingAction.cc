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
#include "G4Ions.hh"

using namespace G4DNAPARSER;

// using MapOfDelayedLists =
// std::map<double, std::map<int, G4TrackList*> >;
// using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(DetectorConstruction *pDetector)
    : G4UserSteppingAction(), fpEventAction(0), fpDetector(pDetector)
{
  fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  Command *commandR(0);
  G4String fileName{"PSfile.bin"};
  if ((command = parser->GetCommandIfActive("-out")) != 0)

    if (command->GetOption().empty() == false)
    {

      if ((commandR = parser->GetCommandIfActive("-outer")) != 0)
      {
        G4int boxes_per_r = pDetector->get_ndiv_Z() * pDetector->get_ndiv_theta();
        G4cout << "n_div_Z: " << pDetector->get_ndiv_Z() << "n_div_theta: " << pDetector->get_ndiv_theta() << G4endl;
        fileName = command->GetOption() + "_r" + commandR->GetOption() + ".bin";
      }
      else
      {
        fileName = command->GetOption() + ".bin";
      }
    }

  PSfile.open(fileName, std::ios::out | std::ios::binary);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{
  PSfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SteppingAction::UserSteppingAction(const G4Step *step)
{

  G4double dE = step->GetTotalEnergyDeposit();

  const G4String &volumeNamePre = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  const G4String &particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4int pdg_enc = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  G4int particleID = -1;

  if (particleName == "alpha" || particleName == "helium" || particleName == "alpha+")
  {
    particleID = 1;
  }
  else if (particleName == "gamma")
  {
    particleID = 2;
  }
  else if (particleName == "e-")
  {
    particleID = 3;
  }
  else if (particleName == "nu_e")
  {
    particleID = 4;
  }
  else if (particleName == "At211")
  {
    particleID = 5;
  }
  else if (G4StrUtil::contains(particleName, "Po211"))
  {
    particleID = 6;
  }
  else if (G4StrUtil::contains(particleName, "Bi207"))
  {
    particleID = 7;
  }
  else if (G4StrUtil::contains(particleName, "Pb207")) //|| (particleName == "Pb207[569.698]") || (particleName == "Pb207[1633.356]") || (particleName == "Pb207[2339.921]"))
  {
    particleID = 8;
  }
  else if (particleName == "proton")
  {
    particleID = 9;
  }
  else if (particleName == "e+")
  {
    particleID = 10;
  }
  else
  {
    G4cout << particleName << " outside  not saved" << G4endl;
    return;
  }
  // G4cout << "DEBUG: " << particleName << ", ID: " << particleID << " vol: " << volumeNamePre << G4endl;
  
  //tracking all particles
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4ThreeVector prePoint = step->GetPreStepPoint()->GetPosition();
  G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *commandR(0);
  G4int copyNo = step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber();
  if ((commandR = parser->GetCommandIfActive("-outer")) != 0)
  {
    /// G4int boxes_per_r = pDetector->get_ndiv_Z() * pDetector->get_ndiv_theta();
    // G4cout << "copyNo: " << copyNo << ", thresh:" << std::stoi(commandR->GetOption()) << ", boxes_per_R: "<<  boxes_per_r << G4endl;
    boxes_per_r = 4000;
    if (copyNo < std::stoi(commandR->GetOption()) * boxes_per_r)
    {
      return;
    }
  }
  G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);
  G4int stepID = step->GetTrack()->GetCurrentStepNumber();
  auto particleEnergy = step->GetPostStepPoint()->GetKineticEnergy();
  G4double steplength = step->GetStepLength();

  const G4VProcess *defprocess = step->GetPostStepPoint()->GetProcessDefinedStep();
  const G4VProcess *creatprocess = step->GetTrack()->GetCreatorProcess();
  
  // if (defprocess)
  // {
  //   G4String defprocessname = defprocess->GetProcessName();

  //   if (((G4StrUtil::contains(defprocessname, "Decay")) || (defprocessname == "RadioactiveDecay")) || (defprocessname == "ionIoni"))
  //   {
  //     G4cout << "Event ID: " << eventID << ", particle: " << particleName << ", step: " << stepID << ", DefProcName: " << defprocessname << G4endl;
  //   }
  // }
  G4int processID = -1;
  if (creatprocess)
  {
    G4String creatprocessname = creatprocess->GetProcessName();
    
    if (creatprocessname == "RadioactiveDecay" || creatprocessname == "Decay" )
    {
      processID = 1;
    }
    else if (creatprocessname == "ionIoni")
    {
      processID = 2;
    }
    else if (creatprocessname == "msc")
    {
      processID = 3;
    }
    else if (creatprocessname == "eIoni")
    {
      processID = 4;
    }
    else if (creatprocessname == "phot")
    {
      processID = 5;
    }
    else if (creatprocessname == "eBrem")
    {
      processID = 6;
    }
    else if (creatprocessname == "compt")
    {
      processID = 7;
    }
    else{G4cout  << creatprocessname << " not saved. " << G4endl;}
    
    }

  if (volumeNamePre == "TrackingVol")
  {
    G4DNAPARSER::CommandLineParser *parser = G4DNAPARSER::CommandLineParser::GetParser();
    G4DNAPARSER::Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
      return;
    fpEventAction->AddPathLength(step->GetStepLength());
  }

  if ((volumeNamePre == "TrackingVol") && (dE > 0))
  {
    G4DNAPARSER::CommandLineParser *parser = G4DNAPARSER::CommandLineParser::GetParser();
    G4DNAPARSER::Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
      return;

    //G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    //G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    // G4ThreeVector prePoint = step->GetPreStepPoint()->GetPosition();
    //G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
    //G4int copyNo = step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber();
    // G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);
    // DEBUG
    // G4cout << "DEBUG: " << particleName << "(ID: " << particleID << ") losing E in chromatinSegment : " << eventID << G4endl;
    
    analysisManager->FillNtupleIColumn(0, 0, eventID);
    analysisManager->FillNtupleDColumn(0, 1, dE / eV);
    analysisManager->FillNtupleDColumn(0, 2, postPoint.x() / nanometer);
    analysisManager->FillNtupleDColumn(0, 3, postPoint.y() / nanometer);
    analysisManager->FillNtupleDColumn(0, 4, postPoint.z() / nanometer);
    analysisManager->FillNtupleIColumn(0, 5, particleID);
    analysisManager->FillNtupleSColumn(0, 6, particleName);
    analysisManager->FillNtupleIColumn(0, 7, copyNo);
    analysisManager->FillNtupleIColumn(0, 8, stepID);
    analysisManager->AddNtupleRow(0);

    fpEventAction->AddEdep(dE);
  }

  // if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-" || step->GetTrack()->GetParticleDefinition()->GetParticleName() == "gamma" || step->GetTrack()->GetParticleDefinition()->GetParticleName() == "alpha") && (volumeName == "waterBox")) // e- from outside the tracking volume
  if ((volumeNamePre == "waterBox") || (volumeNamePre == "bloodVessel"))
  {
    // if (particleName == "At211")
    // {
    //   G4double postEnergy = step->GetPostStepPoint()->GetKineticEnergy();
    //   G4double preEnergy = step->GetPostStepPoint()->GetKineticEnergy();
    //   G4cout << "At: pre E: " << preEnergy << ", post E: " << postEnergy << G4endl;
    // }
    if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "TrackingVol") // step ends on boundary to tracking volume
    {
      //G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
      //G4int stepID = step->GetTrack()->GetCurrentStepNumber();
      G4int trackID = step->GetTrack()->GetTrackID();
      // DEBUG
      //G4cout << "DEBUG: evt:" << eventID << ", " << particleName << "(ID: " << particleID << ") entered TrackingVol. trackID:  " << trackID << " stepID: " << stepID << G4endl;

      G4TouchableHandle touchable = step->GetPostStepPoint()->GetTouchableHandle();

      G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();

      G4ThreeVector localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

      G4ThreeVector worldMomentum = step->GetPostStepPoint()->GetMomentumDirection();

      G4ThreeVector localMomentum = (*(touchable->GetHistory()->GetTopVolume()->GetRotation())) * worldMomentum; // rotate momentum direction by volume rotation
      // G4cout << "world pos: " << worldPos.x() << ", " << worldPos.y() << ", " << worldPos.z() << G4endl;
      // G4cout << "local pos: " << localPos.x() << ", " << localPos.y() << ", " << localPos.z() << G4endl;

      // G4cout << "world mom: " << worldMomentum.x() << ", " << worldMomentum.y() << ", " << worldMomentum.z() << G4endl;
      // G4cout << "local mom: " << localMomentum.x() << ", " << localMomentum.y() << ", " << localMomentum.z() << G4endl;
      G4double time = step->GetPostStepPoint()->GetGlobalTime();

      float output[12];
      output[0] = localPos.x() / mm;
      output[1] = localPos.y() / mm;
      output[2] = localPos.z() / mm;
      output[3] = localMomentum.x();
      output[4] = localMomentum.y();
      output[5] = localMomentum.z();
      output[6] = particleEnergy / MeV;
      output[7] = eventID;
      output[8] = particleID;
      output[9] = copyNo;
      output[10] = time;
      output[11] = processID;
      //output[12] = ((const G4Ions*)( step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy();

      PSfile.write((char *)&output, sizeof(output));
      // step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

      fpEventAction->AddSecondary();
    }
  }

  if ((volumeNamePre == "TrackingVol") && (step->IsFirstStepInVolume())) //(step->GetTrack()->GetCurrentStepNumber() == 1)) // e- from interactions inside the tracking volume
  {
    if (step->GetPreStepPoint()->GetProcessDefinedStep() != nullptr)
      return; // if prestep process is nullptr this is the first step of particle created by interaction in the cell - only save those created by processes in cell

    if (step->GetTrack()->GetCreatorProcess()->GetProcessName() != "RadioactiveDecay")
      return; // only save products of radioactive decay other products are from parents which are saved on entering the cell and will be tracked in DNA simulation.

    G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();

    G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    G4ThreeVector worldMomentum = step->GetPreStepPoint()->GetMomentumDirection();
    G4ThreeVector localMomentum = (*(touchable->GetHistory()->GetTopVolume()->GetRotation())) * worldMomentum;

    G4double time = step->GetPreStepPoint()->GetGlobalTime();

    auto particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    // DEBUG
    //G4cout << "DEBUG: evt:" << eventID << ", " << particleID << "( " << particleName << ") inside tracking vol: " << G4endl;
    G4String process = step->GetTrack()->GetCreatorProcess()->GetProcessTypeName(step->GetTrack()->GetCreatorProcess()->GetProcessType());
    
    float output[12];
    output[0] = localPos.x() / mm;
    output[1] = localPos.y() / mm;
    output[2] = localPos.z() / mm;
    output[3] = localMomentum.x();
    output[4] = localMomentum.y();
    output[5] = localMomentum.z();
    output[6] = particleEnergy / MeV;
    output[7] = eventID;
    output[8] = particleID;
    output[9] = copyNo;
    output[10] = time;
    output[11] = processID;
    //output[12] = ((const G4Ions*)( step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy();
    // G4cout << "DEBUG: evt:" << eventID << ", ID: " << particleID << "( " << particleName << ")" << G4endl;
    PSfile.write((char *)&output, sizeof(output));
    // if (stepID >= 1)
    //{
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    //}
    fpEventAction->AddSecondary();
  }
}
