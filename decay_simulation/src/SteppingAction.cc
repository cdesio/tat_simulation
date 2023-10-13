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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SteppingAction::UserSteppingAction(const G4Step *step)
{
  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  if (particleName == "nu_e")
  {
    return;
  }
  if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "world")
  {
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    return;
  }

  G4double dE = step->GetTotalEnergyDeposit();

  const G4String &volumeNamePre = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  // G4int pdg_enc = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  // G4int particleID = particleMap[particleName];
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
  //G4cout << "DEBUG: " << particleName << ", ID: " << particleID << " vol: " << volumeNamePre << G4endl;
  
  //tracking all particles
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4ThreeVector prePoint = step->GetPreStepPoint()->GetPosition();
  G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
  
  G4int copyNo = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
  G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);
  G4int stepID = step->GetTrack()->GetCurrentStepNumber();
  auto particleEnergy = step->GetPostStepPoint()->GetKineticEnergy();
  G4double steplength = step->GetStepLength();
  const G4VProcess *defprocess = step->GetPostStepPoint()->GetProcessDefinedStep();
  const G4VProcess *creatprocess = step->GetTrack()->GetCreatorProcess();
  G4int parentID = step->GetTrack()->GetParentID();
  G4int trackID = step->GetTrack()->GetTrackID();
  G4double time = step->GetPostStepPoint()->GetGlobalTime();
  G4int processID = -1;
  if (creatprocess)
  {
    
    G4String creatprocessname = creatprocess->GetProcessName();

    if ((fpEventAction->parentParticle.find(trackID) == fpEventAction->parentParticle.end())) // if track not in map keys
    {
      if (creatprocessname == "RadioactiveDecay")
      {
        fpEventAction->parentParticle.insert(std::pair<G4int, G4int>(trackID, particleMap[particleName]));
      }
      else
      {
        G4int parentParticleID = fpEventAction->parentParticle[parentID];
        fpEventAction->parentParticle.insert(std::pair<G4int, G4int>(trackID, parentParticleID));
      }
    }
  }

    if (step->GetPreStepPoint() == nullptr)
  {
    return; // primary particle has no pre-step point
  }


  // PS file and output root file

  //ROOT file
    CommandLineParser *parser = G4DNAPARSER::CommandLineParser::GetParser();
    Command *command(0);
      if ((command = parser->GetCommandIfActive("-out")) == 0)
        return;

    if ((volumeNamePre == "shell") && (dE > 0))
    {      
      G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
      G4int mapped_PID = fpEventAction->parentParticle[trackID];
    

      analysisManager->FillNtupleIColumn(0, 0, eventID);
      analysisManager->FillNtupleDColumn(0, 1, dE / eV);
      analysisManager->FillNtupleDColumn(0, 2, postPoint.x() / nanometer);
      analysisManager->FillNtupleDColumn(0, 3, postPoint.y() / nanometer);
      analysisManager->FillNtupleDColumn(0, 4, postPoint.z() / nanometer);
      analysisManager->FillNtupleIColumn(0, 5, particleID);
      analysisManager->FillNtupleSColumn(0, 6, particleName);
      analysisManager->FillNtupleIColumn(0, 7, copyNo);
      analysisManager->FillNtupleIColumn(0, 8, stepID);
      analysisManager->FillNtupleIColumn(0, 9, mapped_PID);
      analysisManager->FillNtupleDColumn(0, 10, particleEnergy/MeV);
      analysisManager->FillNtupleDColumn(0, 11, steplength);
      analysisManager->AddNtupleRow(0);

      fpEventAction->AddEdep(dE);
    }
    
    else{
      if ((volumeNamePre == "water") && (dE > 0))
    {
      G4DNAPARSER::CommandLineParser *parser = G4DNAPARSER::CommandLineParser::GetParser();
      G4DNAPARSER::Command *command(0);
      if ((command = parser->GetCommandIfActive("-out")) == 0)
        return;

      G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

      G4int mapped_PID = fpEventAction->parentParticle[trackID];

      analysisManager->FillNtupleIColumn(0, 0, eventID);
      analysisManager->FillNtupleDColumn(0, 1, dE / eV);
      analysisManager->FillNtupleDColumn(0, 2, postPoint.x() / nanometer);
      analysisManager->FillNtupleDColumn(0, 3, postPoint.y() / nanometer);
      analysisManager->FillNtupleDColumn(0, 4, postPoint.z() / nanometer);
      analysisManager->FillNtupleIColumn(0, 5, particleID);
      analysisManager->FillNtupleSColumn(0, 6, particleName);
      analysisManager->FillNtupleIColumn(0, 7, -99);
      analysisManager->FillNtupleIColumn(0, 8, stepID);
      analysisManager->FillNtupleIColumn(0, 9, mapped_PID);
      analysisManager->FillNtupleDColumn(0, 10, particleEnergy/MeV);
      analysisManager->FillNtupleDColumn(0, 11, steplength);
      analysisManager->AddNtupleRow(0);

      fpEventAction->AddEdep(dE);
    }
    }

    // PS file
    
    if ((volumeNamePre == "water") && (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "shell") && (step->GetPreStepPoint()->GetKineticEnergy() > 0)) 
    {
          G4int trackID = step->GetTrack()->GetTrackID();
         // G4TouchableHandle touchable = step->GetPostStepPoint()->GetTouchableHandle();

          G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();
          //G4cout << "worldPos: " << worldPos.x() << ", " <<  worldPos.y() << ", " <<  worldPos.z() << G4endl;
          G4double newX = (G4UniformRand() * .0003) - .00015;
          G4double newZ = (G4UniformRand() * .0003) - .00015;

          G4double radius = std::pow(worldPos.x() * worldPos.x() + worldPos.y() * worldPos.y(), 0.5);

          G4double newY = radius - fpDetector->get_start_R() + (step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo()) * fpDetector->get_spacing();
          //G4cout << "newY: " << newY << ", cNo: " << step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo() << ", radius: " << radius << ", s*cno: " << fpDetector->get_spacing() << G4endl;
          //G4ThreeVector localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
          G4ThreeVector newPos = G4ThreeVector(newX, newY, newZ);
          G4ThreeVector worldMomentum = step->GetPostStepPoint()->GetMomentumDirection();
          G4ThreeVector newMomentum = transformDirection(worldPos, worldMomentum);
        
          G4int parentID = step->GetTrack()->GetParentID();
          G4int mapped_PID = fpEventAction->parentParticle[trackID];
          write_to_PS(step->GetTrack(), newPos, newMomentum, particleEnergy, eventID, particleID, copyNo, time, mapped_PID);
        
    }
    else if ((volumeNamePre == "shell") && (step->IsFirstStepInVolume())  && (particleName != "gamma") && (step->GetPreStepPoint()->GetProcessDefinedStep() == nullptr)) //(step->GetTrack()->GetCurrentStepNumber() == 1)) // e- from interactions inside the tracking volume
    {
      if ((step->GetTrack()->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay") && (((const G4Ions *)(step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy() < 1e-15))
      {
        G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();

        G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();
        //G4ThreeVector localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
        G4ThreeVector worldMomentum = step->GetPreStepPoint()->GetMomentumDirection();
        //G4ThreeVector localMomentum = (*(touchable->GetHistory()->GetTopVolume()->GetRotation())) * worldMomentum;

        G4double time = step->GetPreStepPoint()->GetGlobalTime();

        auto particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
        G4int parentID = step->GetTrack()->GetParentID();

        G4int stepID = step->GetTrack()->GetCurrentStepNumber();
        G4int trackID = step->GetTrack()->GetTrackID();
        G4int mapped_PID = fpEventAction->parentParticle[trackID];
       if (fpEventAction->decayPos.find(parentID) == fpEventAction->decayPos.end())
      {
        G4double newX = (G4UniformRand() * .0003) - .00015;
        G4double newZ = (G4UniformRand() * .0003) - .00015;

        G4double radius = std::pow(worldPos.x() * worldPos.x() + worldPos.y() * worldPos.y(), 0.5);

        G4double newY = radius - fpDetector->get_start_R() + (step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo()) * fpDetector->get_spacing();
        
        G4ThreeVector newPos = G4ThreeVector(newX, newY, newZ);

        fpEventAction->decayPos.insert(std::pair<int, G4ThreeVector>(parentID, newPos));

        G4ThreeVector worldMomentum = step->GetPostStepPoint()->GetMomentumDirection();
        G4ThreeVector newMomentum = transformDirection(worldPos, worldMomentum);
        write_to_PS(step->GetTrack(), newPos, newMomentum, particleEnergy, eventID, particleID, copyNo, time, mapped_PID);
      }
      else
      {
        G4ThreeVector newPos = fpEventAction->decayPos[parentID];
        G4ThreeVector newMomentum = transformDirection(worldPos, worldMomentum);
        write_to_PS(step->GetTrack(), newPos, newMomentum, particleEnergy, eventID, particleID, copyNo, time, mapped_PID);
      }
    }
    else
    {
    G4int parentID = step->GetTrack()->GetParentID();
    G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();
    G4double newX = (G4UniformRand() * .0003) - .00015;
    G4double newZ = (G4UniformRand() * .0003) - .00015;

    G4double radius = std::pow(worldPos.x() * worldPos.x() + worldPos.y() * worldPos.y(), 0.5);

    G4double newY = radius - fpDetector->get_start_R() + (step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo()) * fpDetector->get_spacing();
    
    G4ThreeVector newPos = G4ThreeVector(newX, newY, newZ);

    fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID());
    fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));

    fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID()); // erase distance travelled for this track in the box
    fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), G4ThreeVector()));
    
    fpEventAction->tracks.push_back(step->GetTrack()->GetTrackID()); 
    }

    }
    else if ((volumeNamePre == "shell") && ((std::find(fpEventAction->tracks.begin(), fpEventAction->tracks.end(), step->GetTrack()->GetTrackID()) != fpEventAction->tracks.end())) && (particleName != "gamma") && (step->GetPreStepPoint()->GetKineticEnergy() > 0))
    {
    G4ThreeVector entryPosition = fpEventAction->particlePos[step->GetTrack()->GetTrackID()]; // look up position in box frame from last step

    G4ThreeVector deltaWorld = step->GetPostStepPoint()->GetPosition() - step->GetPreStepPoint()->GetPosition(); // change in position in world frame
    G4ThreeVector boxMomentumPre = transformDirection(step->GetPreStepPoint()->GetPosition(), step->GetPreStepPoint()->GetMomentumDirection()); // particle momentum in box frame

    G4ThreeVector delta = deltaWorld.mag() * boxMomentumPre; // change in position in box frame

    G4ThreeVector postStepBox = entryPosition + (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()] + delta;              // post step position in box frame
    if ((std::abs(postStepBox.x()) >= 0.00015) || (std::abs(postStepBox.y()) >= 0.00015) || (std::abs(postStepBox.z()) >= 0.00015)) // if >=0.00015 has crossed the boundary
    {
      // save particle, new position and distance saved
      G4ThreeVector preStepBox = entryPosition + (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()]; // pre step point position in box frame

      if (std::abs(preStepBox.y()) >= 0.00015)
      {
        // Particles which scatter back into the box are not added to the phase space file as scattering is included in the DNA simulation
        return;
      }

      G4double distanceToExit = calculateDistanceToExitBox(preStepBox, boxMomentumPre);

      if (distanceToExit == DBL_MAX)
      // exit y
      // update start position to y exit point and zero distance travelled, in case scattering changes direction
      {
        G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
        G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
        tYneg = tYneg > 1e-15 ? tYneg : DBL_MAX;
        tYpos = tYpos > 1e-15 ? tYpos : DBL_MAX;

        G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

        G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

        fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
        fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

        fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
        fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));

        return;
      }
      else // crosses x or z
      {
        G4double stepDistance = step->GetStepLength();

        G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

        // check which side of the box was crossed and change sign as particle is entering adjacent box
        if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-15))
          newPos.setX(-0.00015);
        else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-15))
          newPos.setX(+0.00015);

        else if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-15))
          newPos.setZ(-0.00015);
        else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-15))
          newPos.setZ(+0.00015);

        G4double percentageOfStep = distanceToExit / stepDistance;

        // calculate KE at point where crossing occurs
        G4double newKE = step->GetPreStepPoint()->GetKineticEnergy() - (step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

        // calculate time at point where crossing occurs
        G4double newTime = step->GetPreStepPoint()->GetGlobalTime() + (step->GetDeltaTime() * percentageOfStep);
        G4int mapped_PID = fpEventAction->parentParticle[trackID];
        //savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);
        write_to_PS(step->GetTrack(), newPos, boxMomentumPre, newKE, eventID, particleID, copyNo, newTime, mapped_PID);
        G4double percentageAccountedFor = percentageOfStep;

        while (percentageAccountedFor < 1)
        {
          G4ThreeVector restOfStep = newPos + (stepDistance - distanceToExit) * boxMomentumPre;

          if ((std::abs(restOfStep.x()) < 0.00015) && (std::abs(restOfStep.y()) < 0.00015) && (std::abs(restOfStep.z()) < 0.00015))
          {
            // remainder of step is contained in the adjacent box
            // save remainder of track travel to next box
            fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
            fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), G4ThreeVector()));

            fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track

            fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos + (stepDistance - distanceToExit) * boxMomentumPre)); // add current box entry position for this track
            percentageAccountedFor = 1;
          }
          else
          {
            // crosses another boundary
            // find which boundary
            G4double distanceToExitRemainder = calculateDistanceToExitBox(newPos, boxMomentumPre);

            if (distanceToExitRemainder == DBL_MAX) // exit y
            // update start position to y exit point and zero distance travelled, in case scattering changes direction
            {
              G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
              G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
              tYneg = tYneg > 1e-15 ? tYneg : DBL_MAX;
              tYpos = tYpos > 1e-15 ? tYpos : DBL_MAX;

              G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

              G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

              fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
              fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

              fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
              fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));

              percentageAccountedFor = 1;

              return;
            }

            newPos += (distanceToExitRemainder * boxMomentumPre);

            // check which side of the box was crossed and change sign as particle is entering adjacent box
            if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-15))
              newPos.setX(-0.00015);
            else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-15))
              newPos.setX(+0.00015);

            else if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-15))
              newPos.setZ(-0.00015);
            else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-15))
              newPos.setZ(+0.00015);

            percentageOfStep = distanceToExitRemainder / stepDistance;

            newKE = newKE - (newKE - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

            newTime += (step->GetDeltaTime() * percentageOfStep);
            distanceToExit += distanceToExitRemainder;

            //savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);
            write_to_PS(step->GetTrack(), newPos, boxMomentumPre, newKE, eventID, particleID, copyNo, newTime, mapped_PID);
            percentageAccountedFor += percentageOfStep;
          }
        }
      }
    }

    else
    {
      //  if position hasn't crossed the box bounday update distance travelled in box
      G4ThreeVector previousDelta = (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()];

      fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
      fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), previousDelta + delta));
    }
  }
  else if ((volumeNamePre == "shel") && ((std::find(fpEventAction->tracks.begin(), fpEventAction->tracks.end(), step->GetTrack()->GetTrackID()) != fpEventAction->tracks.end())) && (particleName == "gamma")) 
{
    G4ThreeVector preStepBox = fpEventAction->particlePos[step->GetTrack()->GetTrackID()]; // look up position in box frame from last step

    G4ThreeVector deltaWorld = step->GetPostStepPoint()->GetPosition() - step->GetPreStepPoint()->GetPosition(); // change in position in world frame

    G4ThreeVector boxMomentumPre = transformDirection(step->GetPreStepPoint()->GetPosition(), step->GetPreStepPoint()->GetMomentumDirection()); // particle momentum in box frame

    G4ThreeVector delta = deltaWorld.mag() * boxMomentumPre; // change in position in box frame

    G4double distanceToExit = calculateDistanceToExitBox(preStepBox, boxMomentumPre);

    if (distanceToExit == DBL_MAX)
    // exit y
    // update start position to y exit point and zero distance travelled, in case scattering changes direction
    {
      G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
      G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
      tYneg = tYneg > 1e-15 ? tYneg : DBL_MAX;
      tYpos = tYpos > 1e-15 ? tYpos : DBL_MAX;

      G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

      G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

      fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
      fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

      fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
      fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));
      
    }
    else if (distanceToExit > deltaWorld.mag())
    {
      G4ThreeVector boxMomentumPost = transformDirection(step->GetPostStepPoint()->GetPosition(), step->GetPostStepPoint()->GetMomentumDirection()); // particle momentum in box frame for next step
      G4int mapped_PID = fpEventAction->parentParticle[trackID];
      // step is contained within box save end point as start of next step
      //savePoint(step->GetTrack(), preStepBox + delta, boxMomentumPost, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), step->GetPreStepPoint()->GetKineticEnergy(), step->GetPreStepPoint()->GetGlobalTime(), fpEventAction->parentParticle[TrackID]);
      write_to_PS(step->GetTrack(), preStepBox+delta, boxMomentumPost, particleEnergy, eventID, particleID, copyNo, time, mapped_PID);
    }
    else
    {
      // crosses x or z boundary
      G4double stepDistance = step->GetStepLength();

      G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

      // check which side of the box was crossed and change sign as particle is entering adjacent box
      if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-15))
        newPos.setX(-0.00015);
      else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-15))
        newPos.setX(+0.00015);

      else if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-15))
        newPos.setZ(-0.00015);
      else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-15))
        newPos.setZ(+0.00015);

      G4double percentageOfStep = distanceToExit / stepDistance;

      // calculate KE at point where crossing occurs
      G4double newKE = step->GetPreStepPoint()->GetKineticEnergy() - (step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

      // calculate time at point where crossing occurs
      G4double newTime = step->GetPreStepPoint()->GetGlobalTime() + (step->GetDeltaTime() * percentageOfStep);
      G4int mapped_PID = fpEventAction->parentParticle[trackID];
      //savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);
      write_to_PS(step->GetTrack(), newPos, boxMomentumPre, newKE, eventID, particleID, copyNo, newTime, mapped_PID);
      G4double percentageAccountedFor = percentageOfStep;

      while (percentageAccountedFor < 1)
      {
        G4ThreeVector restOfStep = newPos + (stepDistance - distanceToExit) * boxMomentumPre;

        if ((std::abs(restOfStep.x()) < 0.00015) && (std::abs(restOfStep.y()) < 0.00015) && (std::abs(restOfStep.z()) < 0.00015))
        {
          // remainder of step is contained in the adjacent box
          // save remainder of track travel to next box
          fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
          fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), G4ThreeVector()));

          fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track

          fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos + (stepDistance - distanceToExit) * boxMomentumPre)); // add current box entry position for this track
          percentageAccountedFor = 1;
        }
        else
        {
          // crosses another boundary
          // find which boundary
          G4double distanceToExitRemainder = calculateDistanceToExitBox(newPos, boxMomentumPre);

          if (distanceToExitRemainder == DBL_MAX) // exit y
          // update start position to y exit point and zero distance travelled, in case scattering changes direction
          {
            G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
            G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
            tYneg = tYneg > 1e-15 ? tYneg : DBL_MAX;
            tYpos = tYpos > 1e-15 ? tYpos : DBL_MAX;

            G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

            newPos = preStepBox + (distanceToExit * boxMomentumPre);

            fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
            fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

            fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
            fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));

            percentageAccountedFor = 1;

            return;
          }

          newPos += (distanceToExitRemainder * boxMomentumPre);

          // check which side of the box was crossed and change sign as particle is entering adjacent box

          // G4cout << newPos.x() << G4endl;
          if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-15))
            newPos.setX(-0.00015);
          else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-15))
            newPos.setX(+0.00015);

          else if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-15))
            newPos.setZ(-0.00015);
          else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-15))
            newPos.setZ(+0.00015);

          percentageOfStep = distanceToExitRemainder / stepDistance;

          newKE = newKE - (newKE - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

          newTime += (step->GetDeltaTime() * percentageOfStep);
          distanceToExit += distanceToExitRemainder;

          //savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);
          write_to_PS(step->GetTrack(), newPos, boxMomentumPre, newKE, eventID, particleID, copyNo, newTime, mapped_PID);
          percentageAccountedFor += percentageOfStep;
        }
      }
    }
  }       
}


G4ThreeVector SteppingAction::transformDirection(const G4ThreeVector &position, const G4ThreeVector &worldMomentum)
{
  G4double theta = std::asin(position.x() / std::pow(position.y() * position.y() + position.x() * position.x(), 0.5));
  if ((position.x() > 0) && (position.y() > 0))
    // # positive-positive quadrant
    theta = theta;
  else if ((position.x() > 0) && (position.y() < 0))
    // # positive-negative quadrant
    theta = 3.14159 - theta;
  else if ((position.x() < 0) && (position.y() < 0))
    // # negative-negative quadrant
    theta = std::abs(theta) + 3.14159;
  else if ((position.x() < 0) && (position.y() > 0))
    // # negative-positive quadrant
    theta = theta;

  G4ThreeVector newMomentum = G4ThreeVector(worldMomentum.x() * std::cos(theta) - worldMomentum.y() * std::sin(theta), worldMomentum.x() * std::sin(theta) + worldMomentum.y() * std::cos(theta), worldMomentum.z());

  return newMomentum;
}




void SteppingAction::write_to_PS(const G4Track *track, const G4ThreeVector &Position, const G4ThreeVector &Momentum, G4double &particleEnergy, G4int &eventID, G4int &PID, G4int &copyNo, G4double &time, G4int &mapped_PID)

{
fpEventAction->particlePos.erase(track->GetTrackID());
fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(track->GetTrackID(), Position));

fpEventAction->particleDist.erase(track->GetTrackID()); 
fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(track->GetTrackID(), G4ThreeVector()));

float output[12];
output[0] = Position.x() / mm;
output[1] = Position.y() / mm;
output[2] = Position.z() / mm;
output[3] = Momentum.x();
output[4] = Momentum.y();
output[5] = Momentum.z();
output[6] = particleEnergy / MeV;
output[7] = eventID;
output[8] = PID;
output[9] = copyNo;
output[10] = time;
output[11] = mapped_PID;

PSfile.write((char *)&output, sizeof(output));
fpEventAction->tracks.push_back(track->GetTrackID());

}



G4double SteppingAction::calculateDistanceToExitBox(const G4ThreeVector &preStepPosition, const G4ThreeVector &preStepMomentumDirection)
{
  // does step exit box in x, y or z?
  G4double tXneg = (-.00015 - preStepPosition.x()) / preStepMomentumDirection.x();
  G4double tXpos = (.00015 - preStepPosition.x()) / preStepMomentumDirection.x();

  G4double tYneg = (-.00015 - preStepPosition.y()) / preStepMomentumDirection.y();
  G4double tYpos = (.00015 - preStepPosition.y()) / preStepMomentumDirection.y();

  G4double tZneg = (-.00015 - preStepPosition.z()) / preStepMomentumDirection.z();
  G4double tZpos = (.00015 - preStepPosition.z()) / preStepMomentumDirection.z();

  // G4cout << tXneg << " " << tXpos << " " << tYneg << " " << tYpos << " " << tZneg << " " << tZpos << " " << G4endl;

  tXneg = tXneg > 1e-15 ? tXneg : DBL_MAX;
  tXpos = tXpos > 1e-15 ? tXpos : DBL_MAX;
  tYneg = tYneg > 1e-15 ? tYneg : DBL_MAX;
  tYpos = tYpos > 1e-15 ? tYpos : DBL_MAX;
  tZneg = tZneg > 1e-15 ? tZneg : DBL_MAX;
  tZpos = tZpos > 1e-15 ? tZpos : DBL_MAX;

  G4double distanceToExit = std::min({tXneg, tXpos, tZneg, tZpos}); // shortest distance travelled to cross box surface

  // G4cout << tXneg << " " << tXpos << " " << tYneg << " " << tYpos << " " << tZneg << " " << tZpos << " " << G4endl;
  if ((tYneg <= distanceToExit) || (tYpos <= distanceToExit)) // exit y
  {
    return DBL_MAX;
  }

  return distanceToExit;
}