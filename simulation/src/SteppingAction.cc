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
#include "G4ITTrackHolder.hh"
#include "G4Track.hh"
#include <map>
#include "globals.hh"
#include "G4Molecule.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4DNAElastic.hh"
#include "G4DNAElectronSolvation.hh"
#include "DetectorConstruction.hh"
#include "CommandLineParser.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hh"

using MapOfDelayedLists =
    std::map<double, std::map<int, G4TrackList *>>;
using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(DetectorConstruction* fpDet)
    : G4UserSteppingAction(), fpEventAction(0)
 , fpDetector(fpDet)
{
  fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SteppingAction::UserSteppingAction(const G4Step *step)
{
  G4double dE = 0;
  G4int flagVolume = 0.;

  const G4String &volumeName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if (volumeName == "sugar0")
  {
    flagVolume = 1;
  }
  else if (volumeName == "sugar1")
  {
    flagVolume = 2;
  }
  else if (volumeName == "histone")
  {
    flagVolume = 3;
  }
  else if (volumeName == "chromatinSegment")
  {
    flagVolume = 4;
  }
  else if (volumeName == "TrackingVol")
  {
    flagVolume = 5;
  }
  if (flagVolume != 0)
  {

    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

    const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    G4String primaryName = generatorAction->primaryName;

    if ((((particleName == "alpha") || (particleName == "alpha+") || (particleName == "helium")) && (primaryName == "alpha")) || ((primaryName=="e-") && (step->GetTrack()->GetTrackID() == 1)) || (primaryName=="gamma")|| (particleName=="neutron"))
    {
      if (false == fpEventAction->GetStartTrackFound())
      {
        fpEventAction->SetStartTrackKE(step->GetPreStepPoint()->GetKineticEnergy());
        fpEventAction->SetStartTrackFound();
        fpEventAction->SetStartTrackPos(step->GetPreStepPoint()->GetPosition());
      }
      if (step->GetPostStepPoint()->GetKineticEnergy() == 0)
      {
        fpEventAction->SetStoppedInBox();
        fpEventAction->SetEndTrackKE(step->GetPostStepPoint()->GetKineticEnergy());
        fpEventAction->SetEndTrackPos(step->GetPostStepPoint()->GetPosition());
        fpEventAction->AddPathLength(step->GetStepLength());
      }
      if ((false == fpEventAction->GetStoppedInBox()))
      {
        fpEventAction->SetEndTrackKE(step->GetPostStepPoint()->GetKineticEnergy());
        fpEventAction->SetEndTrackPos(step->GetPostStepPoint()->GetPosition());
        fpEventAction->AddPathLength(step->GetStepLength());
      }
    }
  }

  dE = step->GetTotalEnergyDeposit();
  if ((volumeName != "world") && (volumeName != "waterBox") && (volumeName != "TrackingVol") && (dE > 0.))
  {

    fpEventAction->AddEdep(dE);
    const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    G4int step1_eventID = generatorAction->step1_eventID;
    G4int step2_eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    G4int step1_copyNo = generatorAction->step1_copyNo;
    G4cout << "stepping - 1evtID: " << step1_eventID << ", 2evtID: " << step2_eventID << ", copyNo: " << step1_copyNo << G4endl;
    // G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    // G4int step2_eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    // const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
    //       G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    
   
    // G4double step1_time = generatorAction->step1_time;
    // G4int step1_PID = generatorAction->step1_PID;
    
    // G4int step1_primaryID = generatorAction->step1_primaryID;
    // G4double step2_time = step1_time + step->GetTrack()->GetGlobalTime();
    // G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    // G4int parentID = step->GetTrack()->GetParentID();
    // G4int trackID = step->GetTrack()->GetTrackID();
    // if (step1_eventID == 9)
    // {
    //   G4cout << "DEBUG: dE" << dE << ", step1: " << step1_eventID << ", step2: " << step2_eventID << ", copyNo: " << step1_copyNo << ", trackID: " << trackID << ", parent: " << parentID << G4endl;
    // }
    // analysisManager->FillNtupleIColumn(5, 0, step2_eventID);
    // analysisManager->FillNtupleIColumn(5, 1, step1_eventID);
    // analysisManager->FillNtupleIColumn(5, 2, step1_copyNo);
    // analysisManager->FillNtupleIColumn(5, 3, step1_PID);
    // analysisManager->FillNtupleDColumn(5, 4, step2_time);
    // analysisManager->FillNtupleDColumn(5, 5, dE/joule);
    // analysisManager->FillNtupleIColumn(5, 6, step1_primaryID);
    // analysisManager->FillNtupleSColumn(5, 7, particleName);
    // analysisManager->AddNtupleRow(5);
    }

  if ((flagVolume == 0) || (flagVolume==5))
  {
    return;
  }

  // remove water molecules created in the sugar volumes or histones
  if ((flagVolume == 1) || (flagVolume == 2) || (flagVolume == 3))
  {
    RemoveTracks(step);
  }

  if ((flagVolume != 3) && (dE > 0))
  {
    CommandLineParser *parser = CommandLineParser::GetParser();
    Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
      return;

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    
    G4String primaryName = generatorAction->primaryName;
    G4int primary_pid = generatorAction->step1_PID;
    G4int step1_copyNo = generatorAction->step1_copyNo;
    G4double step1_time = generatorAction->step1_time;
    G4int step1_primaryID = generatorAction->step1_primaryID;
    G4int step1_PID = generatorAction->step1_PID;
    G4int step1_eventID = generatorAction->step1_eventID;
    G4int step2_eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    G4ThreeVector prePoint = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);
    G4double step2_time = step1_time + step->GetTrack()->GetGlobalTime();
    //G4cout << "dE: " << dE << ", eV: " << eV << ", dE*eV: " << dE/eV << G4endl;

    auto fPositions0 = ((DetectorConstruction *)fpDetector)->fPositions0;
    auto fPositionsBase0 = ((DetectorConstruction *)fpDetector)->fPositionsBase0;
    auto fPositions1 = ((DetectorConstruction *)fpDetector)->fPositions1;
    auto fPositionsBase1 = ((DetectorConstruction *)fpDetector)->fPositionsBase1;
    auto result0 = fPositions0->NearestInRange(point, 1 * nm);
    auto result1 = fPositions1->NearestInRange(point, 1 * nm);
    if ((result0->GetSize() == 0) && (result1->GetSize() == 0))
      return;
    
    
    analysisManager->FillNtupleIColumn(1, 0, step2_eventID);
    analysisManager->FillNtupleIColumn(1, 1, step1_eventID);
    analysisManager->FillNtupleDColumn(1, 2, dE / eV);
    analysisManager->FillNtupleDColumn(1, 3, point.x() / nanometer);
    analysisManager->FillNtupleDColumn(1, 4, point.y() / nanometer);
    analysisManager->FillNtupleDColumn(1, 5, point.z() / nanometer);
    analysisManager->FillNtupleSColumn(1, 6, step->GetTrack()->GetParticleDefinition()->GetParticleName());
    analysisManager->FillNtupleDColumn(1, 7, step->GetPostStepPoint()->GetKineticEnergy());
    analysisManager->FillNtupleIColumn(1, 8, step1_copyNo);
    analysisManager->FillNtupleIColumn(1, 9, step1_PID);
    analysisManager->FillNtupleDColumn(1, 10, step1_time + step->GetTrack()->GetGlobalTime());
    analysisManager->FillNtupleIColumn(1, 11, step1_primaryID);
    analysisManager->AddNtupleRow(1);
  }
}

void SteppingAction::RemoveTracks(const G4Step *step)
{
  MapOfDelayedLists delayList = G4ITTrackHolder::Instance()->GetDelayedLists();
  for (auto &delayedmap_it : delayList)
  {
    for (auto &trackList : delayedmap_it.second)
    {
      if (nullptr == trackList.second)
      {
        continue;
      }
      G4TrackList::iterator itt = trackList.second->begin();
      G4TrackList::iterator endd = trackList.second->end();
      for (; itt != endd; ++itt)
      {
        G4Track *track = *itt;

        if ((track->GetParentID() !=
             step->GetTrack()->GetTrackID()) ||
            (track->GetPosition() !=
             step->GetPostStepPoint()->GetPosition()))
        {
          continue;
        }
        track->SetTrackStatus(fKillTrackAndSecondaries);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
