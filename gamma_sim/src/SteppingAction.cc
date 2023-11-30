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

  const G4String &volumeName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

    if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "gamma") && (volumeName == "chromatinSegment"))
  {
    G4DNAPARSER::CommandLineParser *parser = G4DNAPARSER::CommandLineParser::GetParser();
    // G4DNAPARSER::Command *command(0);
    if (parser->GetCommandIfActive("-out") == 0)
      return;

    fpEventAction->AddPathLength(step->GetStepLength());
  }

  if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() != "e-") && (volumeName == "chromatinSegment") && (dE > 0))
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

    analysisManager->FillNtupleIColumn(0, 0, eventID);
    analysisManager->FillNtupleDColumn(0, 1, dE / eV);
    analysisManager->FillNtupleDColumn(0, 2, postPoint.x() / nanometer);
    analysisManager->FillNtupleDColumn(0, 3, postPoint.y() / nanometer);
    analysisManager->FillNtupleDColumn(0, 4, postPoint.z() / nanometer);
    analysisManager->AddNtupleRow(0);

    fpEventAction->AddEdep(dE);
  }

  if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-") && (volumeName == "waterBox")) // e- from outside the tracking volume
  {
    if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "TrackingVol") // step ends on boundary to tracking volume
    {
      double positionX = step->GetPostStepPoint()->GetPosition().x();
      double positionY = step->GetPostStepPoint()->GetPosition().y();
      double positionZ = step->GetPostStepPoint()->GetPosition().z();
      double momentumX = step->GetPostStepPoint()->GetMomentumDirection().x();
      double momentumY = step->GetPostStepPoint()->GetMomentumDirection().y();
      double momentumZ = step->GetPostStepPoint()->GetMomentumDirection().z();
      double particleEnergy = step->GetPostStepPoint()->GetKineticEnergy();
      G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

      double output[8];
      output[0] = positionX / mm;
      output[1] = positionY / mm;
      output[2] = positionZ / mm;
      output[3] = momentumX;
      output[4] = momentumY;
      output[5] = momentumZ;
      output[6] = particleEnergy / MeV;
      output[7] = eventID;

      PSfile.write((char *)&output, sizeof(output));
      step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

      fpEventAction->AddSecondary();
    }
  }

  if ((step->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-") && ((volumeName == "TrackingVol") || (volumeName == "chromatinSegment"))) // e- from interactions inside the tracking volume
  {

    double positionX = step->GetPreStepPoint()->GetPosition().x();
    double positionY = step->GetPreStepPoint()->GetPosition().y();
    double positionZ = step->GetPreStepPoint()->GetPosition().z();
    double momentumX = step->GetPreStepPoint()->GetMomentumDirection().x();
    double momentumY = step->GetPreStepPoint()->GetMomentumDirection().y();
    double momentumZ = step->GetPreStepPoint()->GetMomentumDirection().z();
    double particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    double output[8];
    output[0] = positionX / mm;
    output[1] = positionY / mm;
    output[2] = positionZ / mm;
    output[3] = momentumX;
    output[4] = momentumY;
    output[5] = momentumZ;
    output[6] = particleEnergy / MeV;
    output[7] = eventID;

    PSfile.write((char *)&output, sizeof(output));
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

    fpEventAction->AddSecondary();
  }
}
