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
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "G4AnalysisManager.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "CommandLineParser.hh"

using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* pDetector) : G4UserEventAction(), dkill(pDetector->Getdkill())
{
  fEdep = 0.;
  fTrackStartKE = 0;
  fTrackEndKE = 0;
  fTrackStartFound = false;
  fTrackStartPos = G4ThreeVector();
  fTrackEndPos = G4ThreeVector();
  fTrackStoppedBox = false;
  fpathLengthTotal = 0;
  fPositions0Event = new G4KDTree();
  fPositions1Event = new G4KDTree();
  fPositionsBase0Event = new G4KDTree();
  fPositionsBase1Event = new G4KDTree();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *event)
{
  fEdep = 0.;
  fTrackStartKE = 0;
  fTrackStartFound = false;
  fTrackStoppedBox = false;
  fTrackEndKE = 0;
  fTrackStartPos = G4ThreeVector();
  fTrackEndPos = G4ThreeVector();
  fpathLengthTotal = 0;
  fPositions0Event->Clear();
  fPositions1Event->Clear();
  fPositionsBase0Event->Clear();
  fPositionsBase1Event->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *)
{
  G4double projectedRange = sqrt(GetStartTrackPos().diff2(GetEndTrackPos()));

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  if (((command = parser->GetCommandIfActive("-PS"))==0) && (projectedRange > 0)) // if the track passes through the chromatin fibre volume
  {
    fTrackMeanKE.push_back((fTrackStartKE + fTrackEndKE) / 2);
    analysisManager->FillNtupleDColumn(0, 0, (fEdep / joule));
    analysisManager->FillNtupleIColumn(0, 1, G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
    analysisManager->FillNtupleDColumn(0, 2, fTrackStartKE / MeV);
    analysisManager->FillNtupleDColumn(0, 3, fTrackEndKE / MeV);
    analysisManager->FillNtupleDColumn(0, 4, projectedRange / nanometer);
    analysisManager->FillNtupleDColumn(0, 5, fpathLengthTotal / nanometer);
    analysisManager->FillNtupleIColumn(0, 6, fTrackStoppedBox);
    analysisManager->AddNtupleRow();
  }
  else if (command = parser->GetCommandIfActive("-PS")) // all secondary electrons in phase space file are in the volume
  {
    analysisManager->FillNtupleDColumn(0, 0, (fEdep / joule));
    analysisManager->FillNtupleDColumn(0, 1, (fEdep / MeV));
    analysisManager->FillNtupleIColumn(0, 2, G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
    analysisManager->AddNtupleRow();
  }
}
