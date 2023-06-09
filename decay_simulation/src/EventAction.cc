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

EventAction::EventAction() : G4UserEventAction()
{
  fEdep = 0.;
  numSecondary = 0;
  fpathLengthTotal = 0;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event * event)
{
  fEdep = 0.;
  fpathLengthTotal = 0;
  while (!parentParticle.empty())
  {
    parentParticle.erase(parentParticle.begin());
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-out")) == 0)
  return;
  auto fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
  if ((fEdep>0) && (fpathLengthTotal>0)) //only save edep>0 to reduce output file size
  {
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(2,0, (fEdep/joule));
  analysisManager->FillNtupleDColumn(2, 1, (fEdep / MeV));
  analysisManager->FillNtupleIColumn(2,2, G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
  analysisManager->FillNtupleDColumn(2, 3, fpathLengthTotal / nanometer);
  analysisManager->FillNtupleIColumn(2, 4, fpEventAction->GetNumSecondaries());

  analysisManager->AddNtupleRow(2);
  }

  if (fpathLengthTotal>0)
  {
      auto fpRunAction = (RunAction *)G4RunManager::GetRunManager()->GetUserRunAction();
      fpRunAction->AddIntersecting();

  }
  // for (auto it = parentParticle.cbegin(); it != parentParticle.cend(); ++it)
  // {
  //     G4int particleID = (int)it->second;
  //     G4cout <<"trackID: " << it->first << ", parentID: " << it->second << ", particleName: "<< particleMapRev[particleID] << G4endl;
  // }
}
