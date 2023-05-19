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

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "globals.hh"
#include <map>
#include "CommandLineParser.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4Event.hh"
#include "DetectorConstruction.hh"
#include "git_version.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
using namespace G4DNAPARSER;

RunAction::RunAction(DetectorConstruction *pDetector)
    : G4UserRunAction(), fpDetector(pDetector)
{
  chromatinVolume = fpDetector->chromatinVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run *)
{
  CreateNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run *run)
{
  WriteNtuple(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::CreateNtuple()
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;

  // Open an output file
  G4String fileName{"output.root"};
  if (command->GetOption().empty() == false)
  {
    fileName = command->GetOption();
  }

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleDirectoryName("ntuple");

  // open output file
  //
  G4bool fileOpen = analysisManager->OpenFile(fileName);
  if (!fileOpen)
  {
    G4cout << "\n---> HistoManager::book(): cannot open " << fileName << G4endl;
    return;
  }

  analysisManager->SetFirstNtupleId(0);

  analysisManager->CreateNtuple("Direct", "Direct");
  analysisManager->CreateNtupleIColumn(0, "EventNo");
  analysisManager->CreateNtupleDColumn(0, "eDep_eV");
  analysisManager->CreateNtupleDColumn(0, "x");
  analysisManager->CreateNtupleDColumn(0, "y");
  analysisManager->CreateNtupleDColumn(0, "z");
  analysisManager->CreateNtupleIColumn(0, "PID");
  analysisManager->CreateNtupleSColumn(0, "particle_name");
  analysisManager->CreateNtupleIColumn(0, "copyNo");
  analysisManager->CreateNtupleIColumn(0, "stepID");

  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("Info", "Info");
  analysisManager->CreateNtupleDColumn(1, "ChromatinVolume_m3");
  analysisManager->CreateNtupleIColumn(1, "NumSecondaries");
  analysisManager->CreateNtupleIColumn(1, "NumEvents");
  analysisManager->CreateNtupleSColumn(1, "GitHash");
  analysisManager->CreateNtupleIColumn(1, "NumIntersecting");
  analysisManager->FinishNtuple(1);

  analysisManager->CreateNtuple("EventData", "EventData");
  analysisManager->CreateNtupleDColumn(2, "Edep_J");
  analysisManager->CreateNtupleDColumn(2, "Edep_MeV");
  analysisManager->CreateNtupleIColumn(2, "EventNo");
  analysisManager->CreateNtupleDColumn(2, "PathLengthChromatin");
  analysisManager->CreateNtupleIColumn(2, "NumSecondaries");

  analysisManager->FinishNtuple(2);

  G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void RunAction::WriteNtuple(const G4Run *run)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  auto fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();

  analysisManager->FillNtupleDColumn(1, 0, chromatinVolume);
  analysisManager->FillNtupleIColumn(1, 1, fpEventAction->GetNumSecondaries());
  analysisManager->FillNtupleIColumn(1, 2, run->GetNumberOfEvent());
  analysisManager->FillNtupleSColumn(1, 3, kGitHash);
  analysisManager->FillNtupleIColumn(1, 4, numIntersecting);

  analysisManager->AddNtupleRow(1);

  analysisManager->Write();
  analysisManager->CloseFile();
  analysisManager->Clear();
  G4cout << "\n----> Histograms are saved" << G4endl;
  G4cout << "\n----> Number of secondaries saved " << fpEventAction->GetNumSecondaries() << G4endl;
}
