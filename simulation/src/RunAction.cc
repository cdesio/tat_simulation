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
#include "G4EventManager.hh"
#include "globals.hh"
#include <map>
#include <vector>
#include <numeric>
#include "CommandLineParser.hh"
#include "DetectorConstruction.hh"
#include "G4ESTARStopping.hh"
#include "G4NistManager.hh"
#include "EventAction.hh"
#include "G4EventManager.hh"
#include "git_version.hh"
#include "G4RunManager.hh"
#include "G4Version.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
using namespace G4DNAPARSER;

RunAction::RunAction(DetectorConstruction *pDetector)
    : G4UserRunAction(), fpDetector(pDetector)
{
  chromatinVolume = fpDetector->chromatinVolume;
  numSugar = fpDetector->numSugar;

  G4RunManager::GetRunManager()->SetPrintProgress(100);

  CreateNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run *)
{

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  G4String fileName{"output.root"};
  if (command->GetOption().empty() == false)
  {
    fileName = command->GetOption();
  }

  G4bool fileOpen = analysisManager->OpenFile(fileName);
  if (!fileOpen)
  {
    G4cout << "\n---> HistoManager::book(): cannot open " << fileName << G4endl;
    return;
  }
  G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run *)
{
  if (isMaster)
  {
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->Write();
    analysisManager->CloseFile();
    return;
  }
  auto fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
  std::vector<G4double> KE = fpEventAction->fTrackMeanKE;
  G4double meanKE = accumulate(KE.begin(), KE.end(), 0.0) / KE.size();

  // ICRUU90 used for alpha/proton stopping power, ESTAR used for electron
  const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
      G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String primaryName = generatorAction->primaryName;

  if (primaryName == "alpha")
  {
    auto fICRU90 = G4NistManager::Instance()->GetICRU90StoppingData();
    fICRU90->Initialise();
    auto material = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    LET = fICRU90->GetElectronicDEDXforAlpha(material, meanKE) * material->GetDensity();
  }
  else if (primaryName == "proton")
  {
    auto fICRU90 = G4NistManager::Instance()->GetICRU90StoppingData();
    fICRU90->Initialise();
    auto material = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    LET = fICRU90->GetElectronicDEDXforProton(material, meanKE) * material->GetDensity();
  }
  else if (primaryName == "e-")
  {
    auto fESTAR = new G4ESTARStopping();
    auto material = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    LET = fESTAR->GetElectronicDEDX(material, meanKE) * material->GetDensity();
  }
  WriteNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::CreateNtuple()
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetNtupleMerging(true);
  // open output file
  //

  analysisManager->SetFirstNtupleId(0);
  analysisManager->CreateNtuple("EventEdep", "EventEdep");
  analysisManager->CreateNtupleDColumn("Edep_J");
  analysisManager->CreateNtupleDColumn("Edep_MeV");
  analysisManager->CreateNtupleIColumn("EventNo");
  if ((command = parser->GetCommandIfActive("-PS")) == 0)
  {
    analysisManager->CreateNtupleDColumn("PrimaryKEEntrance");
    analysisManager->CreateNtupleDColumn("PrimaryKEExit");
    analysisManager->CreateNtupleDColumn("ProjectedRangeChromatin");
    analysisManager->CreateNtupleDColumn("PathLengthChromatin");
    analysisManager->CreateNtupleIColumn("StoppedInTarget");
  }

  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("Direct", "Direct");
  analysisManager->CreateNtupleIColumn(1, "EventNo");
  analysisManager->CreateNtupleDColumn(1, "eDep_eV");
  analysisManager->CreateNtupleDColumn(1, "x");
  analysisManager->CreateNtupleDColumn(1, "y");
  analysisManager->CreateNtupleDColumn(1, "z");
  analysisManager->CreateNtupleSColumn(1, "Particle");
  analysisManager->CreateNtupleDColumn(1, "KE");
  analysisManager->CreateNtupleIColumn(1, "copyNo");
  analysisManager->CreateNtupleDColumn(1, "time");
  analysisManager->FinishNtuple(1);

  // For chemistry

  analysisManager->CreateNtuple("Indirect", "Indirect");

  analysisManager->CreateNtupleIColumn(2, "EventNo");
  analysisManager->CreateNtupleDColumn(2, "x");
  analysisManager->CreateNtupleDColumn(2, "y");
  analysisManager->CreateNtupleDColumn(2, "z");
  analysisManager->CreateNtupleSColumn(2, "DNAmolecule");
  analysisManager->CreateNtupleSColumn(2, "radical");

  analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("Info", "Info");
  analysisManager->CreateNtupleDColumn("ChromatinVolume_m3");
  analysisManager->CreateNtupleDColumn("NumBasepairs");
  analysisManager->CreateNtupleSColumn("GitHash");
  analysisManager->CreateNtupleSColumn("G4Version");
  analysisManager->CreateNtupleSColumn("IRTmodel");
  if ((command = parser->GetCommandIfActive("-PS")) == 0)
  {
    analysisManager->CreateNtupleDColumn("MeanLET");
  }
  analysisManager->FinishNtuple(3);

  if ((command = parser->GetCommandIfActive("-PS")))
  {
    analysisManager->CreateNtuple("Events", "Events");

    analysisManager->CreateNtupleIColumn(4, "EventNo");
    analysisManager->CreateNtupleIColumn(4, "PhotonEventID");
    analysisManager->CreateNtupleIColumn(4, "copyNo");
    analysisManager->CreateNtupleIColumn(4, "particleID");

    analysisManager->FinishNtuple(4);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void RunAction::WriteNtuple()
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillNtupleDColumn(3, 0, chromatinVolume);
  analysisManager->FillNtupleDColumn(3, 1, numSugar);
  analysisManager->FillNtupleSColumn(3, 2, kGitHash);
  analysisManager->FillNtupleSColumn(3, 3, G4Version);
  analysisManager->FillNtupleSColumn(3, 4, "IRT");
  if ((command = parser->GetCommandIfActive("-PS")) == 0)
  {
    analysisManager->FillNtupleDColumn(3, 5, LET);
  }

  analysisManager->AddNtupleRow(3);

  analysisManager->Write();
  analysisManager->CloseFile();
  G4cout << "\n----> Histograms are saved" << G4endl;
}
