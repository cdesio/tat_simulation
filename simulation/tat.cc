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

#include "G4Types.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4ParallelWorldPhysics.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "ParallelWorld.hh"
#include "PhysicsList.hh"
#include "CommandLineParser.hh"
#include "G4DNAChemistryManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using namespace G4DNAPARSER;
CommandLineParser *parser(0);

void Parse(int &argc, char **argv);

int main(int argc, char **argv)
{
  Parse(argc, argv);

  Command *commandLine(0);

  int mySeed{1};
  if ((commandLine = parser->GetCommandIfActive("-seed")))
  {
    mySeed = strtol(commandLine->GetOption(), NULL, 10);
  }

  G4Random::setTheSeed(mySeed);
  G4Random::showEngineStatus();

  //  read in PS file for e- from photon simulation
  std::ifstream ps_file;
  G4String ps_file_name;
  std::vector<std::vector<float>> PS_data;

  if ((commandLine = parser->GetCommandIfActive("-PS")))
  {
    ps_file_name = commandLine->GetOption();

    ps_file.open(ps_file_name); // open file

    if (ps_file.is_open())
    {
      G4cout << "Single file opened OK: " << ps_file_name << G4endl;
    }
    else
    {
      G4cout << "*******FILE NOT FOUND OR OPENED CORRECTY********" << G4endl;
      return 1;
    }

    float line[12];

    // int i = 0;
    while (ps_file.read((char *)&line, sizeof line))
    {
      PS_data.push_back({line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11]});

      // ++i;
    }
  }
  std::unique_ptr<G4RunManager> pRunManager(G4RunManagerFactory::CreateRunManager());
  pRunManager->SetNumberOfThreads(1); // by default

  DetectorConstruction *pDetector = new DetectorConstruction();
  pDetector->RegisterParallelWorld(new ParallelWorld("ChemistryWorld"));
  pRunManager->SetUserInitialization(pDetector);

  PhysicsList *pPhysList = new PhysicsList;
  pRunManager->SetUserInitialization(pPhysList);
  pRunManager->SetUserInitialization(new ActionInitialization(pDetector, PS_data));
  // if (CommandLineParser::GetParser()->GetCommandIfActive("-chemOFF") == 0)
  // {
  //   G4DNAChemistryManager::Instance()->Initialize();
  // }

  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  if ((commandLine = parser->GetCommandIfActive("-mac")))
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + commandLine->GetOption());
  }
  // else
  // {
  //   UImanager->ApplyCommand("/control/execute rbe.in");
  // }
  if ((commandLine = parser->GetCommandIfActive("-PS")))
  {
    pRunManager->Initialize();
    G4cout << "number of beams = " << PS_data.size() << G4endl;

    pRunManager->BeamOn(PS_data.size());
  }
  if ((commandLine = parser->GetCommandIfActive("-gui")))
  {
    // initialize visualization
    G4UIExecutive *ui = new G4UIExecutive(argc, argv);
    G4VisManager *visManager = nullptr;

    visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute vis.mac");

    ui->SessionStart();

    delete ui;

    delete visManager;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GetNameAndPathOfExecutable(char **argv,
                                G4String &executable,
                                G4String &path)
{
  // Get the last position of '/'
  std::string aux(argv[0]);

  // get '/' or '\\' depending on unix/mac or windows.
#if defined(_WIN32) || defined(WIN32)
  int pos = aux.rfind('\\');
#else
  int pos = aux.rfind('/');
#endif

  // Get the path and the name
  path = aux.substr(0, pos + 1);
  executable = aux.substr(pos + 1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Parse(int &argc, char **argv)
{
  //////////
  // Parse options given in commandLine
  //
  parser = CommandLineParser::GetParser();

  parser->AddCommand("-gui",
                     Command::WithoutOption,
                     "Select geant4 UI or just launch a geant4 terminal session",
                     "qt");

  parser->AddCommand("-mac",
                     Command::WithOption,
                     "Give a mac file to execute",
                     "macFile.mac");

  // You cann your own command, as for instance:
  parser->AddCommand("-seed",
                     Command::WithOption,
                     "Give a seed value in argument to be tested", "seed");
  // it is then up to you to manage this option

  G4String exec;
  G4String path;
  GetNameAndPathOfExecutable(argv, exec, path);

  parser->AddCommand("-out",
                     Command::OptionNotCompulsory,
                     "Output files (ROOT is used by default)",
                     exec);
  parser->AddCommand("-chemOFF",
                     Command::WithoutOption,
                     "Deactivate chemistry");
  parser->AddCommand("-sugar",
                     Command::WithOption,
                     "Deoxyribose position file",
                     exec);
  parser->AddCommand("-histone",
                     Command::WithOption,
                     "Histone position file",
                     exec);
  parser->AddCommand("-PS",
                     Command::WithOption,
                     "PS file for simulaiton");
  //////////
  // If -h or --help is given in option : print help and exit
  //
  if (parser->Parse(argc, argv) != 0) // help is being printed
  {
    // if you are using ROOT, create a TApplication in this condition in order
    // to print the help from ROOT as well
    CommandLineParser::DeleteInstance();
    std::exit(0);
  }

  ///////////
  // Kill application if wrong argument in command line
  //
  if (parser->CheckIfNotHandledOptionsExists(argc, argv))
  {
    // if you are using ROOT, you should initialise your TApplication
    // before this condition
    abort();
  }
}