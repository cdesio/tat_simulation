#include "run.hh"
#include "G4AnalysisManager.hh"
MyRunAction::MyRunAction()
{ // if created in the construction, the output file will be created when the software launches

    G4AnalysisManager *man = G4AnalysisManager::Instance();
    // create mc ntuple
    man->CreateNtuple("mc", "mc");      // rows
    man->CreateNtupleIColumn("evt_id"); // columns
    man->CreateNtupleDColumn("posX");   // columns
    man->CreateNtupleDColumn("posY");   // columns
    man->CreateNtupleDColumn("posZ");   // columns
    man->CreateNtupleIColumn("pid");
    man->CreateNtupleSColumn("particle_name");
    man->CreateNtupleDColumn("energy");
    man->FinishNtuple(0);
    // store info - ntuple creation moved to here
    man->CreateNtuple("Hits", "Hits");      // rows
    man->CreateNtupleIColumn("evt_id");     // columns
    man->CreateNtupleDColumn("local_posX"); // columns
    man->CreateNtupleDColumn("local_posY"); // columns
    man->CreateNtupleDColumn("local_posZ"); // columns
    man->CreateNtupleIColumn("copyNo");
    man->FinishNtuple(1); // second ntuple is no. 1

    // creating ntuple to score energy deposition

    man->CreateNtuple("Scoring", "Scoring");  // rows
    man->CreateNtupleDColumn("EdepTarget");   // columns
    man->CreateNtupleDColumn("EdepDetector"); // columns
    man->CreateNtupleDColumn("Etot");         // columns
    man->FinishNtuple(2);                     // third ntuple is no. 2

    man->CreateNtuple("Detector", "Detector"); // rows
    man->CreateNtupleDColumn("Edep");          // columns
    man->CreateNtupleIColumn("copyNo");        // columns
    man->CreateNtupleDColumn("posX");          // columns
    man->CreateNtupleDColumn("posY");          // columns
    man->CreateNtupleDColumn("posZ");
    // man->CreateNtupleDColumn("posDetX"); // columns
    // man->CreateNtupleDColumn("posDetY"); // columns
    // man->CreateNtupleDColumn("posDetZ"); // columns
    man->CreateNtupleIColumn("StepID");  // columns
    man->CreateNtupleIColumn("time_ns"); // columns
    man->CreateNtupleSColumn("process"); // columns
    man->FinishNtuple(3);                // third ntuple is no. 2

    man->CreateNtuple("Target", "Taget"); // rows
    man->CreateNtupleDColumn("Edep");
    man->CreateNtupleDColumn("posX");   // columns
    man->CreateNtupleDColumn("posY");   // columns
    man->CreateNtupleDColumn("posZ");   // columns     // columns
    man->CreateNtupleIColumn("copyNo"); // columns
    man->FinishNtuple(4);               // third ntuple is no. 2
}

MyRunAction::~MyRunAction()
{
}

void MyRunAction::BeginOfRunAction(const G4Run *run)
{ // if we create the file here, it will be created new for every run
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    G4int runID = run->GetRunID();

    std::stringstream strRunID;
    strRunID << runID;

    // open file
    man->OpenFile("output" + strRunID.str() + ".root");
    // ntuple creation was here.
}

void MyRunAction::EndOfRunAction(const G4Run *)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();
}
