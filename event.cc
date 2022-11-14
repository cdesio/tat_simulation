#include "event.hh"
#include "G4AnalysisManager.hh"

MyEventAction::MyEventAction(MyRunAction *)
{
    // start value
    fEdep = 0.;
}

MyEventAction::~MyEventAction()
{
}

void MyEventAction::BeginOfEventAction(const G4Event *)
{ // whenever a new event starts the accumulated edep has to be reset to 0.
    fEdep = 0.;
}

void MyEventAction::EndOfEventAction(const G4Event *)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleDColumn(2, 0, fEdep);
    man->AddNtupleRow(2);
}

void MyEventAction::AddEdep(G4double edep)
{
    fEdep += edep;
};
