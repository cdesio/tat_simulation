#include "event.hh"
#include "G4AnalysisManager.hh"

MyEventAction::MyEventAction(MyRunAction *)
{
    // start value
    fEdep1 = fEdep2 = 0.;
}

MyEventAction::~MyEventAction()
{
}

void MyEventAction::BeginOfEventAction(const G4Event *)
{ // whenever a new event starts the accumulated edep has to be reset to 0.
    fEdep1 = fEdep2 = 0.;
}

void MyEventAction::EndOfEventAction(const G4Event *)
{
    G4double Etot = fEdep1 + fEdep2;
    // if (fEdep > 0.)
    // {
    //     G4cout << "Energy deposition: " << fEdep << G4endl;
    // }
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleDColumn(2, 0, fEdep1);
    man->FillNtupleDColumn(2, 1, fEdep2);
    man->FillNtupleDColumn(2, 2, Etot);
    man->AddNtupleRow(2);
}

void MyEventAction::AddEdep(G4int iVol, G4double edep)
{
    if (iVol == 1)
    {
        fEdep1 += edep;
    }
    if (iVol == 2)
    {
        fEdep2 += edep;
    }
};
