#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"

#include "run.hh"

class MyEventAction : public G4UserEventAction
{
public:
    MyEventAction(MyRunAction *);
    ~MyEventAction();

    virtual void BeginOfEventAction(const G4Event *);
    virtual void EndOfEventAction(const G4Event *);
    // function to accumulate edep

    void AddEdep(G4int iVol, G4double edep);

private:
    // create value to store edep
    G4double fEdep1, fEdep2;
};

#endif