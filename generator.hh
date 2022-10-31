#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
    MyPrimaryGenerator();  // constructor
    ~MyPrimaryGenerator(); // destructor

    virtual void GeneratePrimaries(G4Event *); // function that generates the primaries

private:
    G4ParticleGun *fParticleGun;
};

#endif