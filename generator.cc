#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
    fParticleGun = new G4ParticleGun(1); // 1 primary vertex per event
    // FROM HERE
    //  what kind of particles we want to create
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    // define particle name
    // G4String particleName = "proton";
    // get particle from table
    G4ParticleDefinition *particle = particleTable->FindParticle("geantino"); // placeholder for particle generation, that can be replaced in macro files

    // position
    G4ThreeVector pos(0., 0., 0.); // in the centre of mother volume for atm instead of 000
    // momentum
    G4ThreeVector mom(0., 0., 1.); // z direction

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleMomentum(0. * GeV);
    fParticleGun->SetParticleDefinition(particle);
    // TO HERE
    // this was in the function at the end - it was moved here to make it run at creation, and then use macros to overwrite
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
    // radioactive dacay
    G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition(); // get definition of the particle created in constructor
    // if particle is not specified as input, e.g. it's still a geantino, do this:
    if (particle == G4Geantino::Geantino())
    { // cobalt 60 Co60 decay
        G4int Z = 27;
        G4int A = 60;

        G4double charge = 0. * eplus;
        G4double energy = 0. * keV;

        G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(charge);
    }
    fParticleGun->GeneratePrimaryVertex(anEvent);
}