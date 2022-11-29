#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
    gpsParticleGun = new G4GeneralParticleSource();
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete gpsParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
    G4int Z = 85;
    G4int A = 211;

    G4int eventNum = anEvent->GetEventID();

    G4double charge = 0. * eplus;
    G4double energy = 0. * keV;
    G4ThreeVector pos(0., 0., 0.);
    G4ThreeVector mom(0., 0., 1.);

    G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
    gpsParticleGun->SetParticleDefinition(ion);
    gpsParticleGun->SetParticleCharge(charge);

    // gpsParticleGun->SetParticlePosition(pos);
    // gpsParticleGun->SetParticleEnergy(energy);
    // gpsParticleGun->SetParticleMomentum(mom);

    gpsParticleGun->GeneratePrimaryVertex(anEvent);
}