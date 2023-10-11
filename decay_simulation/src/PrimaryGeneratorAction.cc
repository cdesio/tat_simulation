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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "CommandLineParser.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4AnalysisManager.hh"
#include <math.h>
#include <G4SystemOfUnits.hh>


// using namespace G4DNAPARSER;
// using CLHEP::nanometer;

/*
Using GPS method, but adapted so that the intersection between the initial ray and the tracking cylinder
can be determined. All initially generated positions and directions are recorded but only those which
intersect are added as a primary vertex.
*/

// namespace
// {
//     G4Mutex messangerInit = G4MUTEX_INITIALIZER;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fpParticleGun(nullptr)
{

    // fpParticleGun = new G4GeneralParticleSource();
    G4int n_particle = 1;
    fpParticleGun = new G4ParticleGun(n_particle);

    fpParticleGun->SetParticleEnergy(0*eV);
    fpParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fpParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    G4double vessel_length = 20*um; 
    G4double vessel_radius = 10*um;
    G4double z_part = (G4UniformRand() - 0.5) * vessel_length;
    G4double angle_part = G4UniformRand() * 2 * M_PI;
    G4double x_part = (vessel_radius)*cos(angle_part);
    G4double y_part = (vessel_radius)*sin(angle_part);

    fpParticleGun->SetParticlePosition(G4ThreeVector(x_part,y_part,z_part)); 
      
    fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(cos(angle_part),sin(angle_part),0));
    fpParticleGun->GeneratePrimaryVertex(anEvent);
}
