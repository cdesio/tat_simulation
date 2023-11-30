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
#include "G4GeneralParticleSource.hh"
#include "CommandLineParser.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4AnalysisManager.hh"

using namespace G4DNAPARSER;
using CLHEP::nanometer;

namespace
{
    G4Mutex messangerInit = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fpParticleGun(nullptr)
{

    fpParticleGun = new G4GeneralParticleSource();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    if (G4UniformRand() < 0.5)
    {
        fpParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(1.1732 * CLHEP::MeV);
    }
    else
    {
        fpParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(1.3325 * CLHEP::MeV);
    }

    fpParticleGun->GeneratePrimaryVertex(anEvent);


    // CommandLineParser *parser = CommandLineParser::GetParser();
    // Command *command(0);


    // if ((command = parser->GetCommandIfActive("-out")) != 0)
    // {

    //     G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    //     G4ThreeVector pos = fpParticleGun->GetParticlePosition();
    //     G4ThreeVector direction = fpParticleGun->GetParticleMomentumDirection();

    //     analysisManager->FillNtupleIColumn(1, 0, anEvent->GetEventID());
    //     analysisManager->FillNtupleDColumn(1, 1, fpParticleGun->GetParticleEnergy());
    //     analysisManager->FillNtupleDColumn(1, 2, pos.x());
    //     analysisManager->FillNtupleDColumn(1, 3, pos.y());
    //     analysisManager->FillNtupleDColumn(1, 4, pos.z());
    //     analysisManager->FillNtupleDColumn(1, 5, direction.x());
    //     analysisManager->FillNtupleDColumn(1, 6, direction.y());
    //     analysisManager->FillNtupleDColumn(1, 7, direction.z());
    //     analysisManager->AddNtupleRow(1);

    // }
}
