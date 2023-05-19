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
#include "G4SystemOfUnits.hh"
// #include "CLHEP/Units/SystemOfUnits.h"
#include "G4AnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

using namespace G4DNAPARSER;
using CLHEP::nanometer;

namespace
{
  G4Mutex messangerInit = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(std::vector<std::vector<float>> PS_data)
    : G4VUserPrimaryGeneratorAction(), fpParticleGun(nullptr), fPS_data(PS_data)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-PS")))
  {
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    fParticleGun->SetParticleEnergy(0);
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  }
  else
  {
    fpParticleGun = new G4GeneralParticleSource();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((command = parser->GetCommandIfActive("-PS")))
  {
    G4int step2_eventID = anEvent->GetEventID();

    G4double positionX = fPS_data[step2_eventID][0];
    G4double positionY = fPS_data[step2_eventID][1];
    G4double positionZ = fPS_data[step2_eventID][2];
    G4double momentumX = fPS_data[step2_eventID][3];
    G4double momentumY = fPS_data[step2_eventID][4];
    G4double momentumZ = fPS_data[step2_eventID][5];
    G4double particleEnergy = fPS_data[step2_eventID][6];
    step1_eventID = (int)fPS_data[step2_eventID][7];
    step1_PID = (int)fPS_data[step2_eventID][8];
    step1_copyNo = (int)fPS_data[step2_eventID][9];
    step1_time = fPS_data[step2_eventID][10];
    step1_processID = (int)fPS_data[step2_eventID][11];

    // // DEBUG

    //G4cout << "step1 Evt: " << step1_eventID << ") PID: " << step1_PID << ". processID: " << step1_processID << " copyNo: " << step1_copyNo << " time: "<< step1_time << G4endl;
    // // DEBUG
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4IonTable *ionTable = G4IonTable::GetIonTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("geantino");

    if (step1_PID == 1)
    {
      particle = particleTable->FindParticle("alpha");
    }
    else if (step1_PID == 2)
    {
      particle = particleTable->FindParticle("gamma");
    }
    else if (step1_PID == 3)
    {
      particle = particleTable->FindParticle("e-");
    }
    else if (step1_PID == 4)
      return;
    //      particle = particleTable->FindParticle("geantino");

    else if (step1_PID == 5)
    {
      // particle = particleTable->FindParticle("alpha");
      //  particle = particleTable->FindParticle("At211");
      particle = ionTable->GetIon(85, 211, 0);
      // return;
    }
    else if (step1_PID == 6)
    {
      particle = ionTable->GetIon(84, 211, 0);
    }
    else if (step1_PID == 7)
    {
      particle = ionTable->GetIon(83, 207, 0);
    }
    else if (step1_PID == 8)
    {
      particle = ionTable->GetIon(82, 207, 0);
    }
    else
    {
      G4cout << "DEBUG: unknown particle. Throwing geantino." << G4endl;
      particle = particleTable->FindParticle("geantino");
    }
    // else
    // {
    // particle = particleTable->FindParticle(pdg_enc);
    // }
    // G4cout << "particle: " << particle->GetParticleName() << G4endl;
    fParticleGun->SetParticleDefinition(particle);
    G4String particleName = particle->GetParticleName();
    primaryName = particleName;

    G4cout << "shooting a(n) " << particleName << G4endl;
    fParticleGun->SetParticlePosition(G4ThreeVector(positionX, positionY, positionZ));
    fParticleGun->SetParticleEnergy(particleEnergy);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momentumX, momentumY, momentumZ));
    //primaryName = fpParticleGun->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    if ((command = parser->GetCommandIfActive("-out")) != 0)
    {

      G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

      analysisManager->FillNtupleIColumn(4, 0, step2_eventID);
      analysisManager->FillNtupleIColumn(4, 1, step1_eventID);
      analysisManager->FillNtupleIColumn(4, 2, step1_copyNo);
      analysisManager->FillNtupleIColumn(4, 3, step1_PID);
      analysisManager->FillNtupleIColumn(4, 4, step1_processID);
      analysisManager->AddNtupleRow(4);
    }
  }
  else
  {
    fpParticleGun->GeneratePrimaryVertex(anEvent);
    primaryName = fpParticleGun->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
  }
}