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
#pragma once
#include "G4UserSteppingAction.hh"
#include "G4String.hh"
#include <fstream>
#include <iostream>
#include <map>
#include "RunAction.hh"
#include "G4Track.hh"
class EventAction;
// class G4ParticleDefinition;
// class G4VPhysicalVolume;
// class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(DetectorConstruction* pDetector);
    ~SteppingAction() override;

    void UserSteppingAction(const G4Step* step) override;
    // void Initialize();
private:
    EventAction* fpEventAction;
    DetectorConstruction *fpDetector;
    std::ofstream PSfile;
    void write_to_PS(const G4Track *track, const G4ThreeVector &Position, const G4ThreeVector &Momentum, G4double &particleEnergy, G4int &eventID, G4int &PID, G4int &copyNo, G4double &time, G4int &mapped_PID);
    G4ThreeVector transformDirection(const G4ThreeVector & worldPos, const G4ThreeVector & worldMomentum);
    G4double calculateDistanceToExitBox(const G4ThreeVector & preStepPosition, const G4ThreeVector & preStepMomentumDirection);
  

    std::map<G4String, G4int> particleMap{
        {"Pb212", 0},
        {"alpha", 1},
        {"gamma", 2},
        {"e-", 3},
        {"e+", 4},
        {"nu_e", 5},
        {"Bi212", 6},
        {"Tl208", 7},
        {"Po212", 8},
        {"Pb208", 9}};
    
    std::map<G4String, G4int> particleOriginMap{
      {"Pb212", 0},
      {"Bi212", 1},
      {"Tl208", 2},
      {"Po212", 3},
      {"Pb208", 4},
      {"e-Pb212", 5},
      {"e+Pb212", 6},
      {"alphaBi212", 7},
      {"e-Bi212", 8},
      {"e+Bi212", 9},
      {"e-Tl208", 10},
      {"e+Tl208", 11},
      {"alphaPo212", 12}};

  std::map<G4int, G4String> reverseParticleOriginMap{
    {0, "Pb212"},
      {1, "Bi212"},
      {2, "Tl208"},
      {3, "Po212"},
      {4, "Pb208"},
      {5, "e-Pb212"},
      {6, "e+Pb212"},
      {7, "alphaBi212"},
      {8, "e-Bi212"},
      {9, "e+Bi212"},
      {10, "e-Tl208"},
      {11, "e+Tl208"},
      {12, "alphaPo212"}};

};
