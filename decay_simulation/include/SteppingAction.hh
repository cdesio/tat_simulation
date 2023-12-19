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
    void write_to_PS(const G4Track *track, const G4ThreeVector &globalPosition, const G4ThreeVector &localPosition, const G4ThreeVector &Momentum, G4double &particleEnergy, G4int &eventID, G4int &PID, G4int &copyNo, G4double &time, G4int &mapped_PID);
    G4ThreeVector transformDirection(const G4ThreeVector & worldPos, const G4ThreeVector & worldMomentum);
    G4double calculateDistanceToExitBox(const G4ThreeVector & preStepPosition, const G4ThreeVector & preStepMomentumDirection);
  

    std::map<G4String, G4int> particleMap{
        {"At211", 0},
        {"alpha", 1},
        {"gamma", 2},
        {"e-", 3},
        {"nu_e", 4},
        {"Po211", 5},
        {"Po211*", 6},
        {"Bi207", 7},
        {"Pb207", 8},
        {"Pb207*", 9}};
    
    std::map<G4String, G4int> particleOriginMap{
      {"At211", 0},
      {"Po211", 1},
      {"Po211*", 2},
      {"Bi207", 3},
      {"Pb207", 4},
      {"Pb207*", 5},
      {"alphaAt211", 6},
      {"alphaPo211", 7},
      {"e-At211", 8},
      {"e-Bi207", 9},
      {"e-Pb207*", 10},
      {"gammaAt211", 11},
      {"gammaBi207", 12},
      {"gammaPb207", 13},
      {"gammaPb207*", 14},
      {"gammaPo211", 15},
      {"gammaPo211*", 16}};

  std::map<G4int, G4String> reverseParticleOriginMap{
      {0, "At211"},
      {1, "Po211"},
      {2, "Po211*"},
      {3, "Bi207"},
      {4, "Pb207"},
      {5, "Pb207*"},
      {6, "alphaAt211"},
      {7, "alphaPo211"},
      {8, "e-At211"},
      {9, "e-Bi207"},
      {10, "e-Pb207*"},
      {11, "gammaAt211"},
      {12, "gammaBi207"},
      {13, "gammaPb207"},
      {14, "gammaPb207*"},
      {15, "gammaPo211"}, 
      {16, "gammaPo2011*"}};

};
