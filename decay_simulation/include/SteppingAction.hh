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
class EventAction;
// class G4ParticleDefinition;
// class G4VPhysicalVolume;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction *pDetector);
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step *step) override;
  // void Initialize();
private:
  EventAction* fpEventAction;
  DetectorConstruction *fpDetector;
  std::ofstream PSfile;
  G4int boxes_per_r;

  std::map<G4String, G4int> particleMap{
      {"alpha", 1},
      {"gamma", 2},
      {"e-", 3},
      {"nu_e", 4},
      {"At211", 5},
      {"Po211", 6},
      {"Bi207", 7},
      {"Pb207", 8},
      {"e+", 9}
};

  std::map<G4String, G4int> processMap
  {
        {"RadioactiveDecay", 1},
            {"ionIoni", 2},
            {"msc", 3},
            {"eIoni", 4},
            {"phot", 5},
            {"eBrem", 6},
            {"compt", 7}
  }; 

};
