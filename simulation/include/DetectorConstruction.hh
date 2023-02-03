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

#pragma once
#include "G4VUserDetectorConstruction.hh"
#include <memory>
#include "G4RotationMatrix.hh"
#include "G4KDTree.hh"

class G4VPhysicalVolume;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorConstruction
    : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    G4VPhysicalVolume *Construct() override;
    void UpdateGeometry();

    void SetSize     (G4double);              
    void SetDkill(G4double val) {dkill = val;}     
    G4double Getdkill() {return dkill;} 

    std::vector<G4ThreeVector> fPositions;
    std::vector<G4RotationMatrix *> fRotations;
    G4KDTree* fPositions0;
    G4KDTree* fPositions1;    
    G4KDTree* fPositionsBase0;
    G4KDTree* fPositionsBase1;
    G4double chromatinVolume;
    G4double numSugar;

private:
     G4double              fBoxSize;
     DetectorMessenger* fDetectorMessenger;
    G4bool placeSugars;
    G4bool placeHistones;
    G4double dkill{9*CLHEP::nm};
};
