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
    void set_size     (G4double);  
    void set_spacing (G4double);
    void set_startR(G4double);
    void set_vessellength(G4double);
    void set_ndiv_R (G4int);
    void set_ndiv_theta(G4int);
    void set_ndiv_Z(G4int);



    G4double chromatinVolume;



private:
         DetectorMessenger* fDetectorMessenger;
         G4double              boxSize;
         G4double spacing;
         G4double start_R;
         G4double vessel_length;

         G4int ndiv_R;
         G4int ndiv_theta;
         G4int ndiv_Z;
};
