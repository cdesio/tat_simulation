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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "globals.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include <G4SystemOfUnits.hh>
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"

#include <fstream>
#include <math.h>
#include "CommandLineParser.hh"

#include "Randomize.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4ProductionCuts.hh"
#include "G4RunManager.hh"

using namespace G4DNAPARSER;

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
// using CLHEP::angstrom;
// using CLHEP::degree;
// using CLHEP::micrometer;
// using CLHEP::mm;
// using CLHEP::nanometer;

static G4VisAttributes visInvBlue(false, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visInvWhite(false, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visInvPink(false, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visInvCyan(false, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visInvRed(false, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visInvGreen(false, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visBlue(true, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visWhite(true, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visPink(true, G4Colour(0.9,  0.6,  0.7, 0.5));
static G4VisAttributes visCyan(true, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visRed(true, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visGreen(true, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visGrey(true, G4Colour(0.839216, 0.839216, 0.839216));
static G4VisAttributes visDarkRed(true, G4Colour(0.8, 0.0, 0.3, 0.9));
static G4VisAttributes visLime(true, G4Colour(0.7, 1.0, 0.0, 0.5));
static G4VisAttributes visBlack(true, G4Colour(0,0,0));

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  boxSize = 300 * nm;
  // start_R = 10.5;
  // spacing = 2.0;
  // ndiv_R = 20;
  // ndiv_theta = 80;
  // ndiv_Z = 10;
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{

  /***************************************************************************/
  //                               World
  /***************************************************************************/

  G4NistManager *man = G4NistManager::Instance();
  G4Material *waterMaterial = man->FindOrBuildMaterial("G4_WATER");
  G4Material *air = man->FindOrBuildMaterial("G4_AIR");

  G4Box *solidWorld{nullptr};
  G4Orb *solidWaterBox{nullptr};
  G4Orb *solidVessel{nullptr};

  // larger world volume to provide build up for CPE for photon beam
  solidWorld = new G4Box("world", 0.21 * mm, 0.21 * mm, 0.21 * mm);

  solidWaterBox = new G4Orb("waterBox", .2 * mm);
  solidVessel = new G4Orb("bloodVessel", 10 * um);

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,
                                                    air,
                                                    "world");

  G4PVPlacement *physiWorld = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                "world",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  G4LogicalVolume *logicWaterBox = new G4LogicalVolume(solidWaterBox,
                                                       waterMaterial,
                                                       "waterBox");
  G4LogicalVolume *logicBloodVessel = new G4LogicalVolume(solidVessel,
                                                          waterMaterial, "bloodVessel");

  G4PVPlacement *physiWaterBox = new G4PVPlacement(0,
                                                   G4ThreeVector(),
                                                   logicWaterBox,
                                                   "waterBox",
                                                   logicWorld,
                                                   0,
                                                   false,
                                                   0);

  G4PVPlacement *physiBloodVessel = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                      logicBloodVessel, "bloodVessel", logicWaterBox, 0, false, 0);

  G4VisAttributes *vesselVisAttr = new G4VisAttributes(G4Colour(0.8, 0.0, 0.4, 0.9));//(G4Colour(0.83, 0.83, 0.83, 0.5));

  vesselVisAttr->SetForceSolid(true);
  vesselVisAttr->SetVisibility(true);
  logicBloodVessel->SetVisAttributes(vesselVisAttr);

  logicWorld->SetVisAttributes(&visWhite);
  logicWaterBox->SetVisAttributes(&visInvWhite);
  

  G4double delta = 10 * nm;

  G4Box *solidTrackingVol = new G4Box("TrackingVol", delta * nm + boxSize / 2, delta * nm + boxSize / 2, delta * nm + boxSize / 2); // volume to track radicals in 9nm larger than chromatin
  G4LogicalVolume *logicTrackingVol = new G4LogicalVolume(solidTrackingVol,
                                                          waterMaterial,
                                                          "TrackingVol");
 
  
  int ibox = 0;
  int i_r = 0;
  G4double shell_rin =0;
  G4double shell_rfin = 0;
  G4cout << "pi: " << M_PI << G4endl;
  for (G4int r = 0; r < ndiv_R; r += 1)
  {
  // G4cout << "start: " << start_R+ r*spacing + boxSize/2 << ", stop: " << start_R + (r+1)* spacing - boxSize/2 << G4endl;
  // G4cout << "start_R: " << start_R << ", r * spacing: " << r* spacing << ", box/2: " << boxSize/2 << G4endl;
  if(r==0){
    //shell_rin = 10*um; //use this to make shells that go from vessel to start of first line
    shell_rin = start_R-boxSize/2; //use this to make shells that are as big as the voxels 
  }
  else{
    //shell_rin = start_R+(r-1)*spacing+boxSize; //use this to make shells that go end of the previous line to the start of the next one
    shell_rin = start_R+(r*spacing)-boxSize/2; //use this to make shells that are as big as the voxels
  }
  //shell_rfin = (start_R+(r)*spacing)-boxSize/2; //use this to make shells that go end of the previous line to the start of the next one
  shell_rfin = shell_rin+boxSize; //use this to make shells that are as big as the voxels
  G4Sphere *solidShell = new G4Sphere("shell", shell_rin, shell_rfin, 0,360*degree, 0, 360 * degree);
  
  G4LogicalVolume *logicShell = new G4LogicalVolume(solidShell,
                                                         waterMaterial,
                                                         "shell");
  //logicShell->SetUserLimits(fStepLimit);
  G4VisAttributes *vesselVisAttr = new G4VisAttributes(G4Colour(0.8, 0.0, 0.4, 0.9));//(G4Colour(0.83, 0.83, 0.83, 0.5));

  logicShell->SetVisAttributes(&visRed);
  G4PVPlacement *physiCell = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     logicShell,
                                                     "shell",
                                                     logicWaterBox,
                                                     0,
                                                     r,
                                                     0);
    

    // ndiv_theta = (int)(std::ceil(2 * M_PI * (start_R / um + r * spacing / um) / 0.5));
    // long double th_incr = M_PI / (ndiv_theta / 2);
    // // G4cout << "th_incr: " << th_incr << G4endl;
    // // G4cout << "DEBUG: ndiv_theta: " << ndiv_theta << ", boxes_per_R: " << ndiv_theta*ndiv_Z << ", start_R/um: " << start_R/um <<", spacing: " << spacing/um << ", r: " << r << G4endl;
    
    // for (G4int k = 0; k < ndiv_Z; k++)
    // {
    //   G4double theta = 0;
    //  // for (long double theta = 0; theta < 2 * M_PI; theta += th_incr)
    //   for (G4int th = 0; th < ndiv_theta; th++)
    //   {
    //     G4RotationMatrix *rot = new G4RotationMatrix(theta,
    //                                                  0,
    //                                                  0);
    //     // G4cout << "n_div_R: " << ndiv_R << ", ndiv_Z: " << ndiv_Z << ", ndiv_Theta: " << ndiv_theta << ", start_R: " << start_R << ", spacing: " << spacing << G4endl;

    //     G4PVPlacement *physiTrackingVol = new G4PVPlacement(rot, G4ThreeVector((((start_R + (r * spacing))) * cos(theta)), ((start_R + (r * spacing))) * sin(theta), ((-vessel_length / um + (k + 1) * (vessel_length / um * 2) / ndiv_Z) - 0.3) * um),
    //                                                         logicTrackingVol,
    //                                                         "TrackingVol",
    //                                                         logicWaterBox,
    //                                                         false,
    //                                                         ibox,
    //                                                         false);
    //     //G4cout << "boxes placed: " << ibox << ", r: " << r << ", z: " << k << ", th: " << theta <<  G4endl;
    //     ibox++;
    //     theta+= th_incr;
          
    //     }
      
    // }
    // i_r++;
    
    }
  // G4cout << "DEBUG: vessel length/ndiv_Z: " << vessel_length / ndiv_Z << "boxes: " << ibox << G4endl;
  // logicTrackingVol->SetVisAttributes(&visGrey);

  chromatinVolume = boxSize / m * boxSize / m * boxSize / m; // in m3

  G4Region *aRegion = new G4Region("Target");

  G4ProductionCuts *cuts = new G4ProductionCuts();

  G4double defCut = 1 * nanometer;
  cuts->SetProductionCut(defCut, "gamma");
  cuts->SetProductionCut(defCut, "e-");
  cuts->SetProductionCut(defCut, "e+");
  cuts->SetProductionCut(defCut, "proton");

  aRegion->SetProductionCuts(cuts);

  aRegion->AddRootLogicalVolume(logicTrackingVol); // all volumes placed within are included in the region

  return physiWorld;
}
void DetectorConstruction::set_size(G4double value)
{
  boxSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::set_ndiv_R(G4int value)
{
  ndiv_R = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::set_ndiv_theta(G4int value)
{
  ndiv_theta = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::set_ndiv_Z(G4int value)
{
  ndiv_Z = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::set_spacing(G4double value)
{
  spacing = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::set_startR(G4double value)
{
  start_R = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::set_vessellength(G4double value)
{
  vessel_length = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
