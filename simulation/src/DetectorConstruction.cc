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
#include "G4GenericMessenger.hh"
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
#include "CommandLineParser.hh"

#include "Randomize.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"

#include "G4ProductionCuts.hh"

using namespace G4DNAPARSER;

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::angstrom;
using CLHEP::degree;
using CLHEP::micrometer;
using CLHEP::mm;
using CLHEP::nanometer;

static G4VisAttributes visInvBlue(false, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visInvWhite(false, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visInvPink(false, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visInvCyan(false, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visInvRed(false, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visInvGreen(false, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visBlue(true, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visWhite(true, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visPink(true, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visCyan(true, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visRed(true, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visGreen(true, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visGrey(true, G4Colour(0.839216, 0.839216, 0.839216));

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  fBoxSize = 300 * nm;

  displ_X = 0 * um;
  displ_Y = 0 * um;
  displ_X = 0 * um;

  fDetectorMessenger = new DetectorMessenger(this);

  placeSugars = false; // only needed for visualisation
  placeHistones = true;

  fPositions0 = new G4KDTree();
  fPositions1 = new G4KDTree();
  fPositionsBase0 = new G4KDTree();
  fPositionsBase1 = new G4KDTree();
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
  G4Tubs *solidWaterCylinder{nullptr};
  G4Box *solidPrism{nullptr};

  G4double len_y = 39 * um;
  G4double prism_base = 3.5 * um;

  solidWorld = new G4Box("world", 10 * mm, 10 * mm, 10 * mm);

  solidWaterBox = new G4Orb("waterBox", 1 * mm);

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

  G4PVPlacement *physiWaterBox = new G4PVPlacement(0,
                                                   G4ThreeVector(),
                                                   logicWaterBox,
                                                   "waterBox",
                                                   logicWorld,
                                                   0,
                                                   false,
                                                   0);
  logicWorld->SetVisAttributes(&visWhite);
  logicWaterBox->SetVisAttributes(&visWhite);

  solidWaterCylinder = new G4Tubs("waterCylinder", 0, 100 * um, 3.5 * um, 0, 360);
  G4LogicalVolume *logicWaterCylinder = new G4LogicalVolume(solidWaterCylinder, waterMaterial, "waterCylinder");
  G4PVPlacement *physiWaterCyl = new G4PVPlacement(0, G4ThreeVector(), logicWaterCylinder, "waterCylinder", logicWorld, 0, false, 0);
  logicWaterCylinder->SetVisAttributes(&visWhite);

  solidPrism = new G4Box("solidPrism", prism_base / 2, len_y, prism_base / 2);
  G4LogicalVolume *logicPrism = new G4LogicalVolume(solidPrism, waterMaterial, "logicPrism");

  for (G4int j = 0; j < 16; j++)
  {
    G4Rotate3D rotZ(j * 360 / 16 * deg, G4ThreeVector(0, 0, 1));
    G4Translate3D transPrism(G4ThreeVector(0., (len_y + 10 * um), 0.));
    G4Transform3D transformPrism = (rotZ) * (transPrism); // first translation then rotation

    G4PVPlacement *physPrism = new G4PVPlacement(transformPrism, logicPrism, "physPrism", logicWaterCylinder, false, j, true);
  }

  G4double chromatin_X = fBoxSize;
  G4double chromatin_Y = fBoxSize;
  G4double chromatin_Z = 300 * nm;

  G4Box *solidTrackingVol = new G4Box("TrackingVol", dkill + chromatin_X / 2, dkill + chromatin_Y / 2, dkill + chromatin_Z / 2); // volume to track radicals in dkill larger than chromatin
  G4LogicalVolume *logicTrackingVol = new G4LogicalVolume(solidTrackingVol,
                                                          waterMaterial,
                                                          "TrackingVol");

  G4Rotate3D rotTrackingVol(90 * deg, G4ThreeVector(1, 0, 0));
  G4Translate3D transTrackingVol(G4ThreeVector(displ_X, displ_Y, displ_Z));
  G4Transform3D transformTrackingVol = (rotTrackingVol) * (transTrackingVol);

  G4PVPlacement *physiTrackingVol = new G4PVPlacement(transformTrackingVol,
                                                      logicTrackingVol,
                                                      "TrackingVol",
                                                      logicWaterBox,
                                                      0,
                                                      false,
                                                      false);
  logicTrackingVol->SetVisAttributes(&visWhite);

  G4Box *solidChromatinSegment = new G4Box("chromatinSegment", chromatin_X / 2, chromatin_Y / 2, chromatin_Z / 2);
  G4LogicalVolume *logicChromatinSegment = new G4LogicalVolume(solidChromatinSegment,
                                                               waterMaterial,
                                                               "chromatinSegment");

  chromatinVolume = chromatin_X / m * chromatin_Y / m * chromatin_Z / m; // in m3

  G4PVPlacement *physiChromatinSegment = new G4PVPlacement(0,
                                                           G4ThreeVector(),
                                                           logicChromatinSegment,
                                                           "chromatinSegment",
                                                           logicTrackingVol,
                                                           0,
                                                           0,
                                                           false);
  logicChromatinSegment->SetVisAttributes(&visPink);

  G4double sugarRadius{2.9 * angstrom};

  G4Orb *solidSugar = new G4Orb("sugar", sugarRadius);
  G4LogicalVolume *logicSugar0 = new G4LogicalVolume(solidSugar,
                                                     waterMaterial,
                                                     "sugar0");

  G4LogicalVolume *logicSugar1 = new G4LogicalVolume(solidSugar,
                                                     waterMaterial,
                                                     "sugar1");

  G4Tubs *solidHistone = new G4Tubs("histone",
                                    0,
                                    2.9 * nanometer,
                                    2.75 * nanometer,
                                    0 * degree,
                                    360 * degree);

  G4LogicalVolume *logicHistone = new G4LogicalVolume(solidHistone,
                                                      waterMaterial,
                                                      "histone");

  Command *command = CommandLineParser::GetParser()->GetCommandIfActive("-sugar");

  ifstream f(command->GetOption(), ios::in); // open file

  if (!f.is_open())
  {
    throw std::runtime_error("*******SUGAR GEOMETRY FILE NOT FOUND OR OPENED CORRECTY********");
  }

  G4int copynum_sugar{0};
  while (!f.eof())
  {
    double x1, y1, z1, x2, y2, z2, base_x1, base_y1, base_z1, base_x2, base_y2, base_z2;
    f >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> base_x1 >> base_y1 >> base_z1 >> base_x2 >> base_y2 >> base_z2;

    G4ThreeVector *pos0 = new G4ThreeVector(x1 + displ_X,
                                            y1 + displ_Y,
                                            z1 + displ_Z);
    // G4cout << "before: " << pos0->getX() << " " << pos0->getY() << " " << pos0->getZ() << G4endl;
    pos0->rotate(90 * deg, G4ThreeVector(1, 0, 0));
    fPositions0->Insert(pos0);

    G4ThreeVector *pos1 = new G4ThreeVector(x2 + displ_X,
                                            y2 + displ_Y,
                                            z2 + displ_Z);
    pos1->rotate(90 * deg, G4ThreeVector(1, 0, 0));
    fPositions1->Insert(pos1);

    G4ThreeVector *basepos0 = new G4ThreeVector(base_x2 + displ_X,
                                                base_y2 + displ_Y,
                                                base_z2 + displ_Z);
    basepos0->rotate(90 * deg, G4ThreeVector(1, 0, 0));
    fPositionsBase0->Insert(basepos0);

    G4ThreeVector *basepos1 = new G4ThreeVector(base_x2 + displ_X,
                                                base_y2 + displ_Y,
                                                base_z2 + displ_Z);
    basepos1->rotate(90 * deg, G4ThreeVector(1, 0, 0));
    fPositionsBase1->Insert(basepos1);

    if (placeSugars)
    {

      new G4PVPlacement(0,
                        G4ThreeVector(x1,
                                      y1,
                                      z1),
                        logicSugar0,
                        "sugar0",
                        logicChromatinSegment,
                        false,
                        copynum_sugar,
                        false);

      new G4PVPlacement(0,
                        G4ThreeVector(x2,
                                      y2,
                                      z2),
                        logicSugar1,
                        "sugar1",
                        logicChromatinSegment,
                        false,
                        copynum_sugar,
                        false);
    }
    ++copynum_sugar;
  }

  // Histones
  command = CommandLineParser::GetParser()->GetCommandIfActive("-histone");

  ifstream f2(command->GetOption(), ios::in);
  if (!f2.is_open())
  {
    throw std::runtime_error("*******HISTONE GEOMETRY FILE NOT FOUND OR OPENED CORRECTY********");
  }
  G4int copynum_histone{0};
  while (!f2.eof())

  {
    double x1, y1, z1, rotx, roty, rotz;
    f2 >> x1 >> y1 >> z1 >> rotx >> roty >> rotz;
    fPositions.push_back(G4ThreeVector(x1 + displ_X,
                                       y1 + displ_Y,
                                       z1 + displ_Z));
    fRotations.push_back(new G4RotationMatrix(roty, // geant convention is y,z,x
                                              rotz,
                                              rotx));

    if (placeHistones)
    {
      // for (size_t i = 0; i < fPositions.size(); ++i)
      //{
      new G4PVPlacement(new G4RotationMatrix(rotx, roty, rotz),
                        G4ThreeVector(x1, y1, z1),
                        logicHistone,
                        "histone",
                        logicChromatinSegment,
                        false,
                        copynum_histone,
                        false);
      // }
    }
    ++copynum_histone;
  }

  numSugar = copynum_sugar;

  logicSugar0->SetVisAttributes(&visInvBlue);
  logicSugar1->SetVisAttributes(&visInvRed);
  logicHistone->SetVisAttributes(&visGrey);

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
void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetDisplX(G4double value)
{
  displ_X = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetDisplY(G4double value)
{
  displ_Y = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetDisplZ(G4double value)
{
  displ_Z = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}