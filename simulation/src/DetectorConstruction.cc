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

  solidWorld = new G4Box("world", 1 * mm, 1 * mm, 1 * mm);

  solidWaterBox = new G4Orb("waterBox", 10 * um);

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

  G4double chromatin_X = fBoxSize;
  G4double chromatin_Y = fBoxSize;
  G4double chromatin_Z = 300 * nm;

  G4Box *solidTrackingVol = new G4Box("TrackingVol", dkill + chromatin_X / 2, dkill + chromatin_Y / 2, dkill + chromatin_Z / 2); // volume to track radicals in dkill larger than chromatin
  G4LogicalVolume *logicTrackingVol = new G4LogicalVolume(solidTrackingVol,
                                                          waterMaterial,
                                                          "TrackingVol");
  G4PVPlacement *physiTrackingVol = new G4PVPlacement(0,
                                                      G4ThreeVector(),
                                                      logicTrackingVol,
                                                      "TrackingVol",
                                                      logicWaterBox,
                                                      0,
                                                      false,
                                                      false);
  logicTrackingVol->SetVisAttributes(&visInvWhite);

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

  G4int copynum{0};
  while (!f.eof())
  {
    double x1, y1, z1, x2, y2, z2, base_x1, base_y1, base_z1, base_x2, base_y2, base_z2;
    f >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> base_x1 >> base_y1 >> base_z1 >> base_x2 >> base_y2 >> base_z2;
    fPositions0->Insert(G4ThreeVector(x1,
                                      y1,
                                      z1));
    fPositions1->Insert(G4ThreeVector(x2,
                                      y2,
                                      z2));
    fPositionsBase0->Insert(G4ThreeVector(base_x1,
                                          base_y1,
                                          base_z1));
    fPositionsBase1->Insert(G4ThreeVector(base_x2,
                                          base_y2,
                                          base_z2));
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
                        copynum,
                        false);

      new G4PVPlacement(0,
                        G4ThreeVector(x2,
                                      y2,
                                      z2),
                        logicSugar1,
                        "sugar1",
                        logicChromatinSegment,
                        false,
                        copynum,
                        false);
    }
    ++copynum;
  }

  // Histones
  command = CommandLineParser::GetParser()->GetCommandIfActive("-histone");

  ifstream f2(command->GetOption(), ios::in);
  if (!f2.is_open())
  {
    throw std::runtime_error("*******HISTONE GEOMETRY FILE NOT FOUND OR OPENED CORRECTY********");
  }

  while (!f2.eof())

  {
    double x1, y1, z1, rotx, roty, rotz;
    f2 >> x1 >> y1 >> z1 >> rotx >> roty >> rotz;
    fPositions.push_back(G4ThreeVector(x1,
                                       y1,
                                       z1));
    fRotations.push_back(new G4RotationMatrix(roty, // geant convention is y,z,x
                                              rotz,
                                              rotx));
  }

  if (placeHistones)
  {
    for (size_t i = 0; i < fPositions.size(); ++i)
    {
      new G4PVPlacement(fRotations[i],
                        fPositions[i],
                        logicHistone,
                        "histone",
                        logicChromatinSegment,
                        false,
                        i,
                        false);
    }
  }

  numSugar = copynum;

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