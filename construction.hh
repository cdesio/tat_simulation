#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "detector.hh"
// messenger
#include "G4GenericMessenger.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4VisAttributes.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include <vector>

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    // create scoring volume to store info only for that
    G4LogicalVolume *GetScoringDetector() const { return fScoringDetector; }
    G4LogicalVolume *GetScoringTarget() const { return fScoringTarget; }

    virtual G4VPhysicalVolume *Construct();
    // to access the logical volume from outside the construction function
private:
    G4LogicalVolume *logicDetector;

    // move this declaration here, so it's not created every time it's modified
    G4Tubs *solidTarget, *solidTatDetector, *solidCyilinderEl;
    G4Box *solidWorld, *solidRadiator, *solidDetector, *solidScintillator, *solidAtmosphere;
    G4LogicalVolume *logicWorld, *logicRadiator, *logicScintillator, *logicTarget, *logicAtmosphere[10], *logicTatDetector;
    G4VPhysicalVolume *physWorld, *physRadiator, *physDetector, *physScintillator, *physTatDetector, *physTarget, *physAtmosphere[10], *physCylinder;

    G4Material *SiO2,
        *H2O, *Aerogel, *worldMat, *NaI, *Air[10], *ScaledH2O;
    G4Element *C, *Na, *I, *N, *O;

    // let's put here the definition of the materials
    void
    DefineMaterials();

    // function to construct cherenkov detector
    void ConstructCherenkov();
    // function to construct scintillator cylinder detector
    void ConstructScintillator();

    void ConstructTatDetector();
    void ConstructTatDetector2();

    void ConstructAtmosphere();

    virtual void ConstructSDandField();

    G4GenericMessenger *fMessenger;
    G4LogicalVolume *fScoringDetector, *fScoringTarget, *fScoringVolume;

    G4int nCols, nRows;
    G4double xWorld, yWorld, zWorld;
    G4bool isCherenkov, isScintillator, isTat, isAtmosphere;
    G4OpticalSurface *mirrorSurface;
    G4int n_div_Theta, n_div_Z, n_div_R;
    G4double total_length, inner_radius, outer_radius;
};

#endif