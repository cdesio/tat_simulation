#include "construction.hh"
#include "parameterisation.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
    // change the number of photosensors in one column
    fMessenger->DeclareProperty("nCols", nCols, "Number of columns");
    // change the number of photosensors in one row
    fMessenger->DeclareProperty("nRows", nRows, "Number of rows");
    // add the possibility to choose between the two types of detector: cherenkov or scintillator
    fMessenger->DeclareProperty("isCherenkov", isCherenkov, "Construct Cherenkov ");
    fMessenger->DeclareProperty("isScintillator", isScintillator, "Construct Scintillator");
    fMessenger->DeclareProperty("isTat", isTat, "Construct TaT");
    fMessenger->DeclareProperty("isAtmosphere", isAtmosphere, "Construct Atmosphere");

    fMessenger->DeclareProperty("n_div_Theta", n_div_Theta, "Number of divisions in Theta");
    fMessenger->DeclareProperty("n_div_Z", n_div_Z, "Number of divisions in Cylinder Axis (Z)");
    fMessenger->DeclareProperty("n_div_R", n_div_R, "Number of divisions in radius");

    fMessenger->DeclareProperty("total_length", total_length, "Total length of Cylinder");
    fMessenger->DeclareProperty("inner_radius", inner_radius, "Inner radius of Cylinder");
    fMessenger->DeclareProperty("outer_radius", outer_radius, "Outer radius of Cylinder");

    DefineMaterials();

    // initial values
    nCols = 100;
    nRows = 100;

    isCherenkov = false;
    isScintillator = false;
    isTat = true;
    isAtmosphere = false;

    if (isTat)
    {
        xWorld = .15 * mm;
        yWorld = .15 * mm;
        zWorld = .15 * mm;

        n_div_Theta = 6;
        n_div_Z = 2;
        n_div_R = 4;
        total_length = 40 * um;
        inner_radius = 10 * um;
        outer_radius = 110 * um;
    }
    else if (isAtmosphere /* condition */)
    {
        xWorld = 40 * km;
        yWorld = 40 * km;
        zWorld = 20 * km;
    }
    else
    {
        xWorld = .5 * m;
        yWorld = .5 * m;
        zWorld = .5 * m;
    }
}

MyDetectorConstruction::~MyDetectorConstruction()
{
}

void MyDetectorConstruction::DefineMaterials()
{
    // create instance of nist to use already definer materials
    G4NistManager *nist = G4NistManager::Instance();

    // define aerogelm as new material: name, density, no. of components
    SiO2 = new G4Material("SiO2", 2.201 * g / cm3, 2);
    SiO2->AddElement(nist->FindOrBuildElement("Si"), 1);
    SiO2->AddElement(nist->FindOrBuildElement("O"), 2);

    // Define water
    H2O = new G4Material("H2O", 1.000 * g / cm3, 2);
    H2O->AddElement(nist->FindOrBuildElement("H"), 2);
    H2O->AddElement(nist->FindOrBuildElement("O"), 1);

    // Define scaled water
    ScaledH2O = new G4Material("ScaledH2O", 1.060 * g / cm3, 2);
    ScaledH2O->AddElement(nist->FindOrBuildElement("H"), 2);
    ScaledH2O->AddElement(nist->FindOrBuildElement("O"), 1);

    // carbon
    C = nist->FindOrBuildElement("C");

    // Aerogel from defined elements
    Aerogel = new G4Material("Aerogel", 0.200 * g / cm3, 3);
    Aerogel->AddMaterial(SiO2, 62.5 * perCent);
    Aerogel->AddMaterial(H2O, 37.4 * perCent);
    Aerogel->AddElement(C, 0.1 * perCent);

    // add refractive index of the material for aerogel

    G4double energy[2] = {1.239841939 * eV / 0.9, 1.239841939 * eV / 0.2};
    G4double rindexAerogel[2] = {1.1, 1.1};
    G4double rindexWorld[2] = {1.0, 1.0};
    G4double rindexNaI[2] = {1.78, 1.78};
    G4double reflectivity[2] = {1.0, 1.0};

    G4MaterialPropertiesTable *mptAerogel = new G4MaterialPropertiesTable();
    mptAerogel->AddProperty("RINDEX", energy, rindexAerogel, 2);

    // define material for world
    worldMat = nist->FindOrBuildMaterial("G4_AIR");
    // SAME FOR WORLD - refractive index
    G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
    mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);

    Aerogel->SetMaterialPropertiesTable(mptAerogel);

    worldMat->SetMaterialPropertiesTable(mptWorld);

    Na = nist->FindOrBuildElement("Na");
    I = nist->FindOrBuildElement("I");
    NaI = new G4Material("NaI", 3.67 * g / cm3, 2);
    NaI->AddElement(Na, 1);
    NaI->AddElement(I, 1);

    G4double fraction[2] = {1.0, 1.0}; // same fraction of photons for blue and red wavelength
    // for scintillation light - add refractive index
    G4MaterialPropertiesTable *mptNaI = new G4MaterialPropertiesTable();
    mptNaI->AddProperty("RINDEX", energy, rindexNaI, 2);
    // how many photons to create per wavelength
    mptNaI->AddProperty("SCINTILLATIONCOMPONENT1", energy, fraction, 2);
    // how many photons per energy loss of the particles
    mptNaI->AddConstProperty("SCINTILLATIONYIELD", 38. / keV); // no fluctuation, no array
    // other needed property
    mptNaI->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // decay time of scintillator, that emits particles according to exponential function
    mptNaI->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 250 * ns);
    // related to distribution of the emitted photons
    mptNaI->AddConstProperty("SCINTILLATIONYIELD1", 1.);

    NaI->SetMaterialPropertiesTable(mptNaI);

    mirrorSurface = new G4OpticalSurface("mirrorSurface");
    mirrorSurface->SetType(dielectric_metal);
    mirrorSurface->SetFinish(ground);
    mirrorSurface->SetModel(unified);

    G4MaterialPropertiesTable *mptMirror = new G4MaterialPropertiesTable();
    mptMirror->AddProperty("REFLECTIVITY", energy, reflectivity, 2);

    mirrorSurface->SetMaterialPropertiesTable(mptMirror);

    G4double density0 = 1.29 * kg / m3;
    G4double aN = 14.01 * g / mole;
    G4double aO = 16.00 * g / mole;

    N = new G4Element("Nitrogen", "N", 7, aN);
    O = new G4Element("Oxygen", "O", 8, aO);

    G4double f = 3;
    G4double R = 8.3144626181532;
    G4double g0 = 9.81;
    G4double kappa = (f + 2) / f;
    G4double T = 293.15;
    G4double M = (0.3 * aO + 0.7 * aN) / 1000.;

    for (G4int i = 0; i < 10; i++)
    {
        std::stringstream stri;
        stri << i;
        G4double h = 40e3 / 10. * i;
        G4double density = density0 * pow((1 - (kappa) / kappa * M * g0 * h / (R * T)), (1 / (kappa - 1)));
        Air[i] = new G4Material("G4_AIR" + stri.str(), density, 2);
        Air[i]->AddElement(N, 70 * perCent);
        Air[i]->AddElement(O, 30 * perCent);
    }
}

void MyDetectorConstruction::ConstructCherenkov()
{
    // Create detector
    solidRadiator = new G4Box("solidRadiator", 0.4 * m, 0.4 * m, 0.01 * m);
    // logical volume
    logicRadiator = new G4LogicalVolume(solidRadiator, Aerogel, "logicRadiator");

    // use radiator as scoring volume
    fScoringVolume = logicRadiator;

    // physical volume
    physRadiator = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.25 * m), logicRadiator, "physRadiator", logicWorld, false, 0, true);

    // define 100 photosensors
    // solid
    solidDetector = new G4Box("solidDetector", xWorld / nRows, yWorld / nCols, 0.01 * m);
    // logical
    logicDetector = new G4LogicalVolume(solidDetector, worldMat, "logicalDetector");
    // physical - in this case we have to create an array, so it's different:
    for (G4int i = 0; i < nRows; i++)
    {
        for (G4int j = 0; j < nCols; j++)
        {
            physDetector = new G4PVPlacement(0, G4ThreeVector(-0.5 * m + (i + 0.5) * m / nRows, -0.5 * m + (j + 0.5) * m / nCols, 0.49 * m), logicDetector, "physDetector", logicWorld, false, j + i * nCols, true);
        }
    }
}

void MyDetectorConstruction::ConstructScintillator()
{

    solidScintillator = new G4Box("solidScintillator", 5 * cm, 5 * cm, 6 * cm);

    logicScintillator = new G4LogicalVolume(solidScintillator, NaI, "logicalScintillator");

    G4LogicalSkinSurface *Skin = new G4LogicalSkinSurface("skin", logicWorld, mirrorSurface); // add reflective coating to world volume

    // add photon detectors on top of scintillators
    solidDetector = new G4Box("solidDetector", 1 * cm, 5 * cm, 6 * cm);

    logicDetector = new G4LogicalVolume(solidDetector, worldMat, "logicDetector");

    // use scintillator as scoring volume
    fScoringVolume = logicScintillator;

    for (G4int i = 0; i < 6; i++)
    {
        for (G4int j = 0; j < 16; j++)
        {
            G4Rotate3D rotZ(j * 22.5 * deg, G4ThreeVector(0, 0, 1));
            G4Translate3D transXScint(G4ThreeVector(5. / tan(22.5 / 2 * deg) * cm + 5. * cm, 0 * cm, -40 * cm + i * 15 * cm));

            G4Translate3D transXDet(G4ThreeVector(5. / tan(22.5 / 2 * deg) * cm + 6. * cm + 5. * cm, 0 * cm, -40 * cm + i * 15 * cm));

            G4Transform3D transformScint = (rotZ) * (transXScint); // first translation then rotation
            G4Transform3D transformDet = (rotZ) * (transXDet);     // first translation then rotation

            physScintillator = new G4PVPlacement(transformScint, logicScintillator, "physScintillator", logicWorld, false, 0, true);
            physDetector = new G4PVPlacement(transformDet, logicDetector, "physDetector", logicWorld, false, 0, true);
        }
    }
}

void MyDetectorConstruction::ConstructTatDetector()
{

    solidTarget = new G4Tubs("solidTarget", 0, 10 * um, 20 * um, 0 * deg, 360 * deg);
    logicTarget = new G4LogicalVolume(solidTarget, H2O, "logicTarget");
    fScoringVolume = logicTarget;
    physTarget = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicTarget, "physTarget", logicWorld, false, 0, true);
    G4LogicalSkinSurface *Skin = new G4LogicalSkinSurface("skin", logicWorld, mirrorSurface); // add reflective coating to world volume

    G4int n_j = 4;
    G4double angle = 360. / n_j;
    G4int n_k = 2;
    G4double length = 200 * um;
    G4int n_i = 2;
    G4double l_n = length / n_i;

    // solidTatDetector = new G4Tubs("solidTatDetector", 10 * um, 110 * um, l_n / 2, 0 * deg, angle * deg);
    // logicDetector = new G4LogicalVolume(solidTatDetector, H2O, "logicDetector");
    //  fScoringVolume = logicDetector;

    for (G4int i = 0; i < n_i; i++)
    {
        for (G4int j = 0; j < n_j; j++)
        {
            for (G4int k = 0; k < n_k; k++)
            {

                solidTatDetector = new G4Tubs("solidScintillator2", 10 * um + (k * 100 / n_k) * um, 10 * um + (100 / n_k) * (k + 1) * um, l_n / 2, 0 * deg, angle * deg);
                logicDetector = new G4LogicalVolume(solidTatDetector, H2O, "logicDetector");
                //  fScoringVolume = logicDetector;

                G4Rotate3D rotZTat(j * angle * deg, G4ThreeVector(0, 0, 1));
                G4Translate3D transXTatDetector(G4ThreeVector(0., 0., ((-length / 2) + (i + 1) * l_n - l_n / 2)));
                G4Transform3D transformTat = (rotZTat) * (transXTatDetector);

                physTatDetector = new G4PVPlacement(transformTat, logicDetector, "physTatDetector", logicWorld, false, i * n_j + j * n_k + k, true);
            }
            //}
        }
    }

    G4VisAttributes *targetVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 255.0));
    targetVisAtt->SetVisibility(true);
    logicTarget->SetVisAttributes(targetVisAtt);
}

void MyDetectorConstruction::ConstructTatDetector2()
{

    G4int n_voxels = n_div_R * n_div_Theta * n_div_Z;

    // G4cout << n_voxels << G4endl;

    solidTatDetector = new G4Tubs("solidTatDetector", inner_radius, outer_radius, total_length, 0 * deg, 360 * deg);
    logicTatDetector = new G4LogicalVolume(solidTatDetector, ScaledH2O, "logicTatDetector");
    physTatDetector = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicTatDetector, "physTatDetector", logicWorld, false, 0, true);

    solidCyilinderEl = new G4Tubs("solidCylinderEl", 0, 10 * um, 10 * um, 0 * deg, 360 * deg);
    logicDetector = new G4LogicalVolume(solidCyilinderEl, ScaledH2O, "logicDetector");
    fScoringDetector = logicDetector;

    G4VPVParameterisation *cylinderParam = new CylinderParameterisation(total_length,
                                                                        inner_radius,
                                                                        outer_radius,
                                                                        n_div_R,
                                                                        n_div_Z,
                                                                        n_div_Theta);

    G4VPhysicalVolume *physCylinder = new G4PVParameterised("Cylinder", logicDetector, logicTatDetector, kZAxis, n_voxels, cylinderParam, true);

    solidTarget = new G4Tubs("solidTarget", 0, inner_radius, 40 * um, 0 * deg, 360 * deg);
    logicTarget = new G4LogicalVolume(solidTarget, ScaledH2O, "logicTarget");
    fScoringTarget = logicTarget;
    physTarget = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicTarget, "physTarget", logicWorld, false, 0, true);

    G4VisAttributes *targetVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 255.0, .5));
    targetVisAtt->SetVisibility(true);
    targetVisAtt->SetForceSolid(true);
    logicTarget->SetVisAttributes(targetVisAtt);

    G4VisAttributes *tatDetVisAtt = new G4VisAttributes(G4Colour(1, 1, 1, 0.2));
    tatDetVisAtt->SetVisibility(true);
    tatDetVisAtt->SetForceSolid(true);
    logicTatDetector->SetVisAttributes(tatDetVisAtt);

    G4VisAttributes *voxelsVisAtt = new G4VisAttributes(G4Colour(255., 255, 255., .5));
    voxelsVisAtt->SetVisibility(true);
    voxelsVisAtt->SetLineWidth(2);
    logicDetector->SetVisAttributes(voxelsVisAtt);
    logicDetector->SetUserLimits(new G4UserLimits(.5 * um));
    logicTarget->SetUserLimits(new G4UserLimits(.5 * um));
}

void MyDetectorConstruction::ConstructAtmosphere()
{
    solidAtmosphere = new G4Box("solidAtmosphere", xWorld, yWorld, zWorld / 10.);

    for (G4int i = 0; i < 10; i++)
    {
        logicAtmosphere[i] = new G4LogicalVolume(solidAtmosphere, Air[i], "logicAtmosphere");
        physAtmosphere[i] = new G4PVPlacement(0, G4ThreeVector(0, 0, zWorld / 10. * 2 * i - zWorld + zWorld / 10.), logicAtmosphere[i], "physAtmosphere", logicWorld, false, i, true);
    }
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{

    // define solid (first step)
    solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);

    // define logical volume (Second step)
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");

    // define physical volume: volume placement
    //(0: no rotation, x, y, z: (0,0,0)placed at axes origin), logic volume, name, mother volume: 0, boolean operations: false, copy number (to insert volume several times), true: check for overlaps)
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

    if (isCherenkov)
        ConstructCherenkov();

    if (isScintillator)
        ConstructScintillator();

    if (isTat)
    {
        ConstructTatDetector2();
    }
    if (isAtmosphere)
    {
        ConstructAtmosphere();
    }

    //  return mother volume
    return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");
    if (logicDetector != NULL)
        logicDetector->SetSensitiveDetector(sensDet);

    MySensitiveDetector *sensTarget = new MySensitiveDetector("SensitiveTarget");
    if (logicTarget != NULL)
        logicTarget->SetSensitiveDetector(sensTarget);
}
