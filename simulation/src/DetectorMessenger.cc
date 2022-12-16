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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction *Det)
    : G4UImessenger(), fDetector(Det), fSizeCmd(0), fdkillCmd(0), fDisplCmdX(0), fDisplCmdY(0), fDisplCmdZ(0)
{
  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/det/setSize", this);
  fSizeCmd->SetGuidance("Set x,y size of the box");
  fSizeCmd->SetParameterName("Size", false);
  fSizeCmd->SetRange("Size>0.");
  fSizeCmd->SetUnitCategory("Length");
  fSizeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeCmd->SetToBeBroadcasted(false);

  fdkillCmd = new G4UIcmdWithADoubleAndUnit("/det/dkill", this);
  fdkillCmd->SetGuidance("Set dkill for chemistry");
  fdkillCmd->SetParameterName("dkill", false);
  fdkillCmd->SetRange("dkill>0.");
  fdkillCmd->SetUnitCategory("Length");
  fdkillCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fdkillCmd->SetToBeBroadcasted(false);

  fDisplCmdX = new G4UIcmdWithADoubleAndUnit("/det/displ_X", this);
  fDisplCmdX->SetGuidance("Set displ_X value for DNA box");
  fDisplCmdX->SetParameterName("displ_X", true);
  fDisplCmdX->SetDefaultValue(0);
  fDisplCmdX->SetUnitCategory("Length");
  fDisplCmdX->SetDefaultUnit("um");
  fDisplCmdX->AvailableForStates(G4State_PreInit, G4State_Idle);
  fDisplCmdX->SetToBeBroadcasted(false);

  fDisplCmdY = new G4UIcmdWithADoubleAndUnit("/det/displ_Y", this);
  fDisplCmdY->SetGuidance("Set displ_Y value for DNA box");
  fDisplCmdY->SetParameterName("displ_Y", true);
  fDisplCmdY->SetDefaultValue(0);
  fDisplCmdY->SetUnitCategory("Length");
  fDisplCmdY->SetDefaultUnit("um");
  fDisplCmdY->AvailableForStates(G4State_PreInit, G4State_Idle);
  fDisplCmdY->SetToBeBroadcasted(false);

  fDisplCmdZ = new G4UIcmdWithADoubleAndUnit("/det/displ_Z", this);
  fDisplCmdZ->SetGuidance("Set displ_Z value for DNA box");
  fDisplCmdZ->SetParameterName("displ_Z", true);
  fDisplCmdZ->SetDefaultValue(0);
  fDisplCmdZ->SetUnitCategory("Length");
  fDisplCmdZ->SetDefaultUnit("um");
  fDisplCmdZ->AvailableForStates(G4State_PreInit, G4State_Idle);
  fDisplCmdZ->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{

  delete fSizeCmd;
  delete fdkillCmd;
  delete fDisplCmdX;
  delete fDisplCmdY;
  delete fDisplCmdZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{

  if (command == fSizeCmd)
  {
    fDetector->SetSize(fSizeCmd->GetNewDoubleValue(newValue));
  }
  if (command == fdkillCmd)
  {
    fDetector->SetDkill(fdkillCmd->GetNewDoubleValue(newValue));
  }

  if (command == fDisplCmdX)
  {
    fDetector->SetDisplX(fDisplCmdX->GetNewDoubleValue(newValue));
  }
  if (command == fDisplCmdY)
  {
    fDetector->SetDisplY(fDisplCmdY->GetNewDoubleValue(newValue));
  }
  if (command == fDisplCmdZ)
  {
    fDetector->SetDisplZ(fDisplCmdZ->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
