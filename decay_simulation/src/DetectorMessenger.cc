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
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det), boxSize(0), spacing(0), start_R(0), ndiv_R(0),  ndiv_theta(0), ndiv_Z(0)
{ 
  boxSize = new G4UIcmdWithADoubleAndUnit("/det/set_size",this);
  boxSize->SetGuidance("Set x,y size of the box");
  boxSize->SetParameterName("Size",false);
  boxSize->SetRange("Size>0.");
  boxSize->SetUnitCategory("Length");
  boxSize->AvailableForStates(G4State_PreInit,G4State_Idle);
  boxSize->SetToBeBroadcasted(false);
  
  spacing = new G4UIcmdWithADoubleAndUnit("/det/set_spacing",this);
  spacing->SetGuidance("Set radial spacing of boxes");
  spacing->SetParameterName("spacing",true);
  spacing->SetRange("spacing>0.");
  spacing->SetDefaultValue(2.0);
  spacing->SetDefaultUnit("um");
  spacing->AvailableForStates(G4State_PreInit,G4State_Idle);
  spacing->SetToBeBroadcasted(false);

  start_R = new G4UIcmdWithADoubleAndUnit("/det/set_startR",this);
  start_R->SetGuidance("Set radial spacing of boxes");
  start_R->SetParameterName("start_R",true);
  start_R->SetRange("Size>0.");
  start_R->SetDefaultValue(10.5);
  start_R->SetDefaultUnit("um");
  start_R->AvailableForStates(G4State_PreInit,G4State_Idle);
  start_R->SetToBeBroadcasted(false);

  ndiv_R = new G4UIcmdWithAnInteger("/det/set_ndiv_R",this);
  ndiv_R->SetGuidance("Set no. divisions in R");
  ndiv_R->SetParameterName("ndiv_R",true);
  ndiv_R->SetDefaultValue(10);
  //ndiv_R->AvailableForStates(G4State_PreInit,G4State_Idle);
  //ndiv_R->SetToBeBroadcasted(false);
  
  ndiv_Z = new G4UIcmdWithAnInteger("/det/set_ndiv_Z",this);
  ndiv_Z->SetGuidance("Set no. divisions in z");
  ndiv_Z->SetParameterName("ndiv_Z",true);
  ndiv_Z->SetDefaultValue(10);
  //ndiv_Z->AvailableForStates(G4State_PreInit,G4State_Idle);
  //ndiv_Z->SetToBeBroadcasted(false);
  
  ndiv_theta = new G4UIcmdWithAnInteger("/det/set_ndiv_theta",this);
  ndiv_theta->SetGuidance("Set no. divisions in theta");
  ndiv_theta->SetParameterName("ndiv_theta",true);
  ndiv_theta->SetDefaultValue(20);
  //ndiv_theta->AvailableForStates(G4State_PreInit,G4State_Idle);
  //ndiv_theta->SetToBeBroadcasted(false);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{

  delete boxSize; 
  delete ndiv_R;
  delete spacing;
  delete start_R;
  delete ndiv_theta;
  delete ndiv_Z;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == boxSize )
  {
     fDetector->set_size(boxSize->GetNewDoubleValue(newValue));
  }

  if( command == spacing )
  {
     fDetector->set_spacing(spacing->GetNewDoubleValue(newValue));
  }
  if( command == ndiv_R )
  {
     fDetector->set_ndiv_R(ndiv_R->GetNewIntValue(newValue));
  }
  if( command == ndiv_theta )
  {
     fDetector->set_ndiv_theta(ndiv_theta->GetNewIntValue(newValue));
  }
  if( command == ndiv_Z )
  {
     fDetector->set_ndiv_Z(ndiv_Z->GetNewIntValue(newValue));
  }
  if( command == start_R )
  {
     fDetector->set_startR(start_R->GetNewDoubleValue(newValue));
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
