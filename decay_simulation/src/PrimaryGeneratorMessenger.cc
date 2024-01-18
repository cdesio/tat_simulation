
#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "globals.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                                   PrimaryGeneratorAction* Gun)
:G4UImessenger(),fAction(Gun), fGunDir(0),  gun_length(0)
{ 
  fGunDir = new G4UIdirectory("/shell/gun/");
  fGunDir->SetGuidance("gun parameters");

  gun_length = new G4UIcmdWithADoubleAndUnit("/shell/gun/gun_length",this);
  gun_length->SetGuidance("set generator half length");
  gun_length->SetParameterName("gun_length",false);
  gun_length->SetDefaultValue(1.5);
  gun_length->SetDefaultUnit("micrometer");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete gun_length;
  delete fGunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{ 
  if (command == gun_length)
   {fAction->SetGunLength(gun_length->GetNewDoubleValue(newValue));}   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

