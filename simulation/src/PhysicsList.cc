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

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmAndDNA.hh"
#include "G4EmDNAChemistry_option3.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "CommandLineParser.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace G4DNAPARSER;

PhysicsList::PhysicsList()
    : G4VModularPhysicsList(), fDNAPhysicsList(nullptr), fChemistryList_option3(nullptr)
{
    SetDefaultCutValue(1.0 * nanometer);
    SetVerboseLevel(1);
    RegisterConstructor("G4EmAndDNA");
    if (CommandLineParser::GetParser()->GetCommandIfActive("-chemOFF") == 0)
    {
        RegisterConstructor("G4EmDNAChemistry_option3");
    }

    // Decay
    RegisterPhysics(new G4DecayPhysics());

    // Radioactive decay
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(5 * eV, 1 * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PhysicsList::ConstructParticle()
{
    if (fDNAPhysicsList != nullptr)
    {
        fDNAPhysicsList->ConstructParticle();
    }

    if (fChemistryList_option3 != nullptr)
    {
        fChemistryList_option3->ConstructParticle();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
    AddTransportation();

    if (fDNAPhysicsList != nullptr)
    {
        fDNAPhysicsList->ConstructProcess();
    }
    if (fChemistryList_option3 != nullptr)
    {
        fChemistryList_option3->ConstructProcess();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RegisterConstructor(const G4String &name)
{
    if (name == fPhysDNAName)
    {
        return;
    }

    if (verboseLevel > 0)
    {
        G4cout << "===== Register constructor ==== " << name << G4endl;
    }

    if (name == "G4EmAndDNA")
    {
        fDNAPhysicsList.reset(new G4EmAndDNA(verboseLevel));
        fPhysDNAName = name;
    }
    else if (name == "G4EmDNAChemistry_option3")
    {
        fChemistryList_option3.reset(new G4EmDNAChemistry_option3());
        fChemistryList_option3->SetVerboseLevel(verboseLevel);
    }
    else
    {
        G4cout << "PhysicsList::RegisterConstructor: <" << name << ">"
               << " fails - name is not defined"
               << G4endl;
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......