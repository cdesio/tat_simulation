// Combination of standard physics option 4 world and DNA physics option 2 in the Target

#include "G4EmAndDNA.hh"

#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4VEmModel.hh"
#include "G4DummyModel.hh"
#include "G4PhysicsListHelper.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4ProcessManager.hh"

// Geant4 standard EM MODELS
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4UrbanMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4Generator2BS.hh"
#include "G4NuclearStopping.hh"
#include "G4ePairProduction.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4eplusAnnihilation.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"

// Geant4-DNA MODELS
#include "G4DNAElastic.hh"
#include "G4DNAExcitation.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNABornExcitationModel.hh"
#include "G4DNAIonisation.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNARuddIonisationModel.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"
#include "G4DNAChargeIncrease.hh"
#include "G4DNADingfelderChargeIncreaseModel.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAMeltonAttachmentModel.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAOneStepThermalizationModel.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAVibExcitation.hh"
// #include "G4DNACPA100ExcitationModel.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmAndDNA);

G4EmAndDNA::G4EmAndDNA(G4int ver, const G4String &)
    : G4VPhysicsConstructor("G4EmAndDNA"), verbose(ver)
{
  // defaultCutValue = 1*micrometer;
  SetVerboseLevel(ver);

    G4EmParameters* param = G4EmParameters::Instance();

  param->SetDefaults();
  param->SetFluo(true);  
  param->SetAuger(true);  
  param->SetDeexcitationIgnoreCut(true);
  param->ActivateDNA();

  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmAndDNA::~G4EmAndDNA()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmAndDNA::ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();

  //  baryons
  G4Proton::ProtonDefinition();
  G4GenericIon::GenericIonDefinition();

  // Geant4 DNA new particles
  G4DNAGenericIonsManager *genericIonsManager;
  genericIonsManager = G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmAndDNA::ConstructProcess()
{

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  G4EmParameters *param = G4EmParameters::Instance();

  // nuclear stopping is enabled if the energy limit above zero
  G4double nielEnergyLimit = param->MaxNIELEnergy();
  G4NuclearStopping *pnuc = nullptr;
  if (nielEnergyLimit > 0.0)
  {
    pnuc = new G4NuclearStopping();
    pnuc->SetMaxKinEnergy(nielEnergyLimit);
  }

  // high energy limit for e+- scattering models and bremsstrahlung
  G4double highEnergyLimit = param->MscEnergyLimit();

  while ((*particleIterator)())
  {

    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();

    // *********************************
    // 1) Processes for the World region
    // *********************************

    if (particleName == "e-")
    {

      // STANDARD msc is active in the world
      // multiple scattering
      G4eMultipleScattering *msc = new G4eMultipleScattering();
      // e-/e+ msc gs with Mott-correction
      // (Mott-correction is set through G4EmParameters)
      G4GoudsmitSaundersonMscModel *msc1 = new G4GoudsmitSaundersonMscModel();
      G4WentzelVIModel *msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->SetEmModel(msc1);
      msc->SetEmModel(msc2);

      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4eIonisation *eion = new G4eIonisation();
      eion->SetEmModel(new G4MollerBhabhaModel());
      ph->RegisterProcess(eion, particle);

      // STANDARD Coulomb scattering is active in the world
      G4eCoulombScatteringModel *ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering *ss = new G4CoulombScattering();
      ss->SetEmModel(ssm);
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);

      // STANDARD bremsstrahlung is active in the world
      G4eBremsstrahlung *brem = new G4eBremsstrahlung();
      G4SeltzerBergerModel *br1 = new G4SeltzerBergerModel();
      G4eBremsstrahlungRelModel *br2 = new G4eBremsstrahlungRelModel();
      br1->SetAngularDistribution(new G4Generator2BS());
      br2->SetAngularDistribution(new G4Generator2BS());
      brem->SetEmModel(br1);
      brem->SetEmModel(br2);
      br1->SetHighEnergyLimit(CLHEP::GeV);

      ph->RegisterProcess(brem, particle);

      // DNA elastic is not active in the world
      G4DNAElastic *theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(theDNAElasticProcess);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("e-_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("e-_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA e- solvation
      G4DNAElectronSolvation *solvation =
          new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
      solvation->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(solvation);

      // vibration excitation
      G4DNAVibExcitation *theDNAeVibExcProcess =
          new G4DNAVibExcitation("e-_G4DNAVibExcitation");
      theDNAeVibExcProcess->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(theDNAeVibExcProcess);

      // attachment
      G4DNAAttachment *theDNAAttachProcess =
          new G4DNAAttachment("e-_G4DNAAttachment");
      theDNAAttachProcess->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(theDNAAttachProcess);
    }
    else if (particleName == "proton")
    {

      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4WentzelVIModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD Coulomb scattering is active in the world
      G4CoulombScattering *pcou = new G4CoulombScattering();
      ph->RegisterProcess(pcou, particle);

      // STANDARD ionisation is active in the world
      G4hIonisation *hion = new G4hIonisation();
      hion->SetEmModel(new G4BraggModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("proton_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("proton_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("proton_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge decrease is not active in the world
      G4DNAChargeDecrease *dnacd = new G4DNAChargeDecrease("proton_G4DNAChargeDecrease");
      dnacd->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnacd);
    }
    else if (particleName == "hydrogen")
    {

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("hydrogen_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("hydrogen_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("hydrogen_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge increase is not active in the world
      G4DNAChargeIncrease *dnaci = new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease");
      dnaci->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaci);
    }
    else if (particleName == "GenericIon")
    {

      // WARNING : THIS IS NEEDED FOR STANDARD ALPHA G4ionIonisation PROCESS

      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4UrbanMscModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4ionIonisation *hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);
    }
    else if (particleName == "alpha")
    {

      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4UrbanMscModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4ionIonisation *hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("alpha_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("alpha_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("alpha_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge decrease is not active in the world
      G4DNAChargeDecrease *dnacd = new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease");
      dnacd->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnacd);
    }
    else if (particleName == "alpha+")
    {

      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4UrbanMscModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4ionIonisation *hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("alpha+_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("alpha+_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("alpha+_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge decrease is not active in the world
      G4DNAChargeDecrease *dnacd = new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease");
      dnacd->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnacd);

      // DNA charge increase is not active in the world
      G4DNAChargeIncrease *dnaci = new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease");
      dnaci->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaci);
    }
    else if (particleName == "helium")
    {

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("helium_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("helium_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("helium_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge increase is not active in the world
      G4DNAChargeIncrease *dnaci = new G4DNAChargeIncrease("helium_G4DNAChargeIncrease");
      dnaci->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaci);
    }
    // Warning : the following particles and processes are needed by EM Physics builders
    // They are taken from the default Livermore Physics list
    // These particles are currently not handled by Geant4-DNA

    // e+

    else if (particleName == "e+")
    {
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
    }
    else if (particleName == "gamma")
    {

      // photoelectric effect - Livermore model only
      G4PhotoElectricEffect *thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
      ph->RegisterProcess(thePhotoElectricEffect, particle);

      // Compton scattering - Livermore model only
      G4ComptonScattering *theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
      ph->RegisterProcess(theComptonScattering, particle);

      // gamma conversion - Livermore model below 80 GeV
      G4GammaConversion *theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
      ph->RegisterProcess(theGammaConversion, particle);

      // default Rayleigh scattering is Livermore
      G4RayleighScattering *theRayleigh = new G4RayleighScattering();
      ph->RegisterProcess(theRayleigh, particle);
    }

    // Warning : end of particles and processes are needed by EM Physics builders
  }

    // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  de->SetFluo(true);
  
  // **************************************
  // 2) Define processes for Target region
  // **************************************

  // STANDARD EM processes should be inactivated when
  // corresponding DNA processes are used
  // - STANDARD EM e- processes are inactivated below 1 MeV
  // - STANDARD EM proton & alpha processes are inactivated below
  //   standEnergyLimit
  G4double standEnergyLimit = 10.0 * MeV;
  //
  G4EmConfigurator *em_config =
      G4LossTableManager::Instance()->EmConfigurator();

  G4VEmModel *mod;

  // *** e-

  // ---> STANDARD EM processes are inactivated below 1 MeV in target, G4WentzelVIModel & Coulomb Emin = 100 MeV so not deactivated
  G4EmParameters* theParameters = G4EmParameters::Instance();
  G4double emax = theParameters->MaxKinEnergy();


  mod = new G4GoudsmitSaundersonMscModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("e-", "msc", mod, "Target");

  mod = new G4MollerBhabhaModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("e-",
                             "eIoni",
                             mod,
                             "Target", 0.0, emax, new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNAChampionElasticModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
                             mod, "Target");

  mod = new G4DNABornIonisationModel();
  ((G4DNABornIonisationModel *)(mod))->SelectFasterComputation(true);
  mod->SetLowEnergyLimit(11. * eV);
  mod->SetHighEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
                             mod, "Target", 11 *eV, 1*MeV);

  mod = new G4DNABornExcitationModel();
  mod->SetLowEnergyLimit(9. * eV);
  mod->SetHighEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
                             mod, "Target");

  mod = G4DNASolvationModelFactory::Create("Meesungnoen2002");
  mod->SetHighEnergyLimit(7.4 * eV);

  em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, "Target", 0 * eV, 11 * eV);

  mod = new G4DNASancheExcitationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation", mod, "Target", 0 * eV, 100 * eV);

  mod = new G4DNAMeltonAttachmentModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment", mod, "Target", 0 * eV, 13 * eV);

  // *** proton

  // ---> STANDARD EM processes inactivated below standEnergyLimit
  //      or below 1 MeV for scattering

  // Inactivate following STANDARD processes

  mod = new G4WentzelVIModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("proton", "msc", mod, "Target");

  mod = new G4eCoulombScatteringModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("proton", "CoulombScat", mod, "Target");

  mod = new G4BraggModel();
  mod->SetActivationLowEnergyLimit(100*MeV);
  em_config->SetExtraEmModel("proton", "hIoni",
                             mod, "Target");

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(100*MeV);
  em_config->SetExtraEmModel("proton", "hIoni",
                             mod, "Target");

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAElastic",
                             mod, "Target", 0 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationExtendedModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation",
                             mod, "Target", 0 * eV, 500 * keV);

  mod = new G4DNABornIonisationModel();
  ((G4DNABornIonisationModel*)(mod))->SelectFasterComputation(true);
  em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation",
                             mod, "Target", 0.5 * MeV, 100 * MeV);
                 
  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
                             mod, "Target", 10 * eV, 0.5 * MeV);

  mod = new G4DNABornExcitationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
                             mod, "Target", 0.5 * MeV, 100 * MeV);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAChargeDecrease",
                             mod, "Target", 100 * eV, 100 * MeV);

  // *** hydrogen

  // ---> NO STANDARD EM processes

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAElastic",
                             mod, "Target", 0. * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAIonisation",
                             mod, "Target", 0 * eV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAExcitation",
                             mod, "Target", 10 * eV, 0.5 * MeV);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAChargeIncrease",
                             mod, "Target", 100 * eV, 100 * MeV);

  // *** alpha

  // ---> STANDARD EM processes inactivated below standEnergyLimit
  //      or below 1 MeV for scattering

  // Inactivate following STANDARD processes

  mod = new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("alpha", "msc", mod, "Target");

  mod = new G4BraggIonModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha", "ionIoni",
                             mod, "Target");

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha", "ionIoni",
                             mod, "Target");

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAElastic",
                             mod, "Target", 100 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAIonisation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAExcitation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAChargeDecrease",
                             mod, "Target", 1 * keV, 400 * MeV);

  // *** alpha+

  // ---> STANDARD EM processes inactivated below standEnergyLimit
  //      or below 1 MeV for scattering

  // Inactivate following STANDARD processes

  mod = new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("alpha+", "msc", mod, "Target");

  mod = new G4BraggIonModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha+", "ionIoni",
                             mod, "Target");

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha+", "ionIoni",
                             mod, "Target");

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAElastic",
                             mod, "Target", 100 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAIonisation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAExcitation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeIncrease",
                             mod, "Target", 1 * keV, 400 * MeV);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeDecrease",
                             mod, "Target", 1 * keV, 400 * MeV);

  // *** helium

  // ---> NO STANDARD EM processes

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAElastic",
                             mod, "Target", 100 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAIonisation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAExcitation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAChargeIncrease",
                             mod, "Target", 1 * keV, 400 * MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
