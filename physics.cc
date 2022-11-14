#include "physics.hh"

MyPhysicsList::MyPhysicsList()
{
    G4int verb = 1;
    RegisterPhysics(new G4EmStandardPhysics()); // EM interactions
    // G4EmParameters *param = G4EmParameters::Instance();
    // param->SetAugerCascade(true);
    // param->SetStepFunction(1., 1 * CLHEP::mm);
    // param->SetStepFunctionMuHad(1., 1 * CLHEP::mm);
    RegisterPhysics(new G4OpticalPhysics()); // optical photons for Cherenkov light
    // for radioactive decay
    RegisterPhysics(new G4DecayPhysics(verb));        // any kind of particle decays
    RegisterPhysics(new G4RadioactiveDecayPhysics()); // radioactive ions decay
                                                      // Radioactive decay

    SetVerboseLevel(verb);
    // RegisterPhysics(new BiasedRDPhysics());
    //  Hadron Elastic scattering
    RegisterPhysics(new G4HadronElasticPhysics(verb));
    // Hadron Inelastic physics
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT(verb));
    ////RegisterPhysics( new G4HadronInelasticQBBC(verb));
    ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));

    // Ion Elastic scattering
    RegisterPhysics(new G4IonElasticPhysics(verb));

    // Ion Inelastic physics
    RegisterPhysics(new G4IonPhysics(verb));
    ////RegisterPhysics( new G4IonINCLXXPhysics(verb));

    // Gamma-Nuclear Physics
    G4EmExtraPhysics *gnuc = new G4EmExtraPhysics(verb);
    gnuc->ElectroNuclear(false);
    gnuc->MuonNuclear(false);
    RegisterPhysics(gnuc);

    RegisterPhysics(new G4StepLimiterPhysics());
}

MyPhysicsList::~MyPhysicsList()
{
}

void MyPhysicsList::ConstructParticle()
{
    G4BosonConstructor pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}

void MyPhysicsList::SetCuts()
{
    SetCutValue(0 * mm, "proton");
    SetCutValue(10 * km, "e-");
    SetCutValue(10 * km, "e+");
    SetCutValue(10 * km, "gamma");
}