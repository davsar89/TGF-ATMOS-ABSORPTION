////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere */
//
// //
// // ********************************************************************
// // * License and Disclaimer                                           *
// // *                                                                  *
// // * The  Geant4 software  is  copyright of the Copyright Holders  of *
// // * the Geant4 Collaboration.  It is provided  under  the terms  and *
// // * conditions of the Geant4 Software License,  included in the file *
// // * LICENSE and available at  http://cern.ch/geant4/license .  These *
// // * include a list of copyright holders.                             *
// // *                                                                  *
// // * Neither the authors of this software system, nor their employing *
// // * institutes,nor the agencies providing financial support for this *
// // * work  make  any representation or  warranty, express or implied, *
// // * regarding  this  software system or assume any liability for its *
// // * use.  Please see the license in the file  LICENSE  and URL above *
// // * for the full disclaimer and the limitation of liability.         *
// // *                                                                  *
// // * This  code  implementation is the result of  the  scientific and *
// // * technical work of the GEANT4 collaboration.                      *
// // * By using,  copying,  modifying or  distributing the software (or *
// // * any work based  on the software)  you  agree  to acknowledge its *
// // * use  in  resulting  scientific  publications,  and indicate your *
// // * acceptance of all terms of the Geant4 Software license.          *
// // ********************************************************************
////////////////////////////////////////////////////////////////////////////////
#include "PhysicsList.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option4.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal StepMax *TGF_PhysicsList::fStepMaxProcess_elec = 0;
G4ThreadLocal StepMax *TGF_PhysicsList::fStepMaxProcess_gamma = 0;

TGF_PhysicsList::TGF_PhysicsList(Settings::Settings sets_int) : G4VUserPhysicsList() {

    this->settings = std::move(sets_int);

//    emPhysicsList = new G4EmStandardPhysics_option4_dr(this->settings);
    emPhysicsList = new G4EmStandardPhysics_option4();
    this->DumpCutValuesTable();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TGF_PhysicsList::~TGF_PhysicsList() {
    delete emPhysicsList;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::ConstructParticle() {
    emPhysicsList->ConstructParticle();
    //   G4GenericIon::GenericIon();
    //   G4NeutrinoE::NeutrinoEDefinition();
    //   G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    //   G4Alpha::AlphaDefinition();
    //
    //   G4Geantino::GeantinoDefinition();
    //   G4ChargedGeantino::ChargedGeantinoDefinition();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::ConstructProcess() {
    AddTransportation();
    emPhysicsList->ConstructProcess();
    //   G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
    //   radioactiveDecay->SetICM(true);                //Internal Conversion
    //   radioactiveDecay->SetARM(true);               //Atomic Rearangement
    // radioactiveDecay->SetFBeta(true);
    // radioactiveDecay->SetBRBias(true);
    //   G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    //
    //   ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

    /// STEP LIMITATION for photons in record area

    Add_global_StepMax();

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::SetCuts() {
    defaultCutValue = 100. * m;
    //
    cutForGamma = defaultCutValue;
    cutForElectron = defaultCutValue;
    cutForPositron = defaultCutValue;
    //
    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");

//    G4double lowlimit = settings.MIN_ENERGY_OUTPUT;
  //  G4ProductionCutsTable *aPCTable = G4ProductionCutsTable::GetProductionCutsTable();
   // aPCTable->SetEnergyRange(10.0 * keV, 100 * CLHEP::GeV);

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "StepMax.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"

void TGF_PhysicsList::Add_global_StepMax() {
    // Step limitation seen as a process
    StepMax *stepMaxProcess = new StepMax();

    stepMaxProcess->SetMaxStep(50*km);

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        G4ParticleDefinition *particle = particleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();

        if (stepMaxProcess->IsApplicable(*particle)) {
            pmanager->AddDiscreteProcess(stepMaxProcess);
        }
    }
}