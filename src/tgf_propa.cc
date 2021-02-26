////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in
// Earth's atmosphere */
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

#include "G4UImanager.hh"
#include <DetectorConstruction.hh>
#include <PhysicsList.hh>

#ifdef G4VIS_USE

#include "G4VisExecutive.hh"

#endif // ifdef G4VIS_USE

#ifdef G4UI_USE

#include "G4UIExecutive.hh"

#endif // ifdef G4UI_USE

#include <RunAction.hh>

#include <EventAction.hh>
#include <SteppingAction.hh>

#include <chrono>

#include <G4Threading.hh>

#ifdef G4MULTITHREADED

#include <thread>
#include "G4MTRunManager.hh"

#else
#include "G4RunManager.hh"
#endif

#include "ActionInitialization.hh"

#include <CoordTopocentric.h>
#include <CoordGeodetic.h>
#include <Observer.h>
#include <SGP4.h>
#include "Settings.hh"
#include <string>
#include "myUtils.hh"

using namespace std;

/////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    double wall0 = myUtils::get_wall_time();

    Settings::Settings settings;

    G4String number_st;
    G4String Mode = "run";

    if (argc >= 3)
    {
        number_st = argv[1];
        settings.ALPHA = std::stod(argv[2]);
        settings.E0 = std::stod(argv[3]);
    }
    else
    {
        // default values can be seen/changed in src/include/Settings.hh
        Mode = "run";
        number_st = "10000000";
    }

    ///
#ifdef G4MULTITHREADED
    int nb_cores = G4Threading::G4GetNumberOfCores();
    G4cout << "Number of cores: " << nb_cores << G4endl;
#else
    int nb_cores = 1;
#endif

    if (nb_cores == 8)
    {
        nb_cores = 4;
        G4cout << "Number of cores to be used: " << nb_cores << G4endl;
    }

    if (nb_cores == 64)
    {
        nb_cores = 32;
        G4cout << "Number of cores to be used: " << nb_cores << G4endl;
    }

    if (nb_cores == 12)
    {
        nb_cores = 4;
        G4cout << "Number of cores to be used: " << nb_cores << G4endl;
    }

    //// choose the Random engine and give seed
    G4Random::setTheEngine(new CLHEP::MixMaxRng);

    // random seeds, different each time code is run
    long start = myUtils::generate_a_unique_ID();
    long seeds[2];
    seeds[0] = start;
    seeds[1] = myUtils::generate_a_unique_ID();

    G4cout << "seed set to: " << start << " ns" << G4endl;

    G4Random::setTheSeeds(seeds,2);

    //    The use of random number generators is a bit different in multi-threaded mode with respect to sequential mode.
    //    To guarantee reproducibility independently of the number of threads some modifications have been made:
    //    The master (main program) can be seeded as it is in sequential mode via G4Random::setTheSeed( … ) (and UI command).
    //    To guarantee that RNG numbers do not depend on the number of threads the "history" of random-seeds is not as in sequential.
    //    Each thread has its own RNG engine and the status of these are independent.
    //    Before the run loop starts the master pre-generates a random seed for each event to be processed.
    //    Then threads are spawned and event loop proceeds, at the beginning of each event the worker thread re-seeds
    //    the event with this pre-generated map of seeds from the master thread.

    // Construct the default run manager

#ifdef G4MULTITHREADED
    auto *runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(nb_cores);
#else
    auto *runManager = new G4RunManager;
#endif

    /// check if MAC, magnetic field model must be WMM or GEOLIB
#if __APPLE__
    if (settings.MAGNETIC_FIELD_MODEL == "IGRF")
    {
        G4cout << "ERROR: on macOS, settings.MAGNETIC_FIELD_MODEL cannot be IGRF. You can change it in include/Settings.hh . You may also set up IGRF using the GEOLIB . Aborting." << G4endl;
        std::abort();
    }
#endif

    //
    runManager->SetUserInitialization(new TGFDetectorConstruction(settings));

    TGF_PhysicsList *physicsList = new TGF_PhysicsList(settings);
    runManager->SetUserInitialization(physicsList);

    //
    runManager->SetUserInitialization(new ActionInitialization(settings));

    // Initialize G4 kernel
    runManager->Initialize();
    G4cout << G4endl << "Initialization OK" << G4endl;

    // get the pointer to the User Interface manager
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (Mode == "visu")
    {
#ifdef G4VIS_USE
        G4VisManager *visManager = new G4VisExecutive;
        visManager->Initialize();
#endif // ifdef G4VIS_USE
#ifdef G4UI_USE
        G4UIExecutive *ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute vis.mac");
#endif // ifdef G4VIS_USE
        ui->SessionStart();
        delete ui;
#endif // ifdef G4UI_USE
#ifdef G4VIS_USE
        delete visManager;
#endif // ifdef G4VIS_USE
    }
    else if (Mode == "run")
    {
        UImanager->ApplyCommand("/run/beamOn " + number_st);
    }

    //
    delete runManager;
    //
    double wall1 = myUtils::get_wall_time();
    G4cout << G4endl << "WALL TIME TAKEN : " << wall1 - wall0 << " seconds "
           << G4endl;
    //
    return 0;
} // main
