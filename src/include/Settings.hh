//////////////////////////////////////////////////////////////////////////////////

//// /* GEANT4 code for propagation of gamma-rays, electron and positrons in
/// Earth's atmosphere */
////
//// //
//// // ********************************************************************
//// // * License and Disclaimer                                           *
//// // *                                                                  *
//// // * The  Geant4 software  is  copyright of the Copyright Holders  of *
//// // * the Geant4 Collaboration.  It is provided  under  the terms  and *
//// // * conditions of the Geant4 Software License,  included in the file *
//// // * LICENSE and available at  http://cern.ch/geant4/license .  These *
//// // * include a list of copyright holders.                             *
//// // *                                                                  *
//// // * Neither the authors of this software system, nor their employing *
//// // * institutes,nor the agencies providing financial support for this *
//// // * work  make  any representation or  warranty, express or implied, *
//// // * regarding  this  software system or assume any liability for its *
//// // * use.  Please see the license in the file  LICENSE  and URL above *
//// // * for the full disclaimer and the limitation of liability.         *
//// // *                                                                  *
//// // * This  code  implementation is the result of  the  scientific and *
//// // * technical work of the GEANT4 collaboration.                      *
//// // * By using,  copying,  modifying or  distributing the software (or *
//// // * any work based  on the software)  you  agree  to acknowledge its *
//// // * use  in  resulting  scientific  publications,  and indicate your *
//// // * acceptance of all terms of the Geant4 Software license.          *
//// // ********************************************************************
//////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <vector>
#include <string>

// there is a a copy of it on class that need it
// should not be used to share values around

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Settings {
    struct Settings {

    public:
        //////////////// Parameters are listed below ////////////////
        // ASIM TEB 190324
        // -6.998426914215088, 55.92176818847656 -> magnetic field line tracing
        // -7.0492, 55.9129 −> VAISALA
        // geodetic coordinates

        double SOURCE_LAT = 26.0 ;  // deg
        double SOURCE_LONG = -77.0 ;  // deg
        double SOURCE_ALT = 15.0 ;  // km
        double record_altitude = 500.0 ; // km

        double SOURCE_OPENING_ANGLE = 10.0 ; // deg
        G4String BEAMING_TYPE = "Gaussian";

        double TILT_ANGLE = 0.0;

        double SOURCE_SIGMA_TIME = 0.; // microsecond

        /// Below variables should stay CONSTANT

        // force max step only for layers where particles are recorded
        G4bool USE_STEP_MAX_for_record = true;
        double STEP_MAX_RECORD_AREA = 100.0 * meter;

        G4bool OUTPUT_TO_ASCII_FILE = true;

        double dr_over_R = 0.4; // stepping parameter for ionization, default is 0.2, that may be too high
        double MIN_ENERGY_OUTPUT = 300.0 * keV;

        // date
        int dt_year = 2019;
        double dt_month = 3;
        double dt_day = 24;
        double dt_hour = 0;
        double dt_minute = 31;
        double dt_second = 53;
        double dt_microsecond = 135444;

        // Earth radius
        double earthRadius = 6378.137 * km;

        double E0 = 7.3; // MeV
        double ALPHA = 1.0;

        //std::string SPEC_TYPE = G4String("dwyer_reverse_beam");
        std::string SPEC_TYPE = G4String("usual");

    };
}
