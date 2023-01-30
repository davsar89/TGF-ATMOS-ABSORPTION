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
#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "Settings.hh"
#include <vector>
#include <string>

#include "geodetic_converter.hh"
#include "histogram_sampler.hh"

class G4ParticleGun;

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:

    double Sample_one_RREA_gammaray_energy(double &MinEner, double &MaxEner, double &cut_ener, double & power);

    double BrokenPL(double &p1, double &p2, double &ec, double &x);

    explicit PrimaryGeneratorAction(Settings::Settings set_in);

    ~PrimaryGeneratorAction();

    // methods
    virtual void GeneratePrimaries(G4Event *);



private:

    Settings::Settings settings;

    // data members
    G4ParticleGun *fParticleGun; // pointer a to G4 service class
    G4int nofParticles = 1;

    double GenerateDwyerReverseBeam(const double& emin, const double& emax);

    std::vector<double> data_RB_x{};
    std::vector<double> data_RB_y{};

    histogram_sampler* RB_sampler;

    std::vector<double> data_RB{
    1.18850e-1,  1.81196e+6,
    1.74287e-1,  1.34214e+6,
    2.53878e-1,  1.05300e+6,
    3.59799e-1,  7.10642e+5,
    4.97083e-1,  5.48230e+5,
    7.13792e-1,  3.76268e+5,
    1.05925e+0,  2.54731e+5,
    1.53338e+0,  1.65593e+5,
    2.20317e+0 , 1.16837e+5,
    3.16552e+0  ,7.40698e+4,
    4.57831e+0 , 4.78661e+4,
    6.53493e+0  ,2.97178e+4,
    9.38941e+0 , 1.82601e+4,
    1.26932e+1 , 1.24788e+4,
    1.76372e+1 , 7.33077e+3,
    2.44116e+1 , 3.94288e+3,
    3.23014e+1 , 2.13132e+3,
    4.00155e+1 , 1.16649e+3,
    4.64107e+1 , 6.41986e+2,
    5.20834e+1 , 3.58540e+2,
    5.65551e+1 , 1.96424e+2,
    6.00765e+1 , 9.25315e+1,
    6.34677e+1 , 4.30767e+1,
    6.68344e+1 , 2.29484e+1,
    6.90694e+1 , 1.16654e+1};


};
