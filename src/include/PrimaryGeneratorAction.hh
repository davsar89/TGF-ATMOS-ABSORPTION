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
            1.01658e-1  , 2.05328e+6,
            1.18850e-1  , 1.81196e+6,
            1.47183e-1  , 1.55563e+6,
            1.74287e-1  , 1.34214e+6,
            2.09618e-1  , 1.16185e+6,
            2.37137e-1  , 1.02489e+6,
            2.71142e-1  , 9.28210e+5,
            3.18839e-1  , 8.20809e+5,
            3.59799e-1  , 7.10642e+5,
            4.28690e-1  , 6.21306e+5,
            4.97083e-1  , 5.48230e+5,
            5.95662e-1  , 4.57397e+5,
            7.13792e-1  , 3.76268e+5,
            8.69532e-1  , 3.00981e+5,
            1.05925e+0  , 2.54731e+5,
            1.33352e+0  , 1.98174e+5,
            1.53338e+0  , 1.65593e+5,
            1.82270e+0  , 1.41806e+5,
            2.18513e+0  , 1.11985e+5,
            2.57463e+0  , 9.06807e+4,
            3.16552e+0  , 7.40698e+4,
            3.88406e+0  , 5.96946e+4,
            4.57831e+0  , 4.78661e+4,
            5.48638e+0  , 3.81728e+4,
            6.53493e+0  , 2.97178e+4,
            7.87823e+0  , 2.44153e+4,
            9.38941e+0  , 1.82601e+4,
            1.09468e+1  , 1.51755e+4,
            1.26932e+1  , 1.24788e+4,
            1.47183e+1  , 9.69837e+3,
            1.76372e+1  , 7.33077e+3,
            2.04510e+1  , 5.69741e+3,
            2.44116e+1  , 3.94288e+3,
            2.74970e+1  , 3.07416e+3,
            3.23014e+1  , 2.13132e+3,
            3.66679e+1  , 1.56758e+3,
            4.00155e+1  , 1.16649e+3,
            4.32229e+1  , 8.93476e+2,
            4.64107e+1  , 6.41986e+2,
            4.93012e+1  , 4.87959e+2,
            5.20834e+1  , 3.58540e+2,
            5.44145e+1  , 2.70168e+2,
            5.65551e+1  , 1.96424e+2,
            5.90784e+1  , 1.35494e+2,
            6.10727e+1  , 9.25502e+1,
            6.25790e+1  , 6.42051e+1,
            6.34677e+1  , 4.30767e+1,
            6.46717e+1  , 3.21804e+1,
            6.68344e+1  , 2.29484e+1,
            6.79427e+1  , 1.63616e+1,
            6.90694e+1  , 1.16654e+1
    };


};
