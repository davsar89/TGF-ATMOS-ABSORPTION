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
#include <PrimaryGeneratorAction.hh>
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(Settings::Settings set_in)
        : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr) {

    this->settings=std::move(set_in);

    fParticleGun = new G4ParticleGun(nofParticles);

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
    // this function is called at the begining of event
    // 1/E
    //     G4double energy=MinEner*(pow(((MaxEner)/(MinEner)),G4UniformRand()))
    ////////////// ENERGY /////////////////
    // 1/E * exp (-E / 7.3MeV) (rejection method)
    G4double MaxEner = 40. * MeV; // Max energy
    G4double MinEner = 500. * keV; // Min energy
    G4double cut_ener = settings.E0 * MeV; // exponential cut-off factor
    G4double energy = Sample_one_RREA_gammaray_energy(MinEner, MaxEner, cut_ener, settings.ALPHA);

if (settings.ALPHA<0) {
    G4cout << "ERROR: ALPHA MUST BE POSITIVE AS THE FORMULA IS E^-ALPHA * EXP(-E/E0). ABORTING" << G4endl;
    std::abort();
}

    ////////////// POSITION / ANGLE /////////////////
    // ! : sampling theta uniformly between 0 and SOURCE_OPENING_ANGLE*degree does not sample uniformly over the area
    //  G4double theta = G4UniformRand()*shared_var::SOURCE_OPENING_ANGLE*degree;
    //     G4double theta = CLHEP::RandGauss::shoot(0., shared_var::SOURCE_OPENING_ANGLE * degree);
    //    G4double theta = settings.OpeningAngle() * degree * G4UniformRand();
    //    photons.alpha=acos(1.-rand*(1-cos(isot/180.*!pi)))/!pi*180.
    //    G4double theta;
    G4double X_try = 0, Y_try = 0;
    G4double R_max;
    G4double R_try;
    G4double sigma_sample_R;

    if ((settings.BEAMING_TYPE == "Uniform") || (settings.BEAMING_TYPE == "uniform")) {
        R_max = std::tan(settings.SOURCE_OPENING_ANGLE * degree);
        R_try = R_max + 10.; // just for initialization

        while (R_try > R_max) {
            X_try = 2. * (G4UniformRand() - 0.5) * R_max;
            Y_try = 2. * (G4UniformRand() - 0.5) * R_max;
            R_try = sqrt(X_try * X_try + Y_try * Y_try);
        }

        //         theta = std::acos(1. - G4UniformRand()*(1. - std::cos(settings.OpeningAngle() * degree))); // uniform over spherical(i.e a part of
        // sphere) area
    } else if ((settings.BEAMING_TYPE == "Gaussian") || (settings.BEAMING_TYPE == "gaussian") ||
               (settings.BEAMING_TYPE == "normal") || (settings.BEAMING_TYPE == "Normal")) {
        R_max = 10000.;      // -> maximum angle is atan(10000) = 89.9943 degrees
        sigma_sample_R = std::tan(settings.SOURCE_OPENING_ANGLE * degree);
        R_try = R_max + 10.; // just for initialization

        while (R_try > R_max) {
            X_try = CLHEP::RandGauss::shoot(0., sigma_sample_R); // gaussian position sample
            Y_try = CLHEP::RandGauss::shoot(0., sigma_sample_R); // gaussian position sample
            R_try = sqrt(X_try * X_try + Y_try * Y_try);
            //  G4cout << X_try << G4endl;
        }
    } else {
        G4cout << "ERROR : Beaming type is not Gaussian or Uniform. Aborting." << G4endl;
        std::abort();
    }

    G4double lat = settings.SOURCE_LAT;
    G4double lon = settings.SOURCE_LONG;
    G4double alt = settings.SOURCE_ALT * 1000.0; // km to m
    G4double ecef_x, ecef_y, ecef_z;

    geod_conv::GeodeticConverter::geodetic2ecef(lat, lon, alt, ecef_x, ecef_y, ecef_z);

    G4ThreeVector position;
    position.setX(ecef_x * m);
    position.setY(ecef_y * m);
    position.setZ(ecef_z * m);

    G4ThreeVector localVertical, localVertical_perp1, localVertical_perp2;

    localVertical = position / position.mag();
    // computing two vector perpendicular to the local vertical and perpendicular to each other
    G4double ux = localVertical[0];
    G4double uy = localVertical[1];
    G4double uz = localVertical[2];
    G4double norme = sqrt(pow(ux, 2) + pow(uy, 2));
    G4double vx = -uy / norme;
    G4double vy = ux / norme;
    G4double vz = 0.;
    G4double wx = -uz * vy;
    G4double wy = uz * vx;
    G4double wz = ux * vy - uy * vx;
    localVertical_perp1 = G4ThreeVector(vx, vy, vz);
    localVertical_perp2 = G4ThreeVector(wx, wy, wz);

    // adding tilt angle and recomputing the two perpendicular vectors
    if (settings.TILT_ANGLE != 0.0) {
        G4ThreeVector tilt_shift =
                localVertical_perp1 * std::sin(settings.TILT_ANGLE * degree); // we could also use localVertical_perp2
        localVertical = localVertical + tilt_shift;
        localVertical = position / position.mag();
        ux = localVertical[0];
        uy = localVertical[1];
        uz = localVertical[2];
        norme = sqrt(pow(ux, 2) + pow(uy, 2));
        vx = -uy / norme;
        vy = ux / norme;
        vz = 0.;
        wx = -uz * vy;
        wy = uz * vx;
        wz = ux * vy - uy * vx;
        localVertical_perp1 = G4ThreeVector(vx, vy, vz);
        localVertical_perp2 = G4ThreeVector(wx, wy, wz);
    }

    G4ThreeVector position_virtual_particle;
    position_virtual_particle = position + localVertical + localVertical_perp1 * X_try + localVertical_perp2 * Y_try;
    G4ThreeVector momentumDirection = (position_virtual_particle - position);
    momentumDirection = momentumDirection / momentumDirection.mag();
    //     theta = std::acos(momentumDirection.dot(localVertical))/degree;
    //     G4cout << theta << G4endl;
    ////////////// TIME /////////////////
    G4double time;

    if (settings.SOURCE_SIGMA_TIME == 0.) {
        time = 0.;
    } else {
        time = CLHEP::RandGauss::shoot(0., settings.SOURCE_SIGMA_TIME) * microsecond;
    }

    ////////////// assignments /////////////////
    fParticleGun->SetParticleTime(time);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleMomentumDirection(momentumDirection);
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::BrokenPL(double &p1, double &p2, double &ec, double &x) {
    double broken_PL = 0;

    if (x < ec) {
        broken_PL = std::pow(x, p1);
    } else if (x >= ec) {
        broken_PL = std::pow(ec, p1 - p2) * std::pow(x, p2);
    }

    return broken_PL;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double
PrimaryGeneratorAction::Sample_one_RREA_gammaray_energy(double &MinEner, double &MaxEner, double &cut_ener, double &power) {
    // random samples the energy of one RREA gamma ray
    // RREEA gamma spectrum is approximately = 1/E * exp (-E / 7.3MeV)
    //     G4double energy=MinEner*(pow(((MaxEner)/(MinEner)),G4UniformRand())); // code for sampling 1/E
    // (rejection method)
    G4double H = cut_ener;
    G4double pMax = std::pow(MinEner,-power) * exp(-MinEner / H);
    G4double pMin = std::pow(MaxEner,-power) * exp(-MaxEner / H);
    G4double pOfeRand = 0.0;
    G4double pRand = 1.0;
    G4double eRand;

    while (pOfeRand < pRand) { // rejection
        pRand = pMin + (pMax - pMin) * G4UniformRand();
        eRand = MinEner + (MaxEner - MinEner) * G4UniformRand();
        pOfeRand = std::pow(eRand,-power) * exp(-eRand / H);
    }

    return eRand;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
