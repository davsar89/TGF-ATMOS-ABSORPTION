//
// Created by dsa030 on 6/28/19.
//

#include "histogram_sampler.hh"

histogram_sampler::histogram_sampler(const std::vector<double> binss, const std::vector<double> valuess) {

    // valuess,binss = tabulated function(=distribution) that we want to sample from
    // f is assumed positive, linear per segment, continuous

    if (binss.size() != valuess.size() + 1) {
        std::cout << "ERROR in histogram_sampler : bins.size() != values.size() + 1  ;;; ABORTING" << std::endl;
        std::abort();
    }

    this->bins = binss;
    this->f = valuess;

    nb_bins = bins.size();

    for (int jj = 0; jj <= nb_bins - 2; ++jj) {
        bin_centers.push_back((bins[jj] + bins[jj + 1]) / 2.0);
    }

    nPoints = bin_centers.size();

    initialize_distrubtion_sampling();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

histogram_sampler::~histogram_sampler() {

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double histogram_sampler::Inverse_Cumul_sampling()
// sampling energy from spectrum using the inverse cumulative distribution function method
// also proceeds to interpolation
{
    // tabulated function
    // f is assumed positive, linear per segment, continuous 
    // --> cumulative function is second order polynomial

    //choose y randomly
    double y_rndm = my_UniformRand() * Fc[nPoints - 1];

    //find bin
    int j = nPoints - 2;
    while ((Fc[j] > y_rndm) && (j > 0)) j--;

    //y_rndm --> x_rndm :  Fc(x) is second order polynomial
    double x_rndm = bin_centers[j];
    double aa = a[j];

    if (aa != 0.) {
        double b = f[j] / aa, c = 2 * (y_rndm - Fc[j]) / aa;
        double delta = b * b + c;
        int sign = 1;
        if (aa < 0.) sign = -1;
        x_rndm += sign * std::sqrt(delta) - b;
    } else if (f[j] > 0.) {
        x_rndm += (y_rndm - Fc[j]) / f[j];
    }

    return x_rndm;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void histogram_sampler::initialize_distrubtion_sampling()
// compute useful quantities for the inverse cumulative distribution sampling
{
    // compute fMax
    fMax = 0.;

    for (int ii = 0; ii < nPoints; ii++) {
        if (fMax < f[ii]) fMax = f[ii];
    }

    //compute slopes
    for (int ii = 0; ii < nPoints - 1; ii++) {
        a.push_back((f[ii + 1] - f[ii]) / (bin_centers[ii + 1] - bin_centers[ii]));
    }

    //compute cumulative function
    Fc.push_back(0.0);

    for (int jj = 1; jj < nPoints; jj++) {
        Fc.push_back(Fc[jj - 1] + 0.5 * (f[jj] + f[jj - 1]) * (bin_centers[jj] - bin_centers[jj - 1]));
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double histogram_sampler::my_UniformRand() {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}