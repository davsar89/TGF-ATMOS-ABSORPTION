//
// Created by dsa030 on 6/28/19.
//

#pragma once

#include <vector>
#include <random>
#include <assert.h>
#include <iostream>

class histogram_sampler {

public:
    histogram_sampler(const std::vector<double> bins, const std::vector<double> values);

    double Inverse_Cumul_sampling();


private:
    ~histogram_sampler();

    double my_UniformRand();

    void initialize_distrubtion_sampling();

    std::vector<double> f; // f(between two bins)
    std::vector<double> bins;
    std::vector<double> bin_centers;
    std::vector<double> Fc;// cumulative of f
    double fMax = 0;        // max(f)
    std::vector<double> a;  // slopes

    std::default_random_engine generator;

    int nb_bins = -10; // initialization
    int nPoints;

};