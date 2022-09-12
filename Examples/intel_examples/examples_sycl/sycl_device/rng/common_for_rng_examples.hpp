/*******************************************************************************
* Copyright 2020 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
*
*  Content:
*       File contains service functionality and checkers for statistics of
*       various rng distributions
*
*******************************************************************************/

#ifndef __COMMON_FOR_RNG_DEVICE_EXAMPLES_HPP__
#define __COMMON_FOR_RNG_DEVICE_EXAMPLES_HPP__

// stl includes
#include <iostream>
#include <vector>
#include <math.h>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl/rng/device.hpp"

// local includes
#include "common_for_examples.hpp"

// function to print rng output from buffer
template<typename Type>
void print_output(sycl::buffer<Type, 1>& r, int n_print) {
    auto r_accessor = r.template get_access<sycl::access::mode::read>();
    std::cout << "First "<< n_print << " numbers of " << r.get_count() << ": " << std::endl;
    for(int i = 0 ; i < n_print; i++) {
        std::cout << r_accessor[i] << " ";
    }
    std::cout << std::endl;
}

// function to compare theoretical moments and sample moments
template<typename Type>
int compare_moments(std::vector<Type>& r, double tM, double tD, double tQ) {
    double tD2;
    double sM, sD;
    double sum, sum2;
    double n, s;
    double DeltaM, DeltaD;

    // sample moments
    sum = 0.0;
    sum2 = 0.0;
    for(int i = 0; i < r.size(); i++) {
        sum += (double)r[i];
        sum2 += (double)r[i] * (double)r[i];
    }
    sM = sum / ((double)r.size());
    sD = sum2 / (double)r.size() - (sM * sM);

    // comparison of theoretical and sample moments
    n = (double)r.size();
    tD2 = tD * tD;
    s = ((tQ-tD2) / n) - ( 2 * (tQ - 2 * tD2) / (n * n))+((tQ - 3 * tD2) /
                                                            (n * n * n));

    DeltaM = (tM - sM) / sqrt(tD / n);
    DeltaD = (tD - sD) / sqrt(s);
    if(fabs(DeltaM) > 3.0 || fabs(DeltaD) > 3.0) {
        std::cout << "Error: sample moments (mean=" << sM << ", variance=" << sD
            << ") disagree with theory (mean=" << tM << ", variance=" << tD <<
            ")" << std:: endl;
        return 1;
    }
    std::cout << "Success: sample moments (mean=" << sM << ", variance=" << sD
        << ") agree with theory (mean=" << tM << ", variance=" << tD <<
        ")" << std:: endl;
    return 0;
}

// structure is used to calculate theoretical moments of particular distribution
// and compare them with sample moments
template <typename Type, typename Distribution>
struct statistics {};

template<typename Type>
struct statistics<Type, oneapi::mkl::rng::device::uniform<Type>> {
    int check(std::vector<Type>& r, const oneapi::mkl::rng::device::uniform<Type>& distr) {
        double tM, tD, tQ;
        Type a = distr.a();
        Type b = distr.b();

        // theoretical moments of uniform real type distribution
        tM = (b + a) / 2.0;
        tD = ((b - a) * (b - a)) / 12.0;
        tQ = ((b - a)*(b - a)*(b - a)*(b - a)) / 80.0;

        return compare_moments(r, tM, tD, tQ);
    }
};

template<>
struct statistics<std::int32_t, oneapi::mkl::rng::device::uniform<std::int32_t>> {
    int check(std::vector<std::int32_t>& r, const oneapi::mkl::rng::device::uniform<std::int32_t>& distr) {
        double tM, tD, tQ;
        std::int32_t a = distr.a();
        std::int32_t b = distr.b();

        // theoretical moments of uniform int distribution
        tM = (a + b - 1.0) / 2.0;
        tD = ((b - a)*(b - a) - 1.0) / 12.0;
        tQ = (((b - a) * (b - a)) * ((1.0 / 80.0)*(b - a)*(b - a) - (1.0 / 24.0))) + (7.0 / 240.0);

        return compare_moments(r, tM, tD, tQ);
    }
};

template<typename Type>
struct statistics<Type, oneapi::mkl::rng::device::gaussian<Type>> {
    int check(std::vector<Type>& r, const oneapi::mkl::rng::device::gaussian<Type>& distr) {
        double tM, tD, tQ;
        Type a = distr.mean();
        Type sigma = distr.stddev();

        // theoretical moments of gaussian distribution
        tM = a;
        tD = sigma * sigma;
        tQ = 720.0 * sigma * sigma * sigma * sigma;

        return compare_moments(r, tM, tD, tQ);
    }
};

template<typename Type>
struct statistics<Type, oneapi::mkl::rng::device::lognormal<Type>> {
    int check(std::vector<Type>& r, const oneapi::mkl::rng::device::lognormal<Type>& distr) {
        double tM, tD, tQ;
        Type a = distr.m();
        Type b = distr.displ();
        Type sigma = distr.s();
        Type beta = distr.scale();

        // theoretical moments of lognormal distribution
        tM = b + beta * exp(a + sigma * sigma * 0.5);
        tD = beta * beta * exp(2.0 * a + sigma * sigma) * (exp(sigma * sigma) - 1.0);
        tQ = beta * beta * beta * beta * exp(4.0 * a + 2.0 * sigma * sigma) *
            (exp(6.0 * sigma * sigma) - 4.0 * exp(3.0 * sigma * sigma) + 6.0 * exp(sigma * sigma) - 3.0);

        return compare_moments(r, tM, tD, tQ);
    }
};

#endif // __COMMON_FOR_RNG_DEVICE_EXAMPLES_HPP__
