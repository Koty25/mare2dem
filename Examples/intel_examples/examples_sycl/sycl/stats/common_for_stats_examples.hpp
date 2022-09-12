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
*      Common functionality for statistics examples
*
*******************************************************************************/

#ifndef __COMMON_FOR_STATS_EXAMPLES_HPP__
#define __COMMON_FOR_STATS_EXAMPLES_HPP__

// stl includes
#include <vector>

// local includes
#include "common_for_examples.hpp"

template<typename Vec>
void random_vector(Vec &v, int n_observations, int n_dims) {
    std::srand(1);

    using fp = typename Vec::value_type;
    v.resize(n_observations * n_dims);

    for(int i = 0; i < n_observations * n_dims; i++) {
        v[i] = rand_scalar<fp>();
    }
}

#endif // __COMMON_FOR_STATS_EXAMPLES_HPP__
