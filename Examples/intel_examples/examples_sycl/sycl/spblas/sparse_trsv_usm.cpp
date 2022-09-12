/*******************************************************************************
* Copyright 2019-2020 Intel Corporation.
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
*       This example demonstrates use of oneAPI Math Kernel Library (oneMKL)
*       DPCPP USM API oneapi::mkl::sparse::trsv to perform
*
*       sparse triangular solve on a SYCL device (Host, CPU, GPU).
*
*       op(A) * y = x
*
*       where op() is defined by one of
*oneapi::mkl::transpose::{nontrans,trans,conjtrans}
*
*
*       The supported floating point data types for gemm matrix data are:
*           float
*           double
*
*
*******************************************************************************/

// stl includes
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <vector>

#include "mkl.h"
#include "oneapi/mkl.hpp"
#include <CL/sycl.hpp>

// local includes
#include "../common/common_for_examples.hpp"
#include "common_for_sparse_examples.hpp"

//
// Main example for Sparse Triangular Solver consisting of
// initialization of A matrix, x and y vectors.
// Then the following system is solved
//
// op(A) * y = x
//
// and finally the results are post processed.
//
template <typename fp, typename intType>
int run_sparse_triangular_solve_example(const cl::sycl::device &dev)
{
    // Initialize data for Sparse Triangular Solve
    oneapi::mkl::transpose transA = oneapi::mkl::transpose::nontrans;
    oneapi::mkl::diag diag_val    = oneapi::mkl::diag::nonunit;

    // Matrix data size
    intType size  = 4;
    intType nrows = size * size * size;

    // Catch asynchronous exceptions
    auto exception_handler = [](cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const &e : exceptions) {
            try {
                std::rethrow_exception(e);
            }
            catch (cl::sycl::exception const &e) {
                std::cout << "Caught asynchronous SYCL "
                             "exception during sparse::trsv:\n"
                          << e.what() << std::endl;
            }
        }
    };

    // create execution queue and buffers of matrix data
    cl::sycl::event optimize_done, trsv_done;

    cl::sycl::queue main_queue(dev, exception_handler);
    cl::sycl::context cxt = main_queue.get_context();

    intType *ia, *ja;
    fp *a, *x, *y, *z;
    intType sizea   = 27 * nrows;
    intType sizeja  = 27 * nrows;
    intType sizeia  = nrows + 1;
    intType sizevec = nrows;

    ia = (intType *)malloc_shared(sizeia * sizeof(intType), dev, cxt);
    ja = (intType *)malloc_shared(sizeja * sizeof(intType), dev, cxt);
    a  = (fp *)malloc_shared(sizea * sizeof(fp), dev, cxt);
    x  = (fp *)malloc_shared(sizevec * sizeof(fp), dev, cxt);
    y  = (fp *)malloc_shared(sizevec * sizeof(fp), dev, cxt);
    z  = (fp *)malloc_shared(sizevec * sizeof(fp), dev, cxt);

    if (!ia || !ja || !a || !x || !y || !z)
        throw std::runtime_error("Failed to allocate USM memory");

    generate_sparse_matrix<fp, intType>(size, ia, ja, a);

    // Init vectors x and y
    for (int i = 0; i < nrows; i++) {
        x[i] = set_fp_value(fp(1.0), fp(0.0));
        y[i] = set_fp_value(fp(0.0), fp(0.0));
        z[i] = set_fp_value(fp(0.0), fp(0.0));
    }

    std::vector<intType *> int_ptr_vec;
    int_ptr_vec.push_back(ia);
    int_ptr_vec.push_back(ja);
    std::vector<fp *> fp_ptr_vec;
    fp_ptr_vec.push_back(a);
    fp_ptr_vec.push_back(x);
    fp_ptr_vec.push_back(y);
    fp_ptr_vec.push_back(z);

    //
    // Execute Triangulat Solve
    //

    // create and initialize handle for a Sparse Matrix in CSR format
    oneapi::mkl::sparse::matrix_handle_t handle;

    try {
        oneapi::mkl::sparse::init_matrix_handle(&handle);

        oneapi::mkl::sparse::set_csr_data(handle, nrows, nrows, oneapi::mkl::index_base::zero, ia,
                                          ja, a);

        optimize_done = oneapi::mkl::sparse::optimize_trsv(main_queue, oneapi::mkl::uplo::lower,
                                                           oneapi::mkl::transpose::nontrans,
                                                           diag_val, handle, {});

        // add oneapi::mkl::sparse::trsv to execution queue
        trsv_done = oneapi::mkl::sparse::trsv(main_queue, oneapi::mkl::uplo::lower,
                                              oneapi::mkl::transpose::nontrans, diag_val, handle, x,
                                              y, {optimize_done});

        oneapi::mkl::sparse::release_matrix_handle(&handle, {optimize_done, trsv_done});
    }
    catch (cl::sycl::exception const &e) {
        std::cout << "\t\tCaught synchronous SYCL exception:\n" << e.what() << std::endl;

        cleanup_arrays<fp, intType>(fp_ptr_vec, int_ptr_vec, cxt);
        oneapi::mkl::sparse::release_matrix_handle(&handle);

        return 1;
    }
    catch (std::exception const &e) {
        std::cout << "\t\tCaught std exception:\n" << e.what() << std::endl;

        cleanup_arrays<fp, intType>(fp_ptr_vec, int_ptr_vec, cxt);
        oneapi::mkl::sparse::release_matrix_handle(&handle);

        return 1;
    }

    //
    // Post Processing
    //

    std::cout << "\n\t\tsparse::trsv parameters:\n";
    std::cout << "\t\t\ttransA = "
              << (transA == oneapi::mkl::transpose::nontrans ?
                          "nontrans" :
                          (transA == oneapi::mkl::transpose::trans ? "trans" : "conjtrans"))
              << std::endl;
    std::cout << "\t\t\tnrows = " << nrows << std::endl;

    fp *res = y;
    for (intType row = 0; row < nrows; row++) {
        fp tmp      = x[row];
        fp diag_val = set_fp_value(fp(0.0), fp(0.0));
        for (intType i = ia[row]; i < ia[row + 1]; i++) {
            if (ja[i] < row) {
                tmp -= a[i] * z[ja[i]];
            }
            else if (ja[i] == row) {
                diag_val = a[i];
            }
        }
        z[row] = tmp / diag_val;
        if(!check_result(res[row], z[row], row)) return 1;
    }

    fp diff = set_fp_value(fp(0.0), fp(0.0));
    for (intType i = 0; i < nrows; i++) {
        diff += (z[i] - res[i]) * (z[i] - res[i]);
    }

    std::cout << "\n\t\t sparse::trsv residual:\n"
              << "\t\t\t" << diff << "\n\tFinished" << std::endl;

    cleanup_arrays<fp, intType>(fp_ptr_vec, int_ptr_vec, cxt);

    return 0;
}

//
// Description of example setup, apis used and supported floating point type
// precisions
//
void print_example_banner()
{

    std::cout << "" << std::endl;
    std::cout << "###############################################################"
                 "#########"
              << std::endl;
    std::cout << "# Sparse Triangular Solve Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# y = op(A)^(-1) * x" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A is a sparse matrix in CSR format, x and y are "
                 "dense vectors"
              << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   sparse::trsv" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;
    std::cout << "#   float" << std::endl;
    std::cout << "#   double" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "###############################################################"
                 "#########"
              << std::endl;
    std::cout << std::endl;
}

//
// Main entry point for example.
//
// Dispatches to appropriate device types as set at build time with flag:
// -DSYCL_DEVICES_host -- only runs host implementation
// -DSYCL_DEVICES_cpu -- only runs SYCL CPU implementation
// -DSYCL_DEVICES_gpu -- only runs SYCL GPU implementation
// -DSYCL_DEVICES_all (default) -- runs on all: host, cpu and gpu devices
//
//  For each device selected and each supported data type, MatrixMultiplyExample
//  runs is with all supported data types
//

int main(int argc, char **argv)
{

    print_example_banner();

    std::list<my_sycl_device_types> list_of_devices;
    set_list_of_devices(list_of_devices);

    int status = 0;
    for (auto it = list_of_devices.begin(); it != list_of_devices.end(); ++it) {

        cl::sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device(my_dev, my_dev_is_found, *it);

        if (my_dev_is_found) {
            std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";

            std::cout << "\tRunning with single precision real data type:" << std::endl;
            status = run_sparse_triangular_solve_example<float, std::int32_t>(my_dev);
            if(status != 0) return status;

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                status = run_sparse_triangular_solve_example<double, std::int32_t>(my_dev);
                if(status != 0) return status;
            }
        }
        else {
#ifdef FAIL_ON_MISSING_DEVICES
            std::cout << "No " << sycl_device_names[*it]
                      << " devices found; Fail on missing devices "
                         "is enabled.\n";
            return 1;
#else
            std::cout << "No " << sycl_device_names[*it] << " devices found; skipping "
                      << sycl_device_names[*it] << " tests.\n";
#endif
        }
    }

    return status;
}
