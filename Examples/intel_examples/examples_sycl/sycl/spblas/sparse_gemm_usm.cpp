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
*       This example demonstrates use of oneAPI Math Kernel Library (oneMKL)
*       DPCPP USM API oneapi::mkl::sparse::gemm to perform general sparse matrix-matrix
*       multiplication on a SYCL device (Host, CPU, GPU).
*
*       c = alpha * op(A) * b + beta * c
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
// Main example for Sparse Matrix-Dense Matrix Multiply consisting of
// initialization of A matrix, x and y vectors as well as
// scalars alpha and beta.  Then the product
//
// c = alpha * op(A) * b + beta * c
//
// is performed and finally the results are post processed.
//
template <typename fp, typename intType>
int run_sparse_matrix_dense_matrix_multiply_example(const cl::sycl::device &dev)
{
    // Initialize data for Sparse Matrix-Vector Multiply
    oneapi::mkl::transpose transpose_val = oneapi::mkl::transpose::nontrans;

    // Matrix data size
    intType nrows        = 64;
    intType ncols        = 64;
    std::int64_t columns = 64;
    std::int64_t ldb     = columns;
    std::int64_t ldc     = columns;

    double density_val = 0.05;

    // Input matrix in CSR format
    std::vector<intType, mkl_allocator<intType, 64>> ia_vec;
    std::vector<intType, mkl_allocator<intType, 64>> ja_vec;
    std::vector<fp, mkl_allocator<fp, 64>> a_vec;

    generate_random_sparse_matrix<fp, intType>(nrows, ncols, density_val, ia_vec, ja_vec, a_vec);

    // Catch asynchronous exceptions
    auto exception_handler = [](cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const &e : exceptions) {
            try {
                std::rethrow_exception(e);
            }
            catch (cl::sycl::exception const &e) {
                std::cout << "Caught asynchronous SYCL "
                             "exception during sparse::gemm:\n"
                          << e.what() << std::endl;
            }
        }
    };

    // create execution queue and buffers of matrix data
    cl::sycl::event gemm_done = cl::sycl::event();

    cl::sycl::queue main_queue(dev, exception_handler);
    cl::sycl::context cxt = main_queue.get_context();

    intType *ia, *ja;
    fp *a, *b, *c, *d;
    intType sizea  = a_vec.size();
    intType sizeja = ja_vec.size();
    intType sizeia = ia_vec.size();

    ia = (intType *)malloc_shared(sizeia * sizeof(intType), dev, cxt);
    ja = (intType *)malloc_shared(sizeja * sizeof(intType), dev, cxt);
    a  = (fp *)malloc_shared(sizea * sizeof(fp), dev, cxt);
    b  = (fp *)malloc_shared(ncols * ldb * sizeof(fp), dev, cxt);
    c  = (fp *)malloc_shared(nrows * ldc * sizeof(fp), dev, cxt);
    d  = (fp *)malloc_shared(nrows * ldc * sizeof(fp), dev, cxt);

    if (!ia || !ja || !a || !b || !c || !d)
        throw std::runtime_error("Failed to allocate USM memory");

    intType nrows_b = ncols;

    rand_matrix<fp>(b, oneapi::mkl::transpose::nontrans, columns, nrows_b, ldb);

    for (int i = 0; i < ia_vec.size(); i++)
        ia[i]  = ia_vec[i];

    for (int i = 0; i < ja_vec.size(); i++) {
        ja[i] = ja_vec[i];
        a[i]  = a_vec[i];
    }

    // Init matrices c and d
    for (int i = 0; i < nrows * ldc; i++) {
        c[i] = set_fp_value(fp(0.0), fp(0.0));
        d[i] = set_fp_value(fp(0.0), fp(0.0));
    }

    // Set scalar fp values
    fp alpha = set_fp_value(fp(1.0), fp(0.0));
    fp beta  = set_fp_value(fp(0.0), fp(0.0));

    std::vector<intType *> int_ptr_vec;
    int_ptr_vec.push_back(ia);
    int_ptr_vec.push_back(ja);
    std::vector<fp *> fp_ptr_vec;
    fp_ptr_vec.push_back(a);
    fp_ptr_vec.push_back(b);
    fp_ptr_vec.push_back(c);
    fp_ptr_vec.push_back(d);

    //
    // Execute Matrix Multiply
    //

    // create and initialize handle for a Sparse Matrix in CSR format
    oneapi::mkl::sparse::matrix_handle_t handle;

    try {
        oneapi::mkl::sparse::init_matrix_handle(&handle);

        oneapi::mkl::sparse::set_csr_data(handle, nrows, ncols, oneapi::mkl::index_base::zero, ia,
                                          ja, a);

        // add oneapi::mkl::sparse::gemm to execution queue
        gemm_done = oneapi::mkl::sparse::gemm(main_queue, transpose_val, alpha, handle, b, columns,
                                              ldb, beta, c, ldc, {});

        oneapi::mkl::sparse::release_matrix_handle(&handle, {gemm_done});
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

    std::cout << "\n\t\tsparse::gemm parameters:\n";
    std::cout << "\t\t\ttranspose_val = "
              << (transpose_val == oneapi::mkl::transpose::nontrans ?
                          "nontrans" :
                          (transpose_val == oneapi::mkl::transpose::trans ? "trans" : "conjtrans"))
              << std::endl;
    std::cout << "\t\t\tnrows = " << nrows << std::endl;
    std::cout << "\t\t\tncols = " << ncols << std::endl;
    std::cout << "\t\t\tcolumns = " << columns << std::endl;
    std::cout << "\t\t\tldb = " << ldb << ", ldc = " << ldc << std::endl;
    std::cout << "\t\t\talpha = " << alpha << ", beta = " << beta << std::endl;

    fp *res = c;
    for (intType row = 0; row < nrows; row++) {
        for (intType col = 0; col < columns; col++) {
            intType index = row * ldc + col;

            if (beta == (fp)0) {
                d[index] = set_fp_value(fp(0.0), fp(0.0));
            }

            fp tmp = set_fp_value(fp(0.0), fp(0.0));
            for (intType i = ia[row]; i < ia[row + 1]; i++) {
                tmp += a[i] * b[col + ja[i] * ldb];
            }
            d[index] = alpha * tmp + beta * d[index];
        }
    }

    fp diff  = set_fp_value(fp(0.0), fp(0.0));
    fp diff2 = set_fp_value(fp(0.0), fp(0.0));
    for (intType i = 0; i < nrows * columns; i++) {
        if(!check_result(res[i], d[i], i)) return 1;
        diff += (d[i] - res[i]) * (d[i] - res[i]);
        diff2 += d[i] * d[i];
    }

    std::cout << "\n\t\t sparse::gemm residual:\n"
              << "\t\t\t" << diff / diff2 << "\n\tFinished" << std::endl;

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
    std::cout << "# Sparse Matrix-Dense Matrix Multiply Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# c = alpha * op(A) * b + beta * c" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A is a sparse matrix in CSR format, b and c are "
                 "dense matrices"
              << std::endl;
    std::cout << "# and alpha, beta are floating point type precision scalars." << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   sparse::gemm" << std::endl;
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
            status = run_sparse_matrix_dense_matrix_multiply_example<float, std::int32_t>(my_dev);
            if(status != 0) return status;

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                status = run_sparse_matrix_dense_matrix_multiply_example<double, std::int32_t>(my_dev);
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
