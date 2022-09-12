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
*       DPCPP SPARSE BLAS API to solve a system of linear equations (Ax=b)
*       by preconditioned CG method with symmetric Gauss-Seidel preconditioner:
*       Compute r_0 = b - Ax_0
*       w_0 = B^{-1}*r_0 and p_0 = w_0
*       while not converged
*           {
*                   alpha_k = (r_k , w_k )/(Ap_k , p_k )
*                   x_{k+1} = x_k + alpha_k*p_k
*                   r_{k+1} = r_k - alpha_k*A*p_k
*                   w_{k+1} = B^{-1}*r_{k+1}
*                   beta_k = (r_{k+1}, w_{k+1})/(r_k , w_k )
*                   p_{k+1} = w_{k+1} + beta_k*p_k
*           }
*
*       where A = -L+D-L^t; B = (D-L)*D^{-1}*(D-L^t).
*
*       The supported floating point data types for gemm matrix data are:
*           float
*           double
*
*
*/

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

template <typename fp, typename intType>
class extractDiagonalClass;
template <typename fp, typename intType>
class diagonalMVClass;

template <typename fp, typename intType>
static void diagonal_mv(cl::sycl::queue main_queue,
                        const intType nrows,
                        cl::sycl::buffer<fp, 1> &d_buffer,
                        cl::sycl::buffer<fp, 1> &t_buffer)
{
    main_queue.submit([&](cl::sycl::handler &cgh) {
        auto d = (d_buffer).template get_access<cl::sycl::access::mode::write>(cgh);
        auto t = (t_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
        auto diagonalMVKernel = [=](cl::sycl::item<1> item) {
            const int row = item.get_id(0);
            t[row] *= d[row];
        };
        cgh.parallel_for<class diagonalMVClass<fp, intType>>(cl::sycl::range<1>(nrows),
                                                             diagonalMVKernel);
    });
}

template <typename fp, typename intType>
int run_sparse_cg_example(const cl::sycl::device &dev)
{
    // Matrix data size
    intType size  = 4;
    intType nrows = size * size * size;

    // Input matrix in CSR format
    std::vector<intType, mkl_allocator<intType, 64>> ia;
    std::vector<intType, mkl_allocator<intType, 64>> ja;
    std::vector<fp, mkl_allocator<fp, 64>> a;

    ia.resize(nrows + 1);
    ja.resize(27 * nrows);
    a.resize(27 * nrows);

    generate_sparse_matrix<fp, intType>(size, ia, ja, a);

    // Vectors x and y
    std::vector<fp, mkl_allocator<fp, 64>> x;
    std::vector<fp, mkl_allocator<fp, 64>> b;
    x.resize(nrows);
    b.resize(nrows);

    // Init right hand side and vector x
    for (int i = 0; i < nrows; i++) {
        b[i] = set_fp_value(fp(1.0), fp(0.0));
        x[i] = set_fp_value(fp(0.0), fp(0.0));
    }

    // Catch asynchronous exceptions
    auto exception_handler = [](cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const &e : exceptions) {
            try {
                std::rethrow_exception(e);
            }
            catch (cl::sycl::exception const &e) {
                std::cout << "Caught asynchronous SYCL "
                             "exception during sparse CG:\n"
                          << e.what() << std::endl;
            }
        }
    };

    //
    // Execute CG
    //

    // create execution queue and buffers of matrix data
    cl::sycl::queue main_queue(dev, exception_handler);

    cl::sycl::buffer<intType, 1> ia_buffer(ia.data(), cl::sycl::range<1>(nrows + 1));
    cl::sycl::buffer<intType, 1> ja_buffer(ja.data(), cl::sycl::range<1>(ia[nrows]));
    cl::sycl::buffer<fp, 1> a_buffer(a.data(), cl::sycl::range<1>(ia[nrows]));
    cl::sycl::buffer<fp, 1> x_buffer(x.data(), cl::sycl::range<1>(nrows));
    cl::sycl::buffer<fp, 1> b_buffer(b.data(), cl::sycl::range<1>(nrows));
    cl::sycl::buffer<fp, 1> r_buffer((cl::sycl::range<1>(nrows)));
    cl::sycl::buffer<fp, 1> w_buffer((cl::sycl::range<1>(nrows)));
    cl::sycl::buffer<fp, 1> p_buffer((cl::sycl::range<1>(nrows)));
    cl::sycl::buffer<fp, 1> t_buffer((cl::sycl::range<1>(nrows)));
    cl::sycl::buffer<fp, 1> y_buffer((cl::sycl::range<1>(nrows)));
    cl::sycl::buffer<fp, 1> d_buffer((cl::sycl::range<1>(nrows)));
    cl::sycl::buffer<fp, 1> temp_buffer((cl::sycl::range<1>(1)));

    // create and initialize handle for a Sparse Matrix in CSR format
    oneapi::mkl::sparse::matrix_handle_t handle;

    try {
        oneapi::mkl::sparse::init_matrix_handle(&handle);

        oneapi::mkl::sparse::set_csr_data(handle, nrows, nrows, oneapi::mkl::index_base::zero,
                                          ia_buffer, ja_buffer, a_buffer);

        oneapi::mkl::sparse::set_matrix_property(handle, oneapi::mkl::sparse::property::symmetric);
        oneapi::mkl::sparse::set_matrix_property(handle, oneapi::mkl::sparse::property::sorted);

        oneapi::mkl::sparse::optimize_trsv(main_queue, oneapi::mkl::uplo::lower,
                                           oneapi::mkl::transpose::nontrans,
                                           oneapi::mkl::diag::nonunit, handle);
        oneapi::mkl::sparse::optimize_trsv(main_queue, oneapi::mkl::uplo::upper,
                                           oneapi::mkl::transpose::nontrans,
                                           oneapi::mkl::diag::nonunit, handle);
        oneapi::mkl::sparse::optimize_gemv(main_queue, oneapi::mkl::transpose::nontrans, handle);

        main_queue.submit([&](cl::sycl::handler &cgh) {
            auto ia = (ia_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ja = (ja_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
            auto a  = (a_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
            auto d  = (d_buffer).template get_access<cl::sycl::access::mode::write>(cgh);
            auto extractDiagonalKernel = [=](cl::sycl::item<1> item) {
                const int row = item.get_id(0);
                for (intType i = ia[row]; i < ia[row + 1]; i++) {
                    if (ja[i] == row) {
                        d[row] = a[i];
                        break;
                    }
                }
            };
            cgh.parallel_for<class extractDiagonalClass<fp, intType>>(cl::sycl::range<1>(nrows),
                                                                      extractDiagonalKernel);
        });

        // initial residual equal to RHS cause of zero initial vector
        oneapi::mkl::blas::copy(main_queue, nrows, b_buffer, 1, r_buffer, 1);

        // Calculation B^{-1}r_0
        {
            oneapi::mkl::sparse::trsv(main_queue, oneapi::mkl::uplo::lower,
                                      oneapi::mkl::transpose::nontrans, oneapi::mkl::diag::nonunit,
                                      handle, r_buffer, t_buffer);
            diagonal_mv<fp, intType>(main_queue, nrows, d_buffer, t_buffer);
            oneapi::mkl::sparse::trsv(main_queue, oneapi::mkl::uplo::upper,
                                      oneapi::mkl::transpose::nontrans, oneapi::mkl::diag::nonunit,
                                      handle, t_buffer, w_buffer);
        }

        oneapi::mkl::blas::copy(main_queue, nrows, w_buffer, 1, p_buffer, 1);

        // Calculate initial norm of correction
        fp initial_norm_of_correction = set_fp_value(fp(0.0), fp(0.0));
        oneapi::mkl::blas::nrm2(main_queue, nrows, w_buffer, 1, temp_buffer);
        {
            auto temp_accessor = temp_buffer.template get_access<cl::sycl::access::mode::read>();
            initial_norm_of_correction = temp_accessor[0];
        }
        fp norm_of_correction = initial_norm_of_correction;

        // Start of main PCG algorithm
        std::int32_t k = 0;
        fp alpha, beta, temp;

        oneapi::mkl::blas::dot(main_queue, nrows, r_buffer, 1, w_buffer, 1, temp_buffer);
        {
            auto temp_accessor = temp_buffer.template get_access<cl::sycl::access::mode::read>();
            temp               = temp_accessor[0];
        }

        while (norm_of_correction / initial_norm_of_correction > 1.e-3 && k < 100) {
            // Calculate A*p
            oneapi::mkl::sparse::gemv(main_queue, oneapi::mkl::transpose::nontrans, 1.0, handle,
                                      p_buffer, 0.0, t_buffer);

            // Calculate alpha_k
            oneapi::mkl::blas::dot(main_queue, nrows, p_buffer, 1, t_buffer, 1, temp_buffer);
            {
                auto temp_accessor =
                        temp_buffer.template get_access<cl::sycl::access::mode::read>();
                alpha = temp / temp_accessor[0];
            }

            // Calculate x_k = x_k + alpha*p_k
            oneapi::mkl::blas::axpy(main_queue, nrows, alpha, p_buffer, 1, x_buffer, 1);
            // Calculate r_k = r_k - alpha*A*p_k
            oneapi::mkl::sparse::gemv(main_queue, oneapi::mkl::transpose::nontrans, -alpha, handle,
                                      p_buffer, 1.0, r_buffer);

            // Calculate w_k = B^{-1}r_k
            {
                oneapi::mkl::sparse::trsv(main_queue, oneapi::mkl::uplo::lower,
                                          oneapi::mkl::transpose::nontrans,
                                          oneapi::mkl::diag::nonunit, handle, r_buffer, t_buffer);
                diagonal_mv<fp, intType>(main_queue, nrows, d_buffer, t_buffer);
                oneapi::mkl::sparse::trsv(main_queue, oneapi::mkl::uplo::upper,
                                          oneapi::mkl::transpose::nontrans,
                                          oneapi::mkl::diag::nonunit, handle, t_buffer, w_buffer);
            }

            // Calculate current norm of correction
            oneapi::mkl::blas::nrm2(main_queue, nrows, w_buffer, 1, temp_buffer);
            {
                auto temp_accessor =
                        temp_buffer.template get_access<cl::sycl::access::mode::read>();
                norm_of_correction = temp_accessor[0];
            }
            std::cout << "\t\trelative norm of residual on " << ++k
                      << " iteration: " << norm_of_correction / initial_norm_of_correction
                      << std::endl;
            if (norm_of_correction <= 1.e-3)
                break;

            // Calculate beta_k
            oneapi::mkl::blas::dot(main_queue, nrows, r_buffer, 1, w_buffer, 1, temp_buffer);
            {
                auto temp_accessor =
                        temp_buffer.template get_access<cl::sycl::access::mode::read>();
                beta = temp_accessor[0] / temp;
                temp = temp_accessor[0];
            }

            // Calculate p_k = w_k+beta*p_k
            oneapi::mkl::blas::axpy(main_queue, nrows, beta, p_buffer, 1, w_buffer, 1);
            oneapi::mkl::blas::copy(main_queue, nrows, w_buffer, 1, p_buffer, 1);
        }

        std::cout << "" << std::endl;
        std::cout << "\t\tPreconditioned CG process has successfully converged and following "
                     "solution has been obtained";
        std::cout << "" << std::endl;

        auto result = x_buffer.template get_access<cl::sycl::access::mode::read>();
        for (std::int32_t i = 0; i < 4; i++) {
            std::cout << "\t\tx[" << i << "] = " << result[i] << std::endl;
        }
        std::cout << "\t\t..." << std::endl;

        oneapi::mkl::sparse::release_matrix_handle(&handle);
    }
    catch (cl::sycl::exception const &e) {
        std::cout << "\t\tCaught synchronous SYCL exception:\n" << e.what() << std::endl;
        oneapi::mkl::sparse::release_matrix_handle(&handle);
        return 1;
    }
    catch (std::exception const &e) {
        std::cout << "\t\tCaught std exception:\n" << e.what() << std::endl;
        oneapi::mkl::sparse::release_matrix_handle(&handle);
        return 1;
    }

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
    std::cout << "# Sparse CG Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# A * x = b" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A is a sparse matrix in CSR format, x and b are "
                 "dense vectors"
              << std::endl;
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
//  For each device selected and each supported data type, CG example
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
            status = run_sparse_cg_example<float, std::int32_t>(my_dev);
            if(status != 0) return status;

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                status = run_sparse_cg_example<double, std::int32_t>(my_dev);
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
