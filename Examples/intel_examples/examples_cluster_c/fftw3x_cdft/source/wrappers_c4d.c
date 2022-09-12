/*******************************************************************************
* Copyright 2010-2020 Intel Corporation.
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
! Content:
!    Intel(R) Math Kernel Library (Intel(R) MKL) Cluster DFT using MPI FFTW3
!    interface (via wrappers) example program (C-interface)
!
!    Forward-Backward 4D complex transform for double precision data.
!
!****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fftw3-mpi.h"

#define EPS (1.0E-12)

#define REAL_PTR(x) ((double *)(x))

const char *messages[] = {"PASSED", "FAILED"};

int main(int argc, char *argv[])
{
    fftw_complex *src = NULL, *dst = NULL;
    double err, g_err, scale;
    int comm_size, comm_rank;
    int i, ret;
    fftw_plan plan = NULL;
    ptrdiff_t n[] = {64, 32, 80, 20};
    ptrdiff_t size, data_size;
    ptrdiff_t local_n0, local_0_start;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    fftw_mpi_init();

    if (0 == comm_rank) {
        printf("Number of MPI processes = %d, data cube dimensions = %d x %d x %d x %d\n", 
               comm_size, (int)n[0], (int)n[1], (int)n[2], (int)n[3]);
        printf("Usage of MPI FFTW for Complex Multi-dimensional Transforms\n");
    }

    size = fftw_mpi_local_size(4, n, MPI_COMM_WORLD, &local_n0, &local_0_start);
    data_size = local_n0*n[1]*n[2]*n[3];

    src = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
    dst = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
    if ((NULL == src) || (NULL == dst)) {
        ret = 1;
        printf("Failed to allocate memory for buffers\n");
        goto cleanup;
    }

    srand(1);
    for (i=0; i<2*data_size; i++)
        REAL_PTR(src)[i]=(double)rand()/RAND_MAX;

    plan = fftw_mpi_plan_dft(4, n, src, dst, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
    if (NULL == plan) {
        ret = 1;
        printf("Failed to create plan for forward FFT\n");
        goto cleanup;
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    plan = fftw_mpi_plan_dft(4, n, dst, src, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (NULL == plan) {
        ret = 1;
        printf("Failed to create plan for forward FFT\n");
        goto cleanup;
    }
    fftw_execute(plan);
    scale = 1.0/n[0]/n[1]/n[2]/n[3];
    for (i=0; i<2*data_size; i++)
        REAL_PTR(src)[i]*=scale;
    fftw_destroy_plan(plan);

    srand(1);
    err = 0.0;
    for (i=0; i<2*data_size; i++) {
        double t = fabs(REAL_PTR(src)[i]-(double)rand()/RAND_MAX);
        if (t>err) err=t;
    }
    MPI_Allreduce(&err, &g_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (0 == comm_rank)
        printf("Error=%e\n", g_err);

    ret = (g_err<EPS) ? 0 : 1;

cleanup:

    if (NULL != src) fftw_free(src);
    if (NULL != dst) fftw_free(dst);
    fftw_mpi_cleanup();

    MPI_Finalize();

    if (0 == comm_rank)
        printf("TEST %s\n", messages[ret]);

    return ret;
}

