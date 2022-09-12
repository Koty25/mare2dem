/*******************************************************************************
* Copyright 2006-2020 Intel Corporation.
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
!    Intel(R) Math Kernel Library (Intel(R) MKL) Cluster DFT using FFTW
!    interface (via wrappers) example program (C-interface)
!
!    Forward-Backward 1D complex transform for double/single precision data.
!
!*****************************************************************************
! Configuration parameters:
!    DFTI_FORWARD_DOMAIN = DFTI_COMPLEX            (obligatory)
!    DFTI_PRECISION      = DFTI_DOUBLE/DFTI_SINGLE (obligatory)
!    DFTI_DIMENSION      = 1                       (obligatory)
!    DFTI_LENGTHS        = {len}                   (obligatory)
!    DFTI_FORWARD_SCALE  = 1.0                     (default)
!    DFTI_BACKWARD_SCALE = 1.0/(len)               (default=1.0)
!
!****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fftw_mpi.h"

#ifdef FFTW_ENABLE_FLOAT
    #define EPS 1.0E-6
#else
    #define EPS 1.0E-12
#endif

#define PRINT(...) do { if (PrR == 0) printf(__VA_ARGS__); } while(0)

int usage_of_mpi_fftw( int PrR, int len, int work_is_null)
{
    fftw_complex *local = NULL;
    fftw_complex *work = NULL;
    fftw_real *preal,t,err,scale;
    fftw_mpi_plan plan;

    int i;
    int fail_status = 1;
    int nx,nx_out,start_x,start_x_out,size;

    if(work_is_null)
        PRINT("Usage of MPI FFTW for Complex One-dimensional Transforms (work==NULL)\n");
    else
        PRINT("Usage of MPI FFTW for Complex One-dimensional Transforms\n");
    plan=fftw_mpi_create_plan(MPI_COMM_WORLD,len,FFTW_FORWARD,FFTW_ESTIMATE);
    if(plan == NULL)
        goto error;
    PRINT("     CreatePlan for Forward....DONE\n");
    fftw_mpi_local_sizes(plan,&nx,&start_x,&nx_out,&start_x_out,&size);
    PRINT("     LocalSizes................DONE\n");
    local=(fftw_complex*)fftw_malloc(size*sizeof(fftw_complex));
    if(local == NULL)
        goto error;
    if(work_is_null){
        work = NULL;
    } else{
        work=(fftw_complex*)fftw_malloc(size*sizeof(fftw_complex));
        if(work == NULL)
            goto error;
    }
    preal=(fftw_real*)local;
    srand(PrR*100);
    for (i=0;i<2*nx;i++)
        preal[i]=(fftw_real)rand()/RAND_MAX;

    fftw_mpi(plan,1,local,work);
    PRINT("     Forward FFTW..............DONE\n");

    fftw_mpi_destroy_plan(plan);
    PRINT("     Destroy plan for Forward..DONE\n");
    plan=fftw_mpi_create_plan(MPI_COMM_WORLD,len,FFTW_BACKWARD,FFTW_ESTIMATE);
    if(plan == NULL)
        goto error;
    PRINT("     CreatePlan for Backward...DONE\n");
    fftw_mpi(plan,1,local,work);
    scale=1.0/len;
    PRINT("     Backward FFTW.............DONE\n");
    for (i=0;i<2*nx;i++)
        preal[i]*=scale;
    fftw_mpi_destroy_plan(plan);
    PRINT("     Destroy plan for Backward.DONE\n");
    srand(PrR*100);
    err=0.0;
    for (i=0;i<2*nx;i++) {
        t=fabs(preal[i]-(fftw_real)rand()/RAND_MAX);
        if (t>err) err=t;
    }

    fail_status = err > EPS;
    PRINT("     Error=%e\n",err);

error:
    PRINT(" TEST %s\n", fail_status ? "FAILED" : "PASSED");
    if (local) fftw_free(local);
    if (work) fftw_free(work);

    return fail_status;
}

int main(int argc,char *argv[])
{
    int PrS,PrR;
    int failure = 0;
    int len;

    len = 256;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&PrS);
    MPI_Comm_rank(MPI_COMM_WORLD,&PrR);

    PRINT("PrS=%d, lengths={%d}\n",PrS,len);

/*--- Usage of MPI FFTW for Complex One-dimensional Transforms (work==NULL) ---*/
    failure += usage_of_mpi_fftw(PrR, len, 1);

/*--- Usage of MPI FFTW for Complex One-dimensional Transforms ---*/
    failure += usage_of_mpi_fftw(PrR, len, 0);

    PRINT(" END OF TEST\n");

    MPI_Finalize();
    return (failure != 0);
}
