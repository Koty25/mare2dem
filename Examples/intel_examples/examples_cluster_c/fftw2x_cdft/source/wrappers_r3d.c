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
!    Forward-Backward 3D real-to-complex transform for double/single
!    precision data.
!
!*****************************************************************************
! Configuration parameters:
!    DFTI_FORWARD_DOMAIN = DFTI_REAL               (obligatory)
!    DFTI_PRECISION      = DFTI_DOUBLE/DFTI_SINGLE (obligatory)
!    DFTI_DIMENSION      = 3                       (obligatory)
!    DFTI_LENGTHS        = {lenx,leny,lenz}        (obligatory)
!    DFTI_FORWARD_SCALE  = 1.0                     (default)
!    DFTI_BACKWARD_SCALE = 1.0/(lenx*leny*lenz)    (default=1.0)
!
!****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "rfftw_mpi.h"

#ifdef FFTW_ENABLE_FLOAT
    #define EPS 1.0E-6
#else
    #define EPS 1.0E-12
#endif

#define PRINT(...) do { if (PrR == 0) printf(__VA_ARGS__); } while(0)

int usage_of_mpi_rfftw3d( int PrR, int lenx, int leny, int lenz, int work_is_null)
{
    fftw_real *local = NULL;
    fftw_real *work = NULL;
    fftw_real *preal,t,err,scale;
    rfftwnd_mpi_plan plan;

    int lenz_padded = 2*(lenz/2+1);
    int i,j,k;
    int fail_status = 1;
    int nx,nx_out,start_x,start_x_out,size;

    if(work_is_null)
        PRINT("Usage of MPI FFTW for Real-to-Complex Multi-dimensional Transforms (work==NULL).\nPlan = rfftw3d_mpi_create_plan(...)\n");
    else
        PRINT("Usage of MPI FFTW for Real-to-Complex Multi-dimensional Transforms.\nPlan = rfftw3d_mpi_create_plan(...)\n");
    plan=rfftw3d_mpi_create_plan(MPI_COMM_WORLD, lenx, leny, lenz, FFTW_FORWARD, FFTW_ESTIMATE);
    if(plan == NULL)
        goto error;
    PRINT("     CreatePlan for Forward....DONE\n");
    rfftwnd_mpi_local_sizes(plan,&nx,&start_x,&nx_out,&start_x_out,&size);
    PRINT("     LocalSizes................DONE\n");
    local=(fftw_real*)fftw_malloc(size*sizeof(fftw_complex));
    if(local == NULL)
        goto error;
    if(work_is_null){
        work = NULL;
    } else{
        work=(fftw_real*)fftw_malloc(size*sizeof(fftw_complex));
        if(work == NULL)
            goto error;
    }
    preal=(fftw_real*)local;
    srand(PrR*100);
    for (i=0;i<nx;i++)
        for (j=0;j<leny;j++)
            for(k=0;k<lenz;k++)
                preal[k+(j+i*leny)*lenz_padded]=(fftw_real)rand()/RAND_MAX;

    rfftwnd_mpi(plan,1,local,work,FFTW_NORMAL_ORDER);
    PRINT("     Forward FFTW..............DONE\n");

    rfftwnd_mpi_destroy_plan(plan);
    PRINT("     Destroy plan for Forward..DONE\n");
    plan=rfftw3d_mpi_create_plan(MPI_COMM_WORLD, lenx, leny, lenz, FFTW_BACKWARD, FFTW_ESTIMATE);
    if(plan == NULL)
        goto error;
    PRINT("     CreatePlan for Backward...DONE\n");
    rfftwnd_mpi(plan,1,local,work,FFTW_NORMAL_ORDER);
    scale=1.0/lenx/leny/lenz;
    PRINT("     Backward FFTW.............DONE\n");
     for (i=0;i<nx;i++)
        for (j=0;j<leny;j++)
            for (k=0;k<lenz;k++)
                preal[k+(j+i*leny)*lenz_padded]*=scale;
    rfftwnd_mpi_destroy_plan(plan);
    PRINT("     Destroy plan for Backward.DONE\n");
    srand(PrR*100);
    err=0.0;
    for (i=0;i<nx;i++)
        for (j=0;j<leny;j++)
            for (k=0;k<lenz;k++) {
                t=fabs(preal[k+(j+i*leny)*lenz_padded]-(fftw_real)rand()/RAND_MAX);
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
    int lenx, leny, lenz;
    int failure = 0;

    lenx=128; leny=64; lenz=32;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&PrS);
    MPI_Comm_rank(MPI_COMM_WORLD,&PrR);
    PRINT("PrS=%d, lengths={%d,%d,%d}\n",PrS,lenx,leny,lenz);

/*--- Usage of MPI FFTW for Real-to-Complex One-dimensional Transforms (work==NULL) ---*/
    failure +=  usage_of_mpi_rfftw3d(PrR, lenx, leny, lenz, 1);

/*--- Usage of MPI FFTW for Real-to-Complex One-dimensional Transforms ---*/
    failure +=  usage_of_mpi_rfftw3d(PrR, lenx, leny, lenz, 0);

    PRINT(" END OF TEST\n");

    MPI_Finalize();
    return (failure != 0);
}
