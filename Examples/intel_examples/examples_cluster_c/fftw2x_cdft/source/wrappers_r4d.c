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
!    Forward-Backward 4D real-to-complex transform for double/single
!    precision data.
!
!*****************************************************************************
! Configuration parameters:
!    DFTI_FORWARD_DOMAIN = DFTI_REAL                         (obligatory)
!    DFTI_PRECISION      = DFTI_DOUBLE/DFTI_SINGLE           (obligatory)
!    DFTI_DIMENSION      = 4                                 (obligatory)
!    DFTI_LENGTHS        = {len[0],len[1],len[2],len[3]}     (obligatory)
!    DFTI_FORWARD_SCALE  = 1.0                               (default)
!    DFTI_BACKWARD_SCALE = 1.0/(len[0]*len[1]*len[2]*len[3]) (default=1.0)
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

int usage_of_mpi_rfftw4d( int PrR, int len[], int work_is_null)
{
    fftw_real *local = NULL;
    fftw_real *work = NULL;
    fftw_real *preal,t,err,scale;
    rfftwnd_mpi_plan plan;

    int len3_padded = 2*(len[3]/2+1);
    int i0, i1, i2, i3;
    int fail_status = 1;
    int nx,nx_out,start_x,start_x_out,size;

    if(work_is_null)
        PRINT("Usage of MPI FFTW for Real-to-Complex Multi-dimensional Transforms (work==NULL) n=4.\nPlan = rfftwnd_mpi_create_plan(...)\n");
    else
        PRINT("Usage of MPI FFTW for Real-to-Complex Multi-dimensional Transforms n=4.\nPlan = rfftwnd_mpi_create_plan(...)\n");
    plan=rfftwnd_mpi_create_plan(MPI_COMM_WORLD, 4, len, FFTW_FORWARD, FFTW_ESTIMATE);
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
    for (i0=0;i0<nx;i0++)
        for (i1=0;i1<len[1];i1++)
            for (i2=0;i2<len[2];i2++)
                for (i3=0;i3<len[3];i3++)
                    preal[i3+(i2+i1*len[2]+i0*len[1]*len[2])*len3_padded]=(fftw_real)rand()/RAND_MAX;

    rfftwnd_mpi(plan,1,local,work,FFTW_NORMAL_ORDER);
    PRINT("     Forward FFTW..............DONE\n");

    rfftwnd_mpi_destroy_plan(plan);
    PRINT("     Destroy plan for Forward..DONE\n");

    plan=rfftwnd_mpi_create_plan(MPI_COMM_WORLD,4,len,FFTW_BACKWARD,FFTW_ESTIMATE);
    if(plan == NULL)
        goto error;
    PRINT("     CreatePlan for Backward...DONE\n");
    rfftwnd_mpi(plan,1,local,work,FFTW_NORMAL_ORDER);
    scale=1.0/len[1]/len[2]/len[3]/len[0];
    PRINT("     Backward FFTW.............DONE\n");
    for (i0=0;i0<nx;i0++)
        for (i1=0;i1<len[1];i1++)
            for (i2=0;i2<len[2];i2++)
                for (i3=0;i3<len[3];i3++)
                    preal[i3+(i2+i1*len[2]+i0*len[1]*len[2])*len3_padded]*=scale;

    rfftwnd_mpi_destroy_plan(plan);
    PRINT("     Destroy plan for Backward.DONE\n");
    srand(PrR*100);
    err=0.0;

    for (i0=0;i0<nx;i0++)
        for (i1=0;i1<len[1];i1++)
            for (i2=0;i2<len[2];i2++)
                for (i3=0;i3<len[3];i3++) {
                    t=fabs(preal[i3+(i2+i1*len[2]+i0*len[1]*len[2])*len3_padded]-(fftw_real)rand()/RAND_MAX);
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
    int len[4];
    int failure = 0;

    len[0]=64;
    len[1]=32;
    len[2]=80;
    len[3]=20;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&PrS);
    MPI_Comm_rank(MPI_COMM_WORLD,&PrR);
    PRINT("PrS=%d, lengths={%d,%d,%d,%d}\n",PrS,len[0],len[1],len[2],len[3]);

/*--- Usage of MPI FFTW for Real-to-Complex One-dimensional Transforms (work==NULL) ---*/
    failure +=  usage_of_mpi_rfftw4d(PrR, len, 1);

/*--- Usage of MPI FFTW for Real-to-Complex One-dimensional Transforms ---*/
    failure +=  usage_of_mpi_rfftw4d(PrR, len, 0);

    PRINT(" END OF TEST\n");

    MPI_Finalize();
    return (failure != 0);
}
