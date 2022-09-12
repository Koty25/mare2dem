/*******************************************************************************
* Copyright 2003-2020 Intel Corporation.
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
!     Intel(R) Math Kernel Library (Intel(R) MKL) Cluster DFT interface
!     example program (C-interface)
!
!     Forward-Backward 2D complex transform for single precision data inplace.
!
!*****************************************************************************
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX      (obligatory)
!           DFTI_PRECISION      = DFTI_SINGLE       (obligatory)
!           DFTI_DIMENSION      = 2                 (obligatory)
!           DFTI_LENGTHS        = {m,n}             (obligatory)
!           DFTI_FORWARD_SCALE  = 1.0               (default)
!           DFTI_BACKWARD_SCALE = 1.0/(m*n)         (default=1.0)
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "mkl_service.h"
#include "mkl_cdft.h"
#include "mkl_cdft_examples.h"

#define PREC DFTI_SINGLE

int main(int argc, char *argv[])  /* DM_COMPLEX_2D_SINGLE_EX1 */
{
    mkl_float_complex *x_in;
    mkl_float_complex *x_exp;
    mkl_float_complex *local,*work;

    DFTI_DESCRIPTOR_DM_HANDLE desc;
    long RootRank,ElementSize;
    MKL_LONG nx,nx_out,start_x,start_x_out,size;

    long    m;
    long    n;

    MKL_LONG Status = DFTI_NO_ERROR;
    double  Scale;
    MKL_LONG lengths[2];

    float  maxerr;
    float  global_maxerr;
    float  eps = SINGLE_EPS;

    int MPI_err, err;
    int MPI_nProc;
    int MPI_Rank;

    int local_failure = 0;
    int global_failure = 0;
    x_in = x_exp = local = work = NULL;

    /*
    **  1. Initiate MPI by calling MPI_Init (Perform MPI initialization)
    */
    MPI_err = MPI_Init(&argc, &argv);

    if (MPI_err != MPI_SUCCESS) {
        printf(" MPI initialization error\n");
        printf(" TEST FAILED\n");
        return 1;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);
    if (MPI_Rank == 0)
        printf( " Program is running on %d processes\n", MPI_nProc);

    /*
    **  Read input data from input file
    **  m - size of transform  along first dimension
    **  n - size of transform  along second dimension
    */
    if (MPI_Rank == 0) err = read_data_file_2d(argc, argv, &m, &n);
    MPI_Bcast(&err,1,MPI_INT,0,MPI_COMM_WORLD);
    if (err != 0) {
        global_failure++;
        goto CLOSE_MPI;
    }

    MPI_Bcast(&m,1,MPI_LONG_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&n,1,MPI_LONG_INT,0,MPI_COMM_WORLD);

    if (LEGEND_PRINT && (MPI_Rank == 0)) {
        printf("\n DM_COMPLEX_2D_SINGLE_EX1\n");
        printf(" Forward-Backward 2D complex transform for single precision data inplace\n\n");
        printf(" Configuration parameters:\n\n");
        printf(" DFTI_FORWARD_DOMAIN = DFTI_COMPLEX\n");
        printf(" DFTI_PRECISION      = DFTI_SINGLE\n");
        printf(" DFTI_DIMENSION      = 2\n");
        printf(" DFTI_LENGTHS        = {%ld,%ld}\n", m, n);
        printf(" DFTI_FORWARD_SCALE  = 1.0\n");
        printf(" DFTI_BACKWARD_SCALE = 1.0/(m*n)\n\n");
    }

    lengths[0] = m;
    lengths[1] = n;

    /*
    **  Allocate memory for dynamic arrays
    */
    x_in = (mkl_float_complex *)mkl_malloc(sizeof(mkl_float_complex)*m*n, 64);
    x_exp = (mkl_float_complex *)mkl_malloc(sizeof(mkl_float_complex)*m*n, 64);
    if (x_in == NULL || x_exp == NULL) {
        printf(" Rank %d: Not enough memory\n", MPI_Rank);
        local_failure++;
    }
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_MEM;

    /*
    **  Put input data and expected result
    */
    init_data_2d_c(x_in, m, n);
    memcpy(x_exp,x_in,sizeof(mkl_float_complex)*m*n);

    /* Print input data assembled together */
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0)) {
        printf("\n INPUT Global vector X, n columns\n");
        print_data_2d_c(x_in, m, n, n);
    }

    /*
    **  2. Allocate memory for the descriptor by calling DftiCreateDescriptorDM
    **     (Create DftiDM descriptor for 2D single precision transform)
    */
    Status = DftiCreateDescriptorDM(MPI_COMM_WORLD,&desc,PREC,DFTI_COMPLEX,2,lengths);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Create=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  3. Obtain some values of configuration parameters by calls to
    **     DftiGetValueDM
    */
    Status = DftiGetValueDM(desc,CDFT_LOCAL_SIZE,&size);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,size=%ld\n",(long)Status,(long)size);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    Status = DftiGetValueDM(desc,CDFT_LOCAL_NX,&nx);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,nx=%ld\n",(long)Status,(long)nx);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    Status = DftiGetValueDM(desc,CDFT_LOCAL_X_START,&start_x);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,start_x=%ld\n",(long)Status,(long)start_x);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    Status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_NX,&nx_out);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,nx_out=%ld\n",(long)Status,(long)nx_out);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    Status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_X_START,&start_x_out);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,start_x_out=%ld\n",(long)Status,(long)start_x_out);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    local=(mkl_float_complex*)mkl_malloc(size*sizeof(mkl_float_complex), 64);
    work=(mkl_float_complex*)mkl_malloc(size*sizeof(mkl_float_complex), 64);
    if (local == NULL || work == NULL) {
        printf("Rank %d: Not enough memory\n", MPI_Rank);
        local_failure++;
    }
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  4. Specify a value(s) of configuration parameters by a call(s) to
    **     DftiSetValueDM
    */
    Status = DftiSetValueDM(desc,CDFT_WORKSPACE,work);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Set=%ld,pointer=%p\n",(long)Status,work);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  5. Perform initialization that facilitates DFT computation by a call to
    **     DftiCommitDescriptorDM
    */
    Status = DftiCommitDescriptorDM(desc);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Commit=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  6. Create arrays for local parts of input and output data
    **     (if it is needed) and fill the local part of input data with
    **     values (for more information, see Distributing Data among Processes)
    **     (Spread data among processors)
    */
    RootRank=0;
    ElementSize=sizeof(mkl_float_complex);
    Status = MKL_CDFT_ScatterData(MPI_COMM_WORLD,RootRank,ElementSize,2,lengths,x_in,nx,start_x,local);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Scatter=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  7. Compute the transform by calling
    **     DftiComputeForwardDM or DftiComputeBackward
    **     (Compute Forward transform)
    */
    Status = DftiComputeForwardDM(desc,local);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("ComputeForward=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  Gather data among processors
    */
    Status = MKL_CDFT_GatherData(MPI_COMM_WORLD,RootRank,ElementSize,2,lengths,x_in,nx,start_x,local);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Gather=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  Set Scale number for Backward transform
    */
    Scale = 1.0/(float)(m*n);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("\n\n DFTI_BACKWARD_SCALE  = 1/(m*n)\n");
    Status = DftiSetValueDM(desc,DFTI_BACKWARD_SCALE,Scale);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Set=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  Commit DftiDM descriptor
    */
    Status = DftiCommitDescriptorDM(desc);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Commit=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  Spread data among processors
    */
    Status = MKL_CDFT_ScatterData(MPI_COMM_WORLD,RootRank,ElementSize,2,lengths,x_in,nx,start_x,local);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Scatter=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  Compute Backward transform
    */
    Status = DftiComputeBackwardDM(desc,local);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("ComputeBackward=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /*
    **  8. Gather local output data into the global array using MPI functions
    **     or use them otherwise
    */
    Status = MKL_CDFT_GatherData(MPI_COMM_WORLD,RootRank,ElementSize,2,lengths,x_in,nx,start_x,local);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Gather=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;

    /* Print data after DftiComputeBackwardDM; data assembled together */
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0)) {
        printf("\n Backward result X, n columns\n");
        print_data_2d_c(x_in, m, n, n);
    }

    /*
    **  Check result
    */
    maxerr = check_result_c(x_in, x_exp, m*n);
    if (ACCURACY_PRINT && (MPI_Rank == 0))
        printf("\n Accuracy = %g\n\n", maxerr);
    MPI_Allreduce(&maxerr, &global_maxerr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    if (global_maxerr < eps) {
        if (MPI_Rank == 0) printf(" TEST PASSED\n");
    } else {
        global_failure++;
    }

FREE_DESCRIPTOR:

    /*
    **  Check status of DFTI functions
    */
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) {
        dfti_example_status_print(MPI_Rank, Status);
    }

    /*
    **  9. Release memory allocated for a descriptor by a call to
    **     DftiFreeDescriptorDM
    **     (Free DftiDM descriptor)
    */
    Status = DftiFreeDescriptorDM(&desc);
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("FreeDescriptor=%ld\n",(long)Status);
    if (! DftiErrorClass(Status, DFTI_NO_ERROR)) {
        dfti_example_status_print(MPI_Rank, Status);
        local_failure++;
    }
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

FREE_MEM:

    /*
    **  Deallocate memory for dynamic arrays
    */
    if (x_in  != NULL) mkl_free(x_in);
    if (x_exp != NULL) mkl_free(x_exp);
    if (local != NULL) mkl_free(local);
    if (work  != NULL) mkl_free(work);

CLOSE_MPI:

    /*
    **  10. Finalize communication through MPI by calling MPI_Finalize
    **      (Finalize MPI)
    */
    MPI_Finalize();

    if (global_failure != 0) {
        if (MPI_Rank == 0) printf(" TEST FAILED\n");
        return 1;
    }

    if (MPI_Rank == 0) printf(" END OF TEST\n");

    return 0;
}
