/*******************************************************************************
* Copyright 2005-2020 Intel Corporation.
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
!     Examples support function set
!
!****************************************************************************/

#define _CRT_SECURE_NO_WARNINGS 
#include <stdio.h>
#include <stdlib.h>
#include "mkl_service.h"
#include "mkl_cdft.h"
#include "mkl_cdft_examples.h"

/*
**  Read data from input file routine
*/
int read_data_file_2d(int argc, char *argv[], long *m, long *n)
{
    FILE *in_file;
    char *in_file_name;
    char scan_line[100];
    int  maxline=100;

    if (argc == 1) {
        printf("\n You must specify in_file data file as 1-st parameter");
        return 1;
    }
    in_file_name = argv[1];

    if ( (in_file = fopen( in_file_name, "r" )) == NULL ) {
        printf("\nERROR on OPEN '%s' with mode=%s\n", in_file_name, "r");
        return 1;
    }

    fgets(scan_line, maxline, in_file);
    fscanf(in_file,"%ld\n", m);
    fgets(scan_line, maxline, in_file);
    fscanf(in_file,"%ld\n", n);
    fclose(in_file);
    return 0;
}

/*
**  Init data routines
*/
void init_data_2d_z(void *x, long m, long n)
{
    double* p = x;
    long i;

    p[0] = 1.0;
    for (i = 1; i < m*n*2; i++) p[i] = 0.0;
}

void init_data_2d_c(void *x, long m, long n)
{
    float* p = x;
    long i;

    p[0] = 1.0;
    for (i = 1; i < m*n*2; i++) p[i] = 0.0;
}

/*
**  Print data routines
*/
void print_data_2d_z(void *x, long m, long n, long c)
{
    double* p = x;
    long i, j;

    for (i = 0; i < m; i++) {
        printf("\n Row %ld:\n ",i);
        for (j = 0; j < n; j++) {
            printf("(%8.3f,%8.3f)", p[2*j], p[2*j+1]);
            if ((j%c == c-1) && (j != n-1)) printf("\n ");
        }
        p += 2*n;
    }
    printf("\n");
}

void print_data_2d_c(void *x, long m, long n, long c)
{
    float* p = x;
    long i, j;

    for (i = 0; i < m; i++) {
        printf("\n Row %ld:\n ",i);
        for (j = 0; j < n; j++) {
            printf("(%8.3f,%8.3f)", p[2*j], p[2*j+1]);
            if ((j%c == c-1) && (j != n-1)) printf("\n ");
        }
        p += 2*n;
    }
    printf("\n");
}


/*
**  Print status routine
*/
void dfti_example_status_print(int rank, long status)
{
    printf("Rank %d: Error message: %s\n", rank, DftiErrorMessage(status));
    printf("Rank %d: Error status = %d\n", rank, (int)status);
}

/*
**  Check error routines
*/
double check_result_z(void* x, void* res_exp, long n)
{
    double* x_in = x;
    double* x_exp = res_exp;
    double maxerr, d;
    long i;

    maxerr = 0.0;
    for (i = 0; i < 2*n; i++) {
        d = x_exp[i] - x_in[i];
        if (d < 0.0) d = -d;
        if (d > maxerr) maxerr = d;
    }
    return maxerr;
}

float check_result_c(void* x, void* res_exp, long n)
{
    float* x_in = x;
    float* x_exp = res_exp;
    float maxerr, d;
    long i;

    maxerr = 0.0;
    for (i = 0; i < 2*n; i++) {
        d = x_exp[i] - x_in[i];
        if (d < 0.0) d = -d;
        if (d > maxerr) maxerr = d;
    }
    return maxerr;
}

int MKL_CDFT_Data(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                  MKL_LONG Lengths[],void *global,long nx,long start_x,void *local,
                  int Flag)
{

    int nProc,nRank,MPI_err,i,fd;
    int tmp[2],*counts,*displs,*buf;
    MPI_Request req;
    MPI_Status status;
    counts = displs = buf = NULL;

    MPI_err=MPI_Comm_rank(Comm,&nRank);
    if (MPI_err != MPI_SUCCESS) return MPI_err;

    if (nRank == RootRank) {
        MPI_err=MPI_Comm_size(Comm,&nProc);
        if (MPI_err != MPI_SUCCESS) return MPI_err;
        counts=(int*)mkl_malloc(nProc*sizeof(int), 64);
        displs=(int*)mkl_malloc(nProc*sizeof(int), 64);
        buf=(int*)mkl_malloc(2*nProc*sizeof(int), 64);
        if (counts == NULL || displs == NULL || buf == NULL) {
            printf("Not enough memory\n");
            return 1;
        }
    }

    fd=ElementSize;
    for (i=1;i<(int)Dim;i++) fd*=(int)Lengths[i];

    tmp[0]=(int)nx*fd;
    tmp[1]=(int)start_x*fd;

    MPI_err=MPI_Gather(tmp,2,MPI_INT,buf,2,MPI_INT,RootRank,Comm);
    if (MPI_err != MPI_SUCCESS) return MPI_err;

    if (nRank == RootRank) {
        for (i=0;i<nProc;i++) {
            counts[i]=buf[2*i];
            displs[i]=buf[2*i+1];
        }
    }

    if (Flag == 0) {
        MPI_err = MPI_Irecv(local,tmp[0],MPI_BYTE, RootRank, 123,Comm,&req);
        if (MPI_err != MPI_SUCCESS) return MPI_err;

        if (nRank == RootRank) {
            for (i=0;i<nProc;i++) {
                MPI_err = MPI_Send(((char*)global)+displs[i],counts[i],MPI_BYTE,i,123,Comm);
                if (MPI_err != MPI_SUCCESS) return MPI_err;
            }
        }
        MPI_Wait(&req,&status);
        if (MPI_err != MPI_SUCCESS) return MPI_err;
    }
    else {
        MPI_err = MPI_Gatherv(local,tmp[0],MPI_BYTE,global,counts,displs,MPI_BYTE,RootRank,Comm);
        if (MPI_err != MPI_SUCCESS) return MPI_err;
    }

    if (nRank == RootRank) {
        mkl_free(buf);
        mkl_free(displs);
        mkl_free(counts);
    }
    return 0;
}

int MKL_CDFT_ScatterData(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                         MKL_LONG Lengths[],void *global_in,long nx,long start_x,
                         void *local_in)
{
    return MKL_CDFT_Data(Comm,RootRank,ElementSize,Dim,Lengths,
                         global_in,nx,start_x,local_in,0);
}

int MKL_CDFT_GatherData(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                        MKL_LONG Lengths[],void *global_out,long nx,long start_x,
                        void *local_out)
{
    return MKL_CDFT_Data(Comm,RootRank,ElementSize,Dim,Lengths,
                         global_out,nx,start_x,local_out,1);
}
