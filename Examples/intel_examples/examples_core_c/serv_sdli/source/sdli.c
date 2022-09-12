/*******************************************************************************
* Copyright 1999-2020 Intel Corporation.
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
!  Content:
!      Intel(R) MKL Single Dynamic Library Interface example
!******************************************************************************/

///
/// This example calls DGEMM('N','N',NN,NN,NN,...) and
///   prints the number of the spent processor cycles.
/// If the first parameter is specified as 'seq' or 'par',
///   the example set Intel(R) MKL threading layer accordingly
///   by the call MKL_Set_Threading_Layer().

#include <stdio.h>
#include <string.h>
#include <mkl.h>

#define NN 1000

static void aux_call_dgemm(MKL_INT N, double* A, double* B, double* C, int PRINT_CLOCKS);

int main(int argc, char * args[]) {

  int mkl_thr_i, mkl_thr_o;
  double *a, *b, *c;

  int n = NN;

    if (argc>1) {
        if (!strcmp(args[1],"seq")) {
            mkl_thr_i = MKL_THREADING_SEQUENTIAL;
        } else if (!strcmp(args[1],"par")) {
            mkl_thr_i = MKL_THREADING_INTEL;
        } else {
            printf("\n*** sdli example: unsupported input value '%s'\n\n",args[1]);
            exit(0);
        }
        printf("\n    sdli example: set '%s' threading layer\n",args[1]); fflush(0);
        mkl_thr_o = MKL_Set_Threading_Layer(mkl_thr_i);
        if (mkl_thr_o != mkl_thr_i) {
            printf("\n*** sdli example: couldn't set '%s' threading layer %d/%d\n\n",args[1],mkl_thr_i,mkl_thr_o);
            exit(0);
        }
    }

    (void) MKL_Set_Interface_Layer(MKL_INTERFACE_LP64); /** Force LP64 interface **/

    a = (double*)MKL_malloc(n*n*sizeof(double),128);
    b = (double*)MKL_malloc(n*n*sizeof(double),128);
    c = (double*)MKL_malloc(n*n*sizeof(double),128);

    aux_call_dgemm(n,a,b,c,1);

    MKL_free(a);
    MKL_free(b);
    MKL_free(c);

    return 0;
}

static void aux_call_dgemm(MKL_INT N, double* A, double* B, double* C, int PRINT_CLOCKS) {

  int      i, j;
  char     transa, transb;
  double   alpha, beta;

  unsigned MKL_INT64  dgemm_s, dgemm_e;
  unsigned MKL_INT64  dgemm_clocks, dgemm_cb, dgemm_ct, dgemm_cu;

    transa = 'N'; transb = 'N';
    alpha = 1.1;  beta = -1.2;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
           A[j*N+i] = (double)(i+1);
           B[j*N+i] = (double)(-(i+1));
           C[j*N+i] = 0.0;
        }
    }

    MKL_Get_Cpu_Clocks(&dgemm_s);
    dgemm(&transa,&transb,&N,&N,&N,&alpha,A,&N,B,&N,&beta,C,&N);
    mkl_get_cpu_clocks(&dgemm_e);

    dgemm_clocks = dgemm_e-dgemm_s;
    if (PRINT_CLOCKS) {
        dgemm_cb = dgemm_clocks/1000000;
        dgemm_ct = (dgemm_clocks-(dgemm_cb*1000000))/1000;
        dgemm_cu = (dgemm_clocks-(dgemm_cb*1000000)-(dgemm_ct*1000));
        printf("    DGEMM (%4d) on %2d thread(s):%8.3d,%3d,%3d clocks\n\n",
            (int)N,MKL_Domain_Get_Max_Threads(MKL_DOMAIN_BLAS),(int)dgemm_cb,(int)dgemm_ct,(int)dgemm_cu);
    }
}
