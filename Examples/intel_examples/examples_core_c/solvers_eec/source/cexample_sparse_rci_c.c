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
!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Extended Eigensolvers
!             C example
!
!*******************************************************************************
!
! Example program for using Intel MKL Extended Eigensolvers (sparse format).
!
! The following routines are used in the example:
!          CGEMM MKL_SPARSE_C_MM MKL_SPARSE_C_ADD CFEAST_HRCI PARDISO.
!
! Consider the matrix A:
!
!                 |  10   1+2i  0    0    0    0    0    0    0    0    |
!                 |  1-2i  9   2+3i  0    0    0    0    0    0    0    |
!                 |  0    2-3i  8   3+4i  0    0    0    0    0    0    |
!                 |  0    0    3-4i  7   4+5i  0    0    0    0    0    |
!                 |  0    0    0    4-5i  6   5+6i  0    0    0    0    |
!    A    =       |  0    0    0    0    5-6i  5   6+7i  0    0    0    |,
!                 |  0    0    0    0    0    6-7i  4   7+8i  0    0    |
!                 |  0    0    0    0    0    0    7-8i  3   8+9i  0    |
!                 |  0    0    0    0    0    0    0    8-9i  2   9+10i |
!                 |  0    0    0    0    0    0    0    0    9-10i 1    |
!
! stored as sparse matrix (SINGLE COMPLEX PRECISION version).
! B is a unit matrix:
!
!                 |  1   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0  |
!                 |  0   0   0   0   0   0   0   0   0   1  |
!
!  In what follows the symbol ' represents a transpose conjugate operation.
!
!  The test performs the following operations :
!
!       1. The code calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!       2. The  code solves  the generalized eigenvalue problem  Ax=eBx using
!          CFEAST_HRCI.
!
!       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!           are the expected eigenvalues  and E(i) are eigenvalues computed
!           with the help of CFEAST_HRCI().
!
!       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I
!          where X is the matrix of eigenvectors computed with the help of
!          CFEAST_HRCI. CGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    //!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const MKL_INT N = 10;
    MKL_Complex8  val[19], nval[19], valb[10], *valz, caux[8*10];
    MKL_INT       rows[11], cols[19], rowsb[11], colsb[10], *rowsz, *colsz;

    //!!!!!!!!!!!!!!! Declaration of Spblas variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!! A, B - Handles containing a sparse matrix in internal data structure!!!!!!
    //!!! descrA - Structure specifying sparse matrix properties              !!!!!!
    //!!! indexing - Indicates how input arrays are indexed                   !!!!!!
    //!!! layout - Describes the storage scheme for the dense matrix in mm    !!!!!!
    //!!! opeartion - Specifies operation op() on input matrix in mm          !!!!!!
    MKL_INT rows_zC, cols_zC, *pdum;
    sparse_status_t status;
    sparse_matrix_t A=NULL, B=NULL, nA=NULL, C=NULL;
    struct matrix_descr descrA, descrB;
    sparse_layout_t layout = SPARSE_LAYOUT_COLUMN_MAJOR;
    sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ONE, indexing_zC;

    //!!!!!!!!!!!!!!! Declaration of FEAST variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
    MKL_INT       fpm[128];
    float         Emin, Emax;
    float         epsout;
    MKL_INT       loop;
    MKL_INT       ijob;
    MKL_Complex8  Ze;
    MKL_Complex8  workc[8*10], work[8*10];
    MKL_Complex8  Aq[8*8], Sq[8*8];
    MKL_INT       L = 8;
    MKL_INT       M0, M, info;
    float         E[10];
    MKL_Complex8  X[100];
    float         res[10];

    //!!!!!!!!!!!!!!! Declaration of PARDISO variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!! Internal solver memory pointer pt,                                     !!!
    //!!! 32-bit: int pt[64]; 64-bit: long int pt[64]                            !!!
    //!!! or void *pt[64] should be OK on both architectures                     !!!
    void *pt[64];
    //!!! Pardiso control parameters
    MKL_INT iparm[64];
    MKL_INT phase, maxfct, mnum, mtype, msglvl, error;
    MKL_INT idum = 0; /* Integer dummy. */

    //!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!! Eig - array for storing exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!
    float         Eig[10];
    float         R[10];
    MKL_Complex8  Y[100];
    MKL_INT       i,j;
    MKL_INT       ldx = 10, ldy = 10;
    float         trace;
    float         smax, eigabs, t1;
    char          CGEMMC = 'C', CGEMMN = 'N';
    MKL_Complex8  one, zero;
    MKL_INT       colsX, imem;

    printf("\n    FEAST CFEAST_HRCI EXAMPLE PROGRAM\n");
    one.real  = (float)1.0;
    one.imag  = (float)0.0;
    zero.real = (float)0.0;
    zero.imag = (float)0.0;

    //!!!!!!!!!!!!!!! Exact eigenvalues in range (2.0, 12.0) !!!!!!!!!!!!!!!!!!!!!!
    for ( i=0; i<N; i++ )
        Eig[i] = (float)0.0;

    Eig[0]=(float)2.231051;
    Eig[1]=(float)6.058517;
    Eig[2]=(float)9.109751;
    Eig[3]=(float)11.703148;

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!     Initialize matrices     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for ( i=0; i<N*N; i++ )
    {
        X[i] = zero;
    }

    //!!!!!!!  Initialize upper part of Hermitian matrix A   !!!!!!!!!!!!!!!!!!!!!!!
    for ( i=0; i<9; i++ )
    {
        val[i*2].real = (float)(N-i);
        val[i*2].imag = (float)0.0;
        val[1+i*2].real = (float)(i+1);
        val[1+i*2].imag = (float)(i+2);
    }
    val[18].real= (float)1.0;
    val[18].imag= (float)0.0;

    for ( i=0; i<9; i++ )
    {
        cols[i*2]=i+1;
        cols[1+i*2]=i+2;
        rows[i]=2*i+1;
    }
    cols[18]=10;
    rows[9]=19;
    rows[10]=20;

    status = mkl_sparse_c_create_csr (&A, indexing, N, N, rows, rows+1, cols, val);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("mkl_sparse_c_create_csr status %d \n", status);
        return 1;
    }

    descrA.type = SPARSE_MATRIX_TYPE_HERMITIAN;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    status = mkl_sparse_set_mm_hint (A, operation, descrA, layout, L, 10);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("mkl_sparse_set_mm_hint status %d \n", status);
        return 1;
    }

    status = mkl_sparse_optimize (A);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("mkl_sparse_optimize status %d \n", status);
        return 1;
    }

    //!!!!!!!!!!!!!!!!!!    Initialize unit matrix B   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for ( i=0; i<N; i++ )
    {
        valb[i].real= (float)1.0;
        valb[i].imag= (float)0.0;
        colsb[i] = i+1;
        rowsb[i] = i+1;
    }
    rowsb[N] = N+1;

    status = mkl_sparse_c_create_csr (&B, indexing, N, N, rowsb, rowsb+1, colsb, valb);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("mkl_sparse_c_create_csr status %d \n", status);
        return 1;
    }

    descrB.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrB.mode = SPARSE_FILL_MODE_UPPER;
    descrB.diag = SPARSE_DIAG_NON_UNIT;
    status = mkl_sparse_set_mm_hint (B, operation, descrB, layout, L, 10);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("mkl_sparse_set_mm_hint status %d \n", status);
        return 1;
    }

    status = mkl_sparse_optimize (B);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("mkl_sparse_optimize status %d \n", status);
        return 1;
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!        Initialize  upper part of matrix -A     !!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for ( i=0; i<19; i++ )
    {
        nval[i].real= -val[i].real;
        nval[i].imag= -val[i].imag;
    }

    status = mkl_sparse_c_create_csr (&nA, indexing, N, N, rows, rows+1, cols, nval);
    if (status != SPARSE_STATUS_SUCCESS)
    {
        printf("mkl_sparse_c_create_csr status %d \n", status);
        return 1;
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!   Print matrix dimension          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    printf("Sparse matrix size  %i \n", N);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!      Testing FEAST sparse format drivers !!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!! Search interval [Emin,Emax] including M eigenpairs !!!!!!!!!!!!
    Emin = (float)2.0;
    Emax = (float)12.0;
    M0 = L;
    M = M0;
    printf(" Search interval [ %.15e, %.15e ]  \n", Emin, Emax);
    epsout = (float)0.0;
    loop = 0;
    info = 0;
    //!!!!!!!!!!!!     Initialize PARDISO                   !!!!!!!!!!!!!!!!
    mtype = -4;   /* Hermitian and non-definite matrix   */
    pardisoinit(pt,&mtype,iparm);
    phase    = 11;/* Analysis                            */
    maxfct   = 1; /* Maximal number of factors in memory */
    mnum     = 1; /* The number of matrix to solve       */
    msglvl   = 0; /* No print statistical information    */
    iparm[5] = 1; /* Solution and rhs are input/output   */
    iparm[27] = 1; /* Pardiso single precision           */
    //
    //         Task 1. Initialize FEAST to define the default values for the input
    //         FEAST parameters.
    //
    ijob = -1;
    info = 0;
    feastinit(fpm);
    fpm[0] = 1; /* Generate runtime messages */
    fpm[5] = 1; /* Second stopping criteria  */
    //
    //         Task 2. The  code solves  the generalized eigenvalue problem  Ax=eBx using
    //         CFEAST_SRCI.
    //
    while ( ijob != 0 )
    {
        cfeast_hrci(&ijob,&N,&Ze,work,workc,Aq,Sq,fpm,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);
        if ( info != 0 )
        {
            printf("CFEAST_HRCI info %i \n", info);
            return 1;
        }
        switch ( ijob )
        {
            case -2:
                //!!!!!!!!!!!!!!  New loop                 !!!!!!!!!!!!!!!!!!!!!!!!!
                break;
            case 0:
                //!!!!!!!!!!!!!!  End                      !!!!!!!!!!!!!!!!!!!!!!!!!
                break;
            case 10:
                //!!!!!!!!!!!!!!  Factorize matrix ZeB-A   !!!!!!!!!!!!!!!!!!!!!!!!!
                if (C)
                {
                    status = mkl_sparse_destroy(C);
                    if (status != SPARSE_STATUS_SUCCESS)
                    {
                        printf("mkl_sparse_destroy status %d \n", status);
                        return 1;
                    }
                }
                status = mkl_sparse_c_add (operation, B, Ze, nA, &C);
                if (status != SPARSE_STATUS_SUCCESS)
                {
                    printf("mkl_sparse_c_add status %d \n", status);
                    return 1;
                }
                status = mkl_sparse_c_export_csr (C, &indexing_zC, &rows_zC, &cols_zC, &rowsz, &pdum, &colsz, &valz);
                if (status != SPARSE_STATUS_SUCCESS)
                {
                    printf("mkl_sparse_c_exposr_csr status %d \n", status);
                    return 1;
                }
                phase = 12; /* Reordering and numerical factorization */
                iparm[34] = (indexing_zC == SPARSE_INDEX_BASE_ZERO ? 1 : 0); /*  columns and rows indexing */
                pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,valz,rowsz,colsz,&idum,&M0,iparm,&msglvl,workc,caux,&error);
                if ( error != 0 )
                {
                    printf("PARDISO error %d \n", error);
                    return 1;
                }
                break;
            case 11:
                //!!!!!!! Solve (ZeB-A) caux = workc[0:N-1][0:M0-1] and put result into workc !!
                phase  = 33; /* Solve */
                pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,valz,rowsz,colsz,&idum,&M0,iparm,&msglvl,workc,caux,&error);
                if ( error != 0 )
                {
                    printf("PARDISO error %d \n", error);
                    return 1;
                }
                break;
            case 20:
                //!!!!!!!!!!!!! Factorize matrix (ZeB-A)^H !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //         We can do nothing here;
                break;
            case 21:
                //!!!!! Solve (ZeB-A)^H caux = workc[0:N-1][0:M0-1] and put result into workc !!
                phase  = 33;   /* Solve */
                iparm[11] = 1; /* Transpose conjugate solve */
                pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,valz,rowsz,colsz,&idum,&M0,iparm,&msglvl,workc,caux,&error);
                if ( error != 0 )
                {
                    printf("PARDISO error %d \n", error);
                    return 1;
                }
                break;
            case 30:
                //!!!!!!!!!!!!! Perform multiplication A x[0:N-1][i:j]      !!!!!!!!
                //!!!!!!!!!!!!! and put result into work[0:N-1][i:j]        !!!!!!!!
                //!!!!!!!!!!!!! where i = fpm[23]-1, j = fpm[23]+fpm[24]-2  !!!!!!!!
                colsX = fpm[24];
                imem = N*(fpm[23]-1);
                status = mkl_sparse_c_mm (operation, one, A, descrA, layout, X+imem, colsX, N, zero, work+imem, N);
                if (status != SPARSE_STATUS_SUCCESS)
                {
                    printf("mkl_sparse_c_mm status %d \n", status);
                    return 1;
                }
                break;
            case 40:
                //!!!!!!!!!!!!! Perform multiplication B x[0:N-1][i:j]      !!!!!!!!
                //!!!!!!!!!!!!! and put result into work[0:N-1][i:j]        !!!!!!!!
                //!!!!!!!!!!!!! where i = fpm[23]-1, j = fpm[23]+fpm[24]-2  !!!!!!!!
                colsX = fpm[24];
                imem = N*(fpm[23]-1);
                status = mkl_sparse_c_mm (operation, one, B, descrB, layout, X+imem, colsX, N, zero, work+imem, N);
                if (status != SPARSE_STATUS_SUCCESS)
                {
                    printf("mkl_sparse_z_mm status %d \n", status);
                    return 1;
                }
                break;
            default:
                printf("Wrong ijob %i", ijob); fflush(0);
                return 1;
        }
    }
    //!!!!!!!!!!!!!!! Release memory                              !!!!!!!!!!!!!!!!!!!!
    phase = -1;
    pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,val,rows,cols,&idum,&M0,iparm,&msglvl,workc,caux,&error);
    if ( error != 0 )
    {
        printf("PARDISO RELEASE MEMORY ERROR %d \n", error);
        return 1;
    }
    //
    //         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
    //         are the expected eigenvalues  and E(i) are eigenvalues computed
    //         with the help of CFEAST_HRCI().
    //
    printf("  Number of eigenvalues found %d \n", M);
    printf("   Computed    |    Expected  \n");
    printf("   Eigenvalues |    Eigenvalues \n");
    eigabs = (float)0.0;
    for ( i=0; i<M; i++ )
    {
        R[i] = fabs(E[i]-Eig[i]);
        eigabs = max(eigabs, R[i]);
        printf("%.15e %.15e \n", E[i], Eig[i]);
    }
    printf(" Max value of | computed eigenvalue -expected eigenvalues | %.15e \n", eigabs);
    //
    //         Task 4.  Compute the maximum absolute value of the matrix
    //         Y=(X')*X-I  where X is the matrix of eigenvectors computed with
    //         the help of CFEAST_HRCI.
    //         Call CGEMM (BLAS Level 3 Routine) to compute (X')*X.
    //
    cgemm(&CGEMMC,&CGEMMN,&M,&M,&N,&one,X,&ldx,X,&ldx,&zero,Y,&ldy);

    //          Compute Y=Y-I.
    for ( i=0; i<M; i++ )
        Y[i*M + i].real=Y[i*M + i].real-(float)1.0;

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("# Search interval [Emin,Emax] %.15e %.15e\n",Emin,Emax);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace = (float)0.0;
    for ( i=0; i<M; i++ )
    {
        trace = trace+E[i];
    }
    printf("TRACE %.15e \n", trace);
    printf("Relative error on the Trace %.15e \n",epsout );
    printf("Eigenvalues/Residuals\n");
    for ( i=0; i<M; i++ )
    {
        printf("   %d  %.15e %.15e \n",i, E[i], res[i]);
    }
    smax = (float)0.0;
    for ( i=0; i<M; i++ )
    {
        for ( j=0; j<M; j++ )
        {
            t1 = sqrt(Y[i*M + j].imag*Y[i*M + j].imag+Y[i*M + j].real*Y[i*M + j].real);
            smax = max(smax, t1);
        }
    }
    printf( "Max of (conjugate transposed of X)*X-I %.15e \n", smax);

    mkl_sparse_destroy(A);
    mkl_sparse_destroy(B);
    mkl_sparse_destroy(nA);
    mkl_sparse_destroy(C);
    return 0;
}
