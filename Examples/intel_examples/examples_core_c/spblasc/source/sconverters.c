/*******************************************************************************
* Copyright 2004-2020 Intel Corporation.
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
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse format
*             converters C example
*
********************************************************************************
*
* Example program for using Intel MKL Sparse format converters
* The following Sparse  Sparse format converters are used in the example:
*
*          MKL_SDNSCSR
*          MKL_SCSRCOO
*          MKL_SCSRBSR
*          MKL_SCSRDIA
*          MKL_SCSRCSC
*          MKL_SCSRSKY
*
*/
#include <stdio.h>
#include "mkl_spblas.h"
#include "mkl_types.h"


int main()
   {
/********************************************************************************
*     Definition arrays for sparse matrix formats
********************************************************************************/
#define M      4
#define N      4
#define LDA    4
#define NZMAX  8
#define NNZ    8
#define MBLK   2
#define NN     2
#define INFO   0
#define MN     16
#define IBASE1 1
#define IBASE2 1
#define LOCAT  2
#define IDIAG 3
#define NDIAG 4
#define INDIA 12
    MKL_INT m = M, n=N, lda=LDA, nzmax=NZMAX, nnz = NNZ, mblk=MBLK, nn=NN, info=INFO;
    MKL_INT ibase1 = IBASE1,ibase2 = IBASE2, locat = LOCAT,idiag = IDIAG,ndiag = NDIAG;
    float   Adns[MN];
    float   Adns_standard[MN];
    float   Absr[NZMAX];
    float   Absr_standard[NZMAX]     = {5.0, 9.0, 8.0, 2.0, 3.0, 1.0, 6.0, 4.0};
    float   Acsr[NZMAX];
    float   Acsr_standard[NZMAX]     = {5.0, 8.0, 9.0, 2.0, 3.0, 6.0, 1.0, 4.0};
    float   Acsc[NZMAX];
    float   Acsc_standard[NZMAX]     = {5.0, 9.0, 8.0, 2.0, 3.0, 1.0, 6.0, 4.0};
    float   Adia[INDIA];
    float   Adia_standard[INDIA]     = {0.0, 9.0, 0.0, 1.0, 5.0, 2.0, 3.0, 4.0, 8.0, 0.0, 6.0, 0.0};
    float   Asky[6];
    float   Askyl_standard[6]        = {5.0, 9.0, 2.0, 3.0, 1.0, 4.0};
    float   Acoo[NZMAX];
    float   Acoo_standard[NZMAX]     = {5.0, 8.0, 9.0, 2.0, 3.0, 6.0, 1.0, 4.0};
    MKL_INT AI[M+1];
    MKL_INT AI1[M+1];
    MKL_INT AI_standard[M+1]         = {1, 3, 5, 7, 9};
    MKL_INT AJ[NZMAX];
    MKL_INT AJ1[NZMAX];
    MKL_INT AJ_standard[NZMAX]       = {1, 2, 1, 2, 3, 4, 3, 4};
    MKL_INT AJB[NN];
    MKL_INT AJB_standard[MBLK]       = {1, 2};
    MKL_INT AIB[NN+1];
    MKL_INT AIB_standard[MBLK+1]     = {1, 2, 3};
    MKL_INT ir[NZMAX];
    MKL_INT ir_standard[NZMAX]       = {1, 1, 2, 2, 3, 3, 4, 4};
    MKL_INT jc[NZMAX];
    MKL_INT jc_standard[NZMAX]       = {1, 2, 1, 2, 3, 4, 3, 4};
    MKL_INT pointers[M+1];
    MKL_INT pointersl_standard[M+1]  = {1, 2, 4, 5, 7};
    MKL_INT distance[IDIAG];
    MKL_INT distance_standard[IDIAG] = {-1, 0, 1};
    MKL_INT AJL_standard[6]          = {1, 1, 2, 3, 3, 4};
//*************************************************************************************************
//*    Declaration of local variables :
//*************************************************************************************************
    MKL_INT ifail=1;
    MKL_INT nr,ldAbsr;
    MKL_INT i,j;
    float  rfail=0.0;
    MKL_INT job[8];
    printf("\n EXAMPLE PROGRAM FOR CONVERTER FROM ONE\n");
    printf("\n    SPARSE FORMAT ROUTINES TO OTHER    \n");
    printf("\n       REAL SINGLE PRECISION           \n");
    locat=2;
    ibase1=1;
    ibase2=1;
    job[1]=ibase1;
    job[2]=ibase2;
    job[3]=locat;
    job[4]=nzmax;
//***************************************************************************************************
//* TASK 1    Obtain compressed sparse row matrix from dense matrix
//***************************************************************************************************
    for ( j=0; j<n; j++)
        for ( i=0; i<m; i++)
            Adns[i + lda*j]=0.0;

    Adns[0]=5.0;
    Adns[1]=9.0;
    Adns[4]=8.0;
    Adns[5]=2.0;
    Adns[10]=3.0;
    Adns[11]=1.0;
    Adns[14]=6.0;
    Adns[15]=4.0;


    for ( j=0; j<n; j++)
        for ( i=0; i<m; i++)
            Adns_standard[i + lda*j]=Adns[i + lda*j];

    job[0]=0;
    job[5]=1;
    mkl_sdnscsr(job, &m, &n, Adns, &lda, Acsr, AJ, AI, &info);
    if (info!=0) goto FAILURE1;
    for ( i=0; i<n+1; i++)
    {
        ifail=AI[i]-AI_standard[i];
        if (ifail!=0) goto FAILURE1;
    }
    for ( i=0; i<nzmax; i++)
    {
        ifail=AJ[i]-AJ_standard[i];
        if (ifail!=0) goto FAILURE1;
    }
    for ( i=0; i<nzmax; i++)
    {
        rfail=Acsr[i]-Acsr_standard[i];
        if (rfail!=0) goto FAILURE1;
    }
//***************************************************************************************************
//* TASK 2    Obtain dense matrix from compressed sparse row matrix
//***************************************************************************************************
    for ( j=0; j<n; j++)
        for ( i=0; i<m; i++)
            Adns[i+j*lda]=0.0;

    job[0]=1;
    mkl_sdnscsr(job, &m, &n, Adns, &lda, Acsr, AJ, AI, &info);
    if (info!=0) goto FAILURE2;
    for ( j=0; j<n; j++)
    {
        for ( i=0; i<m; i++)
        {
            rfail = Adns_standard[i + lda*j]-Adns[i + lda*j];
            if (rfail!=0) goto FAILURE2;
        }
    }
//***************************************************************************************************
//* TASK 3    Obtain sparse coordinate matrix from compressed sparse row matrix
//***************************************************************************************************
    job[0]=0;
    job[5]=3;
    mkl_scsrcoo(job, &n, Acsr, AJ, AI,&nnz, Acoo, ir, jc, &info);
    if (info!=0) goto FAILURE3;
    for ( i=0; i<nzmax; i++)
    {
        ifail=ir[i]-ir_standard[i];
        if (ifail!=0) goto FAILURE3;
    }
    for ( i=0; i<nzmax; i++)
    {
        ifail=jc[i]-jc_standard[i];
        if (ifail!=0) goto FAILURE3;
    }

    for ( i=0; i<nzmax; i++)
    {
        rfail=Acoo[i]-Acoo_standard[i];
        if (rfail!=0) goto FAILURE3;
    }
//***************************************************************************************************
//* TASK 4    Obtain compressed sparse row matrix from sparse coordinate matrix
//***************************************************************************************************
    job[0]=1;
    job[5]=2;
    mkl_scsrcoo(job,&n, Acsr, AJ,AI,&nnz,Acoo, ir,jc,&info);
    if (info!=0) goto FAILURE4;
    for ( i=0; i<n+1; i++)
    {
        ifail=AI[i]-AI_standard[i];
        if (ifail!=0) goto FAILURE4;
    }
    for ( i=0; i<nzmax; i++)
    {
        ifail=AJ[i]-AJ_standard[i];
        if (ifail!=0) goto FAILURE4;
    }
    for ( i=0; i<nzmax; i++)
    {
        rfail=Acsr[i]-Acsr_standard[i];
        if (rfail!=0) goto FAILURE4;
    }
//***************************************************************************************************
//* TASK 5    Obtain block sparse row matrix from compressed sparse row matrix
//***************************************************************************************************
    ldAbsr=mblk*mblk;
    for ( i=0; i < nzmax; i++)
        Absr[i]=0.0;
    job[0] = 0;
    job[2] = 1;
    job[5] = 1;
    mkl_scsrbsr(job, &m, &mblk, &ldAbsr, Acsr, AJ, AI, Absr, AJB, AIB, &info);

    nr = 1 + (m-1) / mblk;
    if (info!=0) goto FAILURE5;
    for ( i=0; i < nr+1; i++)
    {
        ifail=AIB[i]-AIB_standard[i];
        if (ifail!=0) goto FAILURE5;
    }
    for ( i=0; i < mblk; i++)
    {
        ifail=AJB[i]-AJB_standard[i];
        if (ifail!=0) goto FAILURE5;
    }

    for ( i=0; i < nzmax; i++)
    {
        rfail=Absr[i]-Absr_standard[i];
        if (rfail!=0) goto FAILURE5;
    }
//***************************************************************************************************
//* TASK 6    Obtain compressed sparse row matrix from block sparse row matrix
//***************************************************************************************************
    ldAbsr=mblk*mblk;
    for ( i=0; i<nzmax; i++)
        Acsr[i]=0;
    job[0] = 1;
    job[5] = 3;
    mkl_scsrbsr(job, &nn, &mblk, &ldAbsr, Acsr, AJ, AI, Absr, AJB, AIB, &info);

    if (info!=0) goto FAILURE6;
    for ( i=0; i< n+1; i++)
    {
        ifail=AI[i]-AI_standard[i];
        if (ifail!=0) goto FAILURE6;
    }
    for ( i=0; i < nzmax; i++)
    {
        ifail=AJ[i]-AJ_standard[i];
        if (ifail!=0) goto FAILURE6;
    }
    for ( i=0; i < nzmax; i++)
    {
        rfail=Acsr[i]-Acsr_standard[i];
        if (rfail!=0) goto FAILURE6;
    }
//***************************************************************************************************
//* TASK 7   Obtain diagonal matrix from compressed sparse row matrix
//***************************************************************************************************
    for ( i=0; i< idiag; i++)
        distance[i]=distance_standard[i];
    job[0] = 0;
    job[5] = 0;
    mkl_scsrdia(job, &n, Acsr, AJ, AI, Adia, &ndiag, distance, &idiag, Acsr, AJ, AI, &info);
    if (info!=0) goto FAILURE7;
    for ( i=0; i< idiag; i++)
    {
        ifail=distance[i]-distance_standard[i];
        if (ifail!=0) goto FAILURE7;
    }
    for ( j=0; j< idiag; j++)
     for ( i=0; i< ndiag; i++)
    {
        rfail=Adia[i]-Adia_standard[i];
        if (rfail!=0) goto FAILURE7;
    }
//***************************************************************************************************
//* TASK 8    Obtain compressed sparse row matrix from  diagonal matrix
//***************************************************************************************************
    job[0] = 1;
    job[5] = 0;

    mkl_scsrdia(job, &n, Acsr, AJ, AI, Adia, &ndiag, distance, &idiag, Acsr, AJ, AI, &info);

    if (info!=0) goto FAILURE8;
    for ( i=0; i< n+1; i++)
    {
        ifail=AI[i]-AI_standard[i];
        if (ifail!=0) goto FAILURE8;
    }
    for ( i=0; i < nzmax; i++)
    {
        ifail=AJ[i]-AJ_standard[i];
        if (ifail!=0) goto FAILURE8;
    }

    for ( i=0; i < nzmax; i++)
    {
        rfail=Acsr[i]-Acsr_standard[i];
        if (rfail!=0) goto FAILURE6;
    };
//***************************************************************************************************
//* TASK 9    Obtain compressed sparse column matrix from compressed sparse row matrix
//***************************************************************************************************
    job[0] = 0;
    job[5] = 1;

    mkl_scsrcsc(job, &n, Acsr, AJ, AI, Acsc, AJ1, AI1, &info);

    if (info!=0) goto FAILURE9;
    for ( i=0; i< n+1; i++)
    {
        ifail=AI1[i]-AI_standard[i];
        if (ifail!=0) goto FAILURE9;
    }
    for ( i=0; i < nzmax; i++)
    {
        ifail=AJ1[i]-AJ_standard[i];
        if (ifail!=0) goto FAILURE9;
    }
    for ( i=0; i < nzmax; i++)
    {
        rfail=Acsc[i]-Acsc_standard[i];
        if (rfail!=0) goto FAILURE9;
    };
//***************************************************************************************************
//* TASK 10    Obtain compressed sparse row matrix from compressed sparse column matrix
//***************************************************************************************************
    job[0] = 1;
    job[5] = 1;
    mkl_scsrcsc(job, &n, Acsr, AJ, AI, Acsc, AJ1, AI1, &info);

    if (info!=0) goto FAILURE10;
    for ( i=0; i< n+1; i++)
    {
        ifail=AI[i]-AI_standard[i];
        if (ifail!=0) goto FAILURE10;
    }
    for ( i=0; i < nzmax; i++)
    {
        ifail=AJ[i]-AJ_standard[i];
        if (ifail!=0) goto FAILURE10;
    }

    for ( i=0; i < nzmax; i++)
    {
        rfail=Acsr[i]-Acsr_standard[i];
        if (rfail!=0) goto FAILURE10;
    };
//***************************************************************************************************
//* TASK 11    Obtain skyline matrix for lower triangle from compressed sparse row matrix
//***************************************************************************************************
    job[0] = 0;
    job[3] = 0;
    job[4] = 6;
    job[5] = 0;

    mkl_scsrsky(job, &n, Acsr, AJ, AI, Asky, pointers, &info);

    if (info!=0) goto FAILURE11;
    for ( i=0; i< n+1; i++)
    {
        ifail=pointers[i]-pointersl_standard[i];
        if (ifail!=0) goto FAILURE11;
    }
    nnz = pointers[n] - pointers[0];
    for ( i=0; i < nnz; i++)
    {
        rfail=Asky[i]-Askyl_standard[i];
        if (rfail!=0) goto FAILURE11;
    };
//***************************************************************************************************
//* TASK 12   Obtain compressed sparse row matrix for lower triangle from skyline matrix
//***************************************************************************************************
    job[0] = 1;
    job[3] = 0;
    job[5] = 0;

    mkl_scsrsky(job, &n, Acsr, AJ, AI, Asky, pointers, &info);

    if (info!=0) goto FAILURE12;
    for ( i=0; i< n+1; i++)
    {
        ifail=AI[i]-pointersl_standard[i];
        if (ifail!=0) goto FAILURE12;
    }
    nnz = pointers[n] - pointers[0];
    for ( i=0; i < nnz; i++)
    {
        ifail=AJ[i]-AJL_standard[i];
        if (ifail!=0) goto FAILURE12;
    }

    for ( i=0; i < nnz; i++)
    {
        rfail=Acsr[i]-Askyl_standard[i];
        if (rfail!=0) goto FAILURE12;
    }
//****************************************************************************************************
    printf("\n           All tests passed  \n");
    return 0;
/* Failure message to print if something went wrong */
FAILURE1: printf("\n Example FAILED to convert from dns to csr...\n");
    return 1;
FAILURE2: printf("\n Example FAILED to convert from csr to dns...\n");
    return 2;
FAILURE3: printf("\n Example FAILED to convert from csr to coo...\n");
    return 3;
FAILURE4: printf("\n Example FAILED to convert from coo to csr...\n");
    return 4;
FAILURE5: printf("\n Example FAILED to convert from csr to bsr...\n");
    return 5;
FAILURE6: printf("\n Example FAILED to convert from bsr to csr...\n");
    return 6;
FAILURE7: printf("\n Example FAILED to convert from csr to dia...\n");
    return 7;
FAILURE8: printf("\n Example FAILED to convert from dia to csr...\n");
    return 8;
FAILURE9: printf("\n Example FAILED to convert from csr to csc...\n");
    return 9;
FAILURE10: printf("\n Example FAILED to convert from csc to csr...\n");
    return 10;
FAILURE11: printf("\n Example FAILED to convert from csr to sky...\n");
    return 11;
FAILURE12: printf("\n Example FAILED to convert from sky to csr...\n");
    return 12;
   }
