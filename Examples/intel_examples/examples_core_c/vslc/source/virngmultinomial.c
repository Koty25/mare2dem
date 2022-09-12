/*******************************************************************************
* Copyright 2018-2020 Intel Corporation.
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
!    viRngMultinomial  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    1
#define BRNG    VSL_BRNG_MCG31
#define METHOD  VSL_RNG_METHOD_MULTINOMIAL_MULTPOISSON
#define N       1000
#define NTRIAL  10
#define K       5
#define NN      10

int main()
{
    int r[N*K];
    int ntrial=NTRIAL;
    int k=K;
    VSLStreamStatePtr stream;
    int i, j, errcode;

    double p[K]={0.2, 0.2, 0.2, 0.2, 0.2};

    double tM,tD,tQ,tD2;
    double sM,sD;
    double sum, sum2;
    double n,s;
    double DeltaM,DeltaD;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, BRNG,  SEED );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = viRngMultinomial( METHOD, stream, N, r, ntrial, k, p );
    CheckVslError( errcode );


    /***** Printing results *****/
    printf("Sample of viRngMultinomial.\n");
    printf("------------------------\n\n");
    printf("Parameters:\n");
    printf("    ntrial=%d\n",ntrial);
    printf("    k=%d\n",k);
    printf("    probability vector:\n");
    for(i=0;i<k;i++) {
        printf("    %.4f",p[i]);
    }


    printf("\n\nResults (first 10 of 1000):\n");
    printf("---------------------------\n");
    for(i=0;i<NN;i++) {
        for(j=0;j<k;j++) {
            printf("r[%d][%d]=%d  ",i,j,r[i*k+j]);
        }
        printf("\n");
    }

    for(i=0;i<k;i++)
    {
        /***** Theoretical moments *****/
        tM=ntrial*p[i];
        tD=ntrial*p[i]*(1.0-p[i]);
        tQ=ntrial*(1.0-p[i])*(4.0*ntrial*p[i]-4.0*ntrial*p[i]*p[i]+4.0*p[i]-3.0);

        /***** Sample moments *****/
        sum=0.0;
        sum2=0.0;
        for(j=0;j<N;j++) {
            sum+=(double)r[j*k+i];
            sum2+=(double)r[j*k+i]*(double)r[j*k+i];
        }
        sM=sum/((double)N);
        sD=sum2/(double)N-(sM*sM);

        /***** Comparison of theoretical and sample moments *****/
        n=(double)N;
        tD2=tD*tD;
        s=((tQ-tD2)/n)-(2*(tQ-2*tD2)/(n*n))+((tQ-3*tD2)/(n*n*n));

        DeltaM=(tM-sM)/sqrt(tD/n);
        DeltaD=(tD-sD)/sqrt(s);

        printf("\n");
        if(fabs(DeltaM)>3.0 || fabs(DeltaD)>3.0) {
            printf("Error: Category=%d sample moments (mean=%.2f, variance=%.2f) disagree with theory (mean=%.2f, variance=%.2f).\n",i,sM,sD,tM,tD);
            return 1;
        }
        else {
            printf("Category=%d Sample moments (mean=%.2f, variance=%.2f) agree with theory (mean=%.2f, variance=%.2f).\n",i,sM,sD,tM,tD);
        }
    }
    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    return 0;
}
