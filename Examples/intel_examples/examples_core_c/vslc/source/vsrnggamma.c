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
!  Content:
!    vsRngGamma  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    7777777
#define BRNG    VSL_BRNG_MT19937
#define METHOD  VSL_RNG_METHOD_GAMMA_GNORM
#define N       1000
#define NN      10

int main()
{
    float r[N];
    VSLStreamStatePtr stream;
    int i, errcode;
    float alpha = 2.0, a=0.0, beta=1.0;

    float tM,tD,tQ,tD2;
    float sM,sD;
    float sum, sum2;
    float n,s;
    float DeltaM,DeltaD;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, BRNG,  SEED );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = vsRngGamma( METHOD, stream, N, r, alpha, a, beta );
    CheckVslError( errcode );

    /***** Theoretical moments *****/
    tM=a+beta * alpha;
    tD=beta * beta * alpha;
    tQ=beta * beta * beta * beta * 3 * alpha * (alpha+2);

    /***** Sample moments *****/
    sum=0.0;
    sum2=0.0;
    for(i=0;i<N;i++) {
        sum+=(float)r[i];
        sum2+=(float)r[i]*(float)r[i];
    }
    sM=sum/((float)N);
    sD=sum2/(float)N-(sM*sM);

    /***** Comparison of theoretical and sample moments *****/
    n=(float)N;
    tD2=tD*tD;
    s=((tQ-tD2)/n)-(2*(tQ-2*tD2)/(n*n))+((tQ-3*tD2)/(n*n*n));

    DeltaM=(tM-sM)/sqrt(tD/n);
    DeltaD=(tD-sD)/sqrt(s);

    /***** Printing results *****/
    printf("Sample of vsRngGamma.\n");
    printf("----------------------\n\n");
    printf("Parameters:\n");
    printf("    alpha=%.4f\n",alpha);
    printf("    a=%.4f\n",a);
    printf("    beta=%.4f\n\n",beta);

    printf("Results (first 10 of 1000):\n");
    printf("---------------------------\n");
    for(i=0;i<NN;i++) {
        printf("r[%d]=%.4f\n",i,r[i]);
    }

    printf("\n");
    if(fabs(DeltaM)>3.0 || fabs(DeltaD)>3.0) {
        printf("Error: sample moments (mean=%.2f, variance=%.2f) disagree with theory (mean=%.2f, variance=%.2f).\n",sM,sD,tM,tD);
        return 1;
    }
    else {
        printf("Sample moments (mean=%.2f, variance=%.2f) agree with theory (mean=%.2f, variance=%.2f).\n",sM,sD,tM,tD);
    }

    printf("\n");

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    return 0;
}
