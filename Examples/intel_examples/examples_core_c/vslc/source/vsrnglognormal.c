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
!    vsRngLognormal  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    1
#define BRNG    VSL_BRNG_MCG31
#define METHOD  VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2
#define N       1000
#define NN      10

int main()
{
    float r[N];
    VSLStreamStatePtr stream;
    int i, errcode;
    float a=0.0,sigma=1.0,b=0.0,beta=1.0;

    double tM,tD,tQ,tD2;
    double sM,sD;
    double sum, sum2;
    double n,s;
    double DeltaM,DeltaD;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, BRNG,  SEED );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = vsRngLognormal( METHOD, stream, N, r, a, sigma, b, beta );
    CheckVslError( errcode );

    /***** Theoretical moments *****/
    tM=b+beta*exp(a+sigma*sigma*0.5);
    tD=beta*beta*exp(2.0*a+sigma*sigma)*(exp(sigma*sigma)-1.0);
    tQ=beta*beta*beta*beta*exp(4.0*a+2.0*sigma*sigma)*(exp(6.0*sigma*sigma)-4.0*exp(3.0*sigma*sigma)+6.0*exp(sigma*sigma)-3.0);

    /***** Sample moments *****/
    sum=0.0;
    sum2=0.0;
    for(i=0;i<N;i++) {
        sum+=(double)r[i];
        sum2+=(double)r[i]*(double)r[i];
    }
    sM=sum/((double)N);
    sD=sum2/(double)N-(sM*sM);

    /***** Comparison of theoretical and sample moments *****/
    n=(double)N;
    tD2=tD*tD;
    s=((tQ-tD2)/n)-(2*(tQ-2*tD2)/(n*n))+((tQ-3*tD2)/(n*n*n));

    DeltaM=(tM-sM)/sqrt(tD/n);
    DeltaD=(tD-sD)/sqrt(s);

    /***** Printing results *****/
    printf("Sample of vsRngLognormal.\n");
    printf("-------------------------\n\n");
    printf("Parameters:\n");
    printf("    a=%.4f\n",a);
    printf("    sigma=%.4f\n",sigma);
    printf("    b=%.4f\n",b);
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

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    return 0;
}
