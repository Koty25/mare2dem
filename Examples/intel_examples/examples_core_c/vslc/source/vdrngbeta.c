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
!    vdRngBeta  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    777
#define BRNG    VSL_BRNG_MCG31
#define METHOD  VSL_RNG_METHOD_BETA_CJA
#define N       1000
#define NN      10

int main()
{
    double r[N];
    VSLStreamStatePtr stream;
    int i, errcode;
    double p=2.0, q=10.0, a=0.0, beta=1.0;

    double tM,tD,tQ,tD2;
    double sM,sD;
    double sum, sum2;
    double b, c, d, e, e2, b2, sum_pq;
    double n,s;
    double DeltaM,DeltaD;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, BRNG,  SEED );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = vdRngBeta( METHOD, stream, N, r, p, q, a, beta );
    CheckVslError( errcode );

    /***** Theoretical moments *****/
    b2 = beta * beta;
    sum_pq = p + q;
    b = (p + 1.0) / (sum_pq + 1.0);
    c = (p + 2.0) / (sum_pq + 2.0);
    d = (p + 3.0) / (sum_pq + 3.0);
    e =  p / sum_pq;
    e2 = e * e;

    tM = a + e * beta;
    tD = b2 * p * q / (sum_pq * sum_pq * (sum_pq + 1.0));
    tQ = b2 * b2 * (e * b * c * d - 4.0 * e2 * b * c + 6.0 * e2 * e * b - 3.0 * e2 * e2);

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
    printf("Sample of vdRngBeta.\n");
    printf("----------------------\n\n");
    printf("Parameters:\n");
    printf("    p=%.4f\n",p);
    printf("    q=%.4f\n",q);
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
