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
!    viRngHypergeometric  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    1
#define BRNG    VSL_BRNG_MCG31
#define METHOD  VSL_RNG_METHOD_HYPERGEOMETRIC_H2PE
#define N       1000
#define NN      10

int main()
{
    int r[N];
    VSLStreamStatePtr stream;
    int i, errcode;
    int l=100,ss=10,m=30;

    double tM,tD,tQ,tD2;
    double sM,sD;
    double sum, sum2;
    double n,s;
    double DeltaM,DeltaD;
    double K, L2, L3, L4, L5, L6, KL, KL4, S2, S3, S4, M2, M3, M4;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, BRNG,  SEED );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = viRngHypergeometric( METHOD, stream, N, r, l, ss, m );
    CheckVslError( errcode );

    /***** Theoretical moments *****/
    K = (l-1)*(l-2)*(l-3);
    L2 = l*l;
    L3 = L2*l;
    L4 = L2*L2;
    L5 = L3*L2;
    L6 = L3*L3;
    KL = K*l;
    KL4 = K*L4;
    S2 = ss*ss;
    S3 = S2*ss;
    S4 = S2*S2;
    M2 = m*m;
    M3 = M2*m;
    M4 = M2*M2;

    tM=(double)m*(double)ss/(double)l;
    tD=(double)(m*ss*(l-m)*(l-ss))/(double)(l*l*(l-1));
    tQ=( (3*l+18)    *S4/KL4 - (6*L2+36*l)  *S3/KL4 + (3*L3+24*L2)   *S2/KL4 - 6        *ss/KL  ) * M4 +
       ( (-6*L2-36*l)*S4/KL4 + (12*L3+72*L2)*S3/KL4 - (6*L4+38*L3)   *S2/KL4 + 12       *ss/K   ) * M3 +
       ( (3*L3+24*L2)*S4/KL4 - (6*L4+48*L3) *S3/KL4 + (31*L4+3*L5+L3)*S2/KL4 - (L4+7*L5)*ss/KL4 ) * M2 +
       ( -6          *S4/KL  + 12           *S3/K   - (4*L4+7*L5)    *S2/KL4 + (L6+L5)  *ss/KL4 ) * m;

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
    printf("Sample of viRngHypergeometric.\n");
    printf("------------------------------\n\n");
    printf("Parameters:\n");
    printf("    l=%d\n",l);
    printf("    s=%d\n",ss);
    printf("    m=%d\n\n",m);

    printf("Results (first 10 of 1000):\n");
    printf("---------------------------\n");
    for(i=0;i<NN;i++) {
        printf("r[%d]=%d\n",i,r[i]);
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
