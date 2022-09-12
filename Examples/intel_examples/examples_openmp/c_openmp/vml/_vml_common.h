/*******************************************************************************
* Copyright 2017-2020 Intel Corporation.
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
#ifndef __VML_COMMON_H__
#define __VML_COMMON_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <omp.h>

#include "mkl.h"
#include "mkl_omp_offload.h"

#ifdef __cplusplus
#include <complex>
using namespace std::complex_literals;
#define VM_COMPLEX8  std::complex<float>
#define VM_COMPLEX16 std::complex<double>
#define crealf(z) (z.real())
#define cimagf(z) (z.imag())
#define creal(z) (z.real())
#define cimag(z) (z.imag())
#define I (1i)
std::complex<float>
operator *(std::complex<double> lhs, float rhs) {
    auto a = lhs * static_cast<double>(rhs);
    return std::complex<float>(a.real(), a.imag());
}

#else
#include <complex.h>
#define VM_COMPLEX8  float complex
#define VM_COMPLEX16 double complex
#endif

#define ACCURACY_HA  0
#define ACCURACY_LA  1
#define ACCURACY_EP  2
#define ACCURACY_NUM 3

#define ARG1_RES1    11
#define ARG2_RES1    21
#define ARG1_RES2    12
#define ARG6_RES1    61
#define ARG1R_RES1C  1011
#define ARG1C_RES1R  1110

// accuracy names
static const char *const accuracy_name[] = { "HA",   
                                             "LA",   
                                             "EP"};
// accuracy modes
static const uint64_t    accuracy_mode[] = { VML_HA, 
                                             VML_LA, 
                                             VML_EP};

// maximum allowed ulps
typedef struct  
{ 
    float  f[3];
    double d[3];
    float  fc[3];
    double dc[3];
}  max_ulp_table_t;

extern max_ulp_table_t max_ulp_table;

//! ===========================================================================
//!
//! @brief Ulp calculation routines
//!
//! ===========================================================================
// float ulp calculation
float ulp_float (float res, float ref)
{
    float resulp = 0.0f;
    const int res_class = fpclassify (res);
    const int ref_class = fpclassify (ref);

    if (isfinite (res) && isfinite (ref))
    {
        int ex;
        double den;

        frexp (ref, &ex);
        den = ldexp (1.0, ex - 24);
        den = (den == 0.0) ? 0x1.p-149 : den;

        resulp = (float)(fabs (((double) (res - ref)) / den));
        if (!(isfinite (resulp))) { resulp = FLT_MAX; }
    }
    else
    {
        if (res_class == ref_class) { resulp = 0.0f; }
        else { resulp = FLT_MAX; }
    }

    return resulp;
}

// double ulp calculation
double ulp_double (double res, double ref)
{
    double resulp = 0.0;
    const int res_class = fpclassify (res);
    const int ref_class = fpclassify (ref);

    if (isfinite (res) && isfinite (ref))
    {
        int ex;
        double den;

        frexp (ref, &ex);
        den = ldexp (1.0, ex - 53);
        den = (den == 0.0) ? 0x1.p-1074 : den;

        resulp = fabs (((double) (res - ref)) / den);
        if (!(isfinite (resulp))) { resulp = DBL_MAX; }
    }
    else
    {
        if (res_class == ref_class) { resulp = 0.0; }
        else { resulp = DBL_MAX; }
    }

    return resulp;
}

//! ===========================================================================
//!
//! @brief Result check routines
//!
//! ===========================================================================
// check float result
int check_result_float (int argtype, float arg1, float arg2, float res1, float res2, float ref1, float ref2, const char* funcname, int accuracy)
{
    int is_error = 0;
    const char* error_message = "";
    float resulp = ulp_float(res1, ref1);
    if(argtype == ARG1_RES2)
    {
	resulp = fmaxf(resulp, ulp_float(res2, ref2));
    }

    is_error = (resulp >= max_ulp_table.f[accuracy]);
    error_message = (is_error)?"FAIL":"PASS";

    if     (argtype == ARG2_RES1) { printf("\t\t%s{%s} (%-12.7g, %-12.7g) = %-12.7g; expected = %-12.7g; ulp = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], arg1, arg2, res1, ref1, resulp, error_message);  }
    else if(argtype == ARG1_RES2) { printf("\t\t%s{%s} (%-12.7g) = %-12.7g, %-12.7g; expected = %-12.7g, %-12.7g; ulp = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], arg1, res1, res2, ref1, ref2, resulp, error_message); }
    else                          { printf("\t\t%s{%s} (%-12.7g) = %-12.7g; expected = %-12.7g; ulp = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], arg1, res1, ref1, resulp, error_message); }
    
    return is_error;
}

// check double result
int check_result_double (int argtype, double arg1, double arg2, double res1, double res2, double ref1, double ref2, const char* funcname, int accuracy)
{
    int is_error = 0;
    const char* error_message = "";
    double resulp = ulp_double(res1, ref1);
    if(argtype == ARG1_RES2)
    {
	resulp = fmax(resulp, ulp_double(res2, ref2));
    }

    is_error = (resulp > max_ulp_table.d[accuracy]);
    error_message = (is_error)?"FAIL":"PASS";

    if     (argtype == ARG2_RES1) { printf("\t\t%s{%s} (%-22.16g, %-22.16g) = %-22.16g; expected = %-22.16g; ulp = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], arg1, arg2, res1, ref1, resulp, error_message); }
    else if(argtype == ARG1_RES2) { printf("\t\t%s{%s} (%-22.16g) = %-22.16g, %-22.16g; expected = %-22.16g, %-22.16g; ulp = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], arg1, res1, res2, ref1, ref2, resulp, error_message); }
    else                          { printf("\t\t%s{%s} (%-22.16g) = %-22.16g; expected = %-22.16g; ulp = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], arg1, res1, ref1, resulp, error_message); }
        
    return is_error;
}

// check complex float result
int check_result_float_complex (int argtype, VM_COMPLEX8 arg1, VM_COMPLEX8 arg2, VM_COMPLEX8 res1, VM_COMPLEX8 res2, VM_COMPLEX8 ref1, VM_COMPLEX8 ref2, const char* funcname, int accuracy)
{
    int is_error = 0;
    const char* error_message = "";
    float resulp = ulp_float(crealf(res1), crealf(ref1));
    resulp = fmaxf(resulp, ulp_float(cimagf(res1), cimagf(ref1)));
    if(argtype == ARG1_RES2)
    {
        resulp = fmaxf(resulp, ulp_float(crealf(res2), crealf(ref2)));
        resulp = fmaxf(resulp, ulp_float(cimagf(res2), cimagf(ref2)));
    }
    
    is_error = (resulp > max_ulp_table.fc[accuracy]);
    error_message = (is_error)?"FAIL":"PASS";

    if     (argtype == ARG2_RES1) { printf("\t\t%s{%s} (%-.7g+i*%-.7g, %-.7g+i*%-15.7g)\n\t\t\t         = %-15.7g+i*%-15.7g;\n\t\t\texpected = %-15.7g+i*%-15.7g;\n\t\t\tulp      = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], crealf(arg1), cimagf(arg1), crealf(arg2), cimagf(arg2), crealf(res1), cimagf(res1), crealf(ref1), cimagf(ref1), resulp, error_message); }
    else if(argtype == ARG1_RES2) { printf("\t\t%s{%s} (%-.7g+i*%-.7g)\n\t\t\t         = %-15.7g+i*%-15.7g, %-15.7g+i*%-15.7g;\n\t\t\texpected = %-15.7g+i*%-15.7g, %-15.7g+i*%-15.7g;\n\t\t\tulp      = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], crealf(arg1), cimagf(arg1), crealf(res1), cimagf(res1), crealf(res2), cimagf(res2), crealf(ref1), cimagf(ref1), crealf(ref2), cimagf(ref2), resulp, error_message); }
    else                          { printf("\t\t%s{%s} (%-.7g+i*%-.7g)\n\t\t\t         = %-15.7g+i*%-15.7g;\n\t\t\texpected = %-15.7g+i*%-15.7g;\n\t\t\tulp      = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], crealf(arg1), cimagf(arg1), crealf(res1), cimagf(res1), crealf(ref1), cimagf(ref1), resulp, error_message); }
       
    return is_error;
}

// check complex double result
int check_result_double_complex (int argtype, VM_COMPLEX16 arg1, VM_COMPLEX16 arg2, VM_COMPLEX16 res1, VM_COMPLEX16 res2, VM_COMPLEX16 ref1, VM_COMPLEX16 ref2, const char* funcname, int accuracy)
{
    int is_error = 0;
    const char* error_message = "";
    double resulp = ulp_double(creal(res1), creal(ref1));
    resulp = fmax(resulp, ulp_double(cimag(res1), cimag(ref1)));
    if(argtype == ARG1_RES2)
    {
        resulp = fmax(resulp, ulp_double(creal(res2), creal(ref2)));
        resulp = fmax(resulp, ulp_double(cimag(res2), cimag(ref2)));
    }
    
    is_error = (resulp > max_ulp_table.dc[accuracy]);
    error_message = (is_error)?"FAIL":"PASS";

    if     (argtype == ARG2_RES1) { printf("\t\t%s{%s} (%-.14g+i*%-.14g, %-.14g+i*%-.14g)\n\t\t\t         = %-22.14g+i*%-22.14g;\n\t\t\texpected = %-22.14g+i*%-22.14g;\n\t\t\tulp      = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], creal(arg1), cimag(arg1), creal(arg2), cimag(arg2), creal(res1), cimag(res1), creal(ref1), cimag(ref1), resulp, error_message); }
    else if(argtype == ARG1_RES2) { printf("\t\t%s{%s} (%-.14g+i*%-.14g)\n\t\t\t         = %-22.14g+i*%-22.14g, %-22.14g+i*%-22.14g;\n\t\t\texpected = %-22.14g+i*%-22.14g, %-22.14g+i*%-22.14g;\n\t\t\tulp      = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], creal(arg1), cimag(arg1), creal(res1), cimag(res1), creal(res2), cimag(res2), creal(ref1), cimag(ref1), creal(ref2), cimag(ref2), resulp, error_message); }
    else                          { printf("\t\t%s{%s} (%-.14g+i*%-.14g)\n\t\t\t         = %-22.14g+i*%-22.14g;\n\t\t\texpected = %-22.14g+i*%-22.14g;\n\t\t\tulp      = %-7.2g (%s)\n", \
                                        funcname, accuracy_name[accuracy], creal(arg1), cimag(arg1), creal(res1), cimag(res1), creal(ref1), cimag(ref1), resulp, error_message); }
       
    return is_error;

}

//! ===========================================================================
//!
//! @brief Data types for tables
//!
//! ===========================================================================
// data types for combined floating point\integer access to floating point values
typedef union float_union_t 
{
    uint32_t u;
    float f;
} float_union_t;

typedef union double_union_t
{
    uint64_t u;
    double f;
} double_union_t;

typedef union complex_float_union_t
{
    struct {uint32_t real; uint32_t imag;} u;
    VM_COMPLEX8 f;
} complex_float_union_t;

typedef union complex_double_union_t
{
    struct {uint64_t real; uint64_t imag;} u;
    VM_COMPLEX16 f;
} complex_double_union_t;

// two float values combined type
typedef struct data_2_f32_t
{
    float_union_t v1;
    float_union_t v2;
} data_2_f32_t;

// two double values combined type
typedef struct data_2_f64_t
{
    double_union_t v1;
    double_union_t v2;
} data_2_f64_t;

// two complex float values combined type
typedef struct data_2_c32_t
{
    complex_float_union_t v1;
    complex_float_union_t v2;
} data_2_c32_t;

// two complex double values combined type
typedef struct data_2_c64_t
{
    complex_double_union_t v1;
    complex_double_union_t v2;
} data_2_c64_t;

// three float values combined type
typedef struct data_3_f32_t
{
    float_union_t v1;
    float_union_t v2;
    float_union_t v3;
} data_3_f32_t;

// three double values combined type
typedef struct data_3_f64_t
{
    double_union_t v1;
    double_union_t v2;
    double_union_t v3;
} data_3_f64_t;

// three complex float values combined type
typedef struct data_3_c32_t
{
    complex_float_union_t v1;
    complex_float_union_t v2;
    complex_float_union_t v3;
} data_3_c32_t;

// three complex double values combined type
typedef struct data_3_c64_t
{
    complex_double_union_t v1;
    complex_double_union_t v2;
    complex_double_union_t v3;
} data_3_c64_t;

// two float\complex float values combined type
typedef struct data_2_cf32_t
{
    complex_float_union_t v1;
    float_union_t         v2;
} data_2_cf32_t;

// two double\complex double values combined type
typedef struct data_2_cf64_t
{
    complex_double_union_t v1;
    double_union_t         v2;
} data_2_cf64_t;

// two complex float\float values combined type
typedef struct data_2_fc32_t
{
    float_union_t         v1;
    complex_float_union_t v2;
} data_2_fc32_t;

// two complex double\double values combined type
typedef struct data_2_fc64_t
{
    double_union_t         v1;
    complex_double_union_t v2;
} data_2_fc64_t;

// seven float values combined type
typedef struct data_7_f32_t
{
    float_union_t v1;
    float_union_t v2;
    float_union_t v3;
    float_union_t v4;
    float_union_t v5;
    float_union_t v6;
    float_union_t v7;
} data_7_f32_t;

// seven double values combined type
typedef struct data_7_f64_t
{
    double_union_t v1;
    double_union_t v2;
    double_union_t v3;
    double_union_t v4;
    double_union_t v5;
    double_union_t v6;
    double_union_t v7;
} data_7_f64_t;

// reference data table structure
// for 1 arg 1 res functions
typedef struct data_2_t
{
    data_2_f32_t data_f32[VLEN];
    data_2_f64_t data_f64[VLEN];
    data_2_c32_t data_c32[VLEN];
    data_2_c64_t data_c64[VLEN];
} data_2_t; // struct data_2_t

// reference data table structure
// for 2 arg 1 res or 1 arg 2 res functions
typedef struct data_3_t
{
    data_3_f32_t data_f32[VLEN];
    data_3_f64_t data_f64[VLEN];
    data_3_c32_t data_c32[VLEN];
    data_3_c64_t data_c64[VLEN];
} data_3_t; // struct data_3_t

// reference data table structure
// for 1 arg 1 res complex to real functions
typedef struct data_2cf_t
{
    data_2_f32_t  data_f32[VLEN];
    data_2_f64_t  data_f64[VLEN];
    data_2_cf32_t data_c32[VLEN];
    data_2_cf64_t data_c64[VLEN];
} data_2cf_t; // struct data_2cf_t

// reference data table structure
// for 1 arg 1 res real to complex functions
typedef struct data_2fc_t
{
    data_2_fc32_t data_c32[VLEN];
    data_2_fc64_t data_c64[VLEN];
} data_2fc_t; // struct data_2fc_t

// reference data table structure
// for 6 arg 1 res functions
typedef struct data_7_t
{
    data_7_f32_t data_f32[VLEN];
    data_7_f64_t data_f64[VLEN];
} data_7_t; // struct data_7_t


#endif // __VML_COMMON_H__
