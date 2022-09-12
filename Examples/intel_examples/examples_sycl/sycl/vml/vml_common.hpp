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
#pragma once

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cfloat>
#include <ctime>
#include <map>

#define ACCURACY_HA  0
#define ACCURACY_LA  1
#define ACCURACY_EP  2
#define ACCURACY_NUM 3

// Accuracy names
static const char *const accuracy_name[] = {"HA",   "LA",   "EP"};
// Accuracy modes
static const oneapi::mkl::vm::mode accuracy_mode[] = {
                                            oneapi::mkl::vm::mode::ha,
                                            oneapi::mkl::vm::mode::la,
                                            oneapi::mkl::vm::mode::ep,
                                            };
using ulp_table_type = std::map<int, double>;

enum {
    MAX_HA_ULP_S = 0,
    MAX_LA_ULP_S,
    MAX_EP_ULP_S,

    MAX_HA_ULP_D,
    MAX_LA_ULP_D,
    MAX_EP_ULP_D,

    MAX_HA_ULP_C,
    MAX_LA_ULP_C,
    MAX_EP_ULP_C,

    MAX_HA_ULP_Z,
    MAX_LA_ULP_Z,
    MAX_EP_ULP_Z,
} ulp_table_index;


//! ===========================================================================
//! @brief Get name of data type.
//!
//! @return  Name of data type
//! ===========================================================================
template <typename T> static std::string get_typename(T arg){ }
template <> std::string get_typename(float arg){  std::string name = "Single"; return name; }
template <> std::string get_typename(double arg){ std::string name = "Double"; return name; }
template <> std::string get_typename(std::complex <float> arg){  std::string name = "Complex Single"; return name; }
template <> std::string get_typename(std::complex <double> arg){ std::string name = "Complex Double"; return name; }

//! ===========================================================================
//! @brief Ulp computation routine.
//!
//! @param[in] res     Computed result
//! @param[in] ref     Reference result
//!
//! @return          Ulp value between computed and reference results
//! ===========================================================================
template <typename T> static T compute_ulp(T res, T ref){ }
// float specialization
template <> float compute_ulp (float res, float ref)
{
  float retval;
  const int res_class = std::fpclassify (res);
  const int ref_class = std::fpclassify (ref);

  if (std::isfinite (res) && std::isfinite (ref))
    {
      int ex;
      double den;

      std::frexp (ref, &ex);
      den = std::ldexp (1.0, ex - 24);
      den = (den == 0.0) ? 0x1.p-149 : den;

      retval = (float)(std::fabs (((double) (res - ref)) / den));
      if (!(std::isfinite (retval))) retval = FLT_MAX;
    }
  else
    {
      if (res_class == ref_class)
    {
      retval = 0.0f;
    }
      else
    {
      retval = FLT_MAX;
    }
    }
  return retval;
}
// double specialization
template <> double compute_ulp (double res, double ref)
{
  double retval;
  const int res_class = std::fpclassify (res);
  const int ref_class = std::fpclassify (ref);

  if (std::isfinite (res) && std::isfinite (ref))
    {
      int ex;
      double den;

      std::frexp (ref, &ex);
      den = std::ldexp (1.0, ex - 53);
      den = (den == 0.0) ? 0x1.p-1074 : den;

      retval = std::fabs (((double) (res - ref)) / den);
      if (!(std::isfinite (retval))) retval = DBL_MAX;
    }
  else
    {
      if (res_class == ref_class)
    {
      retval = 0.0;
    }
      else
    {
      retval = DBL_MAX;
    }
    }
  return retval;
}
// complex float specialization
template <> std::complex <float> compute_ulp (std::complex <float>res, std::complex <float>ref)
{

  std::complex <float> retval;
  const int res_class_re = std::fpclassify (res.real ());
  const int ref_class_re = std::fpclassify (ref.real ());
  const int res_class_im = std::fpclassify (res.imag ());
  const int ref_class_im = std::fpclassify (ref.imag ());

  if (std::isfinite (res.real ()) && std::isfinite (ref.real ())
      && std::isfinite (res.imag ()) && std::isfinite (ref.imag ()))

    {
      int ex;
      double den;
      float ulp_re, ulp_im;

      std::frexp (ref.real (), &ex);
      den = std::ldexp (1.0, ex - 24);
      den = (den == 0.0) ? 0x1.p-149 : den;
      ulp_re = std::fabs (((double) (res.real () - ref.real ())) / den);
      if (!(std::isfinite (ulp_re))) ulp_re = FLT_MAX;

      std::frexp (ref.imag (), &ex);
      den = std::ldexp (1.0, ex - 24);
      den = (den == 0.0) ? 0x1.p-149 : den;
      ulp_im = std::fabs (((double) (res.imag () - ref.imag ())) / den);
      if (!(std::isfinite (ulp_im))) ulp_im = FLT_MAX;

      retval = { ulp_re, ulp_im };
    }
  else
    {
      if ((res_class_re == ref_class_re) && (res_class_im == ref_class_im))

    {
      retval = { 0.0f, 0.0f };
    }
      else
    {
      retval = { FLT_MAX, FLT_MAX };
    }
    }

  return retval;
}
// complex double specialization
template <> std::complex <double> compute_ulp (std::complex <double>res, std::complex <double>ref)
{
  std::complex <double> retval;
  const int res_class_re = std::fpclassify (res.real ());
  const int ref_class_re = std::fpclassify (ref.real ());
  const int res_class_im = std::fpclassify (res.imag ());
  const int ref_class_im = std::fpclassify (ref.imag ());

  if (std::isfinite (res.real ()) && std::isfinite (ref.real ())
      && std::isfinite (res.imag ()) && std::isfinite (ref.imag ()))

    {
      int ex;
      double den;
      double ulp_re, ulp_im;

      std::frexp (ref.real (), &ex);
      den = std::ldexp (1.0, ex - 53);
      den = (den == 0.0) ? 0x1.p-1074 : den;
      ulp_re = std::fabs (((double) (res.real () - ref.real ())) / den);
      if (!(std::isfinite (ulp_re))) ulp_re = DBL_MAX;

      std::frexp (ref.imag (), &ex);
      den = std::ldexp (1.0, ex - 53);
      den = (den == 0.0) ? 0x1.p-1074 : den;
      ulp_im = std::fabs (((double) (res.imag () - ref.imag ())) / den);
      if (!(std::isfinite (ulp_im))) ulp_im = DBL_MAX;

      retval = { ulp_re, ulp_im };
    }
  else
    {
      if ((res_class_re == ref_class_re) && (res_class_im == ref_class_im))

    {
      retval = { 0.0, 0.0 };
    }
      else
    {
      retval = { DBL_MAX, DBL_MAX };
    }
    }

  return retval;
}


//! ===========================================================================
//! @brief Check ulp routine.
//!
//! @param[in] ulp     Computed ulp
//! @param[in] exact   Check for exact match
//! @param[in] acc     Accuracy mode
//!
//! @return            Returns false or true if computed ulp in allowed range
//! ===========================================================================
template <typename T> static bool check_ulp( T ulp,  bool exact, oneapi::mkl::vm::mode acc,  ulp_table_type const & ulp_table) { return true; }
// float specialization
template <> bool check_ulp( float ulp,  bool exact, oneapi::mkl::vm::mode acc, ulp_table_type const & ulp_table)
{
    float max_ulp = (acc == accuracy_mode[ACCURACY_EP])?ulp_table.at(MAX_EP_ULP_S):
                    (acc == accuracy_mode[ACCURACY_LA])?ulp_table.at(MAX_LA_ULP_S):
                                                        ulp_table.at(MAX_HA_ULP_S);
    if(exact == true) max_ulp = 0.0f;
    return (ulp <= max_ulp);
}
// double specialization
template <> bool check_ulp( double ulp,  bool exact, oneapi::mkl::vm::mode acc, ulp_table_type const & ulp_table)
{
    double max_ulp = (acc == accuracy_mode[ACCURACY_EP])?ulp_table.at(MAX_EP_ULP_D):
                     (acc == accuracy_mode[ACCURACY_LA])?ulp_table.at(MAX_LA_ULP_D):
                                                         ulp_table.at(MAX_HA_ULP_D);
    if(exact == true) max_ulp = 0.0;
    return (ulp <= max_ulp);
}
// complex float specialization
template <> bool check_ulp( std::complex <float> ulp,  bool exact, oneapi::mkl::vm::mode acc, ulp_table_type const & ulp_table)
{
    float max_ulp = (acc == accuracy_mode[ACCURACY_EP])?ulp_table.at(MAX_EP_ULP_C):
                    (acc == accuracy_mode[ACCURACY_LA])?ulp_table.at(MAX_LA_ULP_C):
                                                        ulp_table.at(MAX_HA_ULP_C);
    if(exact == true) max_ulp = 0.0f;
    return (std::fmax (ulp.real(), ulp.imag()) <= max_ulp);
}
// complex double specialization
template <> bool check_ulp( std::complex <double> ulp,  bool exact, oneapi::mkl::vm::mode acc, ulp_table_type const & ulp_table)
{
    double max_ulp = (acc == accuracy_mode[ACCURACY_EP])?ulp_table.at(MAX_EP_ULP_Z):
                     (acc == accuracy_mode[ACCURACY_LA])?ulp_table.at(MAX_LA_ULP_Z):
                                                         ulp_table.at(MAX_HA_ULP_Z);
    if(exact == true) max_ulp = 0.0;
    return (std::fmax (ulp.real(), ulp.imag()) <= max_ulp);
}

//! ===========================================================================
//!
//! @brief Data types for tables
//!
//! ===========================================================================
// Data types for combined floating point\integer access to floating point values
using float_union_t = union
        {
            uint32_t u;
            float f;
        };

using double_union_t = union
        {
            uint64_t u;
            double f;
        };

using complex_float_union_t = union
        {
            struct {uint32_t real; uint32_t imag;} u;
            struct {float real; float imag;}       f;
        };

using complex_double_union_t = union
        {
            struct {uint64_t real; uint64_t imag;} u;
            struct {double real; double imag;}     f;
        };

// Two float values combined type
struct data_2_f32_t
{
  float_union_t v1;
  float_union_t v2;
};

// Two double values combined type
struct data_2_f64_t
{
  double_union_t v1;
  double_union_t v2;
};

// Two complex float values combined type
struct data_2_c32_t
{
  complex_float_union_t v1;
  complex_float_union_t v2;
};

// Two complex double values combined type
struct data_2_c64_t
{
  complex_double_union_t v1;
  complex_double_union_t v2;
};

// Three float values combined type
struct data_3_f32_t
{
  float_union_t v1;
  float_union_t v2;
  float_union_t v3;
};

// Three double values combined type
struct data_3_f64_t
{
  double_union_t v1;
  double_union_t v2;
  double_union_t v3;
};

// Three complex float values combined type
struct data_3_c32_t
{
  complex_float_union_t v1;
  complex_float_union_t v2;
  complex_float_union_t v3;
};

// Three complex double values combined type
struct data_3_c64_t
{
  complex_double_union_t v1;
  complex_double_union_t v2;
  complex_double_union_t v3;
};

// Two float\complex float values combined type
struct data_2_cf32_t
{
  complex_float_union_t v1;
  float_union_t         v2;
};

// Two double\complex double values combined type
struct data_2_cf64_t
{
  complex_double_union_t v1;
  double_union_t         v2;
};

// Two complex float\float values combined type
struct data_2_fc32_t
{
  float_union_t         v1;
  complex_float_union_t v2;
};

// Two complex double\double values combined type
struct data_2_fc64_t
{
  double_union_t         v1;
  complex_double_union_t v2;
};

// Seven float values combined type
struct data_7_f32_t
{
  float_union_t v1;
  float_union_t v2;
  float_union_t v3;
  float_union_t v4;
  float_union_t v5;
  float_union_t v6;
  float_union_t v7;
};

// Seven double values combined type
struct data_7_f64_t
{
  double_union_t v1;
  double_union_t v2;
  double_union_t v3;
  double_union_t v4;
  double_union_t v5;
  double_union_t v6;
  double_union_t v7;
};

// Reference data table structure
// for 1 arg 1 res functions
struct data_2_t
{
    int i;
    std::vector<data_2_f32_t> data_f32;
    std::vector<data_2_f64_t> data_f64;
    std::vector<data_2_c32_t> data_c32;
    std::vector<data_2_c64_t> data_c64;

    void get_values(float & v1, float & v2)
    {
      v1 = data_f32[i].v1.f;
      v2 = data_f32[i].v2.f;
      i = (i + 1) % data_f32.size();
    }

    void get_values(double & v1, double & v2)
    {
      v1 = data_f64[i].v1.f;
      v2 = data_f64[i].v2.f;
      i = (i + 1) % data_f64.size();
    }

    void get_values(std::complex<float> & v1, std::complex<float> & v2 )
    {
      v1 = { data_c32[i].v1.f.real, data_c32[i].v1.f.imag};
      v2 = { data_c32[i].v2.f.real, data_c32[i].v2.f.imag};
      i = (i + 1) % data_c32.size();
    }

    void get_values(std::complex<double> & v1, std::complex<double> & v2 )
    {
      v1 = { data_c64[i].v1.f.real, data_c64[i].v1.f.imag};
      v2 = { data_c64[i].v2.f.real, data_c64[i].v2.f.imag};
      i = (i + 1) % data_c64.size();
    }
}; // struct data_2_t

// Reference data table structure
// for 2 arg 1 res or 1 arg 2 res functions
struct data_3_t
{
    int i;
    std::vector<data_3_f32_t> data_f32;
    std::vector<data_3_f64_t> data_f64;
    std::vector<data_3_c32_t> data_c32;
    std::vector<data_3_c64_t> data_c64;

    void get_values(float & v1, float & v2, float & v3)
    {
      v1 = data_f32[i].v1.f;
      v2 = data_f32[i].v2.f;
      v3 = data_f32[i].v3.f;
      i = (i + 1) % data_f32.size();
    }

    void get_values(double & v1, double & v2, double & v3)
    {
      v1 = data_f64[i].v1.f;
      v2 = data_f64[i].v2.f;
      v3 = data_f64[i].v3.f;
      i = (i + 1) % data_f64.size();
    }

    void get_values(std::complex<float> & v1, std::complex<float> & v2, std::complex<float> & v3)
    {
      v1 = { data_c32[i].v1.f.real, data_c32[i].v1.f.imag};
      v2 = { data_c32[i].v2.f.real, data_c32[i].v2.f.imag};
      v3 = { data_c32[i].v3.f.real, data_c32[i].v3.f.imag};
      i = (i + 1) % data_c32.size();
    }

    void get_values(std::complex<double> & v1, std::complex<double> & v2, std::complex<double> & v3)
    {
      v1 = { data_c64[i].v1.f.real, data_c64[i].v1.f.imag};
      v2 = { data_c64[i].v2.f.real, data_c64[i].v2.f.imag};
      v3 = { data_c64[i].v3.f.real, data_c64[i].v3.f.imag};
      i = (i + 1) % data_c64.size();
    }
}; // struct data_3_t

// Reference data table structure
// for 1 arg 1 res complex to real functions
struct data_2cf_t
{
    int i;
    std::vector<data_2_f32_t>  data_f32;
    std::vector<data_2_f64_t>  data_f64;
    std::vector<data_2_cf32_t> data_c32;
    std::vector<data_2_cf64_t> data_c64;

    void get_values(float & v1, float & v2)
    {
      v1 = data_f32[i].v1.f;
      v2 = data_f32[i].v2.f;
      i = (i + 1) % data_f32.size();
    }

    void get_values(double & v1, double & v2)
    {
      v1 = data_f64[i].v1.f;
      v2 = data_f64[i].v2.f;
      i = (i + 1) % data_f64.size();
    }

    void get_values(std::complex<float> & v1, float & v2 )
    {
      v1 = { data_c32[i].v1.f.real, data_c32[i].v1.f.imag};
      v2 = data_c32[i].v2.f;
      i = (i + 1) % data_c32.size();
    }

    void get_values(std::complex<double> & v1, double & v2 )
    {
      v1 = { data_c64[i].v1.f.real, data_c64[i].v1.f.imag};
      v2 = data_c64[i].v2.f;
      i = (i + 1) % data_c64.size();
    }
}; // struct data_2cf_t

// Reference data table structure
// for 1 arg 1 res real to complex functions
struct data_2fc_t
{
    int i;
    std::vector<data_2_fc32_t> data_c32;
    std::vector<data_2_fc64_t> data_c64;

    void get_values(float & v1, std::complex<float> & v2)
    {
      v1 = data_c32[i].v1.f;
      v2 = { data_c32[i].v2.f.real, data_c32[i].v2.f.imag};
      i = (i + 1) % data_c32.size();
    }

    void get_values(double & v1, std::complex<double> & v2 )
    {
      v1 = data_c64[i].v1.f;
      v2 = { data_c64[i].v2.f.real, data_c64[i].v2.f.imag};
      i = (i + 1) % data_c64.size();
    }
}; // struct data_2fc_t

// Reference data table structure
// for 6 arg 1 res functions
struct data_7_t
{
    int i;
    std::vector<data_7_f32_t> data_f32;
    std::vector<data_7_f64_t> data_f64;

    void get_values(float & v1, float & v2, float & v3, float & v4, float & v5, float & v6, float & v7)
    {
      v1 = data_f32[i].v1.f;
      v2 = data_f32[i].v2.f;
      v3 = data_f32[i].v3.f;
      v4 = data_f32[i].v4.f;
      v5 = data_f32[i].v5.f;
      v6 = data_f32[i].v6.f;
      v7 = data_f32[i].v7.f;
      i = (i + 1) % data_f32.size();
    }

    void get_values(double & v1, double & v2, double & v3, double & v4, double & v5, double & v6, double & v7)
    {
      v1 = data_f64[i].v1.f;
      v2 = data_f64[i].v2.f;
      v3 = data_f64[i].v3.f;
      v4 = data_f64[i].v4.f;
      v5 = data_f64[i].v5.f;
      v6 = data_f64[i].v6.f;
      v7 = data_f64[i].v7.f;
      i = (i + 1) % data_f64.size();
    }
}; // struct data_7_t

template <typename T> T * own_malloc_device( cl::sycl::queue & q, size_t s)
{
    T* ptr = static_cast<T *>(cl::sycl::malloc_device(s, q));

    if(ptr == nullptr)
    {
        throw std::runtime_error("Memory allocation error");
    }

    return ptr;
}



template <typename T> T* own_malloc( cl::sycl::queue & q, size_t s)
{
    if(q.get_device().is_gpu())
    {
        T* ptr = static_cast<T *>(cl::sycl::malloc_shared(s, q));

        if(ptr == nullptr)
        {
            throw std::runtime_error("Memory allocation error");
        }

        return ptr;
    }
    else
    {
        return static_cast<T *>(::malloc(s));
    }
}

template <typename T> void own_free_device( cl::sycl::queue & q, T* obj)
{
        cl::sycl::free(obj, q);
} 

template <typename T> void own_free( cl::sycl::queue & q, T* obj)
{
    if(q.get_device().is_gpu())
    {
        cl::sycl::free(obj, q);
    }
    else
    {
        ::free(obj);
    }
    return;
}

template <typename T> T * safe_copy_from_device(T * dev_begin, T * dev_end, T * dst, cl::sycl::queue & q) {
    if (dev_end < dev_begin) { throw std::runtime_error("Safe copy received end < begin"); }
    if (nullptr == dst) { throw std::runtime_error("Safe copy received dst == nullptr"); }

    try
    {
        size_t n = (dev_end - dev_begin);
        size_t s =  n * sizeof(T);
        T * host_ptr = static_cast<T *>(cl::sycl::malloc_host(s, q));
        if (nullptr == host_ptr) { throw std::runtime_error("USM memory allocation error"); }
        q.memcpy(host_ptr, dev_begin, s);
        q.wait_and_throw();
        std::copy(host_ptr, host_ptr + n, dst);
        cl::sycl::free(host_ptr, q);
    }
    catch (cl::sycl::exception & e)
    {
        std::cerr << "SYCL exception during memcpy from device\n"
            << e.what() << std::endl << "OpenCl status: " << e.get_cl_code() << std::endl;
        throw;
    }

    return dst;
}

template <typename T> T * safe_copy_to_device(T * begin, T * end, T * dev_dst, cl::sycl::queue & q) {
    if (end < begin) { throw std::runtime_error("Safe copy received end < begin"); }
    if (nullptr == dev_dst) { throw std::runtime_error("Safe copy received dev_dst == nullptr"); }

    try
    {
        size_t n = (end - begin);
        size_t s =  n * sizeof(T);
        T * host_ptr = static_cast<T *>(cl::sycl::malloc_host(s, q));
        if (nullptr == host_ptr) { throw std::runtime_error("USM memory allocation error"); }
        std::copy(begin, end, host_ptr);
        q.memcpy(dev_dst, host_ptr, s);
        q.wait_and_throw();
        cl::sycl::free(host_ptr, q);
    }
    catch (cl::sycl::exception & e)
    {
        std::cerr << "SYCL exception during memcpy from device\n"
            << e.what() << std::endl << "OpenCl status: " << e.get_cl_code() << std::endl;
        throw;
    }

    return dev_dst;
}

//!
//! @brief Error string output utility
//!
static std::string GetErrString(oneapi::mkl::vm::status st)
{
    std::string retstr = "";
    if(st == oneapi::mkl::vm::status::success)
    {
        retstr += "Success";
    }
    else
    {
        int numstr = 0;
        if(has_any(st, oneapi::mkl::vm::status::errdom))
        {
            retstr += "Error_Domain";
            numstr++;
        }
        if(has_any(st, oneapi::mkl::vm::status::sing))
        {
            if (numstr) retstr += "|";
            retstr += "Singularity";
            numstr++;
        }
        if(has_any(st, oneapi::mkl::vm::status::overflow))
        {
            if (numstr) retstr += "|";
            retstr += "Overflow";
            numstr++;
        }
        if(has_any(st, oneapi::mkl::vm::status::underflow))
        {
            if (numstr) retstr += "|";
            retstr += "Underflow";
            numstr++;
        }
        if(numstr == 0)
        {
            retstr += "Non-identified_status_code";
        }
    }

    retstr = "( " + std::to_string(static_cast<int>(st)) + " )" + retstr;

    return retstr;
}
