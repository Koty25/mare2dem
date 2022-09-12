/*******************************************************************************
* Copyright 2001-2020 Intel Corporation.
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
!    routines for relative error calculation
!******************************************************************************/

#define FABS(x)     ((x)<0.0?(-(x)):(x))

float srelerr(float a, float b)
{
  return FABS((a-b)/a);
}

double drelerr(double a, double b)
{
  return FABS((a-b)/a);
}

float crelerr(MKL_Complex8 a, MKL_Complex8 b)
{
  float re,im;

  re=FABS((a.real-b.real)/a.real);
  im=FABS((a.imag-b.imag)/a.imag);

  if(re>im) return re;
  return im;
}

double zrelerr(MKL_Complex16 a, MKL_Complex16 b)
{
  double re,im;

  re=FABS((a.real-b.real)/a.real);
  im=FABS((a.imag-b.imag)/a.imag);

  if(re>im) return re;
  return im;
}
