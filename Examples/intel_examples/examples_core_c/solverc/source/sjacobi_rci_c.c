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
*   Content : DJACOBI RCI example
*
******************************************************************************** 
*/
/* 
  The program computes the Jacobi matrix of the function on the basis of RCI
  using the central difference.
*/

#include <stdlib.h>
#include <stdio.h>
#include "mkl_service.h"

#include "mkl_rci.h"
#include "mkl_types.h"

int main ()
{
  /* user’s objective function */
  extern void extended_powell (MKL_INT *, MKL_INT *, float *, float *);
  /* n - number of function variables
     m - dimension of function value */
  MKL_INT n = 4, m = 4;
  /* jacobi matrix */
  float *a;
  /* solution vector. contains values x for f(x)
     temporary arrays f1 & f2 that contain f1 = f(x+eps) | f2 = f(x-eps) */
  float *x, *f1, *f2;
  /* precisions for jacobi_matrix calculation */
  float eps = 0.0001;
  /* jacobi-matrix solver handle */
  _JACOBIMATRIX_HANDLE_t handle;
  /* controls of rci cycle */
  MKL_INT successful, rci_request, i, error;

  if ((a = (float *) mkl_malloc(sizeof (float) * n * m, 64)) == NULL)
    {
      printf ("\n#fail: jacobi matrix allocation failed\n");
    }
  if ((x = (float *) mkl_malloc(sizeof (float) * n, 64)) == NULL)
    {
      printf ("\n#fail: vector X allocation failed\n");
      mkl_free (a);
      return 1;
    }
  if ((f1 = (float *) mkl_malloc(sizeof (float) * n, 64)) == NULL)
    {
      printf ("\n#fail: vector F1 allocation failed\n");
      mkl_free (x);
      mkl_free (a);
      return 1;
    }
  if ((f2 = (float *) mkl_malloc(sizeof (float) * n, 64)) == NULL)
    {
      printf ("\n#fail: vector F2 allocation failed\n");
      mkl_free (f1);
      mkl_free (x);
      mkl_free (a);
      return 1;
    }
  error = 0;
  /* set the x values */
  for (i = 0; i < n; i++)
    x[i] = 10.0;
  /* initialize solver (allocate memory, set initial values) */
  if (sjacobi_init (&handle, &n, &m, x, a, &eps) != TR_SUCCESS)
    {
      /* if function does not complete successfully then print error message */
      printf ("\n#fail: error in djacobi_init\n");
      fflush (0);
      error = 1;
      goto end;
    }
  /* set initial rci cycle variables */
  rci_request = 0;
  successful = 0;
  /* rci cycle */
  while (successful == 0)
    {
      /* call solver */
      if (sjacobi_solve (&handle, f1, f2, &rci_request) != TR_SUCCESS)
        {
          /* if function does not complete successfully then print error message */
          printf ("\n#fail: error in djacobi_solve\n");
          fflush (0);
          error = 1;
          goto end;
        }
      if (rci_request == 1)
        {
          /* calculate function value f1 = f(x+eps) */
          extended_powell (&m, &n, x, f1);
        }
      else if (rci_request == 2)
        {
          /* calculate function value f1 = f(x-eps) */
          extended_powell (&m, &n, x, f2);
        }
      else if (rci_request == 0)
        /* exit rci cycle */
        successful = 1;
    } /* rci cycle */
  /* free handle memory */
  if (sjacobi_delete (&handle) != TR_SUCCESS)
    {
      /* if function does not complete successfully then print error message */
      printf ("\n#fail: error in djacobi_delete\n");
      fflush (0);
      error = 1;
      goto end;
    }
end:
  MKL_Free_Buffers ();
  mkl_free (f2);
  mkl_free (f1);
  mkl_free (x);
  mkl_free (a);

  if (error != 0) return 1;

  printf ("#pass\n");
  fflush (0);

  return 0;
}

/* 
    routine for extended Powell function calculation
    m in : dimension of function value
    n in : number of function variables
    x in : vector for function calculation
    f out: function value f(x)
*/
void extended_powell (MKL_INT * m, MKL_INT * n, float *x, float *f)
{
  MKL_INT i;
  for (i = 0; i < (*n) / 4; i++)
    {
      f[4 * i] = x[4 * i] + 10.0 * x[4 * i + 1];
      f[4 * i + 1] = 2.2360679774 * (x[4 * i + 2] - x[4 * i + 3]);
      f[4 * i + 2] = (x[4 * i + 1] - 2.0 * x[4 * i + 2]) * 
                     (x[4 * i + 1] - 2.0 * x[4 * i + 2]);
      f[4 * i + 3] = 3.1622776601 * (x[4 * i] - x[4 * i + 3]) * 
                                    (x[4 * i] - x[4 * i + 3]);
    }
  return;
}
