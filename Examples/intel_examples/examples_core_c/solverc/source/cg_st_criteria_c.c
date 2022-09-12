/*******************************************************************************
* Copyright 2005-2020 Intel Corporation.
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
!
!  Content: Intel(R) Math Kernel Library (Intel(R) MKL) RCI (P)CG C example
!
!*******************************************************************************/

#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

/*---------------------------------------------------------------------------   */
/*  Example program for solving symmetric positive definite system of equations.*/
/*  Medium case: make use of the user-defined stopping test.                    */
/*---------------------------------------------------------------------------   */
int main (void)
{
  /*---------------------------------------------------------------------------   */
  /* Define arrays for the upper triangle of the coefficient matrix and rhs vector */
  /* Compressed sparse row storage is used for sparse representation              */
  /*---------------------------------------------------------------------------   */
  const MKL_INT n = 8;
  MKL_INT rci_request, itercount, expected_itercount = 8, i;
  double rhs[8] = {0.0};
  /* Fill all arrays containing matrix data. */
  MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
  MKL_INT ja[18] = 
  { 1,    3,       6, 7,
       2, 3,    5,
          3,             8,
             4,       7,
                5, 6, 7,
                   6,    8,
                      7,
                         8
  };
  double a[18] = 
  { 7.E0,       1.E0,             2.E0, 7.E0,
         -4.E0, 8.E0,       2.E0,
                1.E0,                         5.E0,
                      7.E0,             9.E0,
                            5.E0, 1.E0, 5.E0,
                                 -1.E0,       5.E0,
                                       11.E0,
                                              5.E0
  };

  /*---------------------------------------------------------------------------*/
  /* Allocate storage for the solver ?par and temporary storage tmp            */
  /*---------------------------------------------------------------------------*/
  MKL_INT length = 128;
  MKL_INT ipar[128];
  double dpar[128], tmp[4 * 8];
  /*---------------------------------------------------------------------------*/
  /* Some additional variables to use with the RCI (P)CG solver                */
  /*---------------------------------------------------------------------------*/
  double solution[8] = {0.0};
  double expected_sol[8] = { 1.E0, 0.E0, 1.E0, 0.E0, 1.E0, 0.E0, 1.E0, 0.E0 };
  double euclidean_norm=0.0, temp[8] = {0.0};
  double eone = -1.E0;
  MKL_INT ione = 1;

  struct matrix_descr descrA;
  // Structure with sparse matrix stored in CSR format
  sparse_matrix_t       csrA;
  sparse_operation_t    transA = SPARSE_OPERATION_NON_TRANSPOSE;
  descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
  descrA.mode = SPARSE_FILL_MODE_UPPER;
  descrA.diag = SPARSE_DIAG_NON_UNIT;
  mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ONE, n, n, ia, ia+1, ja, a );
  /*---------------------------------------------------------------------------*/
  /* Initialize the right hand side through matrix-vector product              */
  /*---------------------------------------------------------------------------*/
  mkl_sparse_d_mv( transA, 1.0, csrA, descrA, expected_sol, 0.0, rhs);
  /*---------------------------------------------------------------------------*/
  /* Initialize the initial guess                                              */
  /*---------------------------------------------------------------------------*/
  for (i = 0; i < n; i++)
    {
      solution[i] = 0.E0;
    }
  /*---------------------------------------------------------------------------*/
  /* Initialize the solver                                                     */
  /*---------------------------------------------------------------------------*/
  dcg_init (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
  if (rci_request != 0)
    goto failure;
  /*---------------------------------------------------------------------------*/
  /* Set the desired parameters:                                               */
  /* LOGICAL parameters:                                                       */
  /* -                                                                         */
  /* INTEGER parameters:                                                       */
  /* set the maximal number of iterations to 100                               */
  /* DOUBLE parameters                                                         */
  /* -                                                                         */
  /*---------------------------------------------------------------------------*/
  ipar[4] = 100;
  /*---------------------------------------------------------------------------*/
  /* Check the correctness and consistency of the newly set parameters         */
  /*---------------------------------------------------------------------------*/
  dcg_check (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
  if (rci_request != 0 && rci_request != -1001)
    goto failure;
  /*---------------------------------------------------------------------------*/
  /* Compute the solution by RCI (P)CG solver                                  */
  /* Reverse Communications starts here                                        */
  /*---------------------------------------------------------------------------*/
rci:dcg (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
  /*---------------------------------------------------------------------------*/
  /* If rci_request=0, then the solution was found according to the requested  */
  /* stopping tests. In this case, this means that it was found after 100      */
  /* iterations.                                                               */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 0)
    goto getsln;
  /*---------------------------------------------------------------------------*/
  /* If rci_request=1, then compute the vector A*TMP[0]                        */
  /* and put the result in vector TMP[n]                                       */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 1)
    {
      mkl_sparse_d_mv( transA, 1.0, csrA, descrA, tmp, 0.0, &tmp[n]);
      goto rci;
    }
  /*---------------------------------------------------------------------------*/
  /* If rci_request=2, then do the user-defined stopping test: compute the     */
  /* Euclidean norm of the actual residual using Intel MKL routines and check if     */
  /* it is less than 1.E-8                                                     */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 2)
    {
	  mkl_sparse_d_mv( transA, 1.0, csrA, descrA, solution, 0.0, temp);
      daxpy (&n, &eone, rhs, &ione, temp, &ione);
      euclidean_norm = dnrm2 (&n, temp, &ione);
      /*---------------------------------------------------------------------------*/
      /* The solution has not been found yet according to the user-defined stopping */
      /* test. Continue RCI (P)CG iterations.                                      */
      /*---------------------------------------------------------------------------*/
      if (euclidean_norm > 1.E-8)
        goto rci;
      /*---------------------------------------------------------------------------*/
      /* The solution has been found according to the user-defined stopping test   */
      /*---------------------------------------------------------------------------*/
      else
        goto getsln;
    }
  /*---------------------------------------------------------------------------*/
  /* If rci_request=anything else, then dcg subroutine failed                  */
  /* to compute the solution vector: solution[n]                               */
  /*---------------------------------------------------------------------------*/
  goto failure;
  /*---------------------------------------------------------------------------*/
  /* Reverse Communication ends here                                           */
  /* Get the current iteration number into itercount                           */
  /*---------------------------------------------------------------------------*/
getsln:dcg_get (&n, solution, rhs, &rci_request, ipar, dpar, tmp, &itercount);
  /*---------------------------------------------------------------------------*/
  /* Print solution vector: solution[n] and number of iterations: itercount    */
  /*---------------------------------------------------------------------------*/
  printf ("The system has been solved\n");
  printf ("The following solution obtained\n");
  for (i = 0; i < n / 2; i++)
    printf ("%6.3f  ", solution[i]);
  printf ("\n");
  for (i = n / 2; i < n; i++)
    printf ("%6.3f  ", solution[i]);
  printf ("\nExpected solution is\n");
  for (i = 0; i < n / 2; i++)
    {
      printf ("%6.3f  ", expected_sol[i]);
      expected_sol[i] -= solution[i];
    }
  printf ("\n");
  for (i = n / 2; i < n; i++)
    {
      printf ("%6.3f  ", expected_sol[i]);
      expected_sol[i] -= solution[i];
    }
  printf ("\nNumber of iterations: %d\n", itercount);
  i = 1;
  euclidean_norm = dnrm2 (&n, expected_sol, &i);

  /*-------------------------------------------------------------------------*/
  /* Release internal Intel MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
  MKL_Free_Buffers ();

  if (itercount == expected_itercount && euclidean_norm < 1.0e-12)
    {
      printf ("This example has successfully PASSED through all steps of computation!\n");
      return 0;
    }
  else
    {
      printf ("This example may have FAILED as either the number of iterations differs\n");
      printf ("from the expected number of iterations %d, or the ", expected_itercount);
      printf ("computed solution\ndiffers much from the expected solution ");
      printf ("(Euclidean norm is %e), or both.\n", euclidean_norm);
      return 1;
    }

  /*-------------------------------------------------------------------------*/
  /* Release internal Intel MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
failure:printf ("This example FAILED as the solver has returned the ERROR code %d", rci_request);
  MKL_Free_Buffers ();
  return 1;
}
