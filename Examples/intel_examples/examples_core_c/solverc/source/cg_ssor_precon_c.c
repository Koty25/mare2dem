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

/*---------------------------------------------------------------------------
* Example program for solving symmetric positive definite system of
* equations A*x =f where A is a symmetric tridiagonal matrix of order n whose
* nonzero elements are defined as follows
*
*         a(i, i) = 1,     i=1, n;
*         a(i, i+1)= 1/2,    i=1, n-1;
*         a(i+1, i) = 1/2,   i=1, n-1;
*
* The system is solved with the help of the conjugate gradient method where
* the Symmetric Successive Over Relaxation (SSOR) stationary iterative solver
* with the given number of iterations is used as the preconditioner. The relaxation
* parameter is set to 1/2, the number of iterations for the SSOR preconditioner
* is equal to 20. For simplicity, we don't check the convergence while solving
* equations with the help of SSOR. Let us recall that the scheme of SSOR is
* the following:
*
*         (D+w*U^{T}) x_{k+1} = (-w*U+(1-w)*D)*x_{k}+w*f
*
* where w is the relaxation parameter, f is the right hand side, D is the diagonal
* of A, U its strict upper part so that
*
*                  A   =  U^{T} + D + U
*
* The compressed sparse row format is used for storing nonzeros of A and
* since D is the identity matrix, we don't need to store diagonal elements. So
* the code given below only uses the sparse representation of U. The Intel
* Intel MKL Sparse BLAS is designed so that it allows the user to perform large
* variety of operations with one sparse representation setting appropriate
* values of the descriptor array (see Intel MKL User Guide for the further details).
*
* Full case: full functionality of RCI (P)CG is used.
*---------------------------------------------------------------------------*/
int main (void)
{
  /*---------------------------------------------------------------------------*/
  /* Define arrays for the upper triangle of the coefficient matrix and        */
  /* preconditioner as well as an array for rhs vector                         */
  /* Compressed sparse row storage is used for sparse representation           */
  /*---------------------------------------------------------------------------*/
  MKL_INT n = 100, rci_request, itercount, lexpected_itercount = 15,
    uexpected_itercount = 19, i;
  double rhs[100] = {0.0};
  MKL_INT ia[100 + 1];
  MKL_INT ja[100 - 1];
  double a[100 - 1], a1[100 - 1];
  /*---------------------------------------------------------------------------*/
  /* Allocate storage for the solver ?par and temporary storage tmp            */
  /*---------------------------------------------------------------------------*/
  MKL_INT length = 128;
  MKL_INT ipar[128];
  double dpar[128], tmp[4 * 100] = {0.0};
  /*---------------------------------------------------------------------------*/
  /* Some additional variables to use with the RCI (P)CG solver                */
  /* OMEGA is the relaxation parameter, NITER_SSOR is the maximum number of    */
  /* iterations for the SSOR preconditioner                                    */
  /*---------------------------------------------------------------------------*/
  double solution[100];
  double expected_sol[100];
  double omega = 0.5E0, one = 1.E0, zero = 0.E0, om = 1.E0 - omega;
  double euclidean_norm, temp[100] = {0.0};
  MKL_INT niter_ssor = 20;
  double eone = -1.E0;
  MKL_INT ione = 1;
  struct matrix_descr descrA;
  // Structure with sparse matrix stored in CSR format
  sparse_matrix_t       csrA, csrA1;
  sparse_operation_t    transA = SPARSE_OPERATION_NON_TRANSPOSE;
  descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
  descrA.mode = SPARSE_FILL_MODE_UPPER;
  descrA.diag = SPARSE_DIAG_UNIT;
  /*---------------------------------------------------------------------------*/
  /* Initialize the coefficient matrix and expected solution                   */
  /*---------------------------------------------------------------------------*/
  for (i = 0; i < n; i++)
    expected_sol[i] = 1.E0;

  for (i = 0; i < n - 1; i++)
    {
      ja[i] = i + 2;
      ia[i] = i + 1;
      a[i] = 0.5E0;
      a1[i] = omega * a[i];
    }
  ia[n - 1] = n;
  ia[n] = ia[n - 1];

  mkl_sparse_d_create_csr ( &csrA , SPARSE_INDEX_BASE_ONE, n, n, ia, ia+1, ja, a );
  mkl_sparse_d_create_csr ( &csrA1, SPARSE_INDEX_BASE_ONE, n, n, ia, ia+1, ja, a1 );

  /*------------------------------------------------------------------------------*/
  /* Initialize vectors rhs, temp, and tmp[n:2*n-1] with zeros as mkl_sparse_d_mv */
  /* routine does not set NAN to zero. Thus, if any of the values in the          */
  /* vectors above accidentally happens to be NAN, the example will fail          */
  /* to complete.                                                                 */
  /* Initialize the right hand side through matrix-vector product                 */
  /*------------------------------------------------------------------------------*/
  for (i = 0; i < n; i++)
    {
      rhs[i] = zero;
      temp[i] = zero;
      tmp[n + i] = zero;
    }
  mkl_sparse_d_mv( transA, 1.0, csrA, descrA, expected_sol, 0.0, rhs);
  /*---------------------------------------------------------------------------*/
  /* Initialize the initial guess                                              */
  /*---------------------------------------------------------------------------*/
  for (i = 0; i < n; i++)
    solution[i] = zero;
  /*---------------------------------------------------------------------------*/
  /* Initialize the solver                                                     */
  /*---------------------------------------------------------------------------*/
  dcg_init (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
  if (rci_request != 0)
    goto failure;
  /*---------------------------------------------------------------------------*/
  /* Set the desired parameters:                                               */
  /* INTEGER parameters:                                                       */
  /* set the maximal number of iterations to 100                               */
  /* LOGICAL parameters:                                                       */
  /* run the Preconditioned version of RCI (P)CG with preconditioner C_inverse */
  /* DOUBLE parameters                                                         */
  /* -                                                                         */
  /*---------------------------------------------------------------------------*/
  ipar[4] = 100;
  ipar[10] = 1;
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
  /* If rci_request=1, then compute the vector A*tmp[0]                        */
  /* and put the result in vector tmp[n]                                       */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 1)
    {
      descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
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
      descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	  mkl_sparse_d_mv( transA, 1.0, csrA, descrA, solution, 0.0, temp);
      daxpy (&n, &eone, rhs, &ione, temp, &ione);
      euclidean_norm = dnrm2 (&n, temp, &ione);
      /*---------------------------------------------------------------------------*/
      /* The solution has not been found yet according to the user-defined stopping */
      /* test. Continue RCI (P)CG iterations.                                      */
      /*---------------------------------------------------------------------------*/
      if (euclidean_norm > 1.E-6)
        goto rci;
      /*---------------------------------------------------------------------------*/
      /* The solution has been found according to the user-defined stopping test   */
      /*---------------------------------------------------------------------------*/
      else
        goto getsln;
    }
  /*---------------------------------------------------------------------------*/
  /* If rci_request=3, then  apply the simplest SSOR preconditioning           */
  /* on vector tmp[2*n] and put the result in vector tmp[3*n]                  */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 3)
    {
      dcopy (&n, &tmp[2 * n], &ione, &tmp[3 * n], &ione);
      descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
      for (i = 1; i <= niter_ssor; i++)
        {
          dcopy (&n, &tmp[2 * n], &ione, temp, &ione);
          descrA.diag = SPARSE_DIAG_NON_UNIT;
	      mkl_sparse_d_mv( transA, -1.0, csrA1, descrA, &tmp[3 * n], omega, temp);
          daxpy (&n, &om, &tmp[3 * n], &ione, temp, &ione);
          descrA.diag = SPARSE_DIAG_UNIT;
	      mkl_sparse_d_trsv( SPARSE_OPERATION_TRANSPOSE, 1.0, csrA1, descrA, temp, &tmp[3 * n]);
        }
      goto rci;
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
  for (i = 0; i < n / 4; i++)
    {
      printf ("%6.3f  %6.3f  %6.3f  %6.3f", solution[4 * i],
            solution[4 * i + 1], solution[4 * i + 2], solution[4 * i + 3]);
      printf ("\n");
    }
  printf ("\nExpected solution is\n");
  for (i = 0; i < n / 4; i++)
    {
      printf ("%6.3f  %6.3f  %6.3f  %6.3f", expected_sol[4 * i],
            expected_sol[4 * i + 1], expected_sol[4 * i + 2],
            expected_sol[4 * i + 3]);
      expected_sol[4 * i] -= solution[4 * i];
      printf ("\n");
    }

  printf ("\nNumber of iterations: %d\n", itercount);
  i = 4;
  n /= 4;
  euclidean_norm = dnrm2 (&n, expected_sol, &i);

  /*-------------------------------------------------------------------------*/
  /* Release internal Intel MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
  MKL_Free_Buffers ();

  if (lexpected_itercount <= itercount <= uexpected_itercount && euclidean_norm < 1.0e-4)
    {
      printf ("This example has successfully PASSED through all steps of computation!");
      printf ("\n");
      return 0;
    }
  else
    {
      printf ("This example may have FAILED as either the number of iterations differs\n");
      printf ("from the expected number of iterations %d-", lexpected_itercount);
      printf ("-%d, or the computed solution\ndiffers much from ", uexpected_itercount);
      printf ("the expected solution (Euclidean norm is %e), or both.\n", euclidean_norm);
      return 1;
    }
  /*-------------------------------------------------------------------------*/
  /* Release internal Intel MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
failure:printf ("This example FAILED as the solver has returned the ERROR code %d", rci_request);
  mkl_sparse_destroy(csrA);
  mkl_sparse_destroy(csrA1);
  MKL_Free_Buffers ();
  return 1;
}
