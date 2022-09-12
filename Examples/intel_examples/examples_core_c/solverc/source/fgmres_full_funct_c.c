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
*  Content:
*  Intel(R) Math Kernel Library (Intel(R) MKL) RCI (P)FGMRES ((Preconditioned)
*  Flexible Generalized Minimal RESidual method) example
********************************************************************************/

/*---------------------------------------------------------------------------
*  Example program for solving non-symmetric indefinite system of equations
*  Fully advanced case: full functionality of RCI FGMRES solver is exploited
*---------------------------------------------------------------------------*/

#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
#define N 5
#define size 128

int main (void)
{

  /*---------------------------------------------------------------------------
  * Define arrays for the upper triangle of the coefficient matrix
  * Compressed sparse row storage is used for sparse representation
  *---------------------------------------------------------------------------*/
  MKL_INT ia[6] = { 1, 3, 6, 9, 12, 14 };
  MKL_INT ja[13] = 
  { 1,    3,
    1, 2,    4,
       2, 3,    5,
          3, 4, 5,
             4, 5
  };
  double A[13] = 
  { 1.0,      -1.0,
   -1.0, 1.0,      -1.0,
         1.0, -2.0,       1.0,
              -1.0, 2.0, -1.0,
                   -1.0, -3.0
  };
  /*---------------------------------------------------------------------------
  * Allocate storage for the ?par parameters and the solution/rhs/residual vectors
  *---------------------------------------------------------------------------*/
  MKL_INT ipar[size];
  double dpar[size], tmp[N * (2 * N + 1) + (N * (N + 9)) / 2 + 1] = {0.0};
  double expected_solution[N] = { -1.0, 1.0, 0.0, 1.0, -1.0 };
  double rhs[N] = {0.0}, b[N] = {0.0};
  double computed_solution[N] = {0.0};
  double residual[N] = {0.0};
  /*---------------------------------------------------------------------------
  * Some additional variables to use with the RCI (P)FGMRES solver
  *---------------------------------------------------------------------------*/
  MKL_INT itercount, expected_itercount = 4;
  MKL_INT RCI_request, i, ivar;
  double dvar;

  // Descriptor of main sparse matrix properties
  struct matrix_descr descrA;
  // Structure with sparse matrix stored in CSR format
  sparse_matrix_t       csrA;
  sparse_operation_t    transA = SPARSE_OPERATION_NON_TRANSPOSE;
  descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
  descrA.mode = SPARSE_FILL_MODE_UPPER;
  descrA.diag = SPARSE_DIAG_NON_UNIT;
  mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ONE, N, N, ia, ia+1, ja, A );
  printf ("--------------------------------------------------------\n");
  printf ("The FULLY ADVANCED example of usage of RCI FGMRES solver\n");
  printf ("   to solve a non-symmetric indefinite non-degenerate\n");
  printf ("          algebraic system of linear equations\n");
  printf ("--------------------------------------------------------\n\n");
  /*---------------------------------------------------------------------------
  * Initialize variables and the right hand side through matrix-vector product
  *---------------------------------------------------------------------------*/
  ivar = N;
  mkl_sparse_d_mv( transA, 1.0, csrA, descrA, expected_solution, 0.0, rhs);
  /*---------------------------------------------------------------------------
  * Save the right-hand side in vector b for future use
  *---------------------------------------------------------------------------*/
  i = 1;
  dcopy (&ivar, rhs, &i, b, &i);
  /*---------------------------------------------------------------------------
  * Initialize the initial guess
  *---------------------------------------------------------------------------*/
  for (i = 0; i < N; i++)
    {
      computed_solution[i] = 1.0;
    }
  /*---------------------------------------------------------------------------
  * Initialize the solver
  *---------------------------------------------------------------------------*/
  dfgmres_init (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
  if (RCI_request != 0)
    goto FAILED;
  /*---------------------------------------------------------------------------
  * Set the desired parameters:
  * do the restart after 2 iterations
  * LOGICAL parameters:
  * do not do the stopping test for the maximal number of iterations
  * do the Preconditioned iterations of FGMRES method
  * DOUBLE PRECISION parameters
  * set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
  *---------------------------------------------------------------------------*/
  ipar[14] = 2;
  ipar[7] = 0;
  ipar[10] = 1;
  dpar[0] = 1.0E-3;
  /*---------------------------------------------------------------------------
  * Check the correctness and consistency of the newly set parameters
  *---------------------------------------------------------------------------*/
  dfgmres_check (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar,
        tmp);
  if (RCI_request != 0 && RCI_request != -1001)
    goto FAILED;
  /*---------------------------------------------------------------------------
  * Print the info about the RCI FGMRES method
  *---------------------------------------------------------------------------*/
  printf ("Some info about the current run of RCI FGMRES method:\n\n");
  if (ipar[7])
    {
      printf ("As ipar[7]=%d, the automatic test for the maximal number of ", ipar[7]);
      printf ("iterations will be\nperformed\n");
    }
  else
    {
      printf ("As ipar[7]=%d, the automatic test for the maximal number of ", ipar[7]);
      printf ("iterations will be\nskipped\n");
    }
  printf ("+++\n");
  if (ipar[8])
    {
      printf ("As ipar[8]=%d, the automatic residual test will be performed\n", ipar[8]);
    }
  else
    {
      printf ("As ipar[8]=%d, the automatic residual test will be skipped\n", ipar[8]);
    }
  printf ("+++\n");
  if (ipar[9])
    {
      printf ("As ipar[9]=%d, the user-defined stopping test will be ", ipar[9]);
      printf ("requested via\nRCI_request=2\n");
    }
  else
    {
      printf ("As ipar[9]=%d, the user-defined stopping test will not be ", ipar[9]);
      printf ("requested, thus,\nRCI_request will not take the value 2\n");
    }
  printf ("+++\n");
  if (ipar[10])
    {
      printf ("As ipar[10]=%d, the Preconditioned FGMRES iterations will be ", ipar[10]);
      printf ("performed, thus,\nthe preconditioner action will be requested via ");
      printf ("RCI_request=3\n");
    }
  else
    {
      printf ("As ipar[10]=%d, the Preconditioned FGMRES iterations will not ", ipar[10]);
      printf ("be performed,\nthus, RCI_request will not take the value 3\n");
    }
  printf ("+++\n");
  if (ipar[11])
    {
      printf ("As ipar[11]=%d, the automatic test for the norm of the next ", ipar[11]);
      printf ("generated vector is\nnot equal to zero up to rounding and ");
      printf ("computational errors will be performed,\nthus, RCI_request will not ");
      printf ("take the value 4\n");
    }
  else
    {
      printf ("As ipar[11]=%d, the automatic test for the norm of the next ", ipar[11]);
      printf ("generated vector is\nnot equal to zero up to rounding and ");
      printf ("computational errors will be skipped,\nthus, the user-defined test ");
      printf ("will be requested via RCI_request=4\n");
    }
  printf ("+++\n\n");
  /*---------------------------------------------------------------------------
  * Compute the solution by RCI (P)FGMRES solver with preconditioning
  * Reverse Communication starts here
  *---------------------------------------------------------------------------*/
ONE:dfgmres (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
  /*---------------------------------------------------------------------------
  * If RCI_request=0, then the solution was found with the required precision
  *---------------------------------------------------------------------------*/
  if (RCI_request == 0)
    goto COMPLETE;
  /*---------------------------------------------------------------------------
  * If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
  * and put the result in vector tmp[ipar[22]-1]
  *---------------------------------------------------------------------------
  * NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
  * therefore, in C code it is required to subtract 1 from them to get C style
  * addresses
  *---------------------------------------------------------------------------*/
  if (RCI_request == 1)
    {
      mkl_sparse_d_mv( transA, 1.0, csrA, descrA, &tmp[ipar[21] - 1], 0.0, &tmp[ipar[22] - 1]);
      goto ONE;
    }
  /*---------------------------------------------------------------------------
  * If RCI_request=2, then do the user-defined stopping test
  * The residual stopping test for the computed solution is performed here
  *---------------------------------------------------------------------------
  * NOTE: from this point vector b[N] is no longer containing the right-hand
  * side of the problem! It contains the current FGMRES approximation to the
  * solution. If you need to keep the right-hand side, save it in some other
  * vector before the call to dfgmres routine. Here we saved it in vector
  * rhs[N]. The vector b is used instead of rhs to preserve the
  * original right-hand side of the problem and guarantee the proper
  * restart of FGMRES method. Vector b will be altered when computing the
  * residual stopping criterion!
  *---------------------------------------------------------------------------*/
  if (RCI_request == 2)
    {
      /* Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]
         --------------------------------------------------------------------------------
         WARNING: beware that the call to dfgmres_get routine with ipar[12]=0 at this
         stage may destroy the convergence of the FGMRES method, therefore, only
         advanced users should exploit this option with care */
      ipar[12] = 1;
      /* Get the current FGMRES solution in the vector b[N] */
      dfgmres_get (&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
      /* Compute the current true residual via Intel MKL (Sparse) BLAS routines */
      mkl_sparse_d_mv( transA, 1.0, csrA, descrA, b, 0.0, residual);
      dvar = -1.0E0;
      i = 1;
      daxpy (&ivar, &dvar, rhs, &i, residual, &i);
      dvar = dnrm2 (&ivar, residual, &i);
      if (dvar < 1.0E-3)
        goto COMPLETE;
      else
        goto ONE;
    }
  /*---------------------------------------------------------------------------
  * If RCI_request=3, then apply the preconditioner on the vector
  * tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]
  *---------------------------------------------------------------------------
  * NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
  * therefore, in C code it is required to subtract 1 from them to get C style
  * addresses
  *---------------------------------------------------------------------------*/
  if (RCI_request == 3)
    {
      if (ipar[3] == 3)
        {
          tmp[ipar[22] - 1 + 0] = -2.0;
          tmp[ipar[22] - 1 + 1] = 0.08519601586107672;
          tmp[ipar[22] - 1 + 2] = -1.1590871369607090;
          tmp[ipar[22] - 1 + 3] = -0.65791939687456980;
          tmp[ipar[22] - 1 + 4] = 0.75660051476696133;
        }
      else
        {
          if (ipar[3] == 4)
            {
              tmp[ipar[22] - 1 + 0] = 0.0;
              tmp[ipar[22] - 1 + 1] = 0.0;
              tmp[ipar[22] - 1 + 2] = 0.0;
              tmp[ipar[22] - 1 + 3] = 1.0;
              tmp[ipar[22] - 1 + 4] = -1.0;
            }
          else
            {
              for (i = 0; i < N; i++)
                {
                  tmp[ipar[22] - 1 + i] = (double) (i) * tmp[ipar[21] - 1 + i];
                }
            }
        }
      goto ONE;
    }
  /*---------------------------------------------------------------------------
  * If RCI_request=4, then check if the norm of the next generated vector is
  * not zero up to rounding and computational errors. The norm is contained
  * in dpar[6] parameter
  *---------------------------------------------------------------------------*/
  if (RCI_request == 4)
    {
      if (dpar[6] < 1.0E-12)
        goto COMPLETE;
      else
        goto ONE;
    }
  /*---------------------------------------------------------------------------
  * If RCI_request=anything else, then dfgmres subroutine failed
  * to compute the solution vector: computed_solution[N]
  *---------------------------------------------------------------------------*/
  else
    {
      goto FAILED;
    }
  /*---------------------------------------------------------------------------
  * Reverse Communication ends here
  * Get the current iteration number and the FGMRES solution (DO NOT FORGET to
  * call dfgmres_get routine as computed_solution is still containing
  * the initial guess!). Request to dfgmres_get to put the solution
  * into vector computed_solution[N] via ipar[12]
  *---------------------------------------------------------------------------*/
COMPLETE:ipar[12] = 0;
  dfgmres_get (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
  /*---------------------------------------------------------------------------
  * Print solution vector: computed_solution[N] and the number of iterations: itercount
  *---------------------------------------------------------------------------*/
  printf (" The system has been solved \n");
  printf ("\n The following solution has been obtained: \n");
  for (i = 0; i < N; i++)
    {
      printf ("computed_solution[%d]=", i);
      printf ("%e\n", computed_solution[i]);
    }
  printf ("\n The expected solution is: \n");
  for (i = 0; i < N; i++)
    {
      printf ("expected_solution[%d]=", i);
      printf ("%e\n", expected_solution[i]);
      expected_solution[i] -= computed_solution[i];
    }
  printf ("\n Number of iterations: %d\n", itercount);
  i = 1;
  dvar = dnrm2 (&ivar, expected_solution, &i);

  /*-------------------------------------------------------------------------*/
  /* Release internal Intel MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
  MKL_Free_Buffers ();

  if (itercount == expected_itercount && dvar < 1.0e-14)
    {
      printf ("\nThis example has successfully PASSED through all steps of ");
      printf ("computation!\n");
      return 0;
    }
  else
    {
      printf ("\nThis example may have FAILED as either the number of iterations differs\n");
      printf ("from the expected number of iterations %d, ", expected_itercount);
      printf ("or the computed solution\n");
      printf ("differs much from the expected solution (Euclidean norm is %e), or both.\n", dvar);
      return 1;
    }
  /*-------------------------------------------------------------------------*/
  /* Release internal Intel MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
FAILED:printf ("\nThis example FAILED as the solver has returned the ERROR code %d", RCI_request);
  mkl_sparse_destroy(csrA);
  MKL_Free_Buffers ();
  return 1;
}
