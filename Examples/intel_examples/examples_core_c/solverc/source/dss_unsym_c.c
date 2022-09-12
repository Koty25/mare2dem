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
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) DSS C example
*
********************************************************************************
*/
/*
**
** Example program to solve real unsymmetric system of equations.
** The example also demonstrates how to solve transposed system ATx=b.
**
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mkl_dss.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
/* Define the array and rhs vectors
*/
#define NROWS       5
#define NCOLS       5
#define NNONZEROS   9
#define NRHS        1
static const MKL_INT nRows = NROWS;
static const MKL_INT nCols = NCOLS;
static const MKL_INT nNonZeros = NNONZEROS;
static const MKL_INT nRhs = NRHS;
static MKL_INT rowIndex[NROWS + 1] = { 1, 3, 5, 7, 9, 10 };
static MKL_INT columns[NNONZEROS] = { 1, 2, 1, 2, 3, 4, 3, 4, 5 };
static double values[NNONZEROS] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
static _DOUBLE_PRECISION_t rhs[NCOLS * 2];
static _DOUBLE_PRECISION_t solValues[NROWS] = { 0, 1, 2, 3, 4 };

MKL_INT
main ()
{
  // Descriptor of main sparse matrix properties
  struct matrix_descr descrA;
  // Structure with sparse matrix stored in CSR format
  sparse_matrix_t       csrA;
  sparse_operation_t    transA;
  MKL_INT i, j;
  /* Allocate storage for the solver handle and the right-hand side. */
  _MKL_DSS_HANDLE_t handle;
  _INTEGER_t error;
  _CHARACTER_t statIn[] = "determinant", *uplo;
  _DOUBLE_PRECISION_t statOut[5], eps = 1e-6;
  MKL_INT opt = MKL_DSS_DEFAULTS, opt1;
  MKL_INT sym = MKL_DSS_NON_SYMMETRIC;
  MKL_INT type = MKL_DSS_INDEFINITE;
/* --------------------- */
/* Initialize the solver */
/* --------------------- */

  error = dss_create (handle, opt);
  if (error != MKL_DSS_SUCCESS)
    goto printError;
/* ------------------------------------------- */
/* Define the non-zero structure of the matrix */
/* ------------------------------------------- */
  error = dss_define_structure (handle, sym, rowIndex, nRows, nCols, columns, nNonZeros);
  if (error != MKL_DSS_SUCCESS)
    goto printError;
/* ------------------ */
/* Reorder the matrix */
/* ------------------ */
  error = dss_reorder (handle, opt, 0);
  if (error != MKL_DSS_SUCCESS)
    goto printError;
/* ------------------ */
/* Factor the matrix  */
/* ------------------ */
  error = dss_factor_real (handle, type, values);
  if (error != MKL_DSS_SUCCESS)
    goto printError;
/* ------------------------ */
/* Get the solution vector for Ax=b and ATx=b and check correctness */
/* ------------------------ */
  descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
  descrA.mode = SPARSE_FILL_MODE_UPPER;
  descrA.diag = SPARSE_DIAG_NON_UNIT;
  mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ONE, nRows, nCols, rowIndex, rowIndex+1, columns, values );
  for (i = 0; i < 3; i++)
    {
      if (i == 0)
        {
          transA = SPARSE_OPERATION_NON_TRANSPOSE;
          opt1 = MKL_DSS_DEFAULTS;
          printf ("\nSolving NON_TRANSPOSE system...\n");
        }
      else if (i == 1)
        {
          transA = SPARSE_OPERATION_TRANSPOSE;
          opt1 = MKL_DSS_TRANSPOSE_SOLVE;
          printf ("\nSolving TRANSPOSE system...\n");
        }
      else if (i == 2)
        {
          transA = SPARSE_OPERATION_CONJUGATE_TRANSPOSE;
          opt1 = MKL_DSS_CONJUGATE_SOLVE;
          printf ("\nSolving CONJUGATE_TRANSPOSE system...\n");
        }

// Compute rhs respectively to uplo to have solution solValue
      mkl_sparse_d_mv( transA, 1.0, csrA, descrA, solValues, 0.0, rhs);

// Nullify solution on entry (for sure)
      for (j = 0; j < nCols; j++)
        solValues[j] = 0.0;

// Apply trans or non-trans option, solve system
      opt |= opt1;
      error = dss_solve_real (handle, opt, rhs, nRhs, solValues);
      if (error != MKL_DSS_SUCCESS)
        goto printError;
      opt &= ~opt1;

// Check solution vector: should be {0,1,2,3,4}
      for (j = 0; j < nCols; j++)
        {
          if ((solValues[j] > j + eps) || (solValues[j] < j - eps))
            {
              printf ("Incorrect solution\n");
              error = 1000 + i;
              goto printError;
            }
        }
      printf ("Print solution array: ");
      for (j = 0; j < nCols; j++)
        printf (" %g", solValues[j]);

      printf ("\n");
    }
/* -------------------------- */
/* Deallocate solver storage  */
/* -------------------------- */
  mkl_sparse_destroy(csrA);
  error = dss_delete (handle, opt);
  if (error != MKL_DSS_SUCCESS)
    goto printError;
/* ---------------------- */
/* Print solution vector  */
/* ---------------------- */
  printf ("\nExample successfully PASSED!\n");
  exit (0);
printError:
  printf ("Solver returned error code %d\n", error);
  exit (1);
}
