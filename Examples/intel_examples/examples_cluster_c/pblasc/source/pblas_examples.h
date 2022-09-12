/*******************************************************************************
* Copyright 2010-2020 Intel Corporation.
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
!      Intel(R) Math Kernel Library PBLAS C example's definitions file
!
!******************************************************************************/

#ifndef _pblas_examples_h
#define _pblas_examples_h

#ifdef _WIN_

/* Definitions for proper work of examples on Windows */
#define blacs_pinfo_ BLACS_PINFO
#define blacs_get_ BLACS_GET
#define blacs_gridinit_ BLACS_GRIDINIT
#define blacs_gridinfo_ BLACS_GRIDINFO
#define blacs_barrier_ BLACS_BARRIER
#define blacs_gridexit_ BLACS_GRIDEXIT
#define blacs_exit_ BLACS_EXIT
#define igebs2d_ IGEBS2D
#define igebr2d_ IGEBR2D
#define sgebs2d_ SGEBS2D
#define sgebr2d_ SGEBR2D
#define dgebs2d_ DGEBS2D
#define dgebr2d_ DGEBR2D
#define sgesd2d_ SGESD2D
#define sgerv2d_ SGERV2D
#define dgesd2d_ DGESD2D
#define dgerv2d_ DGERV2D
#define numroc_ NUMROC
#define descinit_ DESCINIT
#define psnrm2_ PSNRM2
#define pdnrm2_ PDNRM2
#define psscal_ PSSCAL
#define pdscal_ PDSCAL
#define psdot_ PSDOT
#define pddot_ PDDOT
#define pslamch_ PSLAMCH
#define pdlamch_ PDLAMCH
#define indxg2l_ INDXG2L
#define pscopy_ PSCOPY
#define pdcopy_ PDCOPY
#define pstrsv_ PSTRSV
#define pdtrsv_ PDTRSV
#define pstrmv_ PSTRMV
#define pdtrmv_ PDTRMV
#define pslange_ PSLANGE
#define pdlange_ PDLANGE
#define psgemm_ PSGEMM
#define pdgemm_ PDGEMM
#define psgeadd_ PSGEADD
#define pdgeadd_ PDGEADD

#endif

/* Pi-number */
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/* Definition of MIN and MAX functions */
#define MAX(a,b)((a)<(b)?(b):(a))
#define MIN(a,b)((a)>(b)?(b):(a))

/* Definition of matrix descriptor */
typedef MKL_INT MDESC[ 9 ];

#endif
