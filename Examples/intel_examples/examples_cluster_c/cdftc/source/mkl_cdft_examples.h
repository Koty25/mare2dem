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
! Content:
!       Intel(R) Math Kernel Library (Intel(R) MKL) Cluster DFT example's
!       definitions file (C-interface)
!
!****************************************************************************/

/*
**  Print level definition
*/
#define ADVANCED_DATA_PRINT 1
#define ACCURACY_PRINT 1
#define LEGEND_PRINT 1

/*
**  Accuracy definitions
*/
#define SINGLE_EPS 1.0E-6
#define DOUBLE_EPS 1.0E-12

/*
**  Intel MKL test _Complex type definition
*/
typedef struct {
    float re;
    float im;
} mkl_float_complex;

typedef struct {
    double re;
    double im;
} mkl_double_complex;

/*
**  Example support function's interfaces
*/
int read_data_file_2d(int, char*[], long*, long*);

void dfti_example_status_print(int, long);

void print_data_2d_z(void*, long, long, long);

void init_data_2d_z(void*, long, long);

double check_result_z(void*, void*, long);

void print_data_2d_c(void*, long, long, long);

void init_data_2d_c(void*, long, long);

float check_result_c(void*, void*, long);

int MKL_CDFT_Data(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                  MKL_LONG Lengths[],void *global,long nx,long start_x,void *local,
                  int Flag);

int MKL_CDFT_ScatterData(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                         MKL_LONG Lengths[],void *global_in,long nx,long start_x,
                         void *local_in);

int MKL_CDFT_GatherData(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                        MKL_LONG Lengths[],void *global_out,long nx,long start_x,
                        void *local_out);
