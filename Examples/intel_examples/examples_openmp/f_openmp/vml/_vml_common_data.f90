!===============================================================================
! Copyright 2020 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content: common data for examples
!*******************************************************************************

    character(len=1), parameter :: TAB = achar(9)
  
    integer (kind=4), parameter :: VML_ARG1_RES1    = 11    
    integer (kind=4), parameter :: VML_ARG2_RES1    = 21    
    integer (kind=4), parameter :: VML_ARG1_RES2    = 12    
    integer (kind=4), parameter :: VML_ARG6_RES1    = 61    
    integer (kind=4), parameter :: VML_ARG1R_RES1C  = 1011  
    integer (kind=4), parameter :: VML_ARG1C_RES1R  = 1110  
    
    integer (kind=4), parameter :: FP_IEEE_NAN        = 1
    integer (kind=4), parameter :: FP_IEEE_INF        = 2
    integer (kind=4), parameter :: FP_IEEE_ZERO       = 3
    integer (kind=4), parameter :: FP_IEEE_SUBNORMAL  = 4
    integer (kind=4), parameter :: FP_IEEE_NORMAL     = 5
    
    integer (kind=4), parameter :: VML_ACCURACY_HA    = 1
    integer (kind=4), parameter :: VML_ACCURACY_LA    = 2
    integer (kind=4), parameter :: VML_ACCURACY_EP    = 3

    real (kind=8), parameter :: FLOAT_MAXULP_HA          = 4.5
    real (kind=8), parameter :: FLOAT_MAXULP_LA          = 5.0
    real (kind=8), parameter :: FLOAT_MAXULP_EP          = 5D3
    real (kind=8), parameter :: FLOAT_COMPLEX_MAXULP_HA  = 4.5
    real (kind=8), parameter :: FLOAT_COMPLEX_MAXULP_LA  = 5.0
    real (kind=8), parameter :: FLOAT_COMPLEX_MAXULP_EP  = 5D3
    real (kind=8), parameter :: DOUBLE_MAXULP_HA         = 4.5
    real (kind=8), parameter :: DOUBLE_MAXULP_LA         = 5.0
    real (kind=8), parameter :: DOUBLE_MAXULP_EP         = 7D7
    real (kind=8), parameter :: DOUBLE_COMPLEX_MAXULP_HA = 4.5
    real (kind=8), parameter :: DOUBLE_COMPLEX_MAXULP_LA = 5.0
    real (kind=8), parameter :: DOUBLE_COMPLEX_MAXULP_EP = 7D7
    
    character (len=4), parameter, dimension(2) :: TEST_RESULT = (/'PASS','FAIL'/)
    character (len=2), parameter, dimension(3) :: VML_ACCURACY_NAME = (/'HA','LA','EP'/)
  
    real (kind=8), dimension(3) :: FLOAT_MAXULP
    real (kind=8), dimension(3) :: COMPLEX_FLOAT_MAXULP
    real (kind=8), dimension(3) :: DOUBLE_MAXULP
    real (kind=8), dimension(3) :: COMPLEX_DOUBLE_MAXULP
   
    common /heap/ FLOAT_MAXULP, COMPLEX_FLOAT_MAXULP, DOUBLE_MAXULP, COMPLEX_DOUBLE_MAXULP
