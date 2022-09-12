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

!*
!
!*  Content:
!*            CIS example program text (OpenMP offload interface)
!*
!*******************************************************************************/

include "mkl_omp_offload.f90"
include "_vml_common_functions.f90"

! @brief Complex single precision function test begin
integer (kind=4) function test_float_complex(funcname)
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    character (len = *) :: funcname
    real      (kind=4)  :: as_float
    integer   (kind=4)  :: check_result_float_complex
    real      (kind=4),allocatable :: varg1(:)
    complex   (kind=4) carg1
    complex   (kind=4),allocatable :: vres1(:), vmres1(:), vref1(:)
    integer   (kind=4) i, a, errs
    integer   (kind=4) VLEN
    parameter (VLEN = 4)
    integer   (kind=4) test_arg1(VLEN)
    integer   (kind=4) test_ref1(2*VLEN)
    integer   (kind=4) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'40D9B85C', & ! 6.80375481     
                     Z'C007309A', & ! -2.1123414     
                     Z'40B52EFA', & ! 5.66198444     
                     Z'40BF006A'  / ! 5.96880054     
    data test_ref1 / Z'3F5E16D8', Z'3EFEA7D8', & ! 0.867536068     + i * 0.497374296    
                     Z'BF03F53A', Z'BF5B5EAB', & ! -0.51546061     + i * -0.856913269   
                     Z'3F502C93', Z'BF14FEBF', & ! 0.813180149     + i * -0.582012117   
                     Z'3F7373DF', Z'BE9E5396'  / ! 0.950986803     + i * -0.30923146    
    ! Arguments and reference results end
    
    errs = 0

    ! Allocate vectors
    allocate(varg1(VLEN))
    allocate(vres1(VLEN))
    allocate(vmres1(VLEN))
    allocate(vref1(VLEN))

    ! Fill vectors
    do i = 1, VLEN
        varg1(i) = as_float(test_arg1(i))
        vref1(i) = CMPLX(as_float(test_ref1(2*i-1)), as_float(test_ref1(2*i)), 4)
        vres1(i) = CMPLX(as_float(nan_value),as_float(nan_value),4)
        vmres1(i) = CMPLX(as_float(nan_value),as_float(nan_value),4)
    enddo

    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vccis(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmccis(VLEN, varg1, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          carg1 = CMPLX(varg1(i), 0.0, 4)
          errs = errs + check_result_float_complex(VML_ARG1R_RES1C, carg1, carg1, & 
                           vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_float_complex(VML_ARG1R_RES1C, carg1, carg1, & 
                           vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo

    test_float_complex = errs

end function
! @brief Complex single precision function test end

! @brief Complex double precision function test begin
integer (kind=4) function test_double_complex(funcname)
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    character (len = *) :: funcname
    real      (kind=8) :: as_double
    integer   (kind=4) :: check_result_double_complex
    real      (kind=8),allocatable :: varg1(:)
    complex   (kind=8) carg1
    complex   (kind=8),allocatable :: vres1(:), vmres1(:), vref1(:)
    integer   (kind=4) i, a, errs
    integer   (kind=8) VLEN
    parameter (VLEN = 4)
    integer   (kind=8) test_arg1(VLEN)
    integer   (kind=8) test_ref1(2*VLEN)
    integer   (kind=8) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFFFFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'401B370B60E66E18', & ! 6.80375434309419092      
                     Z'C000E6134801CC26', & ! -2.11234146361813924     
                     Z'4016A5DF421D4BBE', & ! 5.66198447517211711      
                     Z'4017E00D485FC01A'  / ! 5.96880066952146571      
    data test_ref1 / Z'3FEBC2DB7AB89950', Z'3FDFD4F93E99B2E0', & ! 0.867536296548488295      + i * 0.497373877652348639     
                     Z'BFE07EA757C4010B', Z'BFEB6BD5549D70BC', & ! -0.51546065465666524      + i * -0.85691324735991925     
                     Z'3FEA05925DBF776B', Z'BFE29FD7C840E7D0', & ! 0.813180144406698502      + i * -0.582012072677793313    
                     Z'3FEE6E7BF8882000', Z'BFD3CA723281D19B'  / ! 0.950986848271895724      + i * -0.309231328318924248    
    ! Arguments and reference results end
    
    errs = 0

    ! Allocate vectors
    allocate(varg1(VLEN))
    allocate(vres1(VLEN))
    allocate(vmres1(VLEN))
    allocate(vref1(VLEN))

    ! Fill vectors
    do i = 1, VLEN
        varg1(i) = as_double(test_arg1(i))
        vref1(i) = CMPLX(as_double(test_ref1(2*i-1)), as_double(test_ref1(2*i)), 8)
        vres1(i) = CMPLX(as_double(nan_value),as_double(nan_value),8)
        vmres1(i) = CMPLX(as_double(nan_value),as_double(nan_value),8)
    enddo

    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vzcis(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmzcis(VLEN, varg1, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          carg1 = CMPLX(varg1(i), 0.0,8)
          errs  = errs + check_result_double_complex(VML_ARG1R_RES1C, carg1, carg1, & 
                            vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs  = errs + check_result_double_complex(VML_ARG1R_RES1C, carg1, carg1, & 
                            vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo
    
    test_double_complex = errs

end function
! @brief Complex double precision function test end

! @brief Main test program begin
program cis_example
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    integer   (kind=4) :: test_float_complex
    integer   (kind=4) :: test_double_complex
    integer   (kind=4) errs, total_errs
    character (len = *), parameter :: funcname = "cis"
    
    total_errs = 0

    data FLOAT_MAXULP /FLOAT_MAXULP_HA,FLOAT_MAXULP_LA,FLOAT_MAXULP_EP/
    data COMPLEX_FLOAT_MAXULP /FLOAT_COMPLEX_MAXULP_HA,FLOAT_COMPLEX_MAXULP_LA,FLOAT_COMPLEX_MAXULP_EP/
    data DOUBLE_MAXULP /DOUBLE_MAXULP_HA,DOUBLE_MAXULP_LA,DOUBLE_MAXULP_EP/
    data COMPLEX_DOUBLE_MAXULP /DOUBLE_COMPLEX_MAXULP_HA,DOUBLE_COMPLEX_MAXULP_LA,DOUBLE_COMPLEX_MAXULP_EP/

    write (*, 111) funcname
    111 format ('Running ', A, ' functions:')

    ! Single precision complex test run begin
    write (*, 115) TAB, funcname
    115 format(A, 'Running ',  A, ' with single precision complex data type:')
    errs = test_float_complex(funcname)    
    write (*, 116) TAB, funcname, TEST_RESULT(blend_int32((errs>0),2,1))
    116 format(A, A, ' single precision complex result: ', A)
    ! Single precision complex test run end

    ! Complex double precision test run begin
    write (*, 119) TAB, funcname
    119 format(A, 'Running ',  A, ' with double precision complex data type:')
    errs = test_double_complex(funcname)    
    write (*, 120) TAB, funcname, TEST_RESULT(blend_int32((errs>0),2,1))
    120 format(A, A, ' double precision complex result: ', A)
    ! Complex double precision  test run end

    write (*, 121) funcname, TEST_RESULT(blend_int32((total_errs>0),2,1))
    121 format(A, ' function result: ', A)

end program
! @brief Main test program end
