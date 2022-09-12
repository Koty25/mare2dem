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
!*            CdfNormInv example program text (OpenMP offload interface)
!*
!*******************************************************************************/

include "mkl_omp_offload.f90"
include "_vml_common_functions.f90"

! @brief Real single precision function test begin
integer (kind=4) function test_float(funcname)
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    character (len = *) :: funcname
    real      (kind=4)  :: as_float
    integer   (kind=4)  :: check_result_float
    real      (kind=4),allocatable :: varg1(:), vres1(:), vmres1(:), vref1(:)
    integer   (kind=4) i, a, errs
    integer   (kind=4) VLEN
    parameter (VLEN = 4)
    integer   (kind=4) test_arg1(VLEN)
    integer   (kind=4) test_ref1(VLEN)
    integer   (kind=4) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'3C82EB10', & ! 0.0159812272   
                     Z'3D780F8E', & ! 0.0605617091   
                     Z'3CB1AF64', & ! 0.0216900781   
                     Z'3CA51E30'  / ! 0.0201559961   
    data test_ref1 / Z'C00945B6', & ! -2.14487982    
                     Z'BFC668DF', & ! -1.55007541    
                     Z'C0014831', & ! -2.02003121    
                     Z'C0033C02'  / ! -2.05053759    
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
        vref1(i) = as_float(test_ref1(i))
        vres1(i) = as_float(nan_value)
        vmres1(i) = as_float(nan_value)
    enddo

    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vscdfnorminv(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmscdfnorminv(VLEN, varg1, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_float(VML_ARG1_RES1, varg1(i), varg1(i), & 
                           vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_float(VML_ARG1_RES1, varg1(i), varg1(i), & 
                           vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo

    test_float = errs

end function
! @brief Real single precision function test end

! @brief Real double precision function test begin
integer (kind=4) function test_double(funcname)
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    character (len = *) :: funcname
    real      (kind=8) :: as_double
    integer   (kind=4) :: check_result_double
    real      (kind=8),allocatable :: varg1(:), vres1(:), vmres1(:), vref1(:)
    integer   (kind=4) i, a, errs
    integer   (kind=8) VLEN
    parameter (VLEN = 4)
    integer   (kind=8) test_arg1(VLEN)
    integer   (kind=8) test_ref1(VLEN)
    integer   (kind=8) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFFFFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'3F905D621353EDF8', & ! 0.0159812282845290532    
                     Z'3FAF01F1B0A46A4A', & ! 0.060561707318090699     
                     Z'3F9635EC782C6BD8', & ! 0.0216900776241394089    
                     Z'3F94A3C609C2E128'  / ! 0.020155996652392677     
    data test_ref1 / Z'C00128B6C129E3D3', & ! -2.14487982663238475     
                     Z'BFF8CD1BE28091B7', & ! -1.5500754211180785      
                     Z'C0002906191BD5FA', & ! -2.02003116241644154     
                     Z'C0006780490CB477'  / ! -2.05053765363714602     
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
        vref1(i) = as_double(test_ref1(i))
        vres1(i) = as_double(nan_value)
        vmres1(i) = as_double(nan_value)
    enddo

    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vdcdfnorminv(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmdcdfnorminv(VLEN, varg1, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_double(VML_ARG1_RES1, varg1(i), varg1(i), & 
                            vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_double(VML_ARG1_RES1, varg1(i), varg1(i), & 
                            vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo
    
    test_double = errs

end function
! @brief Real double precision function test end

! @brief Main test program begin
program cdfnorminv_example
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    integer   (kind=4) :: test_float
    integer   (kind=4) :: test_double
    integer   (kind=4) errs, total_errs
    character (len = *), parameter :: funcname = "cdfnorminv"
    
    total_errs = 0

    data FLOAT_MAXULP /HUGE(0.0),HUGE(0.0),HUGE(0.0)/
    data COMPLEX_FLOAT_MAXULP /FLOAT_COMPLEX_MAXULP_HA,FLOAT_COMPLEX_MAXULP_LA,FLOAT_COMPLEX_MAXULP_EP/
    data DOUBLE_MAXULP /HUGE(0D0),HUGE(0D0),HUGE(0D0)/
    data COMPLEX_DOUBLE_MAXULP /DOUBLE_COMPLEX_MAXULP_HA,DOUBLE_COMPLEX_MAXULP_LA,DOUBLE_COMPLEX_MAXULP_EP/

    write (*, 111) funcname
    111 format ('Running ', A, ' functions:')

    ! Single precision test run begin
    write (*, 112) TAB, funcname
    112 format(A, 'Running ',  A, ' with single precision real data type:')
    errs = test_float(funcname)    
    write (*, 113) TAB, funcname, TEST_RESULT(blend_int32((errs>0),2,1))
    113 format(A, A, ' single precision real result: ', A)
    ! Single precision test run end

    ! Real double precision test run begin
    write (*, 117) TAB, funcname
    117 format(A, 'Running ',  A, ' with double precision real data type:')
    errs = test_double(funcname)    
    write (*, 118) TAB, funcname, TEST_RESULT(blend_int32((errs>0),2,1))
    118 format(A, A, ' double precision real result: ', A)
    ! Real double precision test run end

    write (*, 121) funcname, TEST_RESULT(blend_int32((total_errs>0),2,1))
    121 format(A, ' function result: ', A)

end program
! @brief Main test program end
