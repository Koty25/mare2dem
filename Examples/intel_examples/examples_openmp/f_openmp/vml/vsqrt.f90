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
!*            Sqrt example program text (OpenMP offload interface)
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
    data test_arg1 / Z'41093E24', & ! 8.57767105     
                     Z'4093852F', & ! 4.61000776     
                     Z'41011D03', & ! 8.06958294     
                     Z'41034C40'  / ! 8.20611572     
    data test_ref1 / Z'403B70E7', & ! 2.92876601     
                     Z'400969F8', & ! 2.14709282     
                     Z'4035CE0C', & ! 2.8407011      
                     Z'40375621'  / ! 2.86463189     
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
        call vssqrt(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmssqrt(VLEN, varg1, vmres1, vml_accuracy_mode(a))
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
    data test_arg1 / Z'402127C473A3E923', & ! 8.57767068267691535      
                     Z'401270A5F32DAE19', & ! 4.6100080486899282       
                     Z'402023A0651C4741', & ! 8.06958309145159269      
                     Z'40206988134D9FDD'  / ! 8.20611629793705255      
    data test_ref1 / Z'40076E1CE786FE27', & ! 2.9287660682746437       
                     Z'40012D3F0ED3A6F0', & ! 2.14709292968188237      
                     Z'4006B9C187E1F54A', & ! 2.84070116194076139      
                     Z'4006EAC429F8399F'  / ! 2.86463196553013644      
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
        call vdsqrt(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmdsqrt(VLEN, varg1, vmres1, vml_accuracy_mode(a))
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

! @brief Complex single precision function test begin
integer (kind=4) function test_float_complex(funcname)
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    character (len = *) :: funcname
    real      (kind=4)  :: as_float
    integer   (kind=4)  :: check_result_float_complex
    complex      (kind=4),allocatable :: varg1(:), vres1(:), vmres1(:), vref1(:)
    integer   (kind=4) i, a, errs
    integer   (kind=4) VLEN
    parameter (VLEN = 4)
    integer   (kind=4) test_arg1(2*VLEN)
    integer   (kind=4) test_ref1(2*VLEN)
    integer   (kind=4) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'4093852F', Z'41093E24', & ! 4.61000776      + i * 8.57767105     
                     Z'41034C40', Z'41011D03', & ! 8.20611572      + i * 8.06958294     
                     Z'4036ECDE', Z'41136B29', & ! 2.85820723      + i * 9.21366215     
                     Z'40FDFDE5', Z'4082ABE3'  / ! 7.93724298      + i * 4.08348227     
    data test_ref1 / Z'402B6B72', Z'3FCCF5B2', & ! 2.67843294      + i * 1.60124803     
                     Z'4048F083', Z'3FA47E0B', & ! 3.13967967      + i * 1.28509653     
                     Z'40200838', Z'3FEBD28B', & ! 2.50050163      + i * 1.84236276     
                     Z'4039D6BB', Z'3F34013F'  / ! 2.90373111      + i * 0.703144014    
    ! Arguments and reference results end
    
    errs = 0

    ! Allocate vectors
    allocate(varg1(VLEN))
    allocate(vres1(VLEN))
    allocate(vmres1(VLEN))
    allocate(vref1(VLEN))

    ! Fill vectors
    do i = 1, VLEN
        varg1(i) = CMPLX(as_float(test_arg1(2*i-1)), as_float(test_arg1(2*i)), 4)
        vref1(i) = CMPLX(as_float(test_ref1(2*i-1)), as_float(test_ref1(2*i)), 4)
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
        call vcsqrt(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmcsqrt(VLEN, varg1, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_float_complex(VML_ARG1_RES1, varg1(i), varg1(i), & 
                           vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_float_complex(VML_ARG1_RES1, varg1(i), varg1(i), & 
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
    complex   (kind=8),allocatable :: varg1(:), vres1(:), vmres1(:), vref1(:)
    integer   (kind=4) i, a, errs
    integer   (kind=8) VLEN
    parameter (VLEN = 4)
    integer   (kind=8) test_arg1(2*VLEN)
    integer   (kind=8) test_ref1(2*VLEN)
    integer   (kind=8) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFFFFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'401270A5F32DAE19', Z'402127C473A3E923', & ! 4.6100080486899282        + i * 8.57767068267691535      
                     Z'40206988134D9FDD', Z'402023A0651C4741', & ! 8.20611629793705255       + i * 8.06958309145159269      
                     Z'4006DD9BBAC0EE6B', Z'40226D6509CA7464', & ! 2.8582071867111174        + i * 9.21366148563738108      
                     Z'401FBFBCBB737F7A', Z'4010557C717977C6'  / ! 7.93724339382594657       + i * 4.0834825258625127       
    data test_ref1 / Z'40056D6E42511161', Z'3FF99EB631575346', & ! 2.67843295869731479       + i * 1.60124797128556073      
                     Z'40091E1072F5F2ED', Z'3FF48FC159BA032B', & ! 3.13967981160236898       + i * 1.28509650277573928      
                     Z'40040106EB21DA6F', Z'3FFD7A515924BC6F', & ! 2.50050147721349658       + i * 1.84236273595504563      
                     Z'40073AD76D49B864', Z'3FE68027E8B0FCB8'  / ! 2.90373120671488216       + i * 0.703144030070595782     
    ! Arguments and reference results end
    
    errs = 0

    ! Allocate vectors
    allocate(varg1(VLEN))
    allocate(vres1(VLEN))
    allocate(vmres1(VLEN))
    allocate(vref1(VLEN))

    ! Fill vectors
    do i = 1, VLEN
        varg1(i) = CMPLX(as_double(test_arg1(2*i-1)), as_double(test_arg1(2*i)), 8)
        vref1(i) = CMPLX(as_double(test_ref1(2*i-1)), as_double(test_ref1(2*i)), 8)
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
        call vzsqrt(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmzsqrt(VLEN, varg1, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_double_complex(VML_ARG1_RES1, varg1(i), varg1(i), & 
                            vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_double_complex(VML_ARG1_RES1, varg1(i), varg1(i), & 
                            vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo
    
    test_double_complex = errs

end function
! @brief Complex double precision function test end

! @brief Main test program begin
program sqrt_example
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    integer   (kind=4) :: test_float
    integer   (kind=4) :: test_float_complex
    integer   (kind=4) :: test_double
    integer   (kind=4) :: test_double_complex
    integer   (kind=4) errs, total_errs
    character (len = *), parameter :: funcname = "sqrt"
    
    total_errs = 0

    data FLOAT_MAXULP /FLOAT_MAXULP_HA,FLOAT_MAXULP_LA,HUGE(0.0)/
    data COMPLEX_FLOAT_MAXULP /HUGE(0.0),HUGE(0.0),HUGE(0.0)/
    data DOUBLE_MAXULP /DOUBLE_MAXULP_HA,DOUBLE_MAXULP_LA,DOUBLE_MAXULP_EP/
    data COMPLEX_DOUBLE_MAXULP /HUGE(0D0),HUGE(0D0),HUGE(0D0)/

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
