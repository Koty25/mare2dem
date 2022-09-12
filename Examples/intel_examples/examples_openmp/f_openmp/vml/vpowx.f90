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
!*            Powx example program text (OpenMP offload interface)
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
    real      (kind=4) arg2
    integer   (kind=4) i, a, errs
    integer   (kind=4) VLEN
    parameter (VLEN = 4)
    integer   (kind=4) test_arg1(VLEN)
    integer   (kind=4) test_arg2(VLEN)
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
    data test_arg2 / Z'4048F5C3', & ! 3.1400001      
                     Z'4048F5C3', & ! 3.1400001      
                     Z'4048F5C3', & ! 3.1400001      
                     Z'4048F5C3'  / ! 3.1400001      
    data test_ref1 / Z'44552AC1', & ! 852.66803      
                     Z'42F2B0D7', & ! 121.34539      
                     Z'442FF9C3', & ! 703.902527     
                     Z'44397EBA'  / ! 741.980103     
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
    arg2 = as_float(test_arg2(1))


    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vspowx(VLEN, varg1, arg2, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmspowx(VLEN, varg1, arg2, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_float(VML_ARG2_RES1, varg1(i), arg2, & 
                           vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_float(VML_ARG2_RES1, varg1(i), arg2, & 
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
    real      (kind=8) arg2
    integer   (kind=4) i, a, errs
    integer   (kind=8) VLEN
    parameter (VLEN = 4)
    integer   (kind=8) test_arg1(VLEN)
    integer   (kind=8) test_arg2(VLEN)
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
    data test_arg2 / Z'40091EB851EB851F', & ! 3.14000000000000012      
                     Z'40091EB851EB851F', & ! 3.14000000000000012      
                     Z'40091EB851EB851F', & ! 3.14000000000000012      
                     Z'40091EB851EB851F'  / ! 3.14000000000000012      
    data test_ref1 / Z'408AA55778FCE709', & ! 852.667711234856256      
                     Z'405E561AE4453C46', & ! 121.345391337979066      
                     Z'4085FF381E850103', & ! 703.902401961415649      
                     Z'40872FD74F3829E4'  / ! 741.980131567743683      
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
    arg2 = as_double(test_arg2(1))

    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vdpowx(VLEN, varg1, arg2, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmdpowx(VLEN, varg1, arg2, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_double(VML_ARG2_RES1, varg1(i), arg2, & 
                            vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_double(VML_ARG2_RES1, varg1(i), arg2, & 
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
    complex   (kind=4),allocatable :: varg1(:), vres1(:), vmres1(:), vref1(:)
    complex   (kind=4) arg2
    integer   (kind=4) i, a, errs
    integer   (kind=4) VLEN
    parameter (VLEN = 4)
    integer   (kind=4) test_arg1(2*VLEN)
    integer   (kind=4) test_arg2(2*VLEN)
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
    data test_arg2 / Z'4048F5C3', Z'4048F5C3', & ! 3.1400001       + i * 3.1400001      
                     Z'4048F5C3', Z'4048F5C3', & ! 3.1400001       + i * 3.1400001      
                     Z'4048F5C3', Z'4048F5C3', & ! 3.1400001       + i * 3.1400001      
                     Z'4048F5C3', Z'4048F5C3'  / ! 3.1400001       + i * 3.1400001      
    data test_ref1 / Z'C19A8841', Z'C21A0136', & ! -19.3165302     + i * -38.5011826    
                     Z'C310B7AE', Z'C2ED2BD0', & ! -144.717499     + i * -118.585571    
                     Z'401FC609', Z'C1B5CAED', & ! 2.49646211      + i * -22.7240849    
                     Z'C2D4B709', Z'433D8524'  / ! -106.357491     + i * 189.520081     
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
    arg2 = CMPLX(as_float(test_arg2(1)), as_float(test_arg2(2)), 4)

    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vcpowx(VLEN, varg1, arg2, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmcpowx(VLEN, varg1, arg2, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_float_complex(VML_ARG2_RES1, varg1(i), arg2, & 
                           vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_float_complex(VML_ARG2_RES1, varg1(i), arg2, & 
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
    complex   (kind=8) arg2
    integer   (kind=4) i, a, errs
    integer   (kind=8) VLEN
    parameter (VLEN = 4)
    integer   (kind=8) test_arg1(2*VLEN)
    integer   (kind=8) test_arg2(2*VLEN)
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
    data test_arg2 / Z'40091EB851EB851F', Z'40091EB851EB851F', & ! 3.14000000000000012       + i * 3.14000000000000012      
                     Z'40091EB851EB851F', Z'40091EB851EB851F', & ! 3.14000000000000012       + i * 3.14000000000000012      
                     Z'40091EB851EB851F', Z'40091EB851EB851F', & ! 3.14000000000000012       + i * 3.14000000000000012      
                     Z'40091EB851EB851F', Z'40091EB851EB851F'  / ! 3.14000000000000012       + i * 3.14000000000000012      
    data test_ref1 / Z'C033510970136A6A', Z'C043402653EE7617', & ! -19.3165502593423426      + i * -38.5011696733819733     
                     Z'C06216F614F1FE9F', Z'C05DA579773A3C38', & ! -144.717539284368257      + i * -118.585538679952947     
                     Z'4003F8B937C10DF9', Z'C036B95D4F1E546A', & ! 2.49644702489763093       + i * -22.7240800332114432     
                     Z'C05A96E0D33FD2E3', Z'4067B0A4892EA43F'  / ! -106.357472240760714      + i * 189.520084944817398      
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
    arg2 = CMPLX(as_double(test_arg2(1)), as_double(test_arg2(2)), 8)

    ! Loop by three accuracy flavors
    do a = 1, 3
        ! Call VML function with specific accuracy flavor
        !$omp target variant dispatch
        tmode = vmlsetmode(vml_accuracy_mode(a))
        !$omp end target variant dispatch

        !$omp target data map(varg1,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,vres1)
        call vzpowx(VLEN, varg1, arg2, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmzpowx(VLEN, varg1, arg2, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_double_complex(VML_ARG2_RES1, varg1(i), arg2, & 
                            vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_double_complex(VML_ARG2_RES1, varg1(i), arg2, & 
                            vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo
    
    test_double_complex = errs

end function
! @brief Complex double precision function test end

! @brief Main test program begin
program powx_example
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    integer   (kind=4) :: test_float
    integer   (kind=4) :: test_float_complex
    integer   (kind=4) :: test_double
    integer   (kind=4) :: test_double_complex
    integer   (kind=4) errs, total_errs
    character (len = *), parameter :: funcname = "powx"
    
    total_errs = 0

    data FLOAT_MAXULP /FLOAT_MAXULP_HA,FLOAT_MAXULP_LA,FLOAT_MAXULP_EP/
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
