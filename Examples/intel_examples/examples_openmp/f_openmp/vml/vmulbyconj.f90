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
!*            MulByConj example program text (OpenMP offload interface)
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
    complex      (kind=4),allocatable :: varg1(:), varg2(:), vres1(:), vmres1(:), vref1(:)
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
    data test_arg1 / Z'C007309A', Z'40D9B85C', & ! -2.1123414      + i * 6.80375481     
                     Z'C0C1912F', Z'4103BA28', & ! -6.04897261     + i * 8.2329483      
                     Z'3F8A29C0', Z'C08E3964', & ! 1.07939911      + i * -4.44450569    
                     Z'3E8939C0', Z'C02D136C'  / ! 0.268018723     + i * -2.70431042    
    data test_arg2 / Z'40BF006A', Z'40B52EFA', & ! 5.96880054      + i * 5.66198444     
                     Z'40ABAABC', Z'C052EA36', & ! 5.3645916       + i * -3.2955451     
                     Z'4024F46C', Z'BEE77440', & ! 2.57741833      + i * -0.452058792   
                     Z'41052EB4', Z'4110B6A8'  / ! 8.32390213      + i * 9.04459381     
    data test_ref1 / Z'41CF511F', Z'425247FD', & ! 25.9146099      + i * 52.5703011     
                     Z'C26E544C', Z'41C1DA9C', & ! -59.5823212     + i * 24.2317429     
                     Z'409951D8', Z'C12F7A77', & ! 4.79124069      + i * -10.9673986    
                     Z'C1B1D3D2', Z'C1C779EE'  / ! -22.2284279     + i * -24.934536     
    ! Arguments and reference results end
    
    errs = 0

    ! Allocate vectors
    allocate(varg1(VLEN))
    allocate(varg2(VLEN))
    allocate(vres1(VLEN))
    allocate(vmres1(VLEN))
    allocate(vref1(VLEN))

    ! Fill vectors
    do i = 1, VLEN
        varg1(i) = CMPLX(as_float(test_arg1(2*i-1)), as_float(test_arg1(2*i)), 4)
        varg2(i) = CMPLX(as_float(test_arg2(2*i-1)), as_float(test_arg2(2*i)), 4)
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

        !$omp target data map(varg1,varg2,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vres1)
        call vcmulbyconj(VLEN, varg1, varg2, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,varg2,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vmres1)
        call vmcmulbyconj(VLEN, varg1, varg2, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_float_complex(VML_ARG2_RES1, varg1(i), varg2(i), & 
                           vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_float_complex(VML_ARG2_RES1, varg1(i), varg2(i), & 
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
    complex   (kind=8),allocatable :: varg1(:), varg2(:), vres1(:), vmres1(:), vref1(:)
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
    data test_arg1 / Z'C000E6134801CC26', Z'401B370B60E66E18', & ! -2.11234146361813924      + i * 6.80375434309419092      
                     Z'C0183225E080644C', Z'40207744D998EE8A', & ! -6.04897261413232101      + i * 8.23294715873568705      
                     Z'3FF1453801E28A70', Z'C011C72C86338E59', & ! 1.07939911590861115       + i * -4.44450578393624429     
                     Z'3FD12735D3224E60', Z'C005A26D910B44DC'  / ! 0.268018203912310682      + i * -2.70431054416313366     
    data test_arg2 / Z'4017E00D485FC01A', Z'4016A5DF421D4BBE', & ! 5.96880066952146571       + i * 5.66198447517211711      
                     Z'4015755793FAEAB0', Z'C00A5D46A314BA8E', & ! 5.36459189623808186       + i * -3.2955448857022196      
                     Z'40049E8D96893D1C', Z'BFDCEE88B739DD20', & ! 2.57741849523848821       + i * -0.452058962756796134    
                     Z'4020A5D666294BAC', Z'402216D5173C2DAA'  / ! 8.32390136007401082       + i * 9.04459450349425609      
    data test_ref1 / Z'4039EA23A3CE3157', Z'404A48FF86CC45BE', & ! 25.9146063211822728       + i * 52.5702980515884377      
                     Z'C04DCA8957C19C5A', Z'40383B535E41D16B', & ! -59.5823163695683462      + i * 24.2317408476532528      
                     Z'40132A3B52620689', Z'C025EF4EF83576E2', & ! 4.79124191973972646       + i * -10.9673993649734633     
                     Z'C0363A7ABA8C5B9C', Z'C038EF3D5EDD987E'  / ! -22.2284351914091616      + i * -24.9345301905636845     
    ! Arguments and reference results end
    
    errs = 0

    ! Allocate vectors
    allocate(varg1(VLEN))
    allocate(varg2(VLEN))
    allocate(vres1(VLEN))
    allocate(vmres1(VLEN))
    allocate(vref1(VLEN))

    ! Fill vectors
    do i = 1, VLEN
        varg1(i) = CMPLX(as_double(test_arg1(2*i-1)), as_double(test_arg1(2*i)), 8)
        varg2(i) = CMPLX(as_double(test_arg2(2*i-1)), as_double(test_arg2(2*i)), 8)
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

        !$omp target data map(varg1,varg2,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vres1)
        call vzmulbyconj(VLEN, varg1, varg2, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,varg2,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vmres1)
        call vmzmulbyconj(VLEN, varg1, varg2, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_double_complex(VML_ARG2_RES1, varg1(i), varg2(i), & 
                            vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_double_complex(VML_ARG2_RES1, varg1(i), varg2(i), & 
                            vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo
    
    test_double_complex = errs

end function
! @brief Complex double precision function test end

! @brief Main test program begin
program mulbyconj_example
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    integer   (kind=4) :: test_float_complex
    integer   (kind=4) :: test_double_complex
    integer   (kind=4) errs, total_errs
    character (len = *), parameter :: funcname = "mulbyconj"
    
    total_errs = 0

    data FLOAT_MAXULP /FLOAT_MAXULP_HA,FLOAT_MAXULP_LA,FLOAT_MAXULP_EP/
    data COMPLEX_FLOAT_MAXULP /16.0,16.0,16.0/
    data DOUBLE_MAXULP /DOUBLE_MAXULP_HA,DOUBLE_MAXULP_LA,DOUBLE_MAXULP_EP/
    data COMPLEX_DOUBLE_MAXULP /16.0,16.0,16.0/

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
