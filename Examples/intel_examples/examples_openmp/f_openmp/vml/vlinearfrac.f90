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
!*            LinearFrac example program text (OpenMP offload interface)
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
    real      (kind=4),allocatable :: varg1(:), varg2(:), vres1(:), vmres1(:), vref1(:)
    real      (kind=4) arg3, arg4, arg5, arg6
    integer   (kind=4) i, a, errs
    integer   (kind=4) VLEN
    parameter (VLEN = 4)
    integer   (kind=4) test_arg1(VLEN)
    integer   (kind=4) test_arg2(VLEN)
    integer   (kind=4) test_arg3(VLEN)
    integer   (kind=4) test_arg4(VLEN)
    integer   (kind=4) test_arg5(VLEN)
    integer   (kind=4) test_arg6(VLEN)
    integer   (kind=4) test_ref1(VLEN)
    integer   (kind=4) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'40D9B85C', & ! 6.80375481     
                     Z'40B52EFA', & ! 5.66198444     
                     Z'4103BA28', & ! 8.2329483      
                     Z'C052EA36'  / ! -3.2955451     
    data test_arg2 / Z'C007309A', & ! -2.1123414     
                     Z'40BF006A', & ! 5.96880054     
                     Z'C0C1912F', & ! -6.04897261    
                     Z'40ABAABC'  / ! 5.3645916      
    data test_arg3 / Z'4048F5C3', & ! 3.1400001      
                     Z'4048F5C3', & ! 3.1400001      
                     Z'4048F5C3', & ! 3.1400001      
                     Z'4048F5C3'  / ! 3.1400001      
    data test_arg4 / Z'40C8F5C3', & ! 6.28000021     
                     Z'40C8F5C3', & ! 6.28000021     
                     Z'40C8F5C3', & ! 6.28000021     
                     Z'40C8F5C3'  / ! 6.28000021     
    data test_arg5 / Z'4116B852', & ! 9.42000008     
                     Z'4116B852', & ! 9.42000008     
                     Z'4116B852', & ! 9.42000008     
                     Z'4116B852'  / ! 9.42000008     
    data test_arg6 / Z'4148F5C3', & ! 12.5600004     
                     Z'4148F5C3', & ! 12.5600004     
                     Z'4148F5C3', & ! 12.5600004     
                     Z'4148F5C3'  / ! 12.5600004     
    data test_ref1 / Z'C07117D3', & ! -3.76707911    
                     Z'3EB313C1', & ! 0.349760085    
                     Z'BF392C6C', & ! -0.723334074   
                     Z'BD840B70'  / ! -0.0644749403  
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
        varg1(i) = as_float(test_arg1(i))
        varg2(i) = as_float(test_arg2(i))
        arg3     = as_float(test_arg3(i))
        arg4     = as_float(test_arg4(i))
        arg5     = as_float(test_arg5(i))
        arg6     = as_float(test_arg6(i))
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

        !$omp target data map(varg1,varg2,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vres1)
        call vslinearfrac(VLEN, varg1, varg2, arg3, arg4, arg5, arg6, vres1)
        !$omp end target variant dispatch
        !$omp end target data
        !$omp target data map(varg1,varg2,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vmres1)
        call vmslinearfrac(VLEN, varg1, varg2, arg3, arg4, arg5, arg6, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_float(VML_ARG1_RES1, varg1(i), varg2(i), & 
                           vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_float(VML_ARG1_RES1, varg1(i), varg2(i), & 
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
    real      (kind=8),allocatable :: varg1(:), varg2(:), vres1(:), vmres1(:), vref1(:)
    real      (kind=8) arg3, arg4, arg5, arg6
    integer   (kind=4) i, a, errs
    integer   (kind=8) VLEN
    parameter (VLEN = 4)
    integer   (kind=8) test_arg1(VLEN)
    integer   (kind=8) test_arg2(VLEN)
    integer   (kind=8) test_arg3(VLEN)
    integer   (kind=8) test_arg4(VLEN)
    integer   (kind=8) test_arg5(VLEN)
    integer   (kind=8) test_arg6(VLEN)
    integer   (kind=8) test_ref1(VLEN)
    integer   (kind=8) nan_value
    integer   (kind=8) vml_accuracy_mode(3)
    data vml_accuracy_mode / VML_HA, VML_LA, VML_EP /
    integer tmode

    ! NaN value to fill result vector
    data  nan_value /Z'FFFFFFFFFFFFFFFF'/
    
    ! Arguments and reference results begin
    data test_arg1 / Z'401B370B60E66E18', & ! 6.80375434309419092      
                     Z'4016A5DF421D4BBE', & ! 5.66198447517211711      
                     Z'40207744D998EE8A', & ! 8.23294715873568705      
                     Z'C00A5D46A314BA8E'  / ! -3.2955448857022196      
    data test_arg2 / Z'C000E6134801CC26', & ! -2.11234146361813924     
                     Z'4017E00D485FC01A', & ! 5.96880066952146571      
                     Z'C0183225E080644C', & ! -6.04897261413232101     
                     Z'4015755793FAEAB0'  / ! 5.36459189623808186      
    data test_arg3 / Z'40091EB851EB851F', & ! 3.14000000000000012      
                     Z'40091EB851EB851F', & ! 3.14000000000000012      
                     Z'40091EB851EB851F', & ! 3.14000000000000012      
                     Z'40091EB851EB851F'  / ! 3.14000000000000012      
    data test_arg4 / Z'40191EB851EB851F', & ! 6.28000000000000025      
                     Z'40191EB851EB851F', & ! 6.28000000000000025      
                     Z'40191EB851EB851F', & ! 6.28000000000000025      
                     Z'40191EB851EB851F'  / ! 6.28000000000000025      
    data test_arg5 / Z'4022D70A3D70A3D7', & ! 9.41999999999999993      
                     Z'4022D70A3D70A3D7', & ! 9.41999999999999993      
                     Z'4022D70A3D70A3D7', & ! 9.41999999999999993      
                     Z'4022D70A3D70A3D7'  / ! 9.41999999999999993      
    data test_arg6 / Z'40291EB851EB851F', & ! 12.5600000000000005      
                     Z'40291EB851EB851F', & ! 12.5600000000000005      
                     Z'40291EB851EB851F', & ! 12.5600000000000005      
                     Z'40291EB851EB851F'  / ! 12.5600000000000005      
    data test_ref1 / Z'C00E22FA0DD4CC16', & ! -3.76707850270896483     
                     Z'3FD6627804818A52', & ! 0.349760059738547402     
                     Z'BFE7258D6AC8E3FB', & ! -0.723334034503863577    
                     Z'BFB0816DEA262188'  / ! -0.0644749352123935582   
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
        varg1(i) = as_double(test_arg1(i))
        varg2(i) = as_double(test_arg2(i))
        arg3     = as_double(test_arg3(i))
        arg4     = as_double(test_arg4(i))
        arg5     = as_double(test_arg5(i))
        arg6     = as_double(test_arg6(i))
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

        !$omp target data map(varg1,varg2,vres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vres1)
        call vdlinearfrac(VLEN, varg1, varg2, arg3, arg4, arg5, arg6, vres1)
        !$omp end target variant dispatch
        !$omp end target data
        !$omp target data map(varg1,varg2,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,varg2,vmres1)
        call vmdlinearfrac(VLEN, varg1, varg2, arg3, arg4, arg5, arg6, vmres1, vml_accuracy_mode(a))
        !$omp end target variant dispatch
        !$omp end target data
        ! Check results
        do i = 1, VLEN
          errs = errs + check_result_double(VML_ARG1_RES1, varg1(i), varg2(i), & 
                            vres1(i), vres1(i), vref1(i), vref1(i), "v"//funcname, a)
          errs = errs + check_result_double(VML_ARG1_RES1, varg1(i), varg2(i), & 
                            vmres1(i), vmres1(i), vref1(i), vref1(i), "vm"//funcname, a)
        enddo
    enddo
    
    test_double = errs

end function
! @brief Real double precision function test end

! @brief Main test program begin
program linearfrac_example
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    integer   (kind=4) :: test_float
    integer   (kind=4) :: test_double
    integer   (kind=4) errs, total_errs
    character (len = *), parameter :: funcname = "linearfrac"
    
    total_errs = 0

    data FLOAT_MAXULP /1.0D6,1.0D6,1.0D6/
    data COMPLEX_FLOAT_MAXULP /FLOAT_COMPLEX_MAXULP_HA,FLOAT_COMPLEX_MAXULP_LA,FLOAT_COMPLEX_MAXULP_EP/
    data DOUBLE_MAXULP /1.0D6,1.0D6,1.0D6/
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
