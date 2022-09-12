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
!*            Atanh example program text (OpenMP offload interface)
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
    data test_ref1 / Z'3C82EDEB', & ! 0.0159825888   
                     Z'3D785D5D', & ! 0.0606359132   
                     Z'3CB1B687', & ! 0.0216934811   
                     Z'3CA523EA'  / ! 0.0201587267   
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
        call vsatanh(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmsatanh(VLEN, varg1, vmres1, vml_accuracy_mode(a))
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
    data test_ref1 / Z'3F905DBD64B240D8', & ! 0.0159825890264649606    
                     Z'3FAF0BAB94755F2D', & ! 0.0606359118198141825    
                     Z'3F9636D0CCDE1024', & ! 0.0216934800187261884    
                     Z'3F94A47D42904F8C'  / ! 0.0201587268712298123    
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
        call vdatanh(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmdatanh(VLEN, varg1, vmres1, vml_accuracy_mode(a))
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
    data test_arg1 / Z'3D780F8E', Z'3C82EB10', & ! 0.0605617091    + i * 0.0159812272   
                     Z'3CA51E30', Z'3CB1AF64', & ! 0.0201559961    + i * 0.0216900781   
                     Z'3DA4576B', Z'3C10C1C8', & ! 0.0802448615    + i * 0.00883526355  
                     Z'3CBDDDC8', Z'3D88257A'  / ! 0.0231770426    + i * 0.0664777309   
    data test_ref1 / Z'3D784D08', Z'3C836386', & ! 0.0606203377    + i * 0.0160386674   
                     Z'3CA51005', Z'3CB1BABB', & ! 0.0201492403    + i * 0.0216954853   
                     Z'3DA4AEBF', Z'3C11B0F3', & ! 0.0804114267    + i * 0.00889228564  
                     Z'3CBD1067', Z'3D8804D6'  / ! 0.0230791103    + i * 0.0664154738   
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
        call vcatanh(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmcatanh(VLEN, varg1, vmres1, vml_accuracy_mode(a))
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
    data test_arg1 / Z'3FAF01F1B0A46A4A', Z'3F905D621353EDF8', & ! 0.060561707318090699      + i * 0.0159812282845290532    
                     Z'3F94A3C609C2E128', Z'3F9635EC782C6BD8', & ! 0.020155996652392677      + i * 0.0216900776241394089    
                     Z'3FB48AED668F7C42', Z'3F821839168A96D8', & ! 0.0802448630706616151     + i * 0.00883526420632156639   
                     Z'3F97BBB8DC2F7770', Z'3FB104AF24553C92'  / ! 0.0231770405188095885     + i * 0.0664777244285111035    
    data test_ref1 / Z'3FAF09A0E2B69FE1', Z'3F906C70CDA5C815', & ! 0.060620334315274034      + i * 0.0160386682050060632    
                     Z'3F94A200AD66E706', Z'3F9637574ED5D99F', & ! 0.020149241050353893      + i * 0.0216954843394546702    
                     Z'3FB495D7D955D65A', Z'3F82361E797C0E88', & ! 0.0804114251712574613     + i * 0.00889228637925688903   
                     Z'3F97A20CB1D99626', Z'3FB1009AB2439241'  / ! 0.023079107623195004      + i * 0.0664154706206057238    
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
        call vzatanh(VLEN, varg1, vres1)
        !$omp end target variant dispatch
        !$omp end target data
		
        !$omp target data map(varg1,vmres1)
        !$omp target variant dispatch use_device_ptr(varg1,vmres1)
        call vmzatanh(VLEN, varg1, vmres1, vml_accuracy_mode(a))
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
program atanh_example
    use onemkl_vml_omp_offload
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    integer   (kind=4) :: test_float
    integer   (kind=4) :: test_float_complex
    integer   (kind=4) :: test_double
    integer   (kind=4) :: test_double_complex
    integer   (kind=4) errs, total_errs
    character (len = *), parameter :: funcname = "atanh"
    
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
