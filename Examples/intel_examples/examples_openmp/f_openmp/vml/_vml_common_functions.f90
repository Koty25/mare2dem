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

!  Content: common routines for examples
!*******************************************************************************

! ===========================================================================
!
! @brief Ternary operator
!
! ===========================================================================
integer (kind=4) function blend_int32(cond, v1, v2)

    implicit none
    logical  cond
    integer (kind=4) v1
    integer (kind=4) v2

    if (cond .eq. .true.) then
        blend_int32 = v1
    else
        blend_int32 = v2
    end if

end function
! ===========================================================================

! ===========================================================================
!
! @brief Cast int32 to float
!
! ===========================================================================
real (kind=4) function as_float(v)

    implicit none
    integer (kind=4) v
    real    (kind=4) f
    integer (kind=4) d

    equivalence (f, d)
    d = v
    as_float = f

end function
! ===========================================================================

! ===========================================================================
!
! @brief Cast int64 to double
!
! ===========================================================================
real (kind=8) function as_double(v)

    implicit none
    integer (kind=8) v
    real    (kind=8) d
    integer (kind=8) q

    equivalence (d, q)
    q = v
    as_double = d

end function
! ===========================================================================

! ===========================================================================
!
! @brief Cast float to int32
!
! ===========================================================================
integer (kind=4) function as_int32(v)

    implicit none
    real    (kind=4) v
    integer (kind=4) d
    real    (kind=4) f

    equivalence (f, d)
    f = v
    as_int32 = d

end function
! ===========================================================================

! ===========================================================================
!
! @brief Cast double to int64
!
! ===========================================================================
integer (kind=8) function as_int64(v)

    implicit none
    real    (kind=8) v
    integer (kind=8) q
    real    (kind=8) d

    equivalence (d, q)
    d = v
    as_int64 = q

end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if float value is NaN
!
! ===========================================================================
integer (kind=4) function isnan_float(v)

    implicit none
    integer (kind=4) :: as_int32
    real    (kind=4) v
    integer (kind=4) iv
    real    (kind=4) res
    iv = as_int32(v)
    iv = iv .and. Z'7FFFFFFF'
    if (iv .gt. Z'7F800000') then
        res = 1
    else
        res = 0
    end if

    isnan_float = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if double value is NaN
!
! ===========================================================================
integer (kind=4) function isnan_double(v)

    implicit none
    integer (kind=8) :: as_int64
    real    (kind=8) v
    integer (kind=8) iv
    real    (kind=4) res
    iv = as_int64(v)
    iv = iv .and. Z'7FFFFFFFFFFFFFFF'
    if (iv .gt. Z'7FF0000000000000') then
        res = 1
    else
        res = 0
    end if

    isnan_double = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if float value is Inf
!
! ===========================================================================
integer (kind=4) function isinf_float(v)

    implicit none
    integer (kind=4) :: as_int32
    real    (kind=4) v
    integer (kind=4) iv
    real    (kind=4) res
    iv = as_int32(v)
    iv = iv .and. Z'7FFFFFFF'
    if (iv .eq. Z'7F800000') then
        res = 1
    else
        res = 0
    end if

    isinf_float = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if double value is Inf
!
! ===========================================================================
integer (kind=4) function isinf_double(v)

    implicit none
    integer (kind=8) :: as_int64
    real    (kind=8) v
    integer (kind=8) iv
    real    (kind=4) res
    iv = as_int64(v)
    iv = iv .and. Z'7FFFFFFFFFFFFFFF'
    if (iv .eq. Z'7FF0000000000000') then
        res = 1
    else
        res = 0
    end if

    isinf_double = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if float value is finite
!
! ===========================================================================
integer (kind=4) function isfinite_float(v)

    implicit none
    integer (kind=4) :: isnan_float
    integer (kind=4) :: isinf_float
    real    (kind=4) v
    real    (kind=4) res

    if ((isnan_float(v) .ne. 0) .or. (isinf_float(v) .ne. 0)) then
        res = 0
    else
        res = 1
    end if

    isfinite_float = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if double value is finite
!
! ===========================================================================
integer (kind=4) function isfinite_double(v)

    implicit none
    integer (kind=4) :: isnan_double
    integer (kind=4) :: isinf_double
    real    (kind=8) v
    real    (kind=4) res

    if ((isnan_double(v) .ne. 0) .or. (isinf_double(v) .ne. 0)) then
        res = 0
    else
        res = 1
    end if

    isfinite_double = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if float value is Zero
!
! ===========================================================================
integer (kind=4) function iszero_float(v)

    implicit none
    integer (kind=4) :: as_int32
    real    (kind=4) v
    integer (kind=4) iv
    real    (kind=4) res
    iv = as_int32(v)
    iv = iv .and. Z'7FFFFFFF'
    if (iv .eq. Z'00000000') then
        res = 1
    else
        res = 0
    end if

    iszero_float = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if double value is Zero
!
! ===========================================================================
integer (kind=4) function iszero_double(v)

    implicit none
    integer (kind=8) :: as_int64
    real    (kind=8) v
    integer (kind=8) iv
    real    (kind=4) res
    iv = as_int64(v)
    iv = iv .and. Z'7FFFFFFFFFFFFFFF'
    if (iv .eq. Z'0000000000000000') then
        res = 1
    else
        res = 0
    end if

    iszero_double = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if float value is subnormal
!
! ===========================================================================
integer (kind=4) function issubnormal_float(v)

    implicit none
    integer (kind=4) :: as_int32
    real    (kind=4) v
    integer (kind=4) iv
    real    (kind=4) res
    iv = as_int32(v)
    iv = iv .and. Z'7FFFFFFF'
    if ((iv .ge. Z'00000001') .and. (iv .le. Z'007FFFFF')) then
        res = 1
    else
        res = 0
    end if

    issubnormal_float = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Checks if double value is subnormal
!
! ===========================================================================
integer (kind=4) function issubnormal_double(v)

    implicit none
    integer (kind=8) :: as_int64
    real    (kind=8) v
    integer (kind=8) iv
    real    (kind=4) res
    iv = as_int64(v)
    iv = iv .and. Z'7FFFFFFFFFFFFFFF'
    if ((iv .ge. Z'0000000000000001').and.(iv .le. Z'000FFFFFFFFFFFFF')) then
        res = 1
    else
        res = 0
    end if

    issubnormal_double = res
    
end function
! ===========================================================================

! ===========================================================================
!
! @brief Classify float value
!
! ===========================================================================
integer (kind=4) function fpclassify_float(v)
    
    implicit none
    include "_vml_common_data.f90"
    integer (kind=4) :: isnan_float
    integer (kind=4) :: isinf_float
    integer (kind=4) :: iszero_float
    integer (kind=4) :: issubnormal_float
    real    (kind=4) v
    integer (kind=4) res
    
    if (isnan_float(v) .ne. 0) then
        res = FP_IEEE_NAN
    else if (isinf_float(v) .ne. 0) then
        res = FP_IEEE_INF
    else if (iszero_float(v) .ne. 0) then
        res = FP_IEEE_ZERO
    else if (issubnormal_float(v) .ne. 0) then
        res = FP_IEEE_SUBNORMAL
    else
        res = FP_IEEE_NORMAL
    end if
    
    fpclassify_float = res

end function
! ===========================================================================

! ===========================================================================
!
! @brief Classify double value
!
! ===========================================================================
integer (kind=4) function fpclassify_double(v)

    
    implicit none
    include "_vml_common_data.f90"
    integer (kind=4) :: isnan_double
    integer (kind=4) :: isinf_double
    integer (kind=4) :: iszero_double
    integer (kind=4) :: issubnormal_double
    real    (kind=8) v
    real    (kind=4) res
    
    if (isnan_double(v) .ne. 0) then
        res = FP_IEEE_NAN
    else if (isinf_double(v) .ne. 0) then
        res = FP_IEEE_INF
    else if (iszero_double(v) .ne. 0) then
        res = FP_IEEE_ZERO
    else if (issubnormal_double(v) .ne. 0) then
        res = FP_IEEE_SUBNORMAL
    else
        res = FP_IEEE_NORMAL
    end if
    
    fpclassify_double = res

end function
! ===========================================================================

! ===========================================================================
!
! @brief Compute ulp between two float values
!
! ===========================================================================
real (kind=4) function ulp_float (res, ref)

    
    implicit none
    include "_vml_common_data.f90"
    integer (kind=4) :: isfinite_float
    integer (kind=4) :: fpclassify_float
    real (kind=4) res
    real (kind=4) ref
    real (kind=4) resulp
    integer (kind=4) res_class
    integer (kind=4) ref_class

    if ((isfinite_float(res) .ne. 0) .and. (isfinite_float(ref) .ne. 0)) then
        resulp =  ABS(res-ref)/(2**((FLOOR(LOG(ABS(ref))/LOG(2.0)))-(24.0-1.0)))
    else
        res_class = fpclassify_float(res)
        ref_class = fpclassify_float(ref)
        if (res_class .eq. ref_class) then
            resulp = 0.0
        else
            resulp = HUGE(0.0)
        end if
    end if

    ulp_float = resulp

end function
! ===========================================================================

! ===========================================================================
!
! @brief Compute ulp between two double values
!
! ===========================================================================
real (kind=8) function ulp_double (res, ref)

    
    implicit none
    include "_vml_common_data.f90"
    integer (kind=4) :: isfinite_double
    integer (kind=4) :: fpclassify_double
    real (kind=8) res
    real (kind=8) ref
    real (kind=8) resulp
    integer (kind=4) res_class
    integer (kind=4) ref_class

    if ((isfinite_double(res) .ne. 0) .and. (isfinite_double(ref) .ne. 0)) then
        resulp =  ABS(res-ref)/(2**((FLOOR(LOG(ABS(ref))/LOG(2.0)))-(53.0-1.0)))
    else
        res_class = fpclassify_double(res)
        ref_class = fpclassify_double(ref)
        if (res_class .eq. ref_class) then
            resulp = 0.0
        else
            resulp = HUGE(0D0)
        end if
    end if

    ulp_double = resulp

end function
! ===========================================================================

! ===========================================================================
!
! @brief Check float result
!
! ===========================================================================
integer (kind=4) function check_result_float(argtype, arg1, arg2, & 
                        res1, res2, ref1, ref2, funcname, accuracy)

    
    implicit none
    include "_vml_common_data.f90"
    integer (kind=4) :: blend_int32
    real    (kind=4) :: ulp_float
    real    (kind=4) arg1, arg2, res1, res2, ref1, ref2
    integer (kind=4) argtype, accuracy, check_result
    character (len=*) funcname
    real    (kind=4) resulp
    logical isfail
    
    resulp = ulp_float(res1, ref1)
    if (argtype .eq. VML_ARG1_RES2) then
        resulp = MAX(resulp, ulp_float(res2, ref2))
    end if
    
    isfail = (resulp > FLOAT_MAXULP(accuracy))

    if (argtype .eq. VML_ARG2_RES1) then
        write (*, 661) TAB, TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, arg2, &
                   res1, ref1, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        661 format (A, A, A12, '{', A, '}(', g16.8, ',',g16.8,') = ', &
             g16.8,'; expected = ', g16.8, '; ulp = ', g11.3, ' (', A, ')')
    else if(argtype .eq. VML_ARG1_RES2) then
        write (*, 662) TAB, TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, & 
        res1, res2, ref1, ref2, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        662 format (A, A, A12, '{', A, '}(', g16.8, ') = ', g16.8,',',& 
        g16.8,'; expected = ', g16.8,',',g16.8,'; ulp = ', g11.3, ' (', A, ')')
    else
        write (*, 663) TAB, TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, res1, ref1, &
                                     resulp, TEST_RESULT(blend_int32(isfail,2,1))
        663 format (A, A, A12, '{', A, '}(', g16.8, ') = ', g16.8 &
              ,'; expected = ', g16.8, '; ulp = ', g11.3, ' (', A, ')')
    end if

    check_result_float = blend_int32(isfail,1,0)

end function
! ===========================================================================

! ===========================================================================
!
! @brief Check double result
!
! ===========================================================================
integer (kind=4) function check_result_double(argtype, arg1, arg2, & 
                         res1, res2, ref1, ref2, funcname, accuracy)

    
    implicit none
    include "_vml_common_data.f90"
    integer (kind=4) :: blend_int32
    real    (kind=8) :: ulp_double
    real    (kind=8) arg1, arg2, res1, res2, ref1, ref2
    integer (kind=8) iarg1, iarg2, ires1, ires2, iref1, iref2
    integer (kind=4) argtype, accuracy
    character (len=*) funcname
    real    (kind=8) resulp
    logical isfail

    resulp = ulp_double(res1, ref1)
    if (argtype .eq. VML_ARG1_RES2) then
        resulp = MAX(resulp, ulp_double(res2, ref2))
    end if
    
    isfail = (resulp > DOUBLE_MAXULP(accuracy))

    if (argtype .eq. VML_ARG2_RES1) then
        write (*, 771) TAB, TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, arg2, & 
                   res1, ref1, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        771 format (A, A, A12, '{', A, '}(', g24.16,',',g24.16, ') = ', & 
             g24.16,'; expected = ', g24.16, '; ulp = ', g11.3, ' (', A, ')')
    else if(argtype .eq. VML_ARG1_RES2) then
        write (*, 772) TAB, TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, res1, res2, & 
                         ref1, ref2, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        772 format (A, A, A12, '{', A, '}(', g24.16, ') = ', g24.16,',', &
        g24.16,'; expected = ', g24.16,',',g24.16, '; ulp = ', g11.3, ' (', A, ')')
    else
        write (*, 773) TAB, TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, res1, ref1, & 
                                     resulp, TEST_RESULT(blend_int32(isfail,2,1))
        773 format (A, A, A12, '{', A, '}(', g24.16, ') = ', & 
        g24.16,'; expected = ', g24.16, '; ulp = ', g11.3, ' (', A, ')')
    end if

    check_result_double = blend_int32(isfail,1,0)

end function
! ===========================================================================

! ===========================================================================
!
! @brief Check complex float result
!
! ===========================================================================
integer (kind=4) function check_result_float_complex(argtype, arg1, arg2, & 
                                res1, res2, ref1, ref2, funcname, accuracy)

    
    implicit none
    include "_vml_common_data.f90"
    integer (kind=4) :: blend_int32
    integer (kind=4) :: as_int32
    real    (kind=4) :: ulp_float
    complex (kind=4) arg1, arg2, res1, res2, ref1, ref2
    integer (kind=4) argtype, accuracy, check_result
    character (len=*) funcname
    real    (kind=4) resulp
    logical isfail
    
    resulp = ulp_float(real(res1), real(ref1))

    if (argtype .ne. VML_ARG1C_RES1R) then
        resulp = MAX(resulp,ulp_float(imag(res1), imag(ref1)))
    end if
    if (argtype .eq. VML_ARG1_RES2) then
        resulp = MAX(resulp, ulp_float(real(res2), real(ref2)))
        resulp = MAX(resulp, ulp_float(imag(res2), imag(ref2)))
    end if
    
    isfail = (resulp > COMPLEX_FLOAT_MAXULP(accuracy))

    if (argtype .eq. VML_ARG2_RES1) then
        write (*, 881) TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, arg2
        881 format (A, A12, '{', A, '}(', g16.8,SP,g16.8,'*i, ', g16.8,SP,g16.8,'*i)' )
        write (*, 882) TAB, TAB, res1
        882 format (A, A, '         = ', g16.8,SP,g16.8,'*i')
        write (*, 883) TAB, TAB, ref1
        883 format (A, A, 'expected = ', g16.8,SP,g16.8,'*i;')
        write (*, 884) TAB, TAB, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        884 format (A, A, '     ulp = ', g11.3, ' (', A, ')')
    else if(argtype .eq. VML_ARG1_RES2) then
        write (*, 885) TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1
        885 format (A, A12, '{', A, '}(', g16.8,SP,g16.8,'*i)' )
        write (*, 886) TAB, TAB, res1, res2
        886 format (A, A, '         = ', g16.8,SP,g16.8,'*i ,', g16.8,SP,g16.8,'*i')
        write (*, 887) TAB, TAB, ref1, ref2
        887 format (A, A, 'expected = ', g16.8,SP,g16.8,'*i ,', g16.8,SP,g16.8,'*i;')
        write (*, 888) TAB, TAB, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        888 format (A, A, 'ulp      = ', g11.3, ' (', A, ')')
    else
        write (*, 889) TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1
        889 format (A, A12, '{', A, '}(', g16.8,SP,g16.8,'*i)' )
        write (*, 890) TAB, TAB, res1
        890 format (A, A, '         = ', g16.8,SP,g16.8,'*i')
        write (*, 891) TAB, TAB, ref1
        891 format (A, A, 'expected = ', g16.8,SP,g16.8,'*i;')
        write (*, 892) TAB, TAB, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        892 format (A, A, '     ulp = ', g11.3, ' (', A, ')')
    end if
    write(* , '(A)', ADVANCE = 'YES') ''

    check_result_float_complex = blend_int32(isfail,1,0)

end function
! ===========================================================================

! ===========================================================================
!
! @brief Check complex double result
!
! ===========================================================================
integer (kind=4) function check_result_double_complex(argtype, arg1, arg2, & 
                                res1, res2, ref1, ref2, funcname, accuracy)

    
    implicit none
    include "_vml_common_data.f90"
    integer   (kind=4) :: blend_int32
    real    (kind=8) :: ulp_double
    complex (kind=8) arg1, arg2, res1, res2, ref1, ref2
    integer (kind=4) argtype, accuracy, check_result
    character (len=*) funcname
    real    (kind=8) resulp
    logical isfail
    
    resulp = ulp_double(real(res1), real(ref1))
    if (argtype .ne. VML_ARG1C_RES1R) then
        resulp = MAX(resulp,ulp_double(imag(res1), imag(ref1)))
    end if
    if (argtype .eq. VML_ARG1_RES2) then
        resulp = MAX(resulp, ulp_double(real(res2), real(ref2)))
        resulp = MAX(resulp, ulp_double(imag(res2), imag(ref2)))
    end if
    
    isfail = (resulp > COMPLEX_DOUBLE_MAXULP(accuracy))

    if (argtype .eq. VML_ARG2_RES1) then
        write (*, 881) TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1, arg2
        881 format (A, A12, '{', A, '}(', g24.16,SP,g24.16,'*i, ', g24.16,SP,g24.16,'*i)' )
        write (*, 882) TAB, TAB, res1
        882 format (A, A, '         = ', g24.16,SP,g24.16,'*i')
        write (*, 883) TAB, TAB, ref1
        883 format (A, A, 'expected = ', g24.16,SP,g24.16,'*i;')
        write (*, 884) TAB, TAB, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        884 format (A, A, 'ulp      = ', g11.3, ' (', A, ')')
    else if(argtype .eq. VML_ARG1_RES2) then
        write (*, 885) TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1
        885 format (A, A12, '{', A, '}(', g24.16,SP,g24.16,'*i)' )
        write (*, 886) TAB, TAB, res1, res2
        886 format (A, A, '         = ', g24.16,SP,g24.16,'*i ,', g24.16,SP,g24.16,'*i')
        write (*, 887) TAB, TAB, ref1, ref2
        887 format (A, A, 'expected = ', g24.16,SP,g24.16,'*i ,', g24.16,SP,g24.16,'*i;')
        write (*, 888) TAB, TAB, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        888 format (A, A, 'ulp      = ', g11.3, ' (', A, ')')
    else
        write (*, 889) TAB, funcname, VML_ACCURACY_NAME(accuracy), arg1
        889 format (A, A12, '{', A, '}(', g24.16,SP,g24.16,'*i)' )
        write (*, 890) TAB, TAB, res1
        890 format (A, A, '         = ', g24.16,SP,g24.16,'*i')
        write (*, 891) TAB, TAB, ref1
        891 format (A, A, 'expected = ', g24.16,SP,g24.16,'*i;')
        write (*, 892) TAB, TAB, resulp, TEST_RESULT(blend_int32(isfail,2,1))
        892 format (A, A, 'ulp      = ', g11.3, ' (', A, ')')
    end if
    write(*, '(A)', ADVANCE = 'YES') ''

    check_result_double_complex = blend_int32(isfail,1,0)

end function
! ===========================================================================
