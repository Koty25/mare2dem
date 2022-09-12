!===============================================================================
! Copyright 2011-2020 Intel Corporation.
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

! Content:
!      Example of printing DFTI descriptor's configuration
!
!*****************************************************************************

program config_dump_descriptor

  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Execution status
  integer :: status = 0, ignored_status

  ! DFTI descriptor handle
  type(DFTI_DESCRIPTOR), POINTER :: hand

  hand => null()

  print *,"Example config_dump_descriptor"

  print *,"Create a DFTI descriptor"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_REAL, 3, [5,6,7])
  if (0 /= status) goto 999

  print *,"Dump the descriptor"
  call dump_descriptor(hand)

100 continue

  print *,"Free the descriptor"
  ignored_status = DftiFreeDescriptor(hand)

  if (status == 0) then
    print *, "TEST PASSED"
    call exit(0)
  else
    print *, "TEST FAILED"
    call exit(1)
  end if

999 continue
  print '("  Error, status = ",I0)', status
  goto 100

contains

  subroutine dump_descriptor(hand)
    integer, parameter :: MR = 10 ! maximum rank
    type(DFTI_DESCRIPTOR), POINTER :: hand
    integer :: prec, domain, rank, lengths(MR)
    integer :: placement
    real*8  :: fwd_scale, bwd_scale
    integer :: nut, istrides(1+MR), ostrides(1+MR), ntr, idist, odist
    integer :: storage, packfmt, wspace, cmtstatus

    ! Execution status
    integer :: status

    print '("  PRECISION              = "$)'
    status = DftiGetValue(hand, DFTI_PRECISION, prec)
    if (0 /= status) goto 999

    if      (prec == DFTI_SINGLE) then ; print '("DFTI_SINGLE")'
    else if (prec == DFTI_DOUBLE) then ; print '("DFTI_DOUBLE")'
    else                               ; goto 999
    end if

    print '("  FORWARD_DOMAIN         = "$)'
    status = DftiGetValue(hand, DFTI_FORWARD_DOMAIN, domain)
    if (0 /= status) goto 999

    if      (domain == DFTI_COMPLEX) then ; print '("DFTI_COMPLEX")'
    else if (domain == DFTI_REAL)    then ; print '("DFTI_REAL")'
    else                                  ; goto 999
    end if

    print '("  DIMENSION              = "$)'
    status = DftiGetValue(hand, DFTI_DIMENSION, rank)
    if (0 /= status) goto 999
    print '(I0)', rank

    print '("  LENGTHS                = "$)'
    status = DftiGetValue(hand, DFTI_LENGTHS, lengths)
    if (0 /= status) goto 999
    print '(10(I0:", "))', lengths(1:rank)

    print '("  PLACEMENT              = "$)'
    status = DftiGetValue(hand, DFTI_PLACEMENT, placement)
    if (0 /= status) goto 999

    if      (placement == DFTI_INPLACE)     then ; print '("DFTI_INPLACE")'
    else if (placement == DFTI_NOT_INPLACE) then ; print '("DFTI_NOT_INPLACE")'
    else                                         ; goto 999
    end if

    print '("  F/B SCALES             = "$)'
    status = DftiGetValue(hand, DFTI_FORWARD_SCALE,  fwd_scale)
    if (0 /= status) goto 999
    status = DftiGetValue(hand, DFTI_BACKWARD_SCALE, bwd_scale)
    if (0 /= status) goto 999
    print '(2(G11.4:", "))', fwd_scale, bwd_scale

    print '("  NO OF USER THREADS     = "$)'
    status = DftiGetValue(hand, DFTI_NUMBER_OF_USER_THREADS, nut)
    if (0 /= status) goto 999
    print '(I0)', nut

    print '("  INPUT STRIDES          = "$)'
    status = DftiGetValue(hand, DFTI_INPUT_STRIDES,  istrides)
    if (0 /= status) goto 999
    print '(I0"; "10(I0:", "))', istrides(1:1+rank)

    print '("  OUTPUT STRIDES         = "$)'
    status = DftiGetValue(hand, DFTI_OUTPUT_STRIDES, ostrides)
    if (0 /= status) goto 999
    print '(I0"; "10(I0:", "))', ostrides(1:1+rank)

    print '("  NO OF TRANSFORMS       = "$)'
    status = DftiGetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, ntr)
    if (0 /= status) goto 999
    print '(I0)', ntr

    print '("  I/O DISTANCES          = "$)'
    status = DftiGetValue(hand, DFTI_INPUT_DISTANCE,  idist)
    if (0 /= status) goto 999
    status = DftiGetValue(hand, DFTI_OUTPUT_DISTANCE, odist)
    if (0 /= status) goto 999
    print '(2(I0:", "))', idist, odist

    if (domain == DFTI_COMPLEX) then
      print '("  COMPLEX STORAGE        = "$)'
      status = DftiGetValue(hand, DFTI_COMPLEX_STORAGE, storage)
      if (0 /= status) goto 999
      if (storage == DFTI_COMPLEX_COMPLEX) then
        print '("DFTI_COMPLEX_COMPLEX")'
      else if (storage == DFTI_REAL_REAL) then
        print '("DFTI_REAL_REAL")'
      else
        goto 999
      end if
    else if (domain == DFTI_REAL) then
      print '("  CONJUGATE EVEN STORAGE = "$)'
      status = DftiGetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, storage)
      if (0 /= status) goto 999
      if (storage == DFTI_COMPLEX_COMPLEX) then
        print '("DFTI_COMPLEX_COMPLEX")'
      else if (storage == DFTI_COMPLEX_REAL) then
        print '("DFTI_COMPLEX_REAL")'
      else
        goto 999
      end if
      if (storage == DFTI_COMPLEX_REAL) then
        print '("    PACKED FORMAT        = "$)'
        status = DftiGetValue(hand, DFTI_PACKED_FORMAT, packfmt)
        if (0 /= status) goto 999
        if (packfmt == DFTI_CCS_FORMAT) then
          print '("DFTI_CCS_FORMAT")'
        else if (packfmt == DFTI_PACK_FORMAT) then
          print '("DFTI_PACK_FORMAT")'
        else if (packfmt == DFTI_PERM_FORMAT) then
          print '("DFTI_PERM_FORMAT")'
        else
          goto 999
        end if
      end if
    end if

    print '("  WORKSPACE              = "$)'
    status = DftiGetValue(hand, DFTI_WORKSPACE, wspace)
    if (0 /= status) goto 999

    if      (wspace == DFTI_ALLOW) then ; print '("DFTI_ALLOW")'
    else if (wspace == DFTI_NONE ) then ; print '("DFTI_AVOID")'
    else                                ; goto 999
    end if

    print '("  COMMIT STATUS          = "$)'
    status = DftiGetValue(hand, DFTI_COMMIT_STATUS, cmtstatus)
    if (0 /= status) goto 999

    if (cmtstatus == DFTI_COMMITTED) then
      print '("DFTI_COMMITTED")'
    else if (cmtstatus == DFTI_UNCOMMITTED) then
      print '("DFTI_UNCOMMITTED")'
    else
      goto 999
    end if

    return

999 print *," Error, status = ", status
    call exit(1)

  end subroutine dump_descriptor

end program config_dump_descriptor
