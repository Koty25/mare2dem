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
! An example of error processing when using DFTI functions.
!
!*****************************************************************************

program error_processing

  use MKL_DFTI

  ! Execution status
  integer :: status = 0

  ! DFTI descriptor handle
  type(DFTI_DESCRIPTOR), POINTER :: hand

  hand => null()

  print *,"Example error_processing"

  print *,"Try to create a DFTI descriptor with misplaced arguments"
  status = DftiCreateDescriptor(hand, DFTI_COMPLEX, DFTI_DOUBLE, 1, 2)
  if (0 /= status) call report_error(status)

  if (0 == status) then
    print *,"TEST FAILED"
    call exit(1)
  else
    print *,"TEST PASSED"
    call exit(0)
  end if

contains

  subroutine report_error(s)
    integer s

    print '(" Nonzero status = "I0)', s
    print '(" Check if the status indicates of error")'
    if (.not. DftiErrorClass( s, DFTI_NO_ERROR)) then
      print '("  Error: "A)', DftiErrorMessage(s)
    else
      print '("  Not an error: "A)', DftiErrorMessage(s)
    end if
  end subroutine report_error

end program error_processing
