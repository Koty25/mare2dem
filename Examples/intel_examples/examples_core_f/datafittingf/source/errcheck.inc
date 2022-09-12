!===============================================================================
! Copyright 2010-2020 Intel Corporation.
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

      SUBROUTINE CheckDfError( num )

        USE MKL_DF_TYPE
        INTEGER(KIND=4) :: num

      if ( num == DF_STATUS_OK ) then
        return
      end if

      if ( num == DF_ERROR_NULL_TASK_DESCRIPTOR ) then
        print 33, "Error: null task descriptor (code ", num,")."
        stop 1
      end if

      if ( num == DF_ERROR_MEM_FAILURE ) then
        print 33, "Error: memory allocation failure in DF functionality  &
     &             (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_NX ) then
        print 33, "Error: the number of breakpoints is invalid (code ",  &
     &             num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_X ) then
        print 33, "Error: array of breakpoints is invalid (code ",       &
     &             num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_X_HINT ) then
        print 33, "Error: invalid flag describing structure of partition &
     &             (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_NY ) then
        print 33, "Error: invalid dimension of vector-valued function y  &
     &             (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_Y ) then
        print 33, "Error: array of function values is invalid (code ",   &
     &             num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_Y_HINT ) then
        print 33, "Error: invalid flag describing structure of function  &
     &             y (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_SPLINE_ORDER ) then
        print 33, "Error: invalid spline order (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_SPLINE_TYPE ) then
        print 33, "Error: invalid spline type (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_IC_TYPE ) then
        print 33, "Error: invalid type of internal conditions used       &
     &             in the spline construction (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_IC ) then
        print 33, "Error: array of internal conditions for spline        &
     &             construction is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_BC_TYPE ) then
        print 33, "Error: invalid type of boundary conditions used       &
     &             in the spline construction (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_BC ) then
        print 33, "Error: array of boundary conditions for spline        &
     &             construction is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_PP_COEFF ) then
        print 33, "Error: array of piece-wise polynomial spline          &
     &             coefficients is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_PP_COEFF_HINT ) then
        print 33, "Error: invalid flag describing structure of the       &
     &             piece-wise polynomial spline coefficients (code ",    &
     &             num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_PERIODIC_VAL ) then
        print 33, "Error: function values at the end points of the       &
     &             interpolation interval are not equal as required      &
     &             in periodic boundary conditions (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_DATA_ATTR ) then
        print 33, "Error: invalid attribute of the pointer to be set or  &
     &             modified in Data Fitting task descriptor with         &
     &             EditIdxPtr editor (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_DATA_IDX ) then
        print 33, "Error: index of pointer to be set or modified in Data &
     &             Fitting task descriptor with EditIdxPtr editor is     &
     &             out of range (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_NSITE ) then
        print 33, "Error: invalid number of interpolation sites (code ", &
     &             num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_SITE ) then
        print 33, "Error: array of interpolation sites is not defined    &
     &             (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_SITE_HINT ) then
        print 33, "Error: invalid flag describing structure of           &
     &             interpolation sites (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_NDORDER ) then
        print 33, "Error: invalid size of array that defines order       &
     &             of the derivatives to be computed at the              &
     &             interpolation sites (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_DORDER ) then
        print 33, "Error: array defining derivative orders to be         &
     &             computed at interpolation sites is not defined        &
     &             (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_DATA_HINT ) then
        print 33, "Error: invalid flag providing a-priori information    &
     &             about partition and/or interpolation sites (code ",   &
     &             num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_INTERP ) then
        print 33, "Error: array of spline based interpolation results    &
     &             is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_INTERP_HINT ) then
        print 33, "Error: invalid flag defining structure of spline      &
     &             based interpolation results (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_CELL_IDX ) then
        print 33, "Error: array of indices of partition cells containing &
     &             interpolation sites is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_NLIM ) then
        print 33, "Error: invalid size of arrays containing integration  &
     &             limits (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_LLIM ) then
        print 33, "Error: array of left integration limits               &
     &             is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_RLIM ) then
        print 33, "Error: array of right integration limits              &
     &             is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_INTEGR ) then
        print 33, "Error: array of spline based integration results      &
     &             is not defined (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_INTEGR_HINT ) then
        print 33, "Error: invalid flag defining structure of spline      &
     &             based inetgration results (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_BAD_LOOKUP_INTERP_SITE ) then
        print 33, "Error: bad site provided for interpolation with       &
     &             look-up interpolator (code ",num,")."
        stop 1
      end if

      if ( num == DF_ERROR_NULL_PTR ) then
        print 33, "Error: bad pointer provided in DF function (code ",   &
     &             num,")."
      end if

33    format(A,I5,A)

      END SUBROUTINE
