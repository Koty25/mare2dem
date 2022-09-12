!===============================================================================
! Copyright 2003-2020 Intel Corporation.
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

      SUBROUTINE CheckVslError( num )

        USE MKL_VSL_TYPE
        INTEGER(KIND=4) :: num

      if ( num == VSL_STATUS_OK ) then
        return
      end if

      if ( num == VSL_ERROR_CPU_NOT_SUPPORTED ) then
        print 33, "Error: CPU version is not supported (code ",          &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_FEATURE_NOT_IMPLEMENTED ) then
        print 33, "Error: this feature not implemented yet (code ",      &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_UNKNOWN ) then
        print 33, "Error: unknown error (code ",num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_BADARGS ) then
        print 33, "Error: bad arguments (code ",num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_MEM_FAILURE ) then
        print 34, "Error: memory failure. Memory allocation problem",    &
     &            " maybe (code ",num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_NULL_PTR ) then
        print 33, "Error: null pointer (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_INVALID_BRNG_INDEX ) then
        print 33, "Error: invalid BRNG index (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED ) then
        print 34, "Error: leapfrog initialization is unsupported",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED ) then
        print 34, "Error: skipahead initialization is unsupported",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_SKIPAHEADEX_UNSUPPORTED ) then
        print 34, "Error: skipahead advanced initialization",              &
     &            "is unsupported (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BRNGS_INCOMPATIBLE ) then
        print 34, "Error: BRNGs are not compatible for the operation",   &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_STREAM ) then
        print 33, "Error: random stream is invalid (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BRNG_TABLE_FULL ) then
        print 33, "Error: table of registered BRNGs is full (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE ) then
        print 34, "Error: value in StreamStateSize field is bad",        &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_WORD_SIZE ) then
        print 33, "Error: value in WordSize field is bad (code ",        &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_NSEEDS ) then
        print 33, "Error: value in NSeeds field is bad (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_NBITS ) then
        print 33, "Error: value in NBits field is bad (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_UPDATE ) then
        print 34, "Error: number of updated entries in buffer is",       &
     &            " invalid (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_NO_NUMBERS ) then
        print 34, "Error: zero number of updated entries in buffer",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM ) then
        print 33, "Error: abstract random stream is invalid (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_NONDETERM_NOT_SUPPORTED ) then
        print 35, "Error: Non-deterministic random number generator is", &
     &            " not supported on the CPU which runs the",            &
     &            " application (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_NONDETERM_NRETRIES_EXCEEDED ) then
        print 35, "Error: Number of retries to generate a random",       &
     &            " number by using non-deterministic random number",    &
     &            " generator exceeds threshold (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_ARS5_NOT_SUPPORTED ) then
        print 35, "Error: ARS-5 random number generator is",             &
     &            " not supported on the CPU which runs the",            &
     &            " application (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_CLOSE ) then
        print 33, "Error: can`t close file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_OPEN ) then
        print 33, "Error: can`t open file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_WRITE ) then
        print 33, "Error: can`t write to file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_READ ) then
        print 33, "Error: can`t read from file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_FILE_FORMAT ) then
        print 33, "Error: file format is unknown (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_UNSUPPORTED_FILE_VER ) then
        print 33, "Error: unsupported file version (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_MEM_FORMAT ) then
        print 33, "Error: stream format is unknown (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_ALLOCATION_FAILURE ) then
        print 33, "Error: memory allocation failure (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_DIMEN ) then
        print 33, "Error: bad dimension value (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OBSERV_N ) then
        print 33, "Error: bad number of observations (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_STORAGE_NOT_SUPPORTED ) then
        print 33, "Error: storage format is not supported (code ",       &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_INDC_ADDR ) then
        print 33, "Error: array of indices is not defined (code ",       &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_WEIGHTS ) then
        print 34, "Error: array of weights contains negative values",    &
     &            "(code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MEAN_ADDR ) then
        print 33, "Error: array of means is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_2R_MOM_ADDR ) then
        print 34, "Error: array of 2nd order raw moments is not",        &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_3R_MOM_ADDR ) then
        print 34, "Error: array of 3rd order raw moments is not",        &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_4R_MOM_ADDR ) then
        print 34, "Error: array of 4th order raw moments is not",        &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_2C_MOM_ADDR ) then
        print 34, "Error: array of 2nd order central moments is not",    &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_3C_MOM_ADDR ) then
        print 34, "Error: array of 3rd order central moments is not",    &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_4C_MOM_ADDR ) then
        print 34, "Error: array of 4th order central moments is not",    &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_KURTOSIS_ADDR ) then
        print 34, "Error: array of kurtosis values is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_SKEWNESS_ADDR ) then
        print 34, "Error: array of skewness values is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MIN_ADDR ) then
        print 34, "Error: array of minimum values is not defined",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MAX_ADDR ) then
        print 34, "Error: array of maximum values is not defined",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_VARIATION_ADDR ) then
        print 34, "Error: array of variation coefficients is not",       &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_COV_ADDR ) then
        print 33, "Error: covariance matrix is not defined (code ",      &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_COR_ADDR ) then
        print 33, "Error: correlation matrix is not defined (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_QUANT_ORDER_ADDR ) then
        print 34, "Error: array of quantile orders is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_QUANT_ORDER ) then
        print 33, "Error: bad value of quantile order (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_QUANT_ADDR ) then
        print 33, "Error: array of quantiles is not defined (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ORDER_STATS_ADDR ) then
        print 34, "Error: array of order statistics is not defined",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_MOMORDER_NOT_SUPPORTED ) then
        print 34, "Error: moment of requested order is not supported",   &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_NOT_FULL_RANK_MATRIX ) then
        print 34, "Warning: correlation matrix is not of full rank",     &
     &            " (code ",num,")."
      end if

      if ( num == VSL_SS_ERROR_ALL_OBSERVS_OUTLIERS ) then
        print 33, "Error: all observations are outliers (code ",         &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_COV_ADDR ) then
        print 34, "Error: robust covariance matrix is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_MEAN_ADDR ) then
        print 34, "Error: array of robust means is not defined",         &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_METHOD_NOT_SUPPORTED ) then
        print 33, "Error: requested method is not supported (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_NULL_TASK_DESCRIPTOR ) then
        print 33, "Error: task descriptor is null (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OBSERV_ADDR ) then
        print 33, "Error: dataset matrix is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_SINGULAR_COV ) then
        print 33, "Error: covariance matrix is singular (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_POOLED_COV_ADDR ) then
        print 34, "Error: pooled covariance matrix is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_POOLED_MEAN_ADDR ) then
        print 33, "Error: array of pooled means is not defined (code ",  &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_COV_ADDR ) then
        print 34, "Error: group covariance matrix is not defined",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_MEAN_ADDR ) then
        print 33, "Error: array of group means is not defined (code ",   &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_INDC_ADDR ) then
        print 33, "Error: array of group indices is not defined (code ", &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_INDC ) then
        print 33, "Error: group indices have improper values (code ",    &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_ADDR ) then
        print 34, "Error: array of parameters for outliers detection",   &
     &            " algorithm is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_N_ADDR ) then
        print 35, "Error: pointer to size of parameter array for",       &
     &            " outlier detection algorithm is not defined",         &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OUTLIERS_WEIGHTS_ADDR ) then
        print 34, "Error: output of the outlier detection algorithm",    &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_ADDR ) then
        print 34, "Error: array of parameters of robust covariance",     &
     &            " estimation algorithm is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_N_ADDR ) then
        print 35, "Error: pointer to number of parameters of",           &
     &            " algorithm for robust covariance is not defined",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STORAGE_ADDR ) then
        print 34, "Error: pointer to variable that holds storage",       &
     &            " format is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COV_IDX_ADDR ) then
        print 35, "Error: array that encodes sub-components of",         &
     &            " random vector for partial covariance algorithm",     &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COV_ADDR ) then
        print 34, "Error: partial covariance matrix is not defined",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COR_ADDR ) then
        print 34, "Error: partial correlation matrix is not defined",    &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PARAMS_ADDR ) then
        print 34, "Error: array of parameters for Multiple Imputation",  &
     &            " method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PARAMS_N_ADDR ) then
        print 35, "Error: pointer to number of parameters for",          &
     &            " Multiple Imputation method is not defined",          &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_BAD_PARAMS_N ) then
        print 34, "Error: bad size of the parameter array of Multiple",  &
     &            " Imputation method (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PARAMS ) then
        print 34, "Error: bad parameters of Multiple Imputation",        &
     &            " method (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_N_ADDR ) then
        print 35, "Error: pointer to number of initial estimates in",    &
     &            " Multiple Imputation method is not defined",          &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_ADDR ) then
        print 34, "Error: array of initial estimates for Multiple",      &
     &            " Imputation method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_SIMUL_VALS_ADDR ) then
        print 35, "Error: array of simulated missing values in",         &
     &            " Multiple Imputation method is not defined",          &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N_ADDR ) then
        print 35, "Error: pointer to size of the array of simulated",    &
     &            " missing values in Multiple Imputation method",       &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_ESTIMATES_N_ADDR ) then
        print 35, "Error: pointer to the number of parameter",           &
     &            " estimates in Multiple Imputation method is not",     &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_ESTIMATES_ADDR ) then
        print 34, "Error: array of parameter estimates in Multiple",     &
     &            " Imputation method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N ) then
        print 34, "Error: bad size of the array of simulated values",    &
     &            " in Multiple Imputation method (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_OUTPUT_PARAMS ) then
        print 34, "Error: array of output parameters in Multiple",       &
     &            " Imputation method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PRIOR_N_ADDR ) then
        print 34, "Error: pointer to the number of prior parameters",    &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PRIOR_ADDR ) then
        print 33, "Error: bad number of missing values (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PRIOR_ADDR ) then
        print 33, "Error: bad number of missing values (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_SEMIDEFINITE_COR ) then
        print 34, "Warning: correlation matrix passed into",             &
     &            " parametrization function is semidefinite (code ",    &
     &             num,")."
      end if

      if ( num == VSL_SS_ERROR_BAD_PARAMTR_COR_ADDR ) then
        print 34, "Error: correlation matrix to be parametrized is",     &
     &            " not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_COR ) then
        print 34, "Error: all eigenvalues of correlation matrix to",     &
     &            " be parametrized are non-positive (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N_ADDR ) then
        print 35, "Error: pointer to the number of parameters for",      &
     &            " quantile computation algorithm for streaming data",  &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_ADDR ) then
        print 35, "Error: array of parameters of quantile computation",  &
     &            " algorithm for streaming data is not defined",        &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N ) then
        print 35, "Error: bad number of parameters of quantile",         &
     &            " computation algorithm for streaming data",           &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS ) then
        print 34, "Error: bad parameters of quantile computation",       &
     &            " algorithm for streaming data (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER_ADDR ) then
        print 34, "Error: array of quantile orders for streaming",       &
     &            " data is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER ) then
        print 33, "Error: bad quantile order for streaming data (code ", &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_ADDR ) then
        print 34, "Error: array of quantiles for streaming data",        &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COV_IDX ) then
        print 34, "Error: partial covariance indices have improper",     &
     &            " values (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_SUM_ADDR ) then
        print 34, "Error: array of sums",                                &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_2R_SUM_ADDR ) then
        print 34, "Error: array of 2nd order raw sums",                  &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_3R_SUM_ADDR ) then
        print 34, "Error: array of 3d order raw sums",                   &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_4R_SUM_ADDR ) then
        print 34, "Error: array of 4th order raw sums",                  &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_2C_SUM_ADDR ) then
        print 34, "Error: array of 2nd order central sums",              &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_3C_SUM_ADDR ) then
        print 34, "Error: array of 3d order central sums",               &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_4C_SUM_ADDR ) then
        print 34, "Error: array of 4th order central sums",              &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_CP_ADDR ) then
        print 34, "Error: cross-product matrix",                         &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MDAD_ADDR ) then
        print 34, "Error: array of median absoltute deviations",         &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MNAD_ADDR ) then
        print 34, "Error: array of mean absolute deviations",            &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_SORTED_OBSERV_ADDR ) then
        print 34, "Error: array for sorting",                            &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_INDICES_NOT_SUPPORTED ) then
        print 34, "Error: array of indices is not supported with",       &
     &            "  in-place sorting (code ",num,")."
        stop 1
      end if

33    format(A,I5,A)
34    format(A,A,I5,A)
35    format(A,A,A,I5,A)

      END SUBROUTINE
