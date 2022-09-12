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

!
!     Computes P-values for double precision mean estimate
!

      subroutine dComputePvalsMean( dim, n, mean, a, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: mean(dim)
        real(kind=8),intent(in)  :: a(dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim)

        real(kind=8) s, sqrts
        real(kind=8) d(dim)
        integer i

        do i = 1, dim
          s = C(i, i) / n
          sqrts = SQRT(s)
          d(i) = ABS(mean(i) - a(i)) / sqrts
        end do

        call vdcdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * (1.0 - pvals(i))
        end do

      end subroutine

!
!     Computes P-values for single precision mean estimate
!
      subroutine sComputePvalsMean( dim, n, mean, a, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: mean(dim)
        real(kind=4),intent(in)  :: a(dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim)

        real(kind=4) s, sqrts
        real(kind=4) d(dim)
        integer i

        do i = 1, dim
          s = C(i, i) / n
          sqrts = sqrt(s)
          d(i) = abs(mean(i) - a(i)) / sqrts
        end do

        call vscdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * (1.0 - pvals(i))
        end do

      end subroutine

!
!     Compute P-value for double precision variance estimate
!
      function dPvalVariance( n, var, d )
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: var
        real(kind=8),intent(in)  :: d

        real(kind=8) s(1)
        real(kind=8) pval(1)
        real(kind=8) sqrts

        s(1) = 2.0 / (n - 1)
        sqrts = SQRT(s(1))
        s(1) = ABS(var - d)/(d * sqrts)
        call vdcdfnorm( 1, s, pval )

        dPvalVariance = 2.0 * (1.0 - pval(1))
      end function

!
!     Compute P-value for single precision variance estimate
!
      function sPvalVariance( n, var, d )
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: var
        real(kind=4),intent(in)  :: d

        real(kind=4) s(1)
        real(kind=4) pval(1)
        real(kind=4) sqrts

        s(1) = 2.0 / (n - 1)
        sqrts = SQRT(s(1))
        s(1) = ABS(var - d)/(d * sqrts)
        call vscdfnorm( 1, s, pval )

        sPvalVariance = 2.0 * (1.0 - pval(1))
      end function

!
!     Compute P-values for double precision variance estimates
!
      subroutine dComputePvalsVariance( dim, n, cov, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: cov(dim,dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim,dim)

        integer i

        do i = 1, dim
          pvals(i,i) = dPvalVariance( n, cov(i,i), C(i,i) )
        end do

        return
      end subroutine

!
!     Compute P-values for single precision variance estimates
!
      subroutine sComputePvalsVariance( dim, n, cov, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: cov(dim,dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim,dim)

        integer i

        do i = 1, dim
          pvals(i,i) = sPvalVariance( n, cov(i,i), C(i,i) )
        end do

        return
      end subroutine

!
!     Compute P-value for double precision covariance estimate
!
      function dPvalCovariance( n, C, c_ii, c_ij, c_jj )
        integer,intent(in)      :: n
        real(kind=8),intent(in) :: C
        real(kind=8),intent(in) :: c_ii
        real(kind=8),intent(in) :: c_ij
        real(kind=8),intent(in) :: c_jj

        real(kind=8) s(1)
        real(kind=8) pval(1)
        real(kind=8) sqrts

        s(1) = (n - 1) / ((c_ii * c_jj) + (c_ij * c_ij))
        sqrts = SQRT(s(1))
        s(1) = ABS(C - c_ij) * sqrts
        call vdcdfnorm( 1, s, pval )

        dPvalCovariance = 2.0 * (1.0 - pval(1))
      end function

!
!     Compute P-value for single precision covariance estimate
!
      function sPvalCovariance( n, C, c_ii, c_ij, c_jj )
        integer,intent(in)      :: n
        real(kind=4),intent(in) :: C
        real(kind=4),intent(in) :: c_ii
        real(kind=4),intent(in) :: c_ij
        real(kind=4),intent(in) :: c_jj

        real(kind=4) s(1)
        real(kind=4) pval(1)
        real(kind=4) sqrts

        s(1) = (n - 1) / ((c_ii * c_jj) + (c_ij * c_ij))
        sqrts = SQRT(s(1))
        s(1) = ABS(C - c_ij) * sqrts
        call vscdfnorm( 1, s, pval )

        sPvalCovariance = 2.0 * (1.0 - pval(1))
      end function

!
!     Compute p-values for double precision covariance matrix estimate
!
      subroutine dComputePvalsCovariance( dim, n, cov, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: cov(dim,dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim,dim)

        integer i, j

        do i = 1, dim
          do j = 1, i-1
            pvals(i,j) = dPvalCovariance( n, cov(i,j), C(i,i), C(i,j),   &
     &                                    C(j,j) )
            pvals(j,i) = pvals(i,j)
          end do
        end do
      end subroutine

!
!     Compute p-values for single precision covariance matrix estimate
!
      subroutine sComputePvalsCovariance( dim, n, cov, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: cov(dim,dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim,dim)

        integer i, j

        do i = 1, dim
          do j = 1, i-1
            pvals(i,j) = sPvalCovariance( n, cov(i,j), C(i,i), C(i,j),   &
     &                                    C(j,j) )
            pvals(j,i) = pvals(i,j)
          end do
        end do
      end subroutine

!
!     Implementation of the routine for computation of double precision
!     raw moments estimates
!
      function dComputeRawMoment( order, mean, sigma )
        integer,intent(in)      :: order
        real(kind=8),intent(in) :: mean
        real(kind=8),intent(in) :: sigma

        real(kind=8) mu
        real(kind=8) mean2, sigma2, meansigma2

        if (order == 1) then
          mu = mean
        end if
        if (order == 2) then
          mu = mean**2 + sigma**2
        end if
        if (order == 3) then
          mu = mean * ( mean**2 + 3.0 * ( sigma**2 ) )
        end if
        if (order == 4) then
          mu = mean**4 + 6.0 * ( mean**2 ) * ( sigma**2 ) +                    &
     &         3.0 * (sigma**4)
        end if
        if (order == 6) then
          mean2  = mean**2
          sigma2 = sigma**2
          mu = mean**6 +                                                       &
     &         15.0 * mean2 * sigma2 * ( mean2 + 3.0 * sigma2 ) +              &
     &         15.0 * sigma**6
        end if
        if (order == 8) then
          meansigma2 = (mean * sigma)**2
          mu = mean**8 + 14.0 * meansigma2 * ( 2.0 * ( mean**4 ) +             &
     &         15.0 * meansigma2 + 30.0 * ( sigma**4 ) ) +                     &
     &         105.0 * sigma**8
        end if

        dComputeRawMoment = mu
      end function

!
!     Implementation of the routine for computation of single precision
!     raw moments estimates
!
      function sComputeRawMoment( order, mean, sigma )
        integer,intent(in)      :: order
        real(kind=4),intent(in) :: mean
        real(kind=4),intent(in) :: sigma

        real(kind=4) mu
        real(kind=4) mean2, sigma2, meansigma2

        if (order == 1) then
          mu = mean
        end if
        if (order == 2) then
          mu = mean**2 + sigma**2
        end if
        if (order == 3) then
          mu = mean * ( mean**2 + 3.0 * ( sigma**2 ) )
        end if
        if (order == 4) then
          mu = mean**4 + 6.0 * ( mean**2 ) * ( sigma**2 ) +                    &
     &         3.0 * (sigma**4)
        end if
        if (order == 6) then
          mean2  = mean**2
          sigma2 = sigma**2
          mu = mean**6 +                                                       &
     &         15.0 * mean2 * sigma2 * ( mean2 + 3.0 * sigma2 ) +              &
     &         15.0 * sigma**6
        end if
        if (order == 8) then
          meansigma2 = (mean * sigma)**2
          mu = mean**8 + 14.0 * meansigma2 * ( 2.0 * ( mean**4 ) +             &
     &         15.0 * meansigma2 + 30.0 * ( sigma**4 ) ) +                     &
     &         105.0 * sigma**8
        end if

        sComputeRawMoment = mu
      end function

!
!     Compute p-values for double precision raw moments estimates
!
      subroutine dComputePvalsRawMoments( dim, n, raw, order, a, C,            &
     &                                    pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: raw(dim)
        integer,intent(in)       :: order
        real(kind=8),intent(in)  :: a(dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim)

        integer i, order2
        real(kind=8) mean, sigma, Es, Es2, V, pval
        real(kind=8) d(dim), pvalues(dim)

        do i = 1, dim
          mean = a(i)
          sigma = sqrt( C(i, i) )
          order2 = 2 * order
          Es = dComputeRawMoment( order, mean, sigma )
          Es2 = dComputeRawMoment( order2, mean, sigma )
          V = ( Es2 - Es**2 ) / n
          d(i) = abs( raw(i) - Es ) / sqrt(V)
        end do

        call vdcdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for single precision raw moments estimates
!
      subroutine sComputePvalsRawMoments( dim, n, raw, order, a, C,            &
     &                                    pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: raw(dim)
        integer,intent(in)       :: order
        real(kind=4),intent(in)  :: a(dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim)

        integer i, order2
        real(kind=4) mean, sigma, Es, Es2, V, pval
        real(kind=4) d(dim), pvalues(dim)

        do i = 1, dim
          mean = a(i)
          sigma = sqrt( C(i, i) )
          order2 = 2 * order
          Es = sComputeRawMoment( order, mean, sigma )
          Es2 = sComputeRawMoment( order2, mean, sigma )
          V = ( Es2 - Es**2 ) / n
          d(i) = abs( raw(i) - Es ) / sqrt(V)
        end do

        call vscdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!
!
      function oddFactorial( n )
        integer,intent(in) :: n

        integer res, prev

        res = n
        prev = res - 2
        do while ( prev > 0 )
          res = res * prev
          prev = prev - 2
        end do

        oddFactorial = res
      end function

!
!     Implementation of the routine for computation of double precision
!     central moments estimates
!
      function dComputeCentralMoment( order, V )
        integer,intent(in)      :: order
        real(kind=8),intent(in) :: V

        integer factorial, order1
        real(kind=8) mu

        mu = 0.0

        if ( mod( order, 2 ) == 0 ) then
          order1 = order - 1
          factorial = oddFactorial( order1 )
          mu = factorial * ( V**( order / 2 ) )
        end if

        dComputeCentralMoment = mu
      end function

!
!     Implementation of the routine for computation of single precision
!     central moments estimates
!
      function sComputeCentralMoment( order, V )
        integer,intent(in)      :: order
        real(kind=4),intent(in) :: V

        integer factorial, order1
        real(kind=4) mu

        mu = 0.0

        if ( mod( order, 2 ) == 0 ) then
          order1 = order - 1
          factorial = oddFactorial( order1 )
          mu = factorial * ( V**( order / 2 ) )
        end if

        sComputeCentralMoment = mu
      end function

!
!     Compute p-values for double precision central moments estimates
!
      subroutine dComputePvalsCentralMoments( dim, n, central, order,          &
     &                                        a, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: central(dim)
        integer,intent(in)       :: order
        real(kind=8),intent(in)  :: a(dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim)

        integer i, order2, orderp1, orderm1
        real(kind=8) V, mu_2, mu_k, mu_kp1, mu_km1, mu_2k
        real(kind=8) d(dim)

        order2  = 2 * order
        orderp1 = order + 1
        orderm1 = order - 1

        do i = 1, dim
          mu_2 = C(i, i)
          mu_2k  = dComputeCentralMoment( order2,  mu_2 )
          mu_kp1 = dComputeCentralMoment( orderp1, mu_2 )
          mu_km1 = dComputeCentralMoment( orderm1, mu_2 )
          mu_k   = dComputeCentralMoment( order,   mu_2 )
          V = ( mu_2k - 2.0 * mu_kp1 * mu_km1 - mu_k * mu_k +                  &
     &        order * order * mu_2 * mu_km1 ) / n
          d(i) = abs( central(i) - mu_k ) / sqrt( V )
        end do

        call vdcdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for single precision central moments estimates
!
      subroutine sComputePvalsCentralMoments( dim, n, central, order,          &
     &                                        a, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: central(dim)
        integer,intent(in)       :: order
        real(kind=4),intent(in)  :: a(dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim)

        integer i, order2, orderp1, orderm1
        real(kind=4) V, mu_2, mu_k, mu_kp1, mu_km1, mu_2k
        real(kind=4) d(dim)

        order2  = 2 * order
        orderp1 = order + 1
        orderm1 = order - 1

        do i = 1, dim
          mu_2 = C(i, i)
          mu_2k  = sComputeCentralMoment( order2,  mu_2 )
          mu_kp1 = sComputeCentralMoment( orderp1, mu_2 )
          mu_km1 = sComputeCentralMoment( orderm1, mu_2 )
          mu_k   = sComputeCentralMoment( order,   mu_2 )
          V = ( mu_2k - 2.0 * mu_kp1 * mu_km1 - mu_k * mu_k +                  &
     &        order * order * mu_2 * mu_km1 ) / n
          d(i) = abs( central(i) - mu_k ) / sqrt( V )
        end do

        call vscdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for double precision kurtosis estimates
!
      subroutine dComputePvalsKurtosis( dim, n, kurt, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: kurt(dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim)

        integer factorial
        real(kind=8) sigma2, Es, V
        real(kind=8) mu2, mu4, mu6, mu8
        real(kind=8) d(dim)

        do i = 1, dim
          sigma2 = C(i, i)
          Es = dComputeCentralMoment( 4, sigma2 ) / ( sigma2**2 ) - 3.0

          mu2 = dComputeCentralMoment( 2, sigma2 )
          mu4 = dComputeCentralMoment( 4, sigma2 )
          mu6 = dComputeCentralMoment( 6, sigma2 )
          mu8 = dComputeCentralMoment( 8, sigma2 )

          V = ( ( mu2**2 ) * mu8 - 4.0 * mu2 * mu4 * mu6 +                     &
     &          4.0 * mu4 * ( mu4**2 ) - ( mu2**2 ) * ( mu4**2 ) ) /           &
     &          ( n * ( mu2**6 ) )

          d(i) = abs( kurt(i) - Es ) / sqrt( V )
        end do

        call vdcdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for single precision kurtosis estimates
!
      subroutine sComputePvalsKurtosis( dim, n, kurt, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: kurt(dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim)

        integer factorial
        real(kind=4) sigma2, Es, V
        real(kind=4) mu2, mu4, mu6, mu8
        real(kind=4) d(dim)

        do i = 1, dim
          sigma2 = C(i, i)
          Es = sComputeCentralMoment( 4, sigma2 ) / ( sigma2**2 ) - 3.0

          mu2 = sComputeCentralMoment( 2, sigma2 )
          mu4 = sComputeCentralMoment( 4, sigma2 )
          mu6 = sComputeCentralMoment( 6, sigma2 )
          mu8 = sComputeCentralMoment( 8, sigma2 )

          V = ( ( mu2**2 ) * mu8 - 4.0 * mu2 * mu4 * mu6 +                     &
     &          4.0 * mu4 * ( mu4**2 ) - ( mu2**2 ) * ( mu4**2 ) ) /           &
     &          ( n * ( mu2**6 ) )

          d(i) = abs( kurt(i) - Es ) / sqrt( V )
        end do

        call vscdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for double precision skewness estimates
!
      subroutine dComputePvalsSkewness( dim, n, skew, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: skew(dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim)

        integer factorial
        real(kind=8) sigma2, Es, V
        real(kind=8) mu2, mu4, mu6
        real(kind=8) d(dim)

        do i = 1, dim
          sigma2 = C(i, i)
          mu2 = dComputeCentralMoment( 2, sigma2 )
          mu4 = dComputeCentralMoment( 4, sigma2 )
          mu6 = dComputeCentralMoment( 6, sigma2 )

          V = ( ( mu2**2 ) * mu6 - 6.0 * ( mu2**3 ) * mu4 +                    &
     &          9.0 * ( mu2**5 ) ) / ( n * ( mu2**5 ) )

          d(i) = abs( skew(i) ) / sqrt( V )
        end do

        call vdcdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for single precision skewness estimates
!
      subroutine sComputePvalsSkewness( dim, n, skew, C, pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: skew(dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim)

        integer factorial
        real(kind=4) sigma2, Es, V
        real(kind=4) mu2, mu4, mu6
        real(kind=4) d(dim)

        do i = 1, dim
          sigma2 = C(i, i)
          mu2 = sComputeCentralMoment( 2, sigma2 )
          mu4 = sComputeCentralMoment( 4, sigma2 )
          mu6 = sComputeCentralMoment( 6, sigma2 )

          V = ( ( mu2**2 ) * mu6 - 6.0 * ( mu2**3 ) * mu4 +                    &
     &          9.0 * ( mu2**5 ) ) / ( n * ( mu2**5 ) )

          d(i) = abs( skew(i) ) / sqrt( V )
        end do

        call vscdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for double precision variation estimates
!
      subroutine dComputePvalsVariation( dim, n, variation, a, C,              &
     &                                   pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(in)  :: variation(dim)
        real(kind=8),intent(in)  :: a(dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: pvals(dim)

        integer i
        real(kind=8) s, sigma, V
        real(kind=8) d(dim)

        do i = 1, dim
          s = 2.0 / ( n - 1 )
          V = sqrt( C(i, i) ) / a(i)
          d(i) = abs( variation(i) - V ) / ( V * sqrt( s ) )
        end do

        call vdcdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Compute p-values for single precision variation estimates
!
      subroutine sComputePvalsVariation( dim, n, variation, a, C,              &
     &                                   pvals )
        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(in)  :: variation(dim)
        real(kind=4),intent(in)  :: a(dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: pvals(dim)

        integer i
        real(kind=4) s, sigma, V
        real(kind=4) d(dim)

        do i = 1, dim
          s = 2.0 / ( n - 1 )
          V = sqrt( C(i, i) ) / a(i)
          d(i) = abs( variation(i) - V ) / ( V * sqrt( s ) )
        end do

        call vscdfnorm( dim, d, pvals )

        do i = 1, dim
          pvals(i) = 2.0 * ( 1.0 - pvals(i) )
        end do

        return
      end subroutine

!
!     Implementation of double precision sweep operatior
!
!     Reference: J. L. Schafer, Analysis of Incomplete Multivariate Data, 1999
!
      subroutine dSweep( k, dim, g )
        integer,intent(in)       :: k
        integer,intent(in)       :: dim
        real(kind=8),intent(out) :: g(dim,dim)

        integer i, j
        real(kind=8) h

        g(k, k) = -1.0 / g(k, k)

        h = g(k, k)
        do i = 1, k - 1
          g(i, k) = -g(i, k) * h
          g(k, i) = g(i, k)
        end do

        do i = k + 1, dim
          g(i, k) = -g(i, k) * h
          g(k, i) = g(i, k)
        end do

        do i = 1, dim
          if (i /= k) then
            do j = 1, i - 1
              if (j /= k) then
                g(i, j) = g(i, j) + g(i, k) * g(k, j) / h
                g(j, i) = g(i, j)
              end if
            end do
            g(i, i) = g(i, i) + g(i, k) * g(k, i) / h
          end if
        end do
      end subroutine

!
!     Implementation of single precision sweep operatior
!
!     Reference: J. L. Schafer, Analysis of Incomplete Multivariate Data, 1999
!
      subroutine sSweep( k, dim, g )
        integer,intent(in)       :: k
        integer,intent(in)       :: dim
        real(kind=4),intent(out) :: g(dim,dim)

        integer i, j
        real(kind=4) h

        g(k, k) = -1.0 / g(k, k)

        h = g(k, k)
        do i = 1, k - 1
          g(i, k) = -g(i, k) * h
          g(k, i) = g(i, k)
        end do

        do i = k + 1, dim
          g(i, k) = -g(i, k) * h
          g(k, i) = g(i, k)
        end do

        do i = 1, dim
          if (i /= k) then
            do j = 1, i - 1
              if (j /= k) then
                g(i, j) = g(i, j) + g(i, k) * g(k, j) / h
                g(j, i) = g(i, j)
              end if
            end do
            g(i, i) = g(i, i) + g(i, k) * g(k, i) / h
          end if
        end do
      end subroutine
