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
!  Generate double precision data set using VSL Uniform RNG
!  with left bound a and right bound b
!

      function dGenerateUniformData( n, x, a, b )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: n
        real(kind=8),intent(out) :: x(n)
        real(kind=8),intent(in)  :: a
        real(kind=8),intent(in)  :: b

        integer(kind=4) :: errcode

        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_UNIFORM_STD

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

!       ***** Call RNG *****
        errcode = vdrnguniform( METHOD, stream, n, x, a, b )
        call CheckVslError( errcode )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )
        call CheckVslError( errcode )

        dGenerateUniformData = errcode
      end function

!
!  Generate single precision data set using VSL Uniform RNG
!  with left bound a and right bound b
!

      function sGenerateUniformData( n, x, a, b )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: n
        real(kind=4),intent(out) :: x(n)
        real(kind=4),intent(in)  :: a
        real(kind=4),intent(in)  :: b

        integer(kind=4) errcode

        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_UNIFORM_STD

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

!       ***** Call RNG *****
        errcode = vsrnguniform( METHOD, stream, n, x, a, b )
        call CheckVslError( errcode )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )
        call CheckVslError( errcode )

        sGenerateUniformData = errcode
      end function

!
!  Generate double precision data set using VSL Gaussian RNG
!  with mean a and standard deviation sigma
!
      function dGenerateGaussianData( dim, n, x, a, sigma )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)         :: dim
        integer,intent(in)         :: n
        real(kind=8),intent(out)   :: x(dim,n)
        real(kind=8),intent(in)    :: a
        real(kind=8),intent(in)    :: sigma

        integer(kind=4) :: errcode
        integer x_len

        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIAN_ICDF

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

!       ***** Call RNG *****
        x_len = dim * n
        errcode = vdrnggaussian( METHOD, stream, x_len, x, a, sigma )
        call CheckVslError( errcode )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )
        call CheckVslError( errcode )

        dGenerateGaussianData = errcode
      end function

!
!  Generate single precision data set using VSL Gaussian RNG
!  with mean a and standard deviation sigma
!
      function sGenerateGaussianData( dim, n, x, a, sigma )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)         :: dim
        integer,intent(in)         :: n
        real(kind=4),intent(out)   :: x(dim,n)
        real(kind=4),intent(in)    :: a
        real(kind=4),intent(in)    :: sigma

        integer :: x_len
        integer(kind=4) :: errcode

        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIAN_ICDF

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

!       ***** Call RNG *****
        x_len = dim * n
        errcode = vsrnggaussian( METHOD, stream, x_len, x, a, sigma )
        call CheckVslError( errcode )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )
        call CheckVslError( errcode )

        sGenerateGaussianData = errcode
      end function

!
!
!
      function dGenerateGaussianMVData( dim, n, x, a, C )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(out) :: x(dim,n)
        real(kind=8),intent(in)  :: a(dim)
        real(kind=8),intent(in)  :: C(dim,dim)

        integer :: x_len
        integer(kind=4) :: errcode

        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIANMV_ICDF

!       Following variables are used in Cholesky factorization routine
        integer      :: info
        real(kind=8) :: T(dim,dim)

        T = C
        call dpotrf('U', dim, T, dim, info)

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

!       ***** Call RNG *****
        x_len = dim * n
        errcode = vdrnggaussianmv( METHOD, stream, n, x, dim,                  &
     &                             VSL_MATRIX_STORAGE_FULL, a, T )
        call CheckVslError( errcode )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )
        call CheckVslError( errcode )

        dGenerateGaussianMVData = errcode
      end function

!
!
!
      function sGenerateGaussianMVData( dim, n, x, a, C )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(out) :: x(dim,n)
        real(kind=4),intent(in)  :: a(dim)
        real(kind=4),intent(in)  :: C(dim,dim)

        integer :: x_len
        integer(kind=4) :: errcode

        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIANMV_ICDF

!       Following variables are used in Cholesky factorization routine
        integer      :: info
        character    :: uplo
        real(kind=4) :: T(dim,dim)

        T = C

!       ***** Intel(R) MKL Choelsky factorization routine call *****
        uplo = 'U'
        call spotrf(uplo, dim, T, dim, info)

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

!       ***** Call RNG *****
        x_len = dim * n
        errcode = vsrnggaussianmv( METHOD, stream, n, x, dim,                  &
     &                             VSL_MATRIX_STORAGE_FULL, a, T )
        call CheckVslError( errcode )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )
        call CheckVslError( errcode )

        sGenerateGaussianMVData = errcode
      end function

!
!     Make Cholesky factorization of covariance matrix C and create new
!     stream
!
      function dInitGaussianMVDataGenerator( dim, C, T, stream )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: dim
        real(kind=8),intent(in)  :: C(dim,dim)
        real(kind=8),intent(out) :: T(dim,dim)
        TYPE(VSL_STREAM_STATE),intent(out) :: stream

        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH

        integer :: x_len
        integer(kind=4) :: errcode

!       Following variable is used in Cholesky factorization routine
        integer info
        character uplo

        T = C

!       ***** Intel(R) MKL Choelsky factorization routine call *****
        uplo = 'U'
        call dpotrf( uplo, dim, T, dim, info )

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

        dInitGaussianMVDataGenerator = errcode
      end function

!
!     Make Cholesky factorization of covariance matrix C and create new
!     stream
!
      function sInitGaussianMVDataGenerator( dim, C, T, stream )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: dim
        real(kind=4),intent(in)  :: C(dim,dim)
        real(kind=4),intent(out) :: T(dim,dim)
        TYPE(VSL_STREAM_STATE),intent(out) :: stream

        integer,parameter :: SEED   = 777
        integer,parameter :: BRNG   = VSL_BRNG_WH

        integer(kind=4) :: errcode

!       Following variable is used in Cholesky factorization routine
        integer info
        character uplo

        T = C

!       ***** Intel(R) MKL Choelsky factorization routine call *****
        uplo = 'U'
        call spotrf( uplo, dim, T, dim, info )

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )
        call CheckVslError( errcode )

        sInitGaussianMVDataGenerator = errcode
      end function

!
!     Generate double precision data block using VSL GaussianMV RNG
!
      function dGenerateGaussianMVDataBlock( dim, n, x, stream, a, T )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)                :: dim
        integer,intent(in)                :: n
        real(kind=8),intent(out)          :: x(dim,n)
        TYPE(VSL_STREAM_STATE),intent(in) :: stream
        real(kind=8),intent(in)           :: a(dim)
        real(kind=8),intent(in)           :: T(dim,dim)

        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIANMV_ICDF
        integer(kind=4) errcode

!       ***** Call RNG *****
        errcode = vdrnggaussianmv( METHOD, stream, n, x, dim,                  &
     &                             VSL_MATRIX_STORAGE_FULL, a, T )
        call CheckVslError( errcode )

        dGenerateGaussianMVDataBlock = errcode
      end function

!
!     Generate double precision data block using VSL GaussianMV RNG
!
      function sGenerateGaussianMVDataBlock( dim, n, x, stream, a, T )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)                :: dim
        integer,intent(in)                :: n
        real(kind=4),intent(out)          :: x(dim,n)
        TYPE(VSL_STREAM_STATE),intent(in) :: stream
        real(kind=4),intent(in)           :: a(dim)
        real(kind=4),intent(in)           :: T(dim,dim)

        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIANMV_ICDF
        integer(kind=4) :: errcode

!       ***** Call RNG *****
        errcode = vsrnggaussianmv( METHOD, stream, n, x, dim,                  &
     &                             VSL_MATRIX_STORAGE_FULL, a, T )
        call CheckVslError( errcode )

        sGenerateGaussianMVDataBlock = errcode
      end function

!
!     Generate double precision data set using GaussianMV RNG with two
!     different sets of parameters for different parts of the set
!
      function dGenerateContaminatedDataset( dim, n, x, a, C,                  &
     &                                       ratio, m, coeff )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in) :: dim
        integer,intent(in) :: n
        real(kind=8),intent(out)   :: x(n,dim)
        real(kind=8),intent(in)    :: a(dim)
        real(kind=8),intent(in)    :: C(dim,dim)
        integer,intent(in)         :: ratio
        real(kind=8),intent(in)    :: m
        real(kind=8),intent(in)    :: coeff

        integer i, j, mc
        integer(kind=4) errcode
        real(kind=8) mean(dim), T(dim,dim)
        real(kind=8),dimension(:,:),allocatable :: r1, r2

!       ***** Compute number of outliers *****
        mc = (n * ratio) / 100

!       ***** Allocate memory *****
        allocate(r1(dim, n-mc))
        allocate(r2(dim, mc))

        do i = 1, dim
          do j = 1, dim
            T(i, j) = C(i, j)
          end do
        end do

!       ***** Generate dataset from major distribution *****
        if ( n - mc /= 0 ) then
          errcode = dGenerateGaussianMVData( dim, n - mc, r1, a, T )
          call CheckVslError( errcode )
        end if

!       ***** Generate outliers *****
        if ( mc /= 0 ) then
          do i = 1, dim
            mean(i) = m
            do j = 1, dim
              T(i, j) = coeff * C(i, j)
            end do
          end do

          errcode = dGenerateGaussianMVData( dim, mc, r2, mean, T )
          call CheckVslError( errcode )
        end if

!       ***** Form final dataset *****
        do i = 1, n - mc
          do j = 1, dim
            x(i, j) = r1(j, i)
          end do
        end do

        do i = 1, mc
          do j = 1, dim
            x(i + n - mc, j) = r2(j, i)
          end do
        end do

        deallocate(r2)
        deallocate(r1)

        dGenerateContaminatedDataset = errcode
      end function

!
!     Generate single precision data set using GaussianMV RNG with two
!     different sets of parameters for different parts of the set
!
      function sGenerateContaminatedDataset( dim, n, x, a, C,                &
     &                                         ratio, m, coeff )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in) :: dim
        integer,intent(in) :: n
        real(kind=4),intent(out)   :: x(n,dim)
        real(kind=4),intent(in)    :: a(dim)
        real(kind=4),intent(in)    :: C(dim,dim)
        integer,intent(in)         :: ratio
        real(kind=4),intent(in)    :: m
        real(kind=4),intent(in)    :: coeff

        integer i, j, mc
        integer(kind=4) errcode
        real(kind=4) mean(dim), T(dim,dim)
        real(kind=4),dimension(:,:),allocatable :: r1, r2

!       ***** Compute number of outliers *****
        mc = (n * ratio) / 100

!       ***** Allocate memory *****
        allocate(r1(dim, n-mc))
        allocate(r2(dim, mc))

        do i = 1, dim
          do j = 1, dim
            T(i, j) = C(i, j)
          end do
        end do

!       ***** Generate dataset from major distribution *****
        if ( n - mc /= 0 ) then
          errcode = sGenerateGaussianMVData( dim, n - mc, r1, a, T )
          call CheckVslError( errcode )
        end if

!       ***** Generate outliers *****
        if ( mc /= 0 ) then
          do i = 1, dim
            mean(i) = m
            do j = 1, dim
              T(i, j) = coeff * C(i, j)
            end do
          end do

          errcode = sGenerateGaussianMVData( dim, mc, r2, mean, T )
          call CheckVslError( errcode )
        end if

!       ***** Form final dataset *****
        do i = 1, n - mc
          do j = 1, dim
            x(i, j) = r1(j, i)
          end do
        end do

        do i = 1, mc
          do j = 1, dim
            x(i + n - mc, j) = r2(j, i)
          end do
        end do

        deallocate(r2)
        deallocate(r1)

        sGenerateContaminatedDataset = errcode
      end function

!
!     Generate double precision outliered dataset using GaussianMV RNG.
!     Use GaussianMV RNG with other parameters.
!
      function dGenerateOutliersInput( dim, n, x, eps, a, a_outl,              &
     &                                 C, outlIndex, cntoutl )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)         :: dim
        integer,intent(in)         :: n
        real(kind=8),intent(out)   :: x(dim,n)
        real(kind=8),intent(in)    :: eps
        real(kind=8),intent(in)    :: a(dim)
        real(kind=8),intent(in)    :: a_outl(dim)
        real(kind=8),intent(in)    :: C(dim,dim)
        real(kind=8),intent(out)   :: outlIndex(n)
        integer,intent(out)        :: cntoutl

        integer(kind=4) errcode
        integer i, j, co
        real(kind=8) l, r
        real(kind=8) x_outl(dim,n), tmp_rand(n)

!       Generating random numbers from multivariate normal distribution
        errcode = dGenerateGaussianMVData( dim, n, x, a, C )
        call CheckVslError( errcode )

!       Generating random numbers from multivariate normal distribution
!       to form outliers
        errcode = dGenerateGaussianMVData( dim, n, x_outl, a_outl, C )
        call CheckVslError( errcode )

!       Generating random numbers from uniform distribution to form outliers
        l = 0.0
        r = 1.0
        errcode = dGenerateUniformData( n, tmp_rand, l, r )
        call CheckVslError( errcode )

        co = 0
        do i = 1, n
          if ( tmp_rand(i) < eps ) then
            co = co + 1
            outlIndex(i) = 0.0
            do j = 1, dim
              x(j, i) = x_outl(j, i)
            end do
          else
            outlIndex(i) = 1.0
          end if
        end do

        cntoutl = co

        dGenerateOutliersInput = errcode
      end function

!
!     Generate single precision outliered dataset using GaussianMV RNG.
!     Use GaussianMV RNG with other parameters.
!
      function sGenerateOutliersInput( dim, n, x, eps, a, a_outl,              &
     &                                 C, outlIndex, cntoutl )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)         :: dim
        integer,intent(in)         :: n
        real(kind=4),intent(out)   :: x(dim,n)
        real(kind=4),intent(in)    :: eps
        real(kind=4),intent(in)    :: a(dim)
        real(kind=4),intent(in)    :: a_outl(dim)
        real(kind=4),intent(in)    :: C(dim,dim)
        real(kind=4),intent(out)   :: outlIndex(n)
        integer,intent(out)        :: cntoutl

        integer(kind=4) errcode
        integer i, j, co
        real(kind=4) l, r
        real(kind=4) x_outl(dim,n), tmp_rand(n)

!       Generating random numbers from multivariate normal distribution
        errcode = sGenerateGaussianMVData( dim, n, x, a, C )
        call CheckVslError( errcode )

!       Generating random numbers from multivariate normal distribution
!       to form outliers
        errcode = sGenerateGaussianMVData( dim, n, x_outl, a_outl, C )
        call CheckVslError( errcode )

!       Generating random numbers from uniform distribution to form outliers
        l = 0.0
        r = 1.0
        errcode = sGenerateUniformData( n, tmp_rand, l, r )
        call CheckVslError( errcode )

        co = 0
        do i = 1, n
          if ( tmp_rand(i) < eps ) then
            co = co + 1
            outlIndex(i) = 0.0
            do j = 1, dim
              x(j, i) = x_outl(j, i)
            end do
          else
            outlIndex(i) = 1.0
          end if
        end do

        cntoutl = co

        sGenerateOutliersInput = errcode
      end function

!
!     Generate double precision data set with missing values using
!     GaussianMV RNG with mean a and covariance C.
!     Missing values are generated using patterns patt
!
      function dGenerateMissingValuesInput( dim, n, x, eps, npatt,             &
     &                                      patt, a, C,                        &
     &                                      miss_vals, nmiss_vals )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=8),intent(out) :: x(n,dim)
        real(kind=8),intent(in)  :: eps
        integer,intent(in)       :: npatt
        integer,intent(in)       :: patt(dim*npatt)
        real(kind=8),intent(in)  :: a(dim)
        real(kind=8),intent(in)  :: C(dim,dim)
        integer,intent(out)      :: miss_vals(n)
        integer,intent(out)      :: nmiss_vals

        real(kind=8),parameter   :: dNAN = VSL_SS_DNAN

        real(kind=8) av_weight, threshold
        real(kind=8) left, right
        integer nmv
        integer i, j, itmp, ipatt
        integer(kind=4) errcode
        real(kind=8) tmp_rand(n)
        real(kind=8) itmp_rand(n)
        real(kind=8) xt(dim,n)

!       Calculate average pattern weight
        av_weight = 0.0
        do i = 0, npatt - 1
          do j = 1, dim
            av_weight = av_weight + patt(i*dim + j)
          end do
        end do
        av_weight = av_weight / ( npatt * dim )

        threshold = eps / av_weight

!       Generating random numbers from multivariate normal distribution
        errcode = dGenerateGaussianMVData( dim, n, xt, a, C )
        call CheckVslError( errcode )

        do i = 1, n
          do j = 1, dim
            x(i, j) = xt(j, i)
          end do
        end do

!       Generating random numbers from uniform distribution to form
!       missing values
        left  = 0.0
        right = npatt
        errcode = dGenerateUniformData( n, tmp_rand, left, right )
        call CheckVslError( errcode )

        call vdfloor( n, tmp_rand, itmp_rand )

!       Integration of missing values into matrix of observations
        nmv = 0
        do i = 1, n
          if (tmp_rand(i) / npatt < threshold) then
!           Here if there are missing values in i-th observation
            itmp = i + 1
            if (itmp == n + 1) then
              itmp = 1
            end if

            ipatt = int(itmp_rand(itmp))
            ipatt = ipatt + 1
            miss_vals(i) = ipatt

            do j = 1, dim
              if (patt((ipatt - 1) * dim + j) == 1) then
                x(i, j) = dNAN
                nmv = nmv + 1
              end if
            end do
          end if
        end do

        nmiss_vals = nmv

        dGenerateMissingValuesInput = errcode
      end function

!
!     Generate single precision data set with missing values using
!     GaussianMV RNG with mean a and covariance C.
!     Missing values are generated using patterns patt
!
      function sGenerateMissingValuesInput( dim, n, x, eps, npatt,             &
     &                                      patt, a, C,                        &
     &                                      miss_vals, nmiss_vals )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer,intent(in)       :: dim
        integer,intent(in)       :: n
        real(kind=4),intent(out) :: x(n,dim)
        real(kind=4),intent(in)  :: eps
        integer,intent(in)       :: npatt
        integer,intent(in)       :: patt(dim*npatt)
        real(kind=4),intent(in)  :: a(dim)
        real(kind=4),intent(in)  :: C(dim,dim)
        integer,intent(out)      :: miss_vals(n)
        integer,intent(out)      :: nmiss_vals

        real(kind=4),parameter   :: sNAN = VSL_SS_SNAN

        real(kind=4) av_weight, threshold
        real(kind=4) left, right
        integer nmv
        integer i, j, itmp, ipatt
        integer(kind=4) errcode

        real(kind=4) tmp_rand(n)
        real(kind=4) itmp_rand(n)
        real(kind=4) xt(dim,n)

!       Calculate average pattern weight
        av_weight = 0.0
        do i = 0, npatt - 1
          do j = 1, dim
            av_weight = av_weight + patt(i*dim + j)
          end do
        end do
        av_weight = av_weight / ( npatt * dim )

        threshold = eps / av_weight

!       Generating random numbers from multivariate normal distribution
        errcode = sGenerateGaussianMVData( dim, n, xt, a, C )
        call CheckVslError( errcode )

        do i = 1, n
          do j = 1, dim
            x(i, j) = xt(j, i)
          end do
        end do

!       Generating random numbers from uniform distribution to form
!       missing values
        left  = 0.0
        right = npatt
        errcode = sGenerateUniformData( n, tmp_rand, left, right )
        call CheckVslError( errcode )

        call vsfloor( n, tmp_rand, itmp_rand )

!       Integration of missing values into matrix of observations
        nmv = 0
        do i = 1, n
          if (tmp_rand(i) / npatt < threshold) then
!           Here if there are missing values in i-th observation
            itmp = i + 1
            if (itmp == n + 1) then
              itmp = 1
            end if

            ipatt = int(itmp_rand(itmp))
            ipatt = ipatt + 1
            miss_vals(i) = ipatt

            do j = 1, dim
              if (patt((ipatt - 1) * dim + j) == 1) then
                x(i, j) = sNAN
                nmv = nmv + 1
              end if
            end do
          end if
        end do

        nmiss_vals = nmv

        sGenerateMissingValuesInput = errcode
      end function
