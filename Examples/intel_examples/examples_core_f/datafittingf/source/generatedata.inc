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

      include 'mkl_vsl.fi'


      MODULE DF_GENERATE_DATA

      REAL(4),PARAMETER ::                                               &
     & DF_GEN_SPI = 3.14159265358979311599
      REAL(8),PARAMETER ::                                               &
     & DF_GEN_DPI = 3.14159265358979311599

        INTERFACE
          INTEGER FUNCTION sUniformData(x,left,right,n)
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: left
            REAL(4),INTENT(IN) :: right
            REAL(4),INTENT(OUT) :: x(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dUniformData(x,left,right,n)
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: left
            REAL(8),INTENT(IN) :: right
            REAL(8),INTENT(OUT) :: x(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sQuasiUniformData(x,left,right,n)
            INTEGER,INTENT(IN)  :: n
            REAL(4),INTENT(IN)  :: left
            REAL(4),INTENT(IN)  :: right
            REAL(4),INTENT(OUT) :: x(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dQuasiUniformData(x,left,right,n)
            INTEGER,INTENT(IN)  :: n
            REAL(8),INTENT(IN)  :: left
            REAL(8),INTENT(IN)  :: right
            REAL(8),INTENT(OUT) :: x(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sSinDataUniformGrid(f,freq,left,right,n)
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: freq,left,right,f(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dSinDataUniformGrid(f,freq,left,right,n)
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: freq,left,right,f(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sSinDerDataNotUniformGrid(f,x,freq,n)
            INTEGER,INTENT(IN)   :: n
            REAL(4),INTENT(IN)   :: freq
            REAL(4),INTENT(IN)   :: f(*)
            REAL(4),INTENT(OUT)  :: x(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dSinDerDataNotUniformGrid(f,x,freq,n)
            INTEGER,INTENT(IN)   :: n
            REAL(8),INTENT(IN)   :: freq
            REAL(8),INTENT(IN)   :: f(*)
            REAL(8),INTENT(OUT)  :: x(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sPeakData( x, left, right, n )
            USE MKL_VSL_TYPE
            USE MKL_VSL
            USE MKL_DF_TYPE
            USE MKL_DF

            real(4),intent(out) :: x(n)
            real(4),intent(in)  :: left
            real(4),intent(in)  :: right
            integer,intent(in)  :: n
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dPeakData( x, left, right, n )
            USE MKL_VSL_TYPE
            USE MKL_VSL
            USE MKL_DF_TYPE
            USE MKL_DF

            real(8),intent(out) :: x(n)
            real(8),intent(in)  :: left
            real(8),intent(in)  :: right
            integer,intent(in)  :: n
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sSinDer2DataNotUniformGrid(f,x,freq,n)
            REAL(4),INTENT(OUT) :: f(*)
            REAL(4),INTENT(IN) :: x(*)
            REAL(4),INTENT(IN) :: freq
            INTEGER,INTENT(IN) :: n
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dSinDer2DataNotUniformGrid(f,x,freq,n)
            REAL(8),INTENT(OUT) :: f(*)
            REAL(8),INTENT(IN) :: x(*)
            REAL(8),INTENT(IN) :: freq
            INTEGER,INTENT(IN) :: n
          END FUNCTION
        END INTERFACE

      END MODULE DF_GENERATE_DATA


!
! "Bubble" sort for single precision data
!
      function sSort( n, a )
        USE MKL_DF_TYPE
        USE MKL_DF

        integer,intent(in)          :: n
        real(kind=4),intent(out)    :: a(n)

        integer :: errcode
        integer :: i
        integer :: j
        real(kind=4) :: tmp

        errcode = 0

        do i = 1, n
          do j = i+1, n
            if ( a(j) < a(i) ) then
              tmp  = a(j)
              a(j) = a(i)
              a(i) = tmp
            end if
          end do
        end do

        sSort = errcode
      end function

!
! "Bubble" sort for double precision data
!
      function dSort( n, a )
        USE MKL_DF_TYPE
        USE MKL_DF

        integer,intent(in)          :: n
        real(kind=8),intent(out)    :: a(n)

        integer :: errcode
        integer :: i
        integer :: j
        real(kind=8) :: tmp

        errcode = 0

        do i = 1, n
          do j = i+1, n
            if ( a(j) < a(i) ) then
              tmp  = a(j)
              a(j) = a(i)
              a(i) = tmp
            end if
          end do
        end do

        dSort = errcode
      end function

!
! Generate single precision array of functions with following
! dependency: f(x) = Sin( 2 * Pi * freq * x ),
! the grid is not uniform
!
      function sSinDataNotUniformGrid( f, x, freq, n)
        USE DF_GENERATE_DATA
        USE MKL_DF_TYPE
        USE MKL_DF

        real(kind=4),intent(out)    :: f(n)
        real(kind=4),intent(in)     :: x(n)
        real(kind=4),intent(in)     :: freq
        integer,intent(in)          :: n

        integer :: errcode
        integer :: i
        real(kind=4) :: freq2

        errcode = 0

        freq2     = 2.0 * freq

        do i = 1, n
          f(i) = SIN( freq2 * DF_GEN_SPI * x(i) )
        end do

        sSinDataNotUniformGrid = errcode
      end function

!
! Generate double precision array of functions with following
! dependency: f(x) = Sin( 2 * Pi * freq * x ),
! the grid is not uniform
!
      function dSinDataNotUniformGrid( f, x, freq, n)
        USE DF_GENERATE_DATA
        USE MKL_DF_TYPE
        USE MKL_DF

        real(kind=8),intent(out)    :: f(n)
        real(kind=8),intent(in)     :: x(n)
        real(kind=8),intent(in)     :: freq
        integer,intent(in)          :: n

        integer :: errcode
        integer :: i
        real(kind=8) :: freq2

        errcode = 0

        freq2 = 2.0 * freq

        do i = 1, n
          f(i) = DSIN( freq2 * DF_GEN_DPI * x(i) )
        end do

        dSinDataNotUniformGrid = errcode
      end function

! Generate single precision array of dependent variables with following
! dependency: f(x) = 2 * Pi * freq * Cos( 2 * Pi * freq * x ),
! the grid is not uniform

      FUNCTION sSinDerDataNotUniformGrid(f,x,freq,n)  RESULT (output_4)
        USE DF_GENERATE_DATA, ignore_me => sSinDerDataNotUniformGrid

        INTEGER :: output_4
        INTEGER :: n
        REAL(4) :: freq,f(*),x(*)

        INTEGER :: i
        REAL(4) :: mult


        mult = 2.0 * DF_GEN_SPI * freq
        DO i = 0,n-1
          f(i+1) = mult * COS( mult * x(i+1) )
        END DO
        output_4 = 0
        RETURN



      END FUNCTION

! Generate double precision array of dependent variables with following
! dependency: f(x) = 2 * Pi * freq * Cos( 2 * Pi * freq * x ),
! the grid is not uniform

      FUNCTION dSinDerDataNotUniformGrid(f,x,freq,n)  RESULT (output_4)
        USE DF_GENERATE_DATA, ignore_me => dSinDerDataNotUniformGrid

        INTEGER :: output_4
        INTEGER :: n
        REAL(8) :: freq,f(*),x(*)

        INTEGER :: i
        REAL(8) :: mult

        mult = 2.0 * DF_GEN_DPI * freq
        DO i = 0,n-1
          f(i+1) = mult * DCOS( mult * x(i+1) )
        END DO
        output_4 = 0
        RETURN

      END FUNCTION

! Generate double precision array of dependent variables with following
! dependency: f(x) = -4 * Pi * Pi * freq * freq * Sin( 2 * Pi * freq * x ),
! the grid is not uniform

      FUNCTION sSinDer2DataNotUniformGrid(f,x,freq,n) RESULT (errcode)
        USE DF_GENERATE_DATA, ignore_me => sSinDer2DataNotUniformGrid

        INTEGER :: errcode
        INTEGER :: n
        REAL(4) :: freq,f(*),x(*)

        INTEGER :: i
        REAL(4) :: mult

        mult = 2.0 * DF_GEN_SPI * freq
        DO i = 0,n-1
          f(i+1) = -mult * mult * SIN( mult * x(i+1) )
        END DO
        errcode = 0
        RETURN

      END FUNCTION

! Generate single precision array of dependent variables with following
! dependency: f(x) = Sin( 2 * Pi * freq * x ),
! the grid is uniform

      FUNCTION dSinDer2DataNotUniformGrid(f,x,freq,n) RESULT (errcode)
        USE DF_GENERATE_DATA, ignore_me => dSinDer2DataNotUniformGrid

        INTEGER :: errcode
        INTEGER :: n
        REAL(8) :: freq,f(*),x(*)

        INTEGER :: i
        REAL(8) :: mult

        mult = 2.0 * DF_GEN_DPI * freq
        DO i = 0,n-1
          f(i+1) = -mult * mult * DSIN( mult * x(i+1) )
        END DO
        errcode = 0
        RETURN

      END FUNCTION

! Generate double precision array of dependent variables with following
! dependency: f(x) = Sin( 2 * Pi * freq * x ),
! the grid is uniform

      FUNCTION sSinDataUniformGrid(f,freq,left,right,n)                  &
     & RESULT (output_4)
        USE DF_GENERATE_DATA, ignore_me => sSinDataUniformGrid

        INTEGER :: output_4
        INTEGER :: n
        REAL(4) :: freq,left,right,f(*)

        INTEGER :: i
        REAL(4) :: delta,cur_x

        delta = (right - left) / (FLOAT(n)  - 1)
        cur_x = left
        DO i = 0,n-1
          f(i+1) = SIN( 2.0 * DF_GEN_SPI * freq * cur_x )
          cur_x = cur_x + delta
        END DO
        output_4 = 0
        RETURN

      END FUNCTION

! Generate double precision array of dependent variables with following
! dependency: f(x) = Sin( 2 * Pi * freq * x ),
! the grid is uniform

      FUNCTION dSinDataUniformGrid(f,freq,left,right,n)                  &
     & RESULT (output_4)
        USE DF_GENERATE_DATA, ignore_me => dSinDataUniformGrid
        INTEGER :: output_4
        INTEGER :: n
        REAL(8) :: tmp,freq,left,right,f(*)

        INTEGER :: i
        REAL(8) :: delta,cur_x

        delta = (right - left) / (DBLE(n)  - 1)
        cur_x = left
        DO i = 0,n-1
          f(i+1) = DSIN(2.0 * DF_GEN_DPI * freq * cur_x)
          cur_x = cur_x + delta
        END DO

        output_4 = 0
        RETURN

      END FUNCTION

! Generate single precision array of independent variables with uniform grid

      FUNCTION sUniformData(x,left,right,n)  RESULT (errcode)

        INTEGER :: errcode
        INTEGER :: n
        REAL(4) :: left,right,x(*)

        INTEGER :: i
        REAL(4) :: delta

        delta = (right - left) / (FLOAT(n)  - 1)

        DO i = 1, n
          x(i) = left + FLOAT(i-1)  * delta
        END DO
        errcode = 0
        RETURN

      END FUNCTION

! Generate single precision array of independent variables with uniform grid

      FUNCTION dUniformData(x,left,right,n)  RESULT (errcode)

        INTEGER :: errcode
        INTEGER :: n
        REAL(8) :: left,right,x(*)

        INTEGER :: i
        REAL(8) :: delta

        delta = (right - left) / (DBLE(n)  - 1)

        DO i = 1, n
          x(i) = left + DBLE(i-1)  * delta
        END DO
        errcode = 0
        RETURN

      END FUNCTION

! Generate single precision array of independent variables
! with quasi-uniform grid

      FUNCTION sQuasiUniformData(x,left,right,n)  RESULT (output_4)
        USE DF_GENERATE_DATA, ignore_me => sQuasiUniformData

        INTEGER :: output_4
        INTEGER :: n
        REAL(4) :: left,right,x(*)

        INTEGER :: errcode
        REAL(4) :: eps

        errcode = 0
        errcode = sUniformData(x, left, right, n)
        eps = ( x(0+1) - x(1+1) ) / 5.0
        x(1+1) = x(1+1) + eps
        output_4 = 0
        RETURN

      END FUNCTION

! Generate single precision array of independent variables
! with quasi-uniform grid

      FUNCTION dQuasiUniformData(x,left,right,n)  RESULT (output_4)
        USE DF_GENERATE_DATA, ignore_me => dQuasiUniformData

        INTEGER :: output_4
        INTEGER :: n
        REAL(8) :: left,right,x(*)

        INTEGER :: errcode
        REAL(8) :: eps

        errcode = 0
        errcode = dUniformData(x, left, right, n)
        eps = ( x(0+1) - x(1+1) ) / 5.0
        x(1+1) = x(1+1) + eps
        output_4 = 0
        RETURN

      END FUNCTION

!
! Generate single precision array of uniformly distributed variables
!
      function sUniformRandData( x, left, right, n )
        USE MKL_VSL_TYPE
        USE MKL_VSL
        USE MKL_DF_TYPE
        USE MKL_DF

        real(kind=4),intent(out)    :: x(n)
        real(kind=4),intent(in)     :: left
        real(kind=4),intent(in)     :: right
        integer,intent(in)          :: n

        integer :: errcode
        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 7777777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_UNIFORM_STD

        errcode = 0

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )

!       ***** Call RNG *****
        errcode = vsrnguniform( method, stream, n, x, left, right )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )

        sUniformRandData = errcode
      end function

!
! Generate double precision array of uniformly distributed variables
!
      function dUniformRandData( x, left, right, n )
        USE MKL_VSL_TYPE
        USE MKL_VSL
        USE MKL_DF_TYPE
        USE MKL_DF

        real(kind=8),intent(out)    :: x(n)
        real(kind=8),intent(in)     :: left
        real(kind=8),intent(in)     :: right
        integer,intent(in)          :: n

        integer :: errcode
        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 7777777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_UNIFORM_STD

        errcode = 0

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )

!       ***** Call RNG *****
        errcode = vdrnguniform( method, stream, n, x, left, right )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )

        dUniformRandData = errcode
      end function

!
! Generate single precision sorted array of uniformly distributed variables
!
      function sUniformRandSortedData( x, left, right, n )
        USE MKL_DF_TYPE
        USE MKL_DF

        real(kind=4),intent(out)    :: x(n)
        real(kind=4),intent(in)     :: left
        real(kind=4),intent(in)     :: right
        integer,intent(in)          :: n

        integer :: errcode

        errcode = 0

        errcode = sUniformRandData( x, left, right, n )

        errcode = sSort( n, x )

        sUniformRandSortedData = errcode
      end function

!
! Generate double precision sorted array of uniformly distributed variables
!
      function dUniformRandSortedData( x, left, right, n )
        USE MKL_DF_TYPE
        USE MKL_DF

        real(kind=8),intent(out)    :: x(n)
        real(kind=8),intent(in)     :: left
        real(kind=8),intent(in)     :: right
        integer,intent(in)          :: n

        integer :: errcode

        errcode = 0

        errcode = dUniformRandData( x, left, right, n )

        errcode = dSort( n, x )

        dUniformRandSortedData = errcode
      end function


!
! Generate single precision array of gauss distributed variables
!
      function sPeakData( x, left, right, n ) RESULT(errcode)
        USE DF_GENERATE_DATA, ignore_me => sPeakData
        USE MKL_VSL_TYPE
        USE MKL_VSL

        real(kind=4),intent(out)    :: x(n)
        real(kind=4),intent(in)     :: left
        real(kind=4),intent(in)     :: right
        integer,intent(in)          :: n

        integer :: errcode
        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 7777777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIAN_ICDF

        errcode = 0

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )

!       ***** Call RNG *****
        errcode = vsrnggaussian( method, stream, n, x, left, right )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )

        RETURN
      end function


!
! Generate double precision array of gaussian distributed variables
!
      function dPeakData( x, left, right, n ) RESULT(errcode)
        USE DF_GENERATE_DATA, ignore_me => dPeakData
        USE MKL_VSL_TYPE
        USE MKL_VSL

        real(kind=8),intent(out)    :: x(n)
        real(kind=8),intent(in)     :: left
        real(kind=8),intent(in)     :: right
        integer,intent(in)          :: n

        integer :: errcode
        TYPE(VSL_STREAM_STATE) :: stream
        integer,parameter :: SEED   = 7777777
        integer,parameter :: BRNG   = VSL_BRNG_WH
        integer,parameter :: METHOD = VSL_RNG_METHOD_GAUSSIAN_ICDF

        errcode = 0

!       ***** Initialize *****
        errcode = vslnewstream( stream, BRNG,  SEED )

!       ***** Call RNG *****
        errcode = vdrnggaussian( method, stream, n, x, left, right )

!       ***** Deinitialize *****
        errcode = vsldeletestream( stream )

        RETURN
      end function
