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

!  Content:
!    Use of Intel(R) MKL Sobol QRNG with direction numbers designed
!    by S. Joe,  F. Y. Kuo and supporting dimension up to 21201
!    Example Program Text
!*******************************************************************************

!     This sample demonsrates registration of user-defined parameters
!     in Sobol QRNG. Set of initial Sobol parameters used here
!     follows to [1],[2].

!     Bibliography
!     1. S. Joe,  F. Y. Kuo. Remark on Algorithm 659: Implementing Sobol_s
!     Quasirandom Sequence Generator. ACM Transactions on Mathematical Software,
!     Vol. 29, No. 1, 50-57, March 2003
!     2. http://web.maths.unsw.edu.au/~fkuo/sobol/

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      TYPE (VSL_STREAM_STATE) :: stream

      integer dm
      integer(kind=4) nbits, npoly, maxdeg, nu
!     Dimension of vectors to generate
      parameter (dm=111)
!     Number of bits used for storing unsigned integer
      parameter (nbits=32)
!     Number of primitive polynomials
      parameter (npoly=110)
!     Maximum degree of primitive polynomials for dimensions up to 111
      parameter (maxdeg=10)
      parameter (nu=3+dm*nbits)

      integer brng, method
      integer(kind=4) nn
      integer n, m, nblocks
!     Number of vectors to generate
      parameter (m=10000)
      parameter (nblocks=100)
      parameter (n=m*dm)

      real(kind=8) a, b

      real(kind=8) r(n)

      integer(kind=4) SobolV(dm, nbits)

      integer(kind=4) poly(2:dm),vinit(2:dm,maxdeg)
      integer(kind=4) polydeg(2:dm)
      integer(kind=4) params(nu)
      integer nparams
      integer(kind=4) errnum, errcode

      real(kind=8) res, res_err

      nn     = 10

      brng   = VSL_BRNG_SOBOL
      method = VSL_RNG_METHOD_UNIFORM_STD
      a      = 0.0
      b      = 1.0

!     ***** Read initial Sobol direction numbers, polynomials and
!           polynomial degrees from file *****
      errcode=ReadSobolData(dm, maxdeg, vinit, poly, polydeg)
      if (errcode /= 0) then
        stop 1
      end if

!     Calculate table of Sobol direction numbers using Joe and Kuo
!     initial direction numbers and primitive polynomials
      call InitSobolMatrix( dm, nbits, SobolV, vinit, poly, polydeg )

!     Pack Sobol direction numbers to pass into the generator
      params(1) = dm
      params(2) = VSL_USER_QRNG_INITIAL_VALUES
      params(3) = VSL_USER_DIRECTION_NUMBERS
      k = 4

      do i = 1, dm
         do j = 1, nbits
             params(k) = SobolV(i,j)
             k = k + 1
         end do
      end do
      nparams = 3+dm*nbits

!     ***** Initialize *****
      errcode=vslnewstreamex( stream, brng,  nparams, params )
      call CheckVslError(errcode)

!     ***** Call RNG *****
      errcode=vdrnguniform( method, stream, n, r, a, b )
      call CheckVslError(errcode)

!     ***** Printing results *****
      print *,"Sobol sequence with user-defined direction numbers:"
      print *,"---------------------------------------------------"
      print *,""
      print *,"Parameters:"
      print 103,"    a=",a
      print 103,"    b=",b

      k = 12
      print *,""
      print *,"Results:"
      print *, k, "component of quasi-random sequence"
      print *,"first ", nn, " of ", m
      print *,"----------------------------------------------"

      do i=1,nn
        write (*,102,advance='no') r((k+1)+(i-1)*dm)
      end do
      print *, ""

!     ***** Test generated values by calculating integral of the model
!           function f(X) = Prod((|4*Xj - 2| + Cj)/(1+Cj))
!           Model function is taken from [1] *****
      call Integrate( dm, m, nblocks, method, stream, a, b, r, res,     &
     &                res_err )

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)

      print *, ""
      print *, "Calculate integral of the function:"
      print *, "  f(X) = Prod((|4*Xj - 2| + Cj)/(1+Cj)),  j = 1,...,dim"
      print *, "  Cj = j^(1/3)"
      print *, "over unit hypercube [0,1]^dim"
      print *, ""
      print 103," Computed result: ", res
      print 103," Expected result: ", 1.0
      print 103," Calculation error: ", res_err
      print *, ""

      if ( (1.0 - 3.0*res_err .LE. res) .AND.                           &
     &     (res .LE. 1.0 + 3.0*res_err) ) then
        print *,"Numerical integration results agree with theory"
      else
        print *,"Numerical integration results do not agree"
        print *,"with theory"
      end if

102   format(F7.4)
103   format(A,F10.7)

      contains

!     Read initial Sobol dimension values, polynomials and polynomial
!     degress from file

      function ReadSobolData( dimen, maxdeg, vinit, poly, polydeg )
        integer, intent(in) :: dimen
        integer(kind=4), intent(in) :: maxdeg
        integer(kind=4), intent(out) :: vinit(2:dimen,maxdeg)
        integer(kind=4), intent(out) :: poly(2:dimen)
        integer(kind=4), intent(out) :: polydeg(2:dimen)

        integer(kind=4) j,i,tmp_dim,errcode
        character*100 filename

        read(*,fmt='(A)') filename

        open(unit=10, file=filename, iostat=errcode, status='OLD')
        if ( errcode /= 0 ) then
          ReadSobolData = errcode
          return
        end if

        read(10,fmt="(A32)") line
        do i = 2,dm
          read(10,100,advance='no') tmp_dim, polydeg(i), poly(i)

!         Modify the coding of polynomials to the format expected
!         by the library
          poly(i) = ishft(poly(i),int(1,kind=4))
          poly(i) = ior(poly(i),int(1,kind=4))
          poly(i) = ior(poly(i),ishft(int(1,kind=4),polydeg(i)))

          do j = 1,polydeg(i)-int(1,kind=4)
            read(10,101,advance='no') vinit(i,j)
          end do
          read(10,101) vinit(i,polydeg(i))
        end do

        close(unit=10)

100     format(I8,I8,I8)
101     format(I8)

        ReadSobolData = 0
      end function ReadSobolData

!     Initialization of Sobol matrix

      subroutine InitSobolMatrix( dimen, nbits, SobolV, vinit, poly,    &
     &      polydeg )
        integer, intent(in)  :: dimen
        integer(kind=4), intent(in)  :: nbits
        integer(kind=4), intent(out) :: SobolV(dimen,nbits)
        integer(kind=4), intent(in) :: vinit(2:dimen,maxdeg)
        integer(kind=4), intent(in) :: poly(2:dimen)
        integer(kind=4), intent(in) :: polydeg(2:dimen)

        integer(kind=4) j,i,d
        integer(kind=4) ipolydeg, ipoly, new

        do j=1,nbits
           SobolV(1,j) = 1
        end do

        do i=2,dimen
           ipolydeg = polydeg(i)
           do j=1,ipolydeg
              SobolV(i,j) = vinit(i,j)
           end do

           do j=ipolydeg+1, nbits
              ipoly = poly(i)
              new = SobolV(i,j-ipolydeg)

              do d=ipolydeg,1,-1
                 if ( mod(ipoly,2) .eq. 1 ) then
                   new = ieor(new,ishft(SobolV(i,j-d),d))
                 end if
                 ipoly = ipoly / 2
              end do

              SobolV(i,j) = new
           end do

        end do

        do j=1,nbits
           do i=1,dimen
              SobolV(i,j) = ishft( SobolV(i,j), nbits-j )
           end do
        end do

      end subroutine InitSobolMatrix

!     Calculate integral for the function
!     f(X) = Prod((|4*Xj - 2| + Cj)/(1+Cj))
!     over unit hypercube [0,1]^dimen.
!     Here Cj = j^(1/3). Integration result should be equal to 1.
!     Model function is taken from [1]

      subroutine Integrate( dimen, m, nblocks, method, stream, a, b, r, &
     &                      res, res_err )
        USE MKL_VSL_TYPE
        USE MKL_VSL

        integer, intent(in) :: dimen
        integer, intent(in) :: m
        integer, intent(in) :: nblocks
        integer, intent(in) :: method
        TYPE (VSL_STREAM_STATE), intent(inout) :: stream
        real(kind=8), intent(in) :: a
        real(kind=8), intent(in) :: b
        real(kind=8), intent(inout) :: r(m*dimen)
        real(kind=8), intent(out) :: res, res_err

        integer(kind=4) i, j, errcode
        integer n, r_size
        real(kind=8) C(dimen), D(dimen)
        real(kind=8) inv_n, prod, f, f2

        r_size = m * dimen
        n = m * nblocks
        inv_n = 1.0D0 / dble(n)

        f = 0.0D0
        f2 = 0.0D0

        do j=1,dimen
          D(j) = dble(j)
        end do
        call vdcbrt(dimen,D,C)

        do iblock=1,nblocks
          errcode=vdrnguniform( method, stream, r_size, r, a, b )
          call CheckVslError(errcode)
          do i=0,m-1
            prod = 1.0D0
            do j=1,dimen
              prod = prod * (dabs(4.0D0*r(i*dimen + j) - 2.0D0) + C(j))
              prod = prod / (1.0D0 + C(j))
            end do
            f = f + prod
            f2 = f2 + prod*prod
          end do
        end do

        f = f * inv_n
        f2 = f2 * inv_n

        res = f
        res_err = sqrt((f2 - f*f)*inv_n)

      end subroutine Integrate

      end program
