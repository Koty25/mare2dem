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
!    User-defined direction numbers for Sobol QRNG  Example Program Text
!    Owen-scrambling strategy
!*******************************************************************************


!     This sample demonsrates registration of user-defined parameters
!     in Sobol QRNG. Owen scrambling strategy with slight changes
!     is also implemented, [1]. Set of initial Sobol parameters used here
!     follows to [2].

!     Bibliography
!     1. H.S. Hong, and F.J. Hickernell. Algorithm 823: Implementing
!     Scrambled Digital Sequences. ACM Transactions on Mathematical Software,
!     Vol.29, No. 2, 95-109, 2003
!     2. Sobol, I.M., and Levitan, Yu.L. The production of points uniformly
!     distributed in a multidimensional cube. Preprint 40, Institute of
!     Applied Mathematics, USSR Academy of Sciences, 1976 (In Russian)

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      TYPE (VSL_STREAM_STATE) :: stream

      integer(kind=4) maxdim, nbits, npoly, maxdeg, s_digits, nu
      parameter (maxdim=40)
      parameter (nbits=32)
      parameter (npoly=39)
      parameter (maxdeg=8)
      parameter (s_digits=32)
      parameter (nu=3+maxdim*nbits)

      integer brng, method, dm
      integer(kind=4) nn, dm1
      integer n
      parameter (dm=20)
      parameter (dm1=20)
      parameter (n=1000*dm)

      real(kind=8) a, b

      real(kind=8) r(n), rt(n)
      integer(kind=4) dimen

      integer(kind=4) SobolV(maxdim, nbits), SobolSV(maxdim, nbits)

      integer(kind=4) poly(2:maxdim),vinit(2:maxdim,maxdeg)
      integer(kind=4) polydeg(2:maxdim)
      integer(kind=4) params(nu)
      integer nparams
      integer(kind=4) errnum, errcode


      nn     = 10

      brng   = VSL_BRNG_SOBOL
      method = VSL_RNG_METHOD_UNIFORM_STD
      a      = 0.0
      b      = 1.0

!     Generation of Sobol sequence using direction
!     numbers registered in VSL

!     ***** Initialize *****
      errcode=vslnewstream( stream, brng,  dm )
      call CheckVslError(errcode)

!     ***** Call RNG *****
      errcode=vdrnguniform( method, stream, n, rt, a, b )
      call CheckVslError(errcode)


!     ***** Printing results *****
      print *,"Sobol sequence with default direction numbers:"
      print *,"----------------------------------------------"
      print *,""
      print *,"Parameters:"
      print 11,"    a=",a
      print 11,"    b=",b

      k = 12
      print *,""
      print *,"Results:"
      print *, k, "component of quasi-random sequence"
      print *,"first ", nn, " of ", n/dm
      print *,"---------------------------"
      do i=1,nn
        print 10,rt((k+1)+(i-1)*dm)
      end do

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)



!     GENERATION OF SOBOL SEQUENCE USING EXTERNAL DIRECTION NUMBERS
!     For a moment they are the same as those registered in the library
      call InitSobolMatrix( dm1 )
      params(1) = dm1
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
      print 11,"    a=",a
      print 11,"    b=",b

      k = 12
      print *,""
      print *,"Results:"
      print *, k, "component of quasi-random sequence"
      print *,"first ", nn, " of ", n/dm
      print *,"----------------------------------------------"

      do i=1,nn
        print 10,r((k+1)+(i-1)*dm)
      end do

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)


      errnum = 0
      do i=1,n
        if ( rt(i) .ne. r(i) ) then
          errnum = errnum + 1
        end if
      end do

      print *, ""
      if (errnum .eq. 0 ) then
         print *, "Sobol sequences obtained using two methods are",     &
     &            " similar"
      else
         print *, "Sobol sequences obtained using two methods are",     &
     &            " not similar"
      end if

      print *, ""

!     Owen-like scrambling strategy with slight changes
!     is used for demonstration
      call InitSobolMatrixOs( dm1 )
      params(1) = dm1
      params(2) = VSL_USER_QRNG_INITIAL_VALUES
      params(3) = VSL_USER_DIRECTION_NUMBERS
      k = 4

      do i = 1, dm
         do j = 1, nbits
             params(k) = SobolSV(i,j)
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
      print *,"Scrambled Sobol sequence:"
      print *,"---------------------------------------------------"
      print *,""
      print *,"Parameters:"
      print 11,"    a=",a
      print 11,"    b=",b

      k = 12
      print *,""
      print *,"Results:"
      print *, k, "component of quasi-random sequence"
      print *,"first ", nn, " of ", n/dm
      print *,"----------------------------------------------"

      do i=1,nn
        print 10,r((k+1)+(i-1)*dm)
      end do

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)




10    format(F7.4)
11    format(A,F5.3)
12    format(A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)


c     Set of initial Sobol parameters used here  follows to [2]
c     Table of initial direction numbers
      DATA poly    /3,7,11,13,19,25,37,59,47,
     +              61,55,41,67,97,91,109,103,115,131,
     +              193,137,145,143,241,157,185,167,229,171,
     +              213,191,253,203,211,239,247,285,369,299/
      DATA (vinit(I,1),I=2,40)  /39*1/
      DATA (vinit(I,2),I=3,40)  /1,3,1,3,1,3,3,1,
     +                           3,1,3,1,3,1,1,3,1,3,
     +                           1,3,1,3,3,1,3,1,3,1,
     +                           3,1,1,3,1,3,1,3,1,3/
      DATA (vinit(I,3),I=4,40)  /7,5,1,3,3,7,5,
     +                           5,7,7,1,3,3,7,5,1,1,
     +                           5,3,3,1,7,5,1,3,3,7,
     +                           5,1,1,5,7,7,5,1,3,3/
      DATA (vinit(I,4),I=6,40)  /1,7,9,13,11,
     +                           1,3,7,9,5,13,13,11,3,15,
     +                           5,3,15,7,9,13,9,1,11,7,
     +                           5,15,1,15,11,5,3,1,7,9/
      DATA (vinit(I,5),I=8,40)  /9,3,27,
     +                           15,29,21,23,19,11,25,7,13,17,
     +                           1,25,29,3,31,11,5,23,27,19,
     +                           21,5,1,17,13,7,15,9,31,9/
      DATA (vinit(I,6),I=14,40) /37,33,7,5,11,39,63,
     +                           27,17,15,23,29,3,21,13,31,25,
     +                           9,49,33,19,29,11,19,27,15,25/
      DATA (vinit(I,7),I=20,40) /13,
     +                           33,115,41,79,17,29,119,75,73,105,
     +                           7,59,65,21,3,113,61,89,45,107/
      DATA (vinit(I,8),I=38,40) /7,23,39/


      DATA polydeg    /1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,
     +              7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8/

      contains
c     Initialization of Sobol matrix

        subroutine InitSobolMatrix( dimen )
          integer(kind=4), intent(in)  :: dimen

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


c       Routines for Owen-like scrambling
c       Generator matrix SobolV should be initialized before
c       its transformation for scrambling

        subroutine InitSobolMatrixOs( dimen )
           integer(kind=4), intent(in) :: dimen

           integer(kind=4) lsm(maxdim,s_digits)
           integer(kind=4) i,j,k,l
           integer(kind=4) shift
           integer(kind=4) column, sum, a, b, mask, tmptmp

           tmptmp = NOT(0)
           if ( s_digits < nbits ) then
               mask = ishft(tmptmp, s_digits )
           else
               mask = 0
           endif

           call GenerateLSM( dimen,lsm )

             do i=1, dimen
                do j = 1, nbits
                     shift = 0;
                     column = 0
                     do l = s_digits, 1,-1
                         sum = 0;
                         a = lsm(i,l)
                         b = SobolV(i,j)
                         do k=1,nbits
                             sum = sum + iand(iand(a,1_4), iand(b,1_4))
                             a = ishft( a, -1_4 )
                             b = ishft( b, -1_4 )
                         end do
                         sum = mod ( sum, 2_4 )
                         sum = ishft(sum, shift )
                         column = ior(column, sum )
                         shift = shift + 1_4
                     end do

                     if ( s_digits < nbits ) then
                        SobolSV(i,j) = ior( iand(mask, SobolV(i,j)),    &
     &                                      column )
                     else
                        SobolSV(i,j) = column
                     end if
                end do
            end do

      end subroutine InitSobolMatrixOs

      subroutine GenerateLSM ( dimen, lsm )
          integer(kind=4), intent(in)  :: dimen
          integer(kind=4), intent(out) :: lsm(maxdim,s_digits)

          TYPE (VSL_STREAM_STATE) :: sstream
          integer(kind=4) errcode
          integer brng, method, seed
          real(kind=8) p
          integer(kind=4) rb(maxdim*s_digits*(s_digits-1)/2)
          integer(kind=4) shift, elem
          integer nrb

          method = VSL_RNG_METHOD_BERNOULLI_ICDF
          brng = VSL_BRNG_MCG31
          seed = 77987391
          p = 0.5
          nrb = dimen*0.5*s_digits*(s_digits-1)

          errcode=vslnewstream( sstream, brng,  seed )
          call CheckVslError(errcode)
          errcode = virngbernoulli( method, sstream, nrb, rb, p )
          call CheckVslError(errcode)
          errcode=vsldeletestream( sstream )
          call CheckVslError(errcode)

          l = 1
          do i=1,dimen
             do j=s_digits,1,-1

                shift = nbits-j
                lsm(i,j) = ishft(1,shift)

                do k=j-1,1,-1
                   shift = shift + 1
                   elem =  ishft( rb(l), shift )
                   lsm(i,j) = ior( lsm(i,j), elem )
                   l = l + 1
                end do
             end do
          end do

      end subroutine GenerateLSM

      end
