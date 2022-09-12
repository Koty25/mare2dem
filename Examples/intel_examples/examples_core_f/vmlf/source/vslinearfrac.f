!===============================================================================
! Copyright 2001-2020 Intel Corporation.
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
!    vsLinearFrac  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=4) fA(10)
      real(kind=4) fB(10)
      real(kind=4) fBha0(10)
      real(kind=4) fBha1(10)
      real(kind=4) fBha2(10)
      real(kind=4) fBla1(10)
      real(kind=4) fBla2(10)
      real(kind=4) fBep1(10)
      real(kind=4) fBep2(10)

      real(kind=4) fScale
      real(kind=4) fShift
      real(kind=4) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      real(kind=4) fA_I(10*inca)
      real(kind=4) fB_I(10*incb)
      real(kind=4) fBha0_I(10*incb)
      real(kind=4) fBha1_I(10*incb)
      real(kind=4) fBha2_I(10*incb)
      real(kind=4) fBla1_I(10*incb)
      real(kind=4) fBla2_I(10*incb)
      real(kind=4) fBep1_I(10*incb)
      real(kind=4) fBep2_I(10*incb)
      real(kind=4) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      fA( 1)=-10000.0000
      fA( 2)=-7777.7778
      fA( 3)=-5555.5557
      fA( 4)=-3333.3333
      fA( 5)=-1111.1111
      fA( 6)=1111.1111
      fA( 7)=3333.3333
      fA( 8)=5555.5557
      fA( 9)=7777.7778
      fA(10)=10000.0000
      fB( 1)=1.0000000000000000e+000
      fB( 2)=1.0000000000000000e+000
      fB( 3)=1.0000000000000000e+000
      fB( 4)=1.0000000000000000e+000
      fB( 5)=1.0000000000000000e+000
      fB( 6)=1.0000000000000000e+000
      fB( 7)=1.0000000000000000e+000
      fB( 8)=1.0000000000000000e+000
      fB( 9)=1.0000000000000000e+000
      fB(10)=1.0000000000000000e+000
      fScale=1.0000000000000000e+000
      fShift=0.0000000000000000e+000

      do i=1,10
          fA_I((i-1)*inca+1)=fA(i)
          fB_I((i-1)*incb+1)=fB(i)
      end do

      call VSLINEARFRAC(vec_len,fA,fA,fScale,fShift,fScale,fShift,fBha0)
      call VSLINEARFRACI(vec_len,fA_I,inca,fA_I,inca,fScale,fShift,     &
     &    fScale,fShift,fBha0_I,incb)

      mode=VML_EP
      call VMSLINEARFRAC(vec_len,fA,fA,fScale,fShift,fScale,fShift,     &
     &    fBep1,mode)
      call VMSLINEARFRACI(vec_len,fA_I,inca,fA_I,inca,fScale,fShift,    &
     &    fScale,fShift,fBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VSLINEARFRAC(vec_len,fA,fA,fScale,fShift,fScale,fShift,fBep2)
      call VSLINEARFRACI(vec_len,fA_I,inca,fA_I,inca,fScale,fShift,     &
     &    fScale,fShift,fBep2_I,incb)

      mode=VML_LA
      call VMSLINEARFRAC(vec_len,fA,fA,fScale,fShift,fScale,fShift,     &
     &    fBla1,mode)
      call VMSLINEARFRACI(vec_len,fA_I,inca,fA_I,inca,fScale,fShift,    &
     &    fScale,fShift,fBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VSLINEARFRAC(vec_len,fA,fA,fScale,fShift,fScale,fShift,fBla2)
      call VSLINEARFRACI(vec_len,fA_I,inca,fA_I,inca,fScale,fShift,     &
     &    fScale,fShift,fBla2_I,incb)

      mode=VML_HA
      call VMSLINEARFRAC(vec_len,fA,fA,fScale,fShift,fScale,fShift,     &
     &    fBha1,mode)
      call VMSLINEARFRACI(vec_len,fA_I,inca,fA_I,inca,fScale,fShift,    &
     &    fScale,fShift,fBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VSLINEARFRAC(vec_len,fA,fA,fScale,fShift,fScale,fShift,fBha2)
      call VSLINEARFRACI(vec_len,fA_I,inca,fA_I,inca,fScale,fShift,     &
     &    fScale,fShift,fBha2_I,incb)

      do i=1,10
          if(fBha0(i) .ne. fBha1(i)) then
              print *,"Error! Difference between VSLINEARFRAC and"
              print *," VMSLINEARFRAC in VML_HA mode detected"
              stop 1
          endif
          if(fBha0_I((i-1)*incb+1) .ne. fBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSLINEARFRACI and"
              print *," VMSLINEARFRACI in VML_HA mode detected"
              stop 1
          endif
          if(fBha1(i) .ne. fBha2(i)) then
              print *,"Error! Difference between VSLINEARFRAC and"
              print *," VMSLINEARFRAC in VML_HA mode detected"
              stop 1
          endif
          if(fBha1_I((i-1)*incb+1) .ne. fBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSLINEARFRACI and"
              print *," VMSLINEARFRACI in VML_HA mode detected"
              stop 1
          endif
          if(fBla1(i) .ne. fBla2(i)) then
              print *,"Error! Difference between VSLINEARFRAC and"
              print *," VMSLINEARFRAC in VML_LA mode detected"
              stop 1
          endif
          if(fBla1_I((i-1)*incb+1) .ne. fBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSLINEARFRACI and"
              print *," VMSLINEARFRACI in VML_LA mode detected"
              stop 1
          endif
          if(fBep1(i) .ne. fBep2(i)) then
              print *,"Error! Difference between VSLINEARFRAC and"
              print *," VMSLINEARFRAC in VML_EP mode detected"
              stop 1
          endif
          if(fBep1_I((i-1)*incb+1) .ne. fBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSLINEARFRACI and"
              print *," VMSLINEARFRACI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vsLinearFrac test/example program"
      print *,""
      print 10,"Scalar Parameters: ScaleA = ScaleB = ", fScale
      print 10,"                   ScaleA = ScaleB = ", fShift
      print *,"                    Arguments                         ", &
     &        "      vsLinearFrac"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 11,fA(i),fA(i),fBha0(i)
            CurRMS=srelerr(fB(i),fBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vsLinearFracI test/example program"
      print 10,"Scalar Parameters: ScaleA = ScaleB = ", fScale
      print 10,"                   ScaleA = ScaleB = ", fShift
      print *,"                    Arguments                         ", &
     &        "      vsLinearFrac"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 11,fA_I((i-1)*inca+1),fA_I((i-1)*inca+1),             &
     &          fBha0_I((i-1)*incb+1)
            CurRMS_I=srelerr(fB_I((i-1)*incb+1),fBha0_I((i-1)*incb+1))
            if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I
      end do

      print *,""
      if(MaxRMS>=1e-3) then
            print 10,"Error! Relative accuracy is ",MaxRMS
            stop 1
      else
            print 10,"Relative accuracy is ",MaxRMS
      endif

      if(MaxRMS_I>=1e-3) then
            print 10,"Error! Relative strided accuracy is ",MaxRMS_I
            stop 1
      else
            print 10,"Relative strided accuracy is ",MaxRMS_I
      endif

10    format(A,F25.16)
11    format(E25.14,E25.14,E25.14)

      end
