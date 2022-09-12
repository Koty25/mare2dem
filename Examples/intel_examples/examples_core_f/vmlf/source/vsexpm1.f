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
!    vsExpm1  Example Program Text
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

      fA( 1)=-17.0000
      fA( 2)=-13.1111
      fA( 3)=-9.2222
      fA( 4)=-5.3333
      fA( 5)=-1.4444
      fA( 6)=2.4444
      fA( 7)=6.3333
      fA( 8)=10.2222
      fA( 9)=14.1111
      fA(10)=18.0000
      fB( 1)=-9.9999995860062285e-001
      fB( 2)=-9.9999797734630202e-001
      fB( 3)=-9.9990117899688591e-001
      fB( 4)=-9.9517188962066971e-001
      fB( 5)=-7.6411242224240106e-001
      fB( 6)=1.0523634174991830e+001
      fB( 7)=5.6201153347619049e+002
      fB( 8)=2.7506126448834904e+004
      fB( 9)=1.3439175517311485e+006
      fB(10)=6.5659968137330510e+007

      do i=1,10
          fA_I((i-1)*inca+1)=fA(i)
          fB_I((i-1)*incb+1)=fB(i)
      end do

      call VSEXPM1(vec_len,fA,fBha0)
      call VSEXPM1I(vec_len,fA_I,inca,fBha0_I,incb)

      mode=VML_EP
      call VMSEXPM1(vec_len,fA,fBep1,mode)
      call VMSEXPM1I(vec_len,fA_I,inca,fBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VSEXPM1(vec_len,fA,fBep2)
      call VSEXPM1I(vec_len,fA_I,inca,fBep2_I,incb)

      mode=VML_LA
      call VMSEXPM1(vec_len,fA,fBla1,mode)
      call VMSEXPM1I(vec_len,fA_I,inca,fBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VSEXPM1(vec_len,fA,fBla2)
      call VSEXPM1I(vec_len,fA_I,inca,fBla2_I,incb)

      mode=VML_HA
      call VMSEXPM1(vec_len,fA,fBha1,mode)
      call VMSEXPM1I(vec_len,fA_I,inca,fBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VSEXPM1(vec_len,fA,fBha2)
      call VSEXPM1I(vec_len,fA_I,inca,fBha2_I,incb)

      do i=1,10
          if(fBha0(i) .ne. fBha1(i)) then
              print *,"Error! Difference between VSEXPM1 and"
              print *," VMSEXPM1 in VML_HA mode detected"
              stop 1
          endif
          if(fBha0_I((i-1)*incb+1) .ne. fBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSEXPM1I and"
              print *," VMSEXPM1I in VML_HA mode detected"
              stop 1
          endif
          if(fBha1(i) .ne. fBha2(i)) then
              print *,"Error! Difference between VSEXPM1 and"
              print *," VMSEXPM1 in VML_HA mode detected"
              stop 1
          endif
          if(fBha1_I((i-1)*incb+1) .ne. fBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSEXPM1I and"
              print *," VMSEXPM1I in VML_HA mode detected"
              stop 1
          endif
          if(fBla1(i) .ne. fBla2(i)) then
              print *,"Error! Difference between VSEXPM1 and"
              print *," VMSEXPM1 in VML_LA mode detected"
              stop 1
          endif
          if(fBla1_I((i-1)*incb+1) .ne. fBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSEXPM1I and"
              print *," VMSEXPM1I in VML_LA mode detected"
              stop 1
          endif
          if(fBep1(i) .ne. fBep2(i)) then
              print *,"Error! Difference between VSEXPM1 and"
              print *," VMSEXPM1 in VML_EP mode detected"
              stop 1
          endif
          if(fBep1_I((i-1)*incb+1) .ne. fBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSEXPM1I and"
              print *," VMSEXPM1I in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vsExpm1 test/example program"
      print *,""
      print *,"           Argument                     vsExpm1"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA(i),fBha0(i)
            CurRMS=srelerr(fB(i),fBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vsExpm1I test/example program"
      print *,"           Argument                     vsExpm1"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA_I((i-1)*inca+1),fBha0_I((i-1)*incb+1)
            CurRMS_I=srelerr(fB_I((i-1)*incb+1),fBha0_I((i-1)*incb+1))
            if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I
      end do

      print *,""
      if(MaxRMS>=1e-5) then
            print 11,"Error! Relative accuracy is ",MaxRMS
            stop 1
      else
            print 11,"Relative accuracy is ",MaxRMS
      endif

      if(MaxRMS_I>=1e-5) then
            print 11,"Error! Relative strided accuracy is ",MaxRMS_I
            stop 1
      else
            print 11,"Relative strided accuracy is ",MaxRMS_I
      endif

10    format(E25.14,E25.14)
11    format(A,F25.16)

      end
