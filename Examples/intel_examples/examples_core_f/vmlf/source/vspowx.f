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
!    vsPowx  Example Program Text
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

      fA( 1)=0.1000
      fA( 2)=0.8666
      fA( 3)=1.6333
      fA( 4)=2.4000
      fA( 5)=3.1666
      fA( 6)=3.9333
      fA( 7)=4.7000
      fA( 8)=5.4666
      fA( 9)=6.2333
      fA(10)=7.0000
      fB( 1)=1.1660039110312152e-004
      fB( 2)=5.6940619267520864e-001
      fB( 3)=6.8873522741537876e+000
      fB( 4)=3.1295720312255717e+001
      fB( 5)=9.3107199918766597e+001
      fB( 6)=2.1845299866537277e+002
      fB( 7)=4.4011122728416541e+002
      fB( 8)=7.9737651355941364e+002
      fB( 9)=1.3361785535769964e+003
      fB(10)=2.1087453418813361e+003

      do i=1,10
          fA_I((i-1)*inca+1)=fA(i)
          fB_I((i-1)*incb+1)=fB(i)
      end do

      call VSPOWX(vec_len,fA,fA(6),fBha0)
      call VSPOWXI(vec_len,fA_I,inca,fA(6),fBha0_I,incb)

      mode=VML_EP
      call VMSPOWX(vec_len,fA,fA(6),fBep1,mode)
      call VMSPOWXI(vec_len,fA_I,inca,fA(6),fBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VSPOWX(vec_len,fA,fA(6),fBep2)
      call VSPOWXI(vec_len,fA_I,inca,fA(6),fBep2_I,incb)

      mode=VML_LA
      call VMSPOWX(vec_len,fA,fA(6),fBla1,mode)
      call VMSPOWXI(vec_len,fA_I,inca,fA(6),fBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VSPOWX(vec_len,fA,fA(6),fBla2)
      call VSPOWXI(vec_len,fA_I,inca,fA(6),fBla2_I,incb)

      mode=VML_HA
      call VMSPOWX(vec_len,fA,fA(6),fBha1,mode)
      call VMSPOWXI(vec_len,fA_I,inca,fA(6),fBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VSPOWX(vec_len,fA,fA(6),fBha2)
      call VSPOWXI(vec_len,fA_I,inca,fA(6),fBha2_I,incb)

      do i=1,10
          if(fBha0(i) .ne. fBha1(i)) then
              print *,"Error! Difference between VSPOWX and"
              print *," VMSPOWX in VML_HA mode detected"
              stop 1
          endif
          if(fBha0_I((i-1)*incb+1) .ne. fBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSPOWXI and"
              print *," VMSPOWXI in VML_HA mode detected"
              stop 1
          endif
          if(fBha1(i) .ne. fBha2(i)) then
              print *,"Error! Difference between VSPOWX and"
              print *," VMSPOWX in VML_HA mode detected"
              stop 1
          endif
          if(fBha1_I((i-1)*incb+1) .ne. fBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSPOWXI and"
              print *," VMSPOWXI in VML_HA mode detected"
              stop 1
          endif
          if(fBla1(i) .ne. fBla2(i)) then
              print *,"Error! Difference between VSPOWX and"
              print *," VMSPOWX in VML_LA mode detected"
              stop 1
          endif
          if(fBla1_I((i-1)*incb+1) .ne. fBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSPOWXI and"
              print *," VMSPOWXI in VML_LA mode detected"
              stop 1
          endif
          if(fBep1(i) .ne. fBep2(i)) then
              print *,"Error! Difference between VSPOWX and"
              print *," VMSPOWX in VML_EP mode detected"
              stop 1
          endif
          if(fBep1_I((i-1)*incb+1) .ne. fBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VSPOWXI and"
              print *," VMSPOWXI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vsPowx test/example program"
      print *,""
      print *,"           Argument                     vsPowx"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA(i),fBha0(i)
            CurRMS=srelerr(fB(i),fBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vsPowxI test/example program"
      print *,"           Argument                     vsPowx"
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
