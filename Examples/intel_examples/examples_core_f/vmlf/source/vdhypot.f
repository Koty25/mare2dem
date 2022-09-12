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
!    vdHypot  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=8) dA(10)
      real(kind=8) dB(10)
      real(kind=8) dBha0(10)
      real(kind=8) dBha1(10)
      real(kind=8) dBha2(10)
      real(kind=8) dBla1(10)
      real(kind=8) dBla2(10)
      real(kind=8) dBep1(10)
      real(kind=8) dBep2(10)

      real(kind=8) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      real(kind=8) dA_I(10*inca)
      real(kind=8) dB_I(10*incb)
      real(kind=8) dBha0_I(10*incb)
      real(kind=8) dBha1_I(10*incb)
      real(kind=8) dBha2_I(10*incb)
      real(kind=8) dBla1_I(10*incb)
      real(kind=8) dBla2_I(10*incb)
      real(kind=8) dBep1_I(10*incb)
      real(kind=8) dBep2_I(10*incb)
      real(kind=8) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      dA( 1)=1.0009d+000
      dA( 2)=1112.0008d+000
      dA( 3)=2223.0007d+000
      dA( 4)=3334.0006d+000
      dA( 5)=4445.0005d+000
      dA( 6)=5556.0004d+000
      dA( 7)=6667.0003d+000
      dA( 8)=7778.0002d+000
      dA( 9)=8889.0001d+000
      dA(10)=10000.0000d+000
      dB( 1)=1.4154863545792307d+000
      dB( 2)=1.5726066127297315d+003
      dB( 3)=3.1437977391048839d+003
      dB( 4)=4.7149888654800361d+003
      dB( 5)=6.2861799918551887d+003
      dB( 6)=7.8573711182303405d+003
      dB( 7)=9.4285622446054931d+003
      dB( 8)=1.0999753370980647d+004
      dB( 9)=1.2570944497355797d+004
      dB(10)=1.4142135623730950d+004

      do i=1,10
          dA_I((i-1)*inca+1)=dA(i)
          dB_I((i-1)*incb+1)=dB(i)
      end do

      call VDHYPOT(vec_len,dA,dA,dBha0)
      call VDHYPOTI(vec_len,dA_I,inca,dA_I,inca,dBha0_I,incb)

      mode=VML_EP
      call VMDHYPOT(vec_len,dA,dA,dBep1,mode)
      call VMDHYPOTI(vec_len,dA_I,inca,dA_I,inca,dBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VDHYPOT(vec_len,dA,dA,dBep2)
      call VDHYPOTI(vec_len,dA_I,inca,dA_I,inca,dBep2_I,incb)

      mode=VML_LA
      call VMDHYPOT(vec_len,dA,dA,dBla1,mode)
      call VMDHYPOTI(vec_len,dA_I,inca,dA_I,inca,dBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VDHYPOT(vec_len,dA,dA,dBla2)
      call VDHYPOTI(vec_len,dA_I,inca,dA_I,inca,dBla2_I,incb)

      mode=VML_HA
      call VMDHYPOT(vec_len,dA,dA,dBha1,mode)
      call VMDHYPOTI(vec_len,dA_I,inca,dA_I,inca,dBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VDHYPOT(vec_len,dA,dA,dBha2)
      call VDHYPOTI(vec_len,dA_I,inca,dA_I,inca,dBha2_I,incb)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VDHYPOT and"
              print *," VMDHYPOT in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDHYPOTI and"
              print *," VMDHYPOTI in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VDHYPOT and"
              print *," VMDHYPOT in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDHYPOTI and"
              print *," VMDHYPOTI in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VDHYPOT and"
              print *," VMDHYPOT in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDHYPOTI and"
              print *," VMDHYPOTI in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VDHYPOT and"
              print *," VMDHYPOT in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDHYPOTI and"
              print *," VMDHYPOTI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vdHypot test/example program"
      print *,""
      print *,"                    Arguments                         ", &
     &        "     vdHypot"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA(i),dA(i),dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vdHypotI test/example program"
      print *,"                    Arguments                         ", &
     &        "     vdHypot"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA_I((i-1)*inca+1),dA_I((i-1)*inca+1),             &
     &          dBha0_I((i-1)*incb+1)
            CurRMS_I=drelerr(dB_I((i-1)*incb+1),dBha0_I((i-1)*incb+1))
            if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I
      end do

      print *,""
      if(MaxRMS>=1e-14) then
            print 11,"Error! Relative accuracy is ",MaxRMS
            stop 1
      else
            print 11,"Relative accuracy is ",MaxRMS
      endif

      if(MaxRMS_I>=1e-14) then
            print 11,"Error! Relative strided accuracy is ",MaxRMS_I
            stop 1
      else
            print 11,"Relative strided accuracy is ",MaxRMS_I
      endif

10    format(E25.14,E25.14,E25.14)
11    format(A,F25.16)

      end
