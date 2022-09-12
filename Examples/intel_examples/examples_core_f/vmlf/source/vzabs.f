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
!    vzAbs  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      complex(kind=8) zA(10)
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
      complex(kind=8) zA_I(10*inca)
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

      zA( 1)=(-1000.0000d+000,1000.0000d+000)
      zA( 2)=(-777.7777d+000,777.7777d+000)
      zA( 3)=(-555.5555d+000,555.5555d+000)
      zA( 4)=(-333.3333d+000,333.3333d+000)
      zA( 5)=(-111.1111d+000,111.1111d+000)
      zA( 6)=(111.1111d+000,-111.1111d+000)
      zA( 7)=(333.3333d+000,-333.3333d+000)
      zA( 8)=(555.5555d+000,-555.5555d+000)
      zA( 9)=(777.7777d+000,-777.7777d+000)
      zA(10)=(1000.0000d+000,-1000.0000)
      dB( 1)=1.4142135623730951d+003
      dB( 2)=1.0999437718513525d+003
      dB( 3)=7.8567412275096603d+002
      dB( 4)=4.7140447365057963d+002
      dB( 5)=1.5713482455019320d+002
      dB( 6)=1.5713482455019320d+002
      dB( 7)=4.7140447365057963d+002
      dB( 8)=7.8567412275096603d+002
      dB( 9)=1.0999437718513525d+003
      dB(10)=1.4142135623730951d+003

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          dB_I((i-1)*incb+1)=dB(i)
      end do

      call VZABS(vec_len,zA,dBha0)
      call VZABSI(vec_len,zA_I,inca,dBha0_I,incb)

      mode=VML_EP
      call VMZABS(vec_len,zA,dBep1,mode)
      call VMZABSI(vec_len,zA_I,inca,dBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZABS(vec_len,zA,dBep2)
      call VZABSI(vec_len,zA_I,inca,dBep2_I,incb)

      mode=VML_LA
      call VMZABS(vec_len,zA,dBla1,mode)
      call VMZABSI(vec_len,zA_I,inca,dBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZABS(vec_len,zA,dBla2)
      call VZABSI(vec_len,zA_I,inca,dBla2_I,incb)

      mode=VML_HA
      call VMZABS(vec_len,zA,dBha1,mode)
      call VMZABSI(vec_len,zA_I,inca,dBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZABS(vec_len,zA,dBha2)
      call VZABSI(vec_len,zA_I,inca,dBha2_I,incb)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VZABS and"
              print *," VMZABS in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZABSI and"
              print *," VMZABSI in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VZABS and"
              print *," VMZABS in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZABSI and"
              print *," VMZABSI in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VZABS and"
              print *," VMZABS in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZABSI and"
              print *," VMZABSI in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VZABS and"
              print *," VMZABS in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZABSI and"
              print *," VMZABSI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzAbs test/example program"
      print *,""
      print *,"           Argument                     vzAbs"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),"    ",dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzAbsI test/example program"
      print *,"           Argument                     vzAbs"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,zA_I((i-1)*inca+1),"    ",dBha0_I((i-1)*incb+1)
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

10    format(E15.7,E15.7,A5,E15.7)
11    format(A,F25.16)

      end
