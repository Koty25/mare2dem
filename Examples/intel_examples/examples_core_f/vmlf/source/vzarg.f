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
!    vzArg  Example Program Text
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

      zA( 1)=(-10000.0000d+000,10000.0000d+000)
      zA( 2)=(-7777.7777d+000,7777.7777d+000)
      zA( 3)=(-5555.5555d+000,5555.5555d+000)
      zA( 4)=(-3333.3333d+000,3333.3333d+000)
      zA( 5)=(-1111.1111d+000,1111.1111d+000)
      zA( 6)=(1111.1111d+000,-1111.1111d+000)
      zA( 7)=(3333.3333d+000,-3333.3333d+000)
      zA( 8)=(5555.5555d+000,-5555.5555d+000)
      zA( 9)=(7777.7777d+000,-7777.7777d+000)
      zA(10)=(10000.0000d+000,-10000.0000d+000)
      dB( 1)=2.3561944901923448d+000
      dB( 2)=2.3561944901923448d+000
      dB( 3)=2.3561944901923448d+000
      dB( 4)=2.3561944901923448d+000
      dB( 5)=2.3561944901923448d+000
      dB( 6)=-7.8539816339744828d-001
      dB( 7)=-7.8539816339744839d-001
      dB( 8)=-7.8539816339744828d-001
      dB( 9)=-7.8539816339744839d-001
      dB(10)=-7.8539816339744828d-001

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          dB_I((i-1)*incb+1)=dB(i)
      end do

      call VZARG(vec_len,zA,dBha0)
      call VZARGI(vec_len,zA_I,inca,dBha0_I,incb)

      mode=VML_EP
      call VMZARG(vec_len,zA,dBep1,mode)
      call VMZARGI(vec_len,zA_I,inca,dBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZARG(vec_len,zA,dBep2)
      call VZARGI(vec_len,zA_I,inca,dBep2_I,incb)

      mode=VML_LA
      call VMZARG(vec_len,zA,dBla1,mode)
      call VMZARGI(vec_len,zA_I,inca,dBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZARG(vec_len,zA,dBla2)
      call VZARGI(vec_len,zA_I,inca,dBla2_I,incb)

      mode=VML_HA
      call VMZARG(vec_len,zA,dBha1,mode)
      call VMZARGI(vec_len,zA_I,inca,dBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZARG(vec_len,zA,dBha2)
      call VZARGI(vec_len,zA_I,inca,dBha2_I,incb)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VZARG and"
              print *," VMZARG in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZARGI and"
              print *," VMZARGI in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VZARG and"
              print *," VMZARG in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZARGI and"
              print *," VMZARGI in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VZARG and"
              print *," VMZARG in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZARGI and"
              print *," VMZARGI in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VZARG and"
              print *," VMZARG in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZARGI and"
              print *," VMZARGI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzArg test/example program"
      print *,""
      print *,"           Argument                     vzArg"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),"    ",dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzArgI test/example program"
      print *,"           Argument                     vzArg"
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
