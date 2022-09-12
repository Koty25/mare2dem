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
!    vzAsin  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      complex(kind=8) zA(10)
      complex(kind=8) zB(10)
      complex(kind=8) zBha0(10)
      complex(kind=8) zBha1(10)
      complex(kind=8) zBha2(10)
      complex(kind=8) zBla1(10)
      complex(kind=8) zBla2(10)
      complex(kind=8) zBep1(10)
      complex(kind=8) zBep2(10)

      real(kind=8) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      complex(kind=8) zA_I(10*inca)
      complex(kind=8) zB_I(10*incb)
      complex(kind=8) zBha0_I(10*incb)
      complex(kind=8) zBha1_I(10*incb)
      complex(kind=8) zBha2_I(10*incb)
      complex(kind=8) zBla1_I(10*incb)
      complex(kind=8) zBla2_I(10*incb)
      complex(kind=8) zBep1_I(10*incb)
      complex(kind=8) zBep2_I(10*incb)
      real(kind=8) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      zA( 1)=(-1.0000d+000,1.0000d+000)
      zA( 2)=(-0.7777d+000,0.7777d+000)
      zA( 3)=(-0.5555d+000,0.5555d+000)
      zA( 4)=(-0.3333d+000,0.3333d+000)
      zA( 5)=(-0.1111d+000,0.1111d+000)
      zA( 6)=(0.1111d+000,-0.1111d+000)
      zA( 7)=(0.3333d+000,-0.3333d+000)
      zA( 8)=(0.5555d+000,-0.5555d+000)
      zA( 9)=(0.7777d+000,-0.7777d+000)
      zA(10)=(1.0000d+000,-1.0000d+000)
      zB( 1)=(-6.6623943249251527d-001,1.0612750619050357d+000)
      zB( 2)=(-6.0135587708509519d-001,8.4060455895425301d-001)
      zB( 3)=(-4.8944080432477105d-001,5.9386557353258917d-001)
      zB( 4)=(-3.1990734647712815d-001,3.4427275439844690d-001)
      zB( 5)=(-1.1063788765091863d-001,1.1155195888225219d-001)
      zB( 6)=(1.1063788765091863d-001,-1.1155195888225219d-001)
      zB( 7)=(3.1990734647712815d-001,-3.4427275439844690d-001)
      zB( 8)=(4.8944080432477105d-001,-5.9386557353258917d-001)
      zB( 9)=(6.0135587708509519d-001,-8.4060455895425301d-001)
      zB(10)=(6.6623943249251527d-001,-1.0612750619050357d+000)

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZASIN(vec_len,zA,zBha0)
      call VZASINI(vec_len,zA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZASIN(vec_len,zA,zBep1,mode)
      call VMZASINI(vec_len,zA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZASIN(vec_len,zA,zBep2)
      call VZASINI(vec_len,zA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZASIN(vec_len,zA,zBla1,mode)
      call VMZASINI(vec_len,zA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZASIN(vec_len,zA,zBla2)
      call VZASINI(vec_len,zA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZASIN(vec_len,zA,zBha1,mode)
      call VMZASINI(vec_len,zA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZASIN(vec_len,zA,zBha2)
      call VZASINI(vec_len,zA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZASIN and"
              print *," VMZASIN in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZASINI and"
              print *," VMZASINI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZASIN and"
              print *," VMZASIN in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZASINI and"
              print *," VMZASINI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZASIN and"
              print *," VMZASIN in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZASINI and"
              print *," VMZASINI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZASIN and"
              print *," VMZASIN in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZASINI and"
              print *," VMZASINI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzAsin test/example program"
      print *,""
      print *,"           Argument                     vzAsin"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzAsinI test/example program"
      print *,"           Argument                     vzAsin"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,zA_I((i-1)*inca+1),"    ",zBha0_I((i-1)*incb+1)
            CurRMS_I=zrelerr(zB_I((i-1)*incb+1),zBha0_I((i-1)*incb+1))
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

10    format(E15.7,E15.7,A5,E15.7,E15.7)
11    format(A,F25.16)

      end
