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
!    vzAcosh  Example Program Text
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

      zA( 1)=(3.5000d+000,10000.0000d+000)
      zA( 2)=(1114.2222d+000,8889.2777d+000)
      zA( 3)=(2224.9444d+000,7778.5555d+000)
      zA( 4)=(3335.6666d+000,6667.8333d+000)
      zA( 5)=(4446.3888d+000,5557.1111d+000)
      zA( 6)=(5557.1111d+000,4446.3888d+000)
      zA( 7)=(6667.8333d+000,3335.6666d+000)
      zA( 8)=(7778.5555d+000,2224.9444d+000)
      zA( 9)=(8889.2777d+000,1114.2222d+000)
      zA(10)=(10000.0000d+000,3.5000d+000)
      zB( 1)=(9.9034876162861227d+000,1.5704463268109383d+000)
      zB( 2)=(9.7935428119074004d+000,1.4461021328656807d+000)
      zB( 3)=(9.6915938598109470d+000,1.2921995583383095d+000)
      zB( 4)=(9.6098741976279989d+000,1.1069387861764741d+000)
      zB( 5)=(9.5633904086913990d+000,8.9597859512086708d-001)
      zB( 6)=(9.5633904065260058d+000,6.7481774130490502d-001)
      zB( 7)=(9.6098741922340434d+000,4.6385754781665739d-001)
      zB( 8)=(9.6915938533276709d+000,2.7859677249597142d-001)
      zB( 9)=(9.7935428058704304d+000,1.2469419546677529d-001)
      zB(10)=(9.9034876112861259d+000,3.4999998745833399d-004)

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZACOSH(vec_len,zA,zBha0)
      call VZACOSHI(vec_len,zA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZACOSH(vec_len,zA,zBep1,mode)
      call VMZACOSHI(vec_len,zA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZACOSH(vec_len,zA,zBep2)
      call VZACOSHI(vec_len,zA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZACOSH(vec_len,zA,zBla1,mode)
      call VMZACOSHI(vec_len,zA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZACOSH(vec_len,zA,zBla2)
      call VZACOSHI(vec_len,zA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZACOSH(vec_len,zA,zBha1,mode)
      call VMZACOSHI(vec_len,zA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZACOSH(vec_len,zA,zBha2)
      call VZACOSHI(vec_len,zA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZACOSH and"
              print *," VMZACOSH in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZACOSHI and"
              print *," VMZACOSHI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZACOSH and"
              print *," VMZACOSH in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZACOSHI and"
              print *," VMZACOSHI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZACOSH and"
              print *," VMZACOSH in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZACOSHI and"
              print *," VMZACOSHI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZACOSH and"
              print *," VMZACOSH in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZACOSHI and"
              print *," VMZACOSHI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzAcosh test/example program"
      print *,""
      print *,"           Argument                     vzAcosh"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzAcoshI test/example program"
      print *,"           Argument                     vzAcosh"
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
