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
!    vzCosh  Example Program Text
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

      zA( 1)=(-7.0000d+000,7.0000d+000)
      zA( 2)=(-5.4444d+000,5.4444d+000)
      zA( 3)=(-3.8888d+000,3.8888d+000)
      zA( 4)=(-2.3333d+000,2.3333d+000)
      zA( 5)=(-0.7777d+000,0.7777d+000)
      zA( 6)=(0.7777d+000,-0.7777d+000)
      zA( 7)=(2.3333d+000,-2.3333d+000)
      zA( 8)=(3.8888d+000,-3.8888d+000)
      zA( 9)=(5.4444d+000,-5.4444d+000)
      zA(10)=(7.0000d+000,-7.0000d+000)
      zB( 1)=(4.1337744889835142d+002,-3.6023634485196249d+002)
      zB( 2)=(7.7350991837190122d+001,8.6081439544691477d+001)
      zB( 3)=(-1.7926251549472390d+001,1.6592854640916766d+001)
      zB( 4)=(-3.5948872906678093d+000,-3.6932568898797853d+000)
      zB( 5)=(9.3908576795773324d-001,-6.0235973079599203d-001)
      zB( 6)=(9.3908576795773324d-001,-6.0235973079599203d-001)
      zB( 7)=(-3.5948872906678093d+000,-3.6932568898797853d+000)
      zB( 8)=(-1.7926251549472390d+001,1.6592854640916766d+001)
      zB( 9)=(7.7350991837190122d+001,8.6081439544691477d+001)
      zB(10)=(4.1337744889835142d+002,-3.6023634485196249d+002)

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZCOSH(vec_len,zA,zBha0)
      call VZCOSHI(vec_len,zA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZCOSH(vec_len,zA,zBep1,mode)
      call VMZCOSHI(vec_len,zA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZCOSH(vec_len,zA,zBep2)
      call VZCOSHI(vec_len,zA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZCOSH(vec_len,zA,zBla1,mode)
      call VMZCOSHI(vec_len,zA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZCOSH(vec_len,zA,zBla2)
      call VZCOSHI(vec_len,zA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZCOSH(vec_len,zA,zBha1,mode)
      call VMZCOSHI(vec_len,zA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZCOSH(vec_len,zA,zBha2)
      call VZCOSHI(vec_len,zA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZCOSH and"
              print *," VMZCOSH in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCOSHI and"
              print *," VMZCOSHI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZCOSH and"
              print *," VMZCOSH in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCOSHI and"
              print *," VMZCOSHI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZCOSH and"
              print *," VMZCOSH in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCOSHI and"
              print *," VMZCOSHI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZCOSH and"
              print *," VMZCOSH in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCOSHI and"
              print *," VMZCOSHI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzCosh test/example program"
      print *,""
      print *,"           Argument                     vzCosh"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzCoshI test/example program"
      print *,"           Argument                     vzCosh"
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
