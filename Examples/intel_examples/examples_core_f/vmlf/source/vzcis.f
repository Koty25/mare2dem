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
!    vzCIS  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=8) dA(10)
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
      real(kind=8) dA_I(10*inca)
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

      dA( 1)=-10000.0000d+000
      dA( 2)=-7777.7777d+000
      dA( 3)=-5555.5555d+000
      dA( 4)=-3333.3333d+000
      dA( 5)=-1111.1111d+000
      dA( 6)=1111.1111d+000
      dA( 7)=3333.3333d+000
      dA( 8)=5555.5555d+000
      dA( 9)=7777.7777d+000
      dA(10)=10000.0000d+000
      zB( 1)=(-9.5215536825901481d-001,3.0561438888825215d-001)
      zB( 2)=(6.9259906263181004d-001,7.2132276994528466d-001)
      zB( 3)=(3.4393830299014322d-001,-9.3899224902885137d-001)
      zB( 4)=(-9.9464921859866051d-001,1.0330988307550266d-001)
      zB( 5)=(5.2957287328011915d-001,8.4826444690664649d-001)
      zB( 6)=(5.2957287328011915d-001,-8.4826444690664649d-001)
      zB( 7)=(-9.9464921859866051d-001,-1.0330988307550266d-001)
      zB( 8)=(3.4393830299014322d-001,9.3899224902885137d-001)
      zB( 9)=(6.9259906263181004d-001,-7.2132276994528466d-001)
      zB(10)=(-9.5215536825901481d-001,-3.0561438888825215d-001)

      do i=1,10
          dA_I((i-1)*inca+1)=dA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZCIS(vec_len,dA,zBha0)
      call VZCISI(vec_len,dA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZCIS(vec_len,dA,zBep1,mode)
      call VMZCISI(vec_len,dA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZCIS(vec_len,dA,zBep2)
      call VZCISI(vec_len,dA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZCIS(vec_len,dA,zBla1,mode)
      call VMZCISI(vec_len,dA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZCIS(vec_len,dA,zBla2)
      call VZCISI(vec_len,dA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZCIS(vec_len,dA,zBha1,mode)
      call VMZCISI(vec_len,dA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZCIS(vec_len,dA,zBha2)
      call VZCISI(vec_len,dA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZCIS and"
              print *," VMCCIS in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCISI and"
              print *," VMZCISI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZCIS and"
              print *," VMCCIS in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCISI and"
              print *," VMZCISI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZCIS and"
              print *," VMCCIS in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCISI and"
              print *," VMZCISI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZCIS and"
              print *," VMCCIS in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZCISI and"
              print *," VMZCISI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzCIS test/example program"
      print *,""
      print *,"   Argument                     vzCIS"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,dA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzCISI test/example program"
      print *,"   Argument                     vzCIS"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA_I((i-1)*inca+1),"    ",zBha0_I((i-1)*incb+1)
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

10    format(E15.7,A5,E20.10,E20.10)
11    format(A,F25.16)

      end
