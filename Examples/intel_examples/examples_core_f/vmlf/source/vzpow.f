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
!    vzPow  Example Program Text
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

      zA( 1)=(0.1000d+000,7.0000d+000)
      zA( 2)=(0.8666d+000,6.2333d+000)
      zA( 3)=(1.6333d+000,5.4666d+000)
      zA( 4)=(2.4000d+000,4.6999d+000)
      zA( 5)=(3.1666d+000,3.9333d+000)
      zA( 6)=(3.9333d+000,3.1666d+000)
      zA( 7)=(4.7000d+000,2.3999d+000)
      zA( 8)=(5.4666d+000,1.6333d+000)
      zA( 9)=(6.2333d+000,0.8666d+000)
      zA(10)=(7.0000d+000,0.0999d+000)
      zB( 1)=(7.9222171309252660d-006,2.1083585790345553d-005)
      zB( 2)=(6.4512100825773902d-004,9.1691995410974377d-005)
      zB( 3)=(9.0503657535304545d-003,-1.2801470190500392d-002)
      zB( 4)=(-1.5959437975218965d-001,-2.6566220521977135d-001)
      zB( 5)=(-4.8997205920402100d+000,1.1363610159141824d+000)
      zB( 6)=(4.1042328112091866d+000,6.8100661769434140d+001)
      zB( 7)=(7.9823060676727789d+002,-5.7795808684580606d+001)
      zB( 8)=(-2.3514338397471715d+003,-8.1467677119527843d+003)
      zB( 9)=(-6.5479995707607908d+004,5.3650124957407497d+004)
      zB(10)=(7.8757355451202882d+005,2.3871476187431620d+005)

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZPOW(vec_len,zA,zA,zBha0)
      call VZPOWI(vec_len,zA_I,inca,zA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZPOW(vec_len,zA,zA,zBep1,mode)
      call VMZPOWI(vec_len,zA_I,inca,zA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZPOW(vec_len,zA,zA,zBep2)
      call VZPOWI(vec_len,zA_I,inca,zA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZPOW(vec_len,zA,zA,zBla1,mode)
      call VMZPOWI(vec_len,zA_I,inca,zA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZPOW(vec_len,zA,zA,zBla2)
      call VZPOWI(vec_len,zA_I,inca,zA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZPOW(vec_len,zA,zA,zBha1,mode)
      call VMZPOWI(vec_len,zA_I,inca,zA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZPOW(vec_len,zA,zA,zBha2)
      call VZPOWI(vec_len,zA_I,inca,zA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZPOW and"
              print *," VMZPOW in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZPOWI and"
              print *," VMZPOWI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZPOW and"
              print *," VMZPOW in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZPOWI and"
              print *," VMZPOWI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZPOW and"
              print *," VMZPOW in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZPOWI and"
              print *," VMZPOWI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZPOW and"
              print *," VMZPOW in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZPOWI and"
              print *," VMZPOWI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzPow test/example program"
      print *,""
      print *,"                    Arguments                         ", &
     &        "      vzPow"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),zA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzPowI test/example program"
      print *,"                    Arguments                         ", &
     &        "      vzPow"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,zA_I((i-1)*inca+1),zA_I((i-1)*inca+1),"    ",      &
     &          zBha0_I((i-1)*incb+1)
            CurRMS_I=zrelerr(zB_I((i-1)*incb+1),zBha0_I((i-1)*incb+1))
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

10    format(E12.3,E12.3,E12.3,E12.3,A5,E12.3,E12.3)
11    format(A,F25.16)

      end
