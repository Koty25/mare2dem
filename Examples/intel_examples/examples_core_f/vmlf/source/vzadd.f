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
!    vzAdd  Example Program Text
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
      zB( 1)=(-2.0000000000000000d+004,2.0000000000000000d+004)
      zB( 2)=(-1.5555555399999999d+004,1.5555555400000001d+004)
      zB( 3)=(-1.1111111000000001d+004,1.1111111000000001d+004)
      zB( 4)=(-6.6666665999999996d+003,6.6666666000000005d+003)
      zB( 5)=(-2.2222222000000002d+003,2.2222222000000002d+003)
      zB( 6)=(2.2222222000000002d+003,-2.2222222000000002d+003)
      zB( 7)=(6.6666665999999996d+003,-6.6666666000000005d+003)
      zB( 8)=(1.1111111000000001d+004,-1.1111111000000001d+004)
      zB( 9)=(1.5555555399999999d+004,-1.5555555400000001d+004)
      zB(10)=(2.0000000000000000d+004,-2.0000000000000000d+004)

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZADD(vec_len,zA,zA,zBha0)
      call VZADDI(vec_len,zA_I,inca,zA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZADD(vec_len,zA,zA,zBep1,mode)
      call VMZADDI(vec_len,zA_I,inca,zA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZADD(vec_len,zA,zA,zBep2)
      call VZADDI(vec_len,zA_I,inca,zA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZADD(vec_len,zA,zA,zBla1,mode)
      call VMZADDI(vec_len,zA_I,inca,zA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZADD(vec_len,zA,zA,zBla2)
      call VZADDI(vec_len,zA_I,inca,zA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZADD(vec_len,zA,zA,zBha1,mode)
      call VMZADDI(vec_len,zA_I,inca,zA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZADD(vec_len,zA,zA,zBha2)
      call VZADDI(vec_len,zA_I,inca,zA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZADD and"
              print *," VMZADD in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZADDI and"
              print *," VMZADDI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZADD and"
              print *," VMZADD in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZADDI and"
              print *," VMZADDI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZADD and"
              print *," VMZADD in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZADDI and"
              print *," VMZADDI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZADD and"
              print *," VMZADD in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZADDI and"
              print *," VMZADDI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzAdd test/example program"
      print *,""
      print *,"                    Arguments                         ", &
     &        "      vzAdd"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),zA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzAddI test/example program"
      print *,"                    Arguments                         ", &
     &        "      vzAdd"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,zA_I((i-1)*inca+1),zA_I((i-1)*inca+1),"    ",      &
     &          zBha0_I((i-1)*incb+1)
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

10    format(E12.3,E12.3,E12.3,E12.3,A5,E12.3,E12.3)
11    format(A,F25.16)

      end
